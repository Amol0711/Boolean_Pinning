#!/usr/bin/env python3
"""Compare a completed J1 replay with the fixed 20,000-draw release.

This read-only scientific comparison does not regenerate networks or substitute
for the independent dynamical verifier. The optional --output writes only its
JSON report. Production source, protocols and archived results are not edited.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import zlib

from representation import load_representation, protocol_reference

ROOT = Path(__file__).resolve().parent
SCIENTIFIC_FILES = (
    'raw/J1_n25600.csv.gz',
    'details/J1_n25600.jsonl.gz',
    'summary.json',
    'summary.csv',
    'rule_counts.json',
)
DRAW_COUNT = 20000


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference_identity():
    protocol = read_json(ROOT / 'Protocol_J1_Replication.json')
    receipt, current = load_representation(ROOT)
    require(protocol['total_draws'] == protocol['replicates_per_cell'] == DRAW_COUNT,
            'Fixed protocol does not describe exactly 20000 draws')
    require(protocol['sizes'] == [25600] and list(protocol['families']) == ['J1'],
            'Fixed family/size mismatch')
    prefix = [protocol['master_seed'], protocol['families']['J1']['id'], 25600]
    require(prefix == [2026101601, 3101, 25600], 'Unexpected replication seed namespace')
    with (ROOT/'seed_manifest.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == DRAW_COUNT, 'Seed manifest row count differs from 20000')
    for i, row in enumerate(rows):
        actual = [int(row[k]) for k in ('master_seed','family_id','n','replicate')]
        require(actual == prefix+[i], 'Seed manifest identity mismatch at '+str(i))
    return prefix, receipt, current


def inspect_run(directory, prefix, current):
    manifest = read_json(directory/'manifest.json')
    require(manifest['status'] == 'completed' and manifest['total_draws'] == DRAW_COUNT,
            'Run is not completed with 20000 draws: '+str(directory))
    expected_links = {
        'protocol_sha256': current['protocol_sha256'],
        'integrity_record_sha256': digest((ROOT/'integrity_record.json').read_bytes()),
        'independent_gate_sha256': current['independent_preflight_sha256'],
        'source_sha256': current['source_sha256'],
    }
    for key, expected in expected_links.items():
        require(manifest.get(key) == expected, 'Manifest provenance mismatch: '+key)
    integrity = read_json(directory/'io_integrity.json')
    require(integrity['status'] == 'PASS', 'Run integrity gate is not PASS')
    streams = integrity['streams']
    expected_streams = set(SCIENTIFIC_FILES[:2])
    require(len(streams) == 2 and {x['path'] for x in streams} == expected_streams,
            'Integrity stream list must contain exactly the raw and detail files')
    integrity_by_path = {x['path']:x for x in streams}
    payloads = {}
    file_info = {}
    for name in SCIENTIFIC_FILES:
        encoded = (directory/name).read_bytes()
        decoded = gzip.decompress(encoded) if name.endswith('.gz') else encoded
        payloads[name] = decoded
        file_info[name] = {'sha256':digest(encoded),'payload_sha256':digest(decoded)}
        if name in expected_streams:
            entry = integrity_by_path[name]
            require(entry['records'] == DRAW_COUNT and entry['sha256'] == digest(encoded),
                    'Integrity count/hash mismatch: '+name)
    raw = list(csv.DictReader(io.StringIO(payloads[SCIENTIFIC_FILES[0]].decode(), newline='')))
    details = [json.loads(line) for line in payloads[SCIENTIFIC_FILES[1]].decode().splitlines()]
    require(len(raw) == len(details) == DRAW_COUNT, 'Raw/detail row count differs from 20000')
    seen_seeds = set()
    identities = hashlib.sha256()
    aggregate = {'1':[0]*8, '3':[0]*8}
    unresolved = 0
    for i, (row, detail) in enumerate(zip(raw, details)):
        seed = json.loads(row['seed_entropy'])
        require(int(row['replicate']) == detail['replicate'] == i,
                'Missing, duplicated or reordered replicate at '+str(i))
        require(seed == detail['seed_entropy'] == prefix+[i], 'Seed identity mismatch at '+str(i))
        require(tuple(seed) not in seen_seeds, 'Duplicate seed at '+str(i))
        seen_seeds.add(tuple(seed))
        require(row['family'] == detail['family'] == 'J1' and int(row['n']) == detail['n'] == 25600,
                'Family/size mismatch at '+str(i))
        network_hash = row['network_sha256']
        require(re.fullmatch('[0-9a-f]{64}', network_hash) is not None and
                network_hash == detail['network_sha256'], 'Raw/detail network hash mismatch at '+str(i))
        lower, upper = int(row['cost_lower']), int(row['cost_upper'])
        require(0 <= lower <= upper and lower == detail['lower'] and upper == detail['upper'],
                'Raw/detail cost bounds mismatch at '+str(i))
        unresolved += lower != upper
        for degree in (1,3):
            counts = json.loads(row[f'rule_counts_d{degree}'])
            require(counts == detail[f'rule_counts_d{degree}'] and len(counts) == 8 and
                    all(isinstance(x,int) and x >= 0 for x in counts) and sum(counts) == 12800,
                    'Raw/detail rule-count mismatch at '+str(i))
            aggregate[str(degree)] = [a+b for a,b in zip(aggregate[str(degree)], counts)]
        identities.update((json.dumps(seed,separators=(',',':'))+' '+network_hash+'\n').encode())
    rule_counts = json.loads(payloads['rule_counts.json'])
    require(rule_counts == {'order':'corner_code0..7','draws':DRAW_COUNT,
                           'rows_per_degree':256000000,'counts':aggregate},
            'Aggregate rule counts do not agree with every raw record')
    summary = json.loads(payloads['summary.json'])
    require(len(summary) == 1 and summary[0]['family'] == 'J1' and summary[0]['n'] == 25600 and
            summary[0]['draws'] == DRAW_COUNT and summary[0]['unresolved'] == unresolved and
            manifest['unresolved'] == unresolved, 'Summary/completion count mismatch')
    return {'files':file_info,'raw_records':len(raw),'detail_records':len(details),
            'unique_seed_identities':len(seen_seeds),'ordered_identity_sha256':identities.hexdigest(),
            'unresolved':unresolved,'completion_integrity_and_provenance_verified':True}


def verify(candidate):
    reference = ROOT/'results/main'
    prefix, receipt, current = reference_identity()
    archive_manifest = read_json(ROOT/'SHA256_MANIFEST.json')['files']
    for name in SCIENTIFIC_FILES:
        rel = 'results/main/'+name
        require(archive_manifest.get(rel) == digest((reference/name).read_bytes()),
                'Recorded scientific file differs from current distribution inventory: '+rel)
    expected = inspect_run(reference, prefix, current)
    actual = inspect_run(candidate, prefix, current)
    comparisons = []
    for name in SCIENTIFIC_FILES:
        a,b = actual['files'][name], expected['files'][name]
        comparisons.append({'path':name,'bytes_identical':a['sha256']==b['sha256'],
                            'decoded_payload_identical':a['payload_sha256']==b['payload_sha256'],
                            'candidate_sha256':a['sha256'],'recorded_sha256':b['sha256']})
    passed = all(x['bytes_identical'] and x['decoded_payload_identical'] for x in comparisons)
    return {'status':'PASS' if passed else 'FAIL','draws':DRAW_COUNT,'scientific_files':comparisons,
            'candidate_identity_check':{k:v for k,v in actual.items() if k!='files'},
            'recorded_identity_check':{k:v for k,v in expected.items() if k!='files'},
            'protocol_sha256':current['protocol_sha256'],
            'protocol_record':protocol_reference(receipt),
            'integrity_note':'Current protocol and source links are checked for both runs.',
            'comparison_scope':'Fixed scientific payloads and full ordered identities; no network regeneration or independent dynamical optimization.',
            'not_byte_compared':['run_started.json','resources.json','manifest.json','io_integrity.json',
                                 'shards/','statistical_summary.csv','rule_frequency_summary.csv'],
            'exclusion_note':'Completion/integrity/provenance fields are checked separately. Statistical CSVs are generated by the independent statistics task, not the production runner.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,required=True,help='Completed fresh J1 results directory')
    parser.add_argument('--output',type=Path,help='Optional new JSON report (must not exist); stdout receives the same report')
    args = parser.parse_args()
    output = None
    try:
        if args.output:
            proposed = args.output.resolve()
            require(not proposed.exists(), 'Report path already exists; refusing overwrite: '+str(proposed))
            require(not proposed.is_relative_to((ROOT/'results/main').resolve()),
                    'Report path may not be inside recorded results/main')
            output = proposed
        result = verify(args.candidate.resolve())
    except (OSError, ValueError, KeyError, TypeError, EOFError, zlib.error) as error:
        result = {'status':'FAIL','error_type':type(error).__name__,'error':str(error)}
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if output:
        output.parent.mkdir(parents=True,exist_ok=True)
        # Exclusive creation also protects against a file appearing during the check.
        with output.open('x') as stream:
            stream.write(text)
    print(text,end='')
    return 0 if result['status']=='PASS' else 1


if __name__=='__main__':
    sys.exit(main())
