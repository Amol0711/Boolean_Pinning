"""Replay the frozen Bernoulli sample and compare every scientific stream.

Run from any working directory. A default temporary output is retained so the
comparison report and regenerated data remain inspectable after completion.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def inspect(path, expected_count):
    compressed = path.read_bytes()
    payload = gzip.decompress(compressed)
    text = payload.decode('utf-8')
    records = (csv.DictReader(io.StringIO(text, newline=''))
               if path.name.endswith('.csv.gz')
               else (json.loads(line) for line in text.splitlines()))
    count = 0
    for count, row in enumerate(records, 1):
        if int(row['replicate']) != count - 1:
            raise ValueError(f'Nonconsecutive record IDs in {path}')
    if count != expected_count:
        raise ValueError(f'{path}: expected {expected_count} records, read {count}')
    return {'records': count, 'gzip_sha256': sha(compressed), 'payload_sha256': sha(payload)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, help='New replay output directory; must not exist')
    parser.add_argument('--pilot', action='store_true', help='First 100 original replicates in all 25 cells')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    if args.out is not None and args.out.exists():
        parser.error('--out must not already exist')
    baseline = ROOT/'results'/('pilot_10m' if args.pilot else 'main')
    count = 100 if args.pilot else 10000
    protocol = json.loads((ROOT/'Protocol_Bernoulli.json').read_text())
    relative = [Path(kind)/f'{mixture}_n{n}.{suffix}'
                for mixture in 'ABCDE' for n in protocol['sizes']
                for kind, suffix in [('raw','csv.gz'),('details','jsonl.gz')]]
    expected = {str(path): inspect(baseline/path, count) for path in relative}
    out = (args.out.resolve() if args.out is not None else
           Path(tempfile.mkdtemp(prefix='boolean_bernoulli_reproduction_'))/'replay')
    command = [sys.executable, str(ROOT/'src/replay_scc.py'), '--workers', str(args.workers),
               '--work-cap', str(protocol['component_caps']['max_work']), '--output', str(out)]
    if args.pilot:
        command.append('--pilot')
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise SystemExit(completed.returncode)
    manifest = json.loads((out/'run_manifest.json').read_text())
    integrity = json.loads((out/'io_integrity.json').read_text())
    if manifest['status'] != 'completed' or integrity['status'] != 'PASS':
        raise RuntimeError('Replay did not pass the completion integrity gate')
    comparisons = []
    for path in relative:
        actual = inspect(out/path, count)
        reference = expected[str(path)]
        comparisons.append({'path': str(path), 'records': count,
                            'gzip_byte_identical': actual['gzip_sha256'] == reference['gzip_sha256'],
                            'decoded_payload_identical': actual['payload_sha256'] == reference['payload_sha256'],
                            'reference': reference, 'replayed': actual})
    summary_equal = all((baseline/name).read_bytes() == (out/name).read_bytes()
                        for name in ['operational_summary.json','operational_summary.csv'])
    passed = summary_equal and all(r['gzip_byte_identical'] and r['decoded_payload_identical'] for r in comparisons)
    report = {'status': 'PASS' if passed else 'FAIL', 'pilot': args.pilot,
              'draws': 25*count, 'scientific_streams': len(comparisons),
              'operational_summaries_byte_identical': summary_equal,
              'timing_and_backend_metadata_compared': False,
              'protocol_sha256': sha((ROOT/'Protocol_Bernoulli.json').read_bytes()),
              'streams': comparisons}
    report_path = out/'reproduction_comparison.json'
    temporary = report_path.with_name(report_path.name+'.crdownload')
    temporary.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    temporary.replace(report_path)
    print(json.dumps({'status':report['status'],'draws':25*count,
                      'scientific_streams':len(comparisons),'report':str(report_path)}),flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
