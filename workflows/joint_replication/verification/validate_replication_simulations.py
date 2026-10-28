#!/usr/bin/env python3
"""Independently replay every J1 replication network and validate physical certificates.

The generator is independently reconstructed from the fixed protocol. No
production module is imported. Support uses synchronous set iteration, SCCs
use Tarjan, and physical pin success uses functional-graph indegree peeling.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from representation import load_representation, protocol_reference

from independent_network_oracle import (
    reconstruct, literal_support, export_rows, graph_certificate, exact_cost,
    witness_success, normalized_key,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_chunk(job):
    protocol, records, full_residual_cap = job
    counts = Counter()
    cache = {}
    complex_sizes = Counter()
    reasons = Counter()
    rule_counts = {degree: np.zeros(8, dtype=np.int64) for degree in (1, 3)}
    incomplete = []
    max_iterations = max_component = 0
    repeated_residuals = []
    noncopy_cases = []
    for record, detail in records:
        ident = (record['family'], int(record['n']), int(record['replicate']))
        family, n, rep = ident
        assert (detail['family'], detail['n'], detail['replicate']) == ident
        entropy = [protocol['master_seed'], protocol['families'][family]['id'], n, rep]
        assert json.loads(record['seed_entropy']) == detail['seed_entropy'] == entropy
        p, b, digest = reconstruct(protocol, *ident)
        assert record['network_sha256'] == detail['network_sha256'] == digest, ident
        expected_degrees = np.repeat([1, 3], n // 2)
        assert np.array_equal(np.bincount(p.ravel(), minlength=n), expected_degrees)
        codes = b[:, 0].astype(np.int64) + 2*b[:, 1] + 4*b[:, 2]
        for degree, sl in [(1, slice(0, n//2)), (3, slice(n//2, n))]:
            values = np.bincount(codes[sl], minlength=8)
            field = f'rule_counts_d{degree}'
            assert values.tolist() == json.loads(record[field]) == detail[field], (ident, field)
            assert sum(values) == n//2
            assert not any(values[j] for j in [1, 2, 3, 7])
            rule_counts[degree] += values
        residual, iterations = literal_support(p, b)
        rows = export_rows(p, b, residual)
        assert rows == detail['residual_rows'], ident
        counts['full_seeded_networks_reconstructed'] += 1
        counts['source_coordinates_reconstructed'] += n
        counts['full_source_support_recomputed'] += 1
        counts['degree_counts_checked'] += 1
        counts['degree_conditioned_rule_count_vectors_checked'] += 2
        counts['residual_rows_compared'] += len(rows)
        counts['physical_repeated_parent_rows_full_network'] += int(np.count_nonzero(p[:, 0] == p[:, 1]))
        repeats = [v for v, a, c, *_ in rows if a == c]
        counts['physical_repeated_parent_rows_residual'] += len(repeats)
        if repeats:
            repeated_residuals.append(dict(replicate=rep, vertices=repeats))
        max_iterations = max(max_iterations, iterations)
        graph = graph_certificate(rows)
        assert len(rows) == int(record['residual_size'])
        for key in ['G', 'G_prime']:
            assert bool(int(record[key])) == detail[key] == graph[key], (ident, key)
        assert graph['double_rows'] == int(record['double_rows'])
        assert graph['cyclic_component_count'] == int(record['cyclic_components'])
        assert graph['largest_complex'] == int(record['largest_complex'])
        sizes = sorted(len(c['nodes']) for c in graph['components'])
        cyclic_sizes = sorted(len(c['nodes']) for c in graph['components'] if c['cyclic'])
        for field, expected in [('component_sizes', sizes), ('cyclic_component_sizes', cyclic_sizes)]:
            assert json.loads(record[field]) == detail[field] == expected, (ident, field)
        cc = graph['cyclic_component_count'] if graph['G_prime'] else None
        assert detail['certificate_cycle_count'] == cc
        assert record['certificate_cycle_count'] == (str(cc) if cc is not None else '')
        counts['G'] += graph['G']
        counts['G_prime'] += graph['G_prime']
        counts['G_false_G_prime_true'] += not graph['G'] and graph['G_prime']
        lower, upper = int(record['cost_lower']), int(record['cost_upper'])
        assert (lower, upper) == (detail['lower'], detail['upper'])
        assert 0 <= lower <= upper <= len(rows)
        pins = detail['pins']
        assert pins == ([] if not record['pins'] else list(map(int, record['pins'].split(';'))))
        assert pins == sorted(set(pins)) and len(pins) == upper
        assert set(pins) <= set(residual)
        counts['exact_draws'] += lower == upper
        counts['interval_draws'] += lower != upper
        expected_reason = ('original_G' if graph['G'] else 'weaker_G_prime' if graph['G_prime']
                           else 'SCC_exact' if lower == upper else 'bounded_unresolved')
        assert record['exact_reason'] == detail['exact_reason'] == expected_reason
        expected_components = graph['components']
        produced = {}
        if detail['components']:
            mapping = detail['local_to_original']
            assert mapping == residual + [None]
            for c in detail['components']:
                nodes = tuple(mapping[v] for v in c['nodes'])
                assert None not in nodes and nodes not in produced
                produced[nodes] = c
            assert set(produced) == {tuple(c['nodes']) for c in expected_components}
        else:
            assert graph['G'] and detail['components'] is None
            assert detail['local_to_original'] is None
            mapping = None
        exact_sum = scc_lower = scc_upper = work = 0
        all_known = True
        stops = Counter()
        for component in expected_components:
            nodes = component['nodes']
            counts['SCC_records_checked'] += 1
            max_component = max(max_component, len(nodes))
            c = produced.get(tuple(nodes))
            if c:
                assert c['cyclic'] == component['cyclic']
                assert c['marked_slot_cycle'] == component['marked_slot_cycle']
                slots = sorted((mapping[x['source']], mapping[x['receiver']], x['slot'])
                               for x in c['internal_slots'])
                assert slots == sorted(component['internal_slots'])
                c_lo, c_hi = c['lower'], c['upper']
                local_pins = [mapping[x] for x in c['pins']]
                reasons[c['reason']] += 1
                work += c['work']
                assert c['exact_cost'] == (c_lo if c_lo == c_hi else None)
                if c['exact_cost'] is None:
                    stops[c['reason']] += 1
                if 'limits' in c:
                    assert c['limits'] == protocol['component_caps']
                    assert 0 <= c['work'] <= c['limits']['max_work']
                    assert c['work'] == c['row_evaluations'] + c['successor_traversals']
                    if c['reason'] == 'component_size_cap':
                        assert len(nodes) > c['limits']['max_component_size']
                    if c['reason'] == 'work_cap_successor_construction':
                        assert len(nodes) * (1 << len(nodes)) > c['limits']['max_work']
            else:
                c_lo = c_hi = int(component['cyclic'])
                local_pins = [min(nodes)] if component['cyclic'] else []
            assert 0 <= c_lo <= c_hi <= len(nodes)
            assert len(set(local_pins)) == len(local_pins) == c_hi
            assert set(local_pins) <= set(nodes)
            scc_lower += c_lo
            scc_upper += c_hi
            if not component['cyclic']:
                assert c_lo == c_hi == 0
                counts['acyclic_components_symbolically_verified'] += 1
            elif component['marked_slot_cycle']:
                assert c_lo == c_hi == 1
                assert len(set(pins) & set(nodes)) >= 1
                exact_sum += 1
                counts['copy_cycle_components_symbolically_verified'] += 1
            else:
                counts['complex_components_checked'] += 1
                complex_sizes[len(nodes)] += 1
                key = normalized_key(rows, nodes)
                if key not in cache:
                    cache[key] = exact_cost(rows, nodes, max_compatible_states=100_000_000)
                    counts['chunk_distinct_complex_oracle_problems'] += 1
                    counts['chunk_distinct_complex_pin_candidates'] += cache[key]['candidates']
                    counts['chunk_distinct_complex_compatible_states'] += cache[key]['states']
                oracle = cache[key]
                if oracle['complete']:
                    assert c_lo <= oracle['cost'] <= c_hi, (ident, nodes, c, oracle)
                    if c and c['exact_cost'] is not None:
                        assert c['exact_cost'] == oracle['cost']
                    exact_sum += oracle['cost']
                    counts['complex_component_optima_exhaustively_verified'] += 1
                else:
                    all_known = False
                    assert c_lo == 0 or c_lo <= oracle.get('lower', 0)
                    incomplete.append(dict(replicate=rep, nodes=nodes, reason=oracle['reason'],
                                           production_interval=[c_lo, c_hi]))
                    counts['complex_component_oracle_guards'] += 1
                noncopy_cases.append(dict(replicate=rep, nodes=nodes,
                                          production_interval=[c_lo, c_hi], oracle=oracle))
                for label, witness in [('component', local_pins), ('global', pins)]:
                    result = witness_success(rows, nodes, witness)
                    assert result['checked'] and result['success'], (ident, nodes, label, result)
                    counts[label + '_complex_witnesses_verified'] += 1
                    counts[label + '_complex_witness_states'] += result['states']
                    counts[label + '_complex_witness_' + result['method']] += 1
        assert (scc_lower, scc_upper) == (lower, upper)
        assert work == int(record['solver_work'])
        assert dict(stops) == json.loads(record['stop_summary'])
        if all_known:
            assert lower <= exact_sum <= upper
            counts['draw_optima_known_from_independent_components'] += 1
        if not graph['G'] and len(rows) <= full_residual_cap:
            oracle = exact_cost(rows, residual, max_size=full_residual_cap)
            assert oracle['complete']
            assert lower <= oracle['cost'] <= upper
            if all_known:
                assert oracle['cost'] == exact_sum
            counts['full_residual_global_optima_enumerated'] += 1
            counts['full_residual_global_pin_candidates'] += oracle['candidates']
            counts['full_residual_global_compatible_states'] += oracle['states']
        counts['final_witnesses_verified_by_components'] += 1
    return dict(first_replicate=int(records[0][0]['replicate']),
                last_replicate=int(records[-1][0]['replicate']), draws=len(records),
                counts=dict(counts), maximum_literal_support_iterations=max_iterations,
                maximum_SCC_size=max_component, complex_component_size_counts=dict(complex_sizes),
                component_reason_counts=dict(reasons), independent_optimum_guards=incomplete,
                degree_conditioned_rule_counts={str(d): v.tolist() for d, v in rule_counts.items()},
                repeated_parent_residual_cases=repeated_residuals, noncopy_component_cases=noncopy_cases)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=root/'Protocol_J1_Replication.json')
    parser.add_argument('--results', type=Path, default=root/'results/main')
    parser.add_argument('--output', type=Path, default=root/'qa/independent_network_replay.json')
    parser.add_argument('--preflight', type=Path, default=root/'qa/independent_preflight.json')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--chunk-size', type=int, default=1000)
    parser.add_argument('--full-residual-cap', type=int, default=12)
    args = parser.parse_args()
    assert args.workers >= 1 and args.chunk_size >= 1
    protocol = json.loads(args.protocol.read_text())
    assert protocol['master_seed'] == 2026101601 and protocol['sizes'] == [25600]
    assert set(protocol['families']) == {'J1'} and protocol['families']['J1']['id'] == 3101
    assert protocol['replicates_per_cell'] == 20000
    assert protocol['component_caps'] == dict(max_component_size=18, max_pins=18, max_work=10000000)
    expected_types = [(1, [30, 0, 0, 0, 15, 45, 10, 0]),
                      (3, [60, 0, 0, 0, 10, 10, 20, 0])]
    for actual, (degree, weights) in zip(protocol['families']['J1']['types'], expected_types):
        assert actual['degree'] == degree and actual['fraction'] == [1, 2]
        assert actual['weights'] == weights and actual['probability_denominator'] == 100
    assert len(protocol['families']['J1']['types']) == 2
    preflight = json.loads(args.preflight.read_text())
    assert preflight['status'] == 'PASS'
    oracle_path = Path(__file__).with_name('independent_network_oracle.py')
    assert preflight['source_sha256'][oracle_path.name] == sha(oracle_path)
    manifest_path = args.results/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    assert manifest['status'] == 'completed'
    assert manifest['protocol_sha256'] == sha(args.protocol)
    assert manifest['total_draws'] == 20000
    assert manifest['independent_gate_sha256'] == sha(args.preflight)
    freeze_path = root/'integrity_record.json'
    freeze, current = load_representation(root)
    assert manifest['integrity_record_sha256'] == sha(freeze_path)
    assert freeze['schema'] == 1
    assert current['protocol_sha256'] == sha(args.protocol)
    assert current['independent_preflight_sha256'] == sha(args.preflight)
    assert current['independent_oracle_sha256'] == sha(oracle_path)
    seed_path = root/'seed_manifest.csv'
    assert current['seed_manifest_sha256'] == sha(seed_path)
    with seed_path.open(newline='') as seed_stream:
        seeds = list(csv.DictReader(seed_stream))
    assert len(seeds) == 20000
    assert [[int(row[k]) for k in ['replicate', 'master_seed', 'family_id', 'n']]
            for row in seeds] == [[r, 2026101601, 3101, 25600] for r in range(20000)]
    source_hashes = manifest.get('source_sha256', {})
    assert source_hashes, 'Production source hashes missing'
    assert source_hashes == current['source_sha256']
    for relative, digest in source_hashes.items():
        assert sha(root/relative) == digest, relative
    raw_path = args.results/'raw/J1_n25600.csv.gz'
    detail_path = args.results/'details/J1_n25600.jsonl.gz'
    integrity_path = args.results/'io_integrity.json'
    integrity = json.loads(integrity_path.read_text())
    assert integrity['status'] == 'PASS'
    stream_hashes = {str(p.relative_to(args.results)): sha(p) for p in [raw_path, detail_path]}
    for relative, digest in stream_hashes.items():
        matches = [r for r in integrity['streams'] if r['path'] == relative]
        assert len(matches) == 1
        assert matches[0]['records'] == 20000 and matches[0]['sha256'] == digest
    records = []
    with gzip.open(raw_path, 'rt', newline='') as raws, gzip.open(detail_path, 'rt') as details:
        for rep, (record, line) in enumerate(itertools.zip_longest(csv.DictReader(raws), details)):
            assert record is not None and line is not None
            detail = json.loads(line)
            assert (record['family'], int(record['n']), int(record['replicate'])) == ('J1', 25600, rep)
            records.append((record, detail))
    assert len(records) == 20000
    chunks = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(validate_chunk, (protocol, records[j:j+args.chunk_size], args.full_residual_cap))
                   for j in range(0, len(records), args.chunk_size)]
        for future in as_completed(futures):
            chunk = future.result()
            chunks.append(chunk)
            print(f"verified replicates {chunk['first_replicate']}--{chunk['last_replicate']}", flush=True)
    chunks.sort(key=lambda c: c['first_replicate'])
    totals, sizes, reasons = Counter(), Counter(), Counter()
    degree_counts = {d: np.zeros(8, dtype=np.int64) for d in ('1', '3')}
    for chunk in chunks:
        totals.update(chunk['counts'])
        sizes.update(chunk['complex_component_size_counts'])
        reasons.update(chunk['component_reason_counts'])
        for d in degree_counts:
            degree_counts[d] += chunk['degree_conditioned_rule_counts'][d]
    assert totals['full_seeded_networks_reconstructed'] == 20000
    assert totals['final_witnesses_verified_by_components'] == 20000
    assert sha(raw_path) == stream_hashes[str(raw_path.relative_to(args.results))]
    assert sha(detail_path) == stream_hashes[str(detail_path.relative_to(args.results))]
    result = dict(schema_version=1, status='PASS', draws=20000, n=25600,
                  protocol_sha256=sha(args.protocol), preflight_sha256=sha(args.preflight),
                  integrity_record_sha256=sha(freeze_path), seed_manifest_sha256=sha(seed_path),
                  production_manifest_sha256=sha(manifest_path),
                  production_io_integrity_sha256=sha(integrity_path), stream_sha256=stream_hashes,
                  validation_source_sha256={p.name: sha(p) for p in [Path(__file__), oracle_path]},
                  protocol_record=protocol_reference(freeze),
                  integrity_note='Validation performed with the distributed file identities.',
                  coverage=dict(full_source_replay=True, full_residual_global_optimum_dimension_cap=args.full_residual_cap,
                                complex_component_oracle_dimension_cap=18,
                                complex_component_oracle_candidate_cap=100000,
                                complex_component_oracle_compatible_state_cap=100000000),
                  independence='No production imports. Independent fixed PCG64 reconstruction, synchronous nominal support iteration, Tarjan SCCs, physical successor construction and indegree-peeling pin success; copy and acyclic components use exact symbolic certificates.',
                  counts=dict(sorted(totals.items())),
                  degree_conditioned_rule_counts={d: c.tolist() for d, c in degree_counts.items()},
                  complex_component_size_counts=dict(sorted(sizes.items())),
                  component_reason_counts=dict(sorted(reasons.items())),
                  independent_optimum_guards=[g for c in chunks for g in c['independent_optimum_guards']],
                  unverified_witnesses=[], chunks=chunks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'output': str(args.output), 'counts': result['counts']}), flush=True)


if __name__ == '__main__':
    main()
