#!/usr/bin/env python3
"""Replay recorded JointRules networks and independently validate every exported certificate.

Run from the archive root:
  python verification/validate_joint_simulations.py --results validation/results/main
The default reads validation/results/main; --results may select a bounded pilot.
Only the independent_network_oracle module is imported for reference answers.
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

from independent_network_oracle import (reconstruct, literal_support, export_rows,
    graph_certificate, exact_cost, witness_success, normalized_key,
    self_check)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with gzip.open(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream))


def worker(args):
    root, result_dir, family, n, expected_rows, full_residual_cap = args
    root, result_dir = Path(root), Path(result_dir)
    protocol = json.loads((root / 'validation/Protocol_JointRules.json').read_text())
    stem = f'{family}_n{n}'
    raw_path = result_dir / 'raw' / (stem + '.csv.gz')
    detail_path = result_dir / 'details' / (stem + '.jsonl.gz')
    new_rows = read_csv(raw_path)
    assert len(new_rows) == expected_rows, (stem, len(new_rows), expected_rows)
    assert [int(r['replicate']) for r in new_rows] == list(range(expected_rows))
    counts = Counter()
    max_iterations = max_component = 0
    oracle_cache = {}
    unvalidated = []
    incomplete_oracles = []
    complex_size_counts = Counter()
    bound_reason_counts = Counter()
    with gzip.open(detail_path, 'rt') as details:
        for expected_replicate, (record, line) in enumerate(itertools.zip_longest(new_rows, details)):
            assert record is not None and line is not None, 'detail/raw length mismatch'
            detail = json.loads(line)
            ident = (family, n, expected_replicate)
            assert (record['family'], int(record['n']), int(record['replicate'])) == ident
            assert (detail['family'], detail['n'], detail['replicate']) == ident
            p, b, digest = reconstruct(protocol, *ident)
            assert record['network_sha256'] == detail['network_sha256'] == digest, ident
            residual, iterations = literal_support(p, b)
            rows = export_rows(p, b, residual)
            assert rows == detail['residual_rows'], ident
            counts['full_seeded_networks_reconstructed'] += 1
            counts['full_source_support_recomputed'] += 1
            counts['source_coordinates_reconstructed'] += n
            counts['residual_rows_compared'] += len(rows)
            max_iterations = max(max_iterations, iterations)
            certificate = graph_certificate(rows)
            assert len(rows) == int(record['residual_size'])
            for key in ['G', 'G_prime']:
                assert bool(int(record[key])) == detail[key] == certificate[key], (ident, key)
            assert certificate['double_rows'] == int(record['double_rows'])
            assert certificate['cyclic_component_count'] == int(record['cyclic_components'])
            assert certificate['largest_complex'] == int(record['largest_complex'])
            cc = certificate['cyclic_component_count'] if certificate['G_prime'] else None
            assert detail['certificate_cycle_count'] == cc
            assert record['certificate_cycle_count'] == (str(cc) if cc is not None else '')
            counts['G'] += certificate['G']
            counts['G_prime'] += certificate['G_prime']
            counts['G_false_G_prime_true'] += not certificate['G'] and certificate['G_prime']
            lower, upper = int(record['cost_lower']), int(record['cost_upper'])
            assert (lower, upper) == (detail['lower'], detail['upper'])
            assert 0 <= lower <= upper <= len(rows), ident
            pins = detail['pins']
            assert pins == ([] if not record['pins'] else list(map(int, record['pins'].split(';'))))
            assert len(pins) == len(set(pins)) == upper
            assert set(pins) <= set(residual)
            counts['exact_draws'] += lower == upper
            counts['interval_draws'] += lower != upper
            # Compare SCC identities and certificates with explicit index mapping.
            expected_components = certificate['components']
            if detail['components']:
                mapping = detail['local_to_original']
                produced = {}
                for c in detail['components']:
                    nodes = tuple(mapping[v] for v in c['nodes'])
                    assert None not in nodes
                    produced[nodes] = c
                assert set(produced) == {tuple(c['nodes']) for c in expected_components}, ident
            else:
                assert certificate['G']
                produced = {}
                mapping = None
            exact_sum = 0
            all_component_optima_known = True
            scc_lower = scc_upper = 0
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
                    bound_reason_counts[c['reason']] += 1
                    assert c['exact_cost'] == (c_lo if c_lo == c_hi else None)
                    if 'limits' in c:
                        assert 0 <= c['work'] <= c['limits']['max_work']
                        assert c['work'] == c['row_evaluations'] + c['successor_traversals']
                        if c['reason'] == 'component_size_cap':
                            assert len(nodes) > c['limits']['max_component_size']
                        if c['reason'] == 'work_cap_successor_construction':
                            assert len(nodes) * (1 << len(nodes)) > c['limits']['max_work']
                else:
                    c_lo = c_hi = int(component['cyclic'])
                    local_pins = [min(nodes)] if component['cyclic'] else []
                scc_lower += c_lo
                scc_upper += c_hi
                assert len(local_pins) == c_hi
                if not component['cyclic']:
                    assert c_lo == c_hi == 0
                    exact_sum += 0
                    counts['acyclic_components_symbolically_verified'] += 1
                elif component['marked_slot_cycle']:
                    assert c_lo == c_hi == 1
                    assert len(set(local_pins) & set(nodes)) == 1
                    assert set(pins) & set(nodes), (ident, 'global witness misses copy cycle', nodes)
                    exact_sum += 1
                    counts['copy_cycle_components_symbolically_verified'] += 1
                else:
                    counts['complex_components_checked'] += 1
                    complex_size_counts[len(nodes)] += 1
                    key = normalized_key(rows, nodes)
                    if key not in oracle_cache:
                        oracle_cache[key] = exact_cost(rows, nodes, max_compatible_states=100_000_000)
                        counts['distinct_complex_oracle_problems'] += 1
                        counts['distinct_complex_pin_candidates'] += oracle_cache[key]['candidates']
                        counts['distinct_complex_compatible_states'] += oracle_cache[key]['states']
                    oracle = oracle_cache[key]
                    if oracle['complete']:
                        assert c_lo <= oracle['cost'] <= c_hi, (ident, nodes, c, oracle)
                        if c and c['exact_cost'] is not None:
                            assert oracle['cost'] == c['exact_cost'], (ident, nodes)
                        exact_sum += oracle['cost']
                        counts['complex_component_optima_exhaustively_verified'] += 1
                    else:
                        all_component_optima_known = False
                        # Incomplete production results may still carry simple valid
                        # bounds; no exhaustive optimum is claimed for this branch.
                        assert c_lo == 0 or c_lo <= oracle.get('lower', 0)
                        incomplete_oracles.append({'id': list(ident), 'nodes': nodes,
                                                   'reason': oracle['reason'],
                                                   'production_interval': [c_lo, c_hi]})
                        counts['complex_component_oracle_guards'] += 1
                    for label, witness in [('component', local_pins), ('global', pins)]:
                        w = witness_success(rows, nodes, witness)
                        if w['checked']:
                            assert w['success'], (ident, nodes, label, witness)
                            counts[label + '_complex_witnesses_verified'] += 1
                            counts[label + '_complex_witness_states'] += w['states']
                            counts[label + '_complex_witness_' + w['method']] += 1
                        else:
                            unvalidated.append({'id': list(ident), 'nodes': nodes,
                                                'witness': label, **w})
            assert (scc_lower, scc_upper) == (lower, upper)
            if all_component_optima_known:
                assert lower <= exact_sum <= upper, (ident, exact_sum, lower, upper)
                if lower == upper:
                    assert lower == exact_sum
                counts['draw_optima_known_from_independent_components'] += 1
            if not certificate['G'] and len(rows) <= full_residual_cap:
                oracle = exact_cost(rows, residual, max_size=full_residual_cap)
                assert oracle['complete']
                assert lower <= oracle['cost'] <= upper, (ident, oracle)
                if all_component_optima_known:
                    assert oracle['cost'] == exact_sum
                counts['full_residual_global_optima_enumerated'] += 1
                counts['full_residual_global_pin_candidates'] += oracle['candidates']
                counts['full_residual_global_compatible_states'] += oracle['states']
            counts['final_witnesses_verified_by_components'] += 1
            assert not unvalidated, unvalidated[-1:]  # A verification gap is a failure.
    return {'family': family, 'n': n, 'draws': expected_rows,
            'raw_sha256': sha(raw_path), 'detail_sha256': sha(detail_path),
            'counts': dict(sorted(counts.items())),
            'maximum_literal_support_iterations': max_iterations,
            'maximum_SCC_size': max_component,
            'complex_component_size_counts': dict(sorted(complex_size_counts.items())),
            'component_reason_counts': dict(sorted(bound_reason_counts.items())),
            'independent_optimum_guards': incomplete_oracles,
            'unverified_witnesses': unvalidated}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--results', type=Path, default=root/'validation/results/main')
    parser.add_argument('--output', type=Path, default=root/'qa/independent_joint_replay.json')
    parser.add_argument('--expected-rows', type=int, default=2000)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--full-residual-cap', type=int, default=12)
    args = parser.parse_args()
    self_check_record = self_check()
    protocol_path = root/'validation/Protocol_JointRules.json'
    protocol = json.loads(protocol_path.read_text())
    assert args.expected_rows <= protocol['replicates_per_cell']
    manifest = json.loads((args.results/'manifest.json').read_text())
    assert manifest['status'] == 'completed'
    assert manifest['protocol_sha256'] == sha(protocol_path)
    gate_path = root/'qa/independent_exact_checks.json'
    assert manifest['independent_gate_sha256'] == sha(gate_path)
    representation=json.loads((root/'REPRESENTATION.json').read_text())
    assert sha(gate_path) == representation['current_gate_sha256']
    for relative,digest in manifest['source_sha256'].items():
        assert sha(root/'validation'/relative) == digest
    integrity = json.loads((args.results/'io_integrity.json').read_text())
    assert integrity['status'] == 'PASS'
    for stream in integrity['streams']:
        assert sha(args.results/stream['path']) == stream['sha256']
        assert stream['records'] == args.expected_rows
    jobs = [(str(root), str(args.results.resolve()), m, n, args.expected_rows,
             args.full_residual_cap) for m in protocol['families'] for n in protocol['sizes']]
    cells = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(worker, job) for job in jobs]
        for future in as_completed(futures):
            cell = future.result()
            cells.append(cell)
            print(f"verified {cell['family']} n={cell['n']} draws={cell['draws']}", flush=True)
    cells.sort(key=lambda c:(c['family'], c['n']))
    totals = Counter()
    for c in cells:
        totals.update(c['counts'])
    record = {'schema_version': 1, 'status': 'PASS',
              'protocol_sha256': sha(protocol_path),
              'validation_gate_sha256': sha(gate_path),
              'production_manifest_sha256': sha(args.results/'manifest.json'),
              'production_io_integrity_sha256': sha(args.results/'io_integrity.json'),
              'oracle_self_check': self_check_record,
              'validation_source_sha256': {p.name:sha(p) for p in
                  [Path(__file__), Path(__file__).with_name('independent_network_oracle.py')]},
              'coverage': {'cells': len(cells), 'draws': sum(c['draws'] for c in cells),
                           'full_source_replay': True,
                           'full_residual_global_optimum_dimension_cap': args.full_residual_cap,
                           'complex_component_oracle_dimension_cap': 18,
                           'complex_component_oracle_candidate_cap': 100000,
                           'complex_component_oracle_compatible_state_cap': 100000000},
              'independence': 'No production functions imported. Fixed PCG64 call sequence independently implemented; synchronous support iteration; Tarjan SCCs; independently constructed physical successors and indegree-peeling pin success. Copy and acyclic components use exact symbolic certificates.',
              'counts': dict(sorted(totals.items())), 'cells': cells}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': 'PASS', 'output': str(args.output), 'counts': record['counts']}), flush=True)


if __name__ == '__main__':
    main()
