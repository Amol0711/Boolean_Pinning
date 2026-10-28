#!/usr/bin/env python3
"""Resolve the one finite oracle guard on an Bernoulli production-exact component.

The full independent audit found one size-17 component whose own 20M-state
search guard stopped during cardinality four. This targeted replay raises only
that component's compatible-state guard to 100M, retaining all other caps.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json

from independent_network_oracle import (reconstruct, literal_support, export_rows,
    graph_certificate, exact_cost, witness_success)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=root/'simulation/results/main')
    parser.add_argument('--output', type=Path, default=root/'results/independent_bernoulli_targeted.json')
    args = parser.parse_args()
    full_path = root/'results/independent_bernoulli_validation.json'
    full = json.loads(full_path.read_text())
    guards = [g for c in full['cells'] for g in c['independent_optimum_guards']
              if len(g['nodes']) <= 18 and g['production_interval'][0] == g['production_interval'][1]]
    assert len(guards) == 1
    guard = guards[0]
    mixture, n, replicate = guard['id']
    details_path = args.results/'details'/f'{mixture}_n{n}.jsonl.gz'
    with gzip.open(details_path, 'rt') as stream:
        detail = next(json.loads(line) for j,line in enumerate(stream) if j == replicate)
    protocol_path = root/'simulation/baseline/Protocol_Baseline.json'
    protocol = json.loads(protocol_path.read_text())
    p, b, digest = reconstruct(protocol, mixture, n, replicate)
    residual, iterations = literal_support(p, b)
    rows = export_rows(p, b, residual)
    assert digest == detail['network_sha256']
    assert rows == detail['residual_rows']
    graph = graph_certificate(rows)
    components = []
    total = 0
    for component in graph['components']:
        nodes = component['nodes']
        if not component['cyclic']:
            cost = 0
        elif component['marked_slot_cycle']:
            cost = 1
        else:
            limit = 100_000_000 if nodes == guard['nodes'] else 20_000_000
            result = exact_cost(rows, nodes, max_compatible_states=limit)
            assert result['complete'], result
            cost = result['cost']
            witness = witness_success(rows, nodes, detail['pins'])
            assert witness['checked'] and witness['success']
            components.append({'nodes': nodes, 'compatible_state_guard': limit,
                               'oracle': result, 'production_witness_check': witness})
        total += cost
    assert total == detail['lower'] == detail['upper']
    selected = next(c for c in components if c['nodes'] == guard['nodes'])
    assert selected['oracle']['cost'] == guard['production_interval'][0]
    record = {'schema_version': 1, 'status': 'PASS', 'id': guard['id'],
              'full_audit_sha256': sha(full_path),
              'detail_file_sha256': sha(details_path),
              'protocol_sha256': sha(protocol_path),
              'network_sha256': digest,
              'source_sha256': {p.name:sha(p) for p in [Path(__file__),
                                  Path(__file__).with_name('independent_network_oracle.py')]},
              'original_finite_guard': guard, 'components': components,
              'independently_verified_draw_optimum': total,
              'combined_full_and_targeted_coverage': {
                  'production_exact_draws_independently_verified':
                      full['counts']['draw_optima_known_from_independent_components']+1,
                  'complex_component_optima_exhaustively_verified':
                      full['counts']['complex_component_optima_exhaustively_verified']+1,
                  'remaining_component_optimum_guards':
                      full['counts']['complex_component_oracle_guards']-1,
                  'remaining_guard_reason': 'All remaining guards have more than 18 coordinates and are retained as intervals, never asserted exact.'}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
