#!/usr/bin/env python3
"""Independent J1 replication preflight, before reading the new production outcomes.

No production module is imported. The exhaustive reference check compares
physical functional-graph peeling with direct literal trajectories. Additional
fixtures exercise repeated physical parents and truthful bounded-oracle exits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from independent_network_oracle import (
    exact_cost, export_rows, graph_certificate, literal_support, self_check,
    witness_success,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixture(name, parents, codes, expected_cost, expected_g, expected_gp):
    p = np.asarray(parents, dtype=np.int64)
    b = np.asarray([[(code >> j) & 1 for j in range(3)] for code in codes],
                   dtype=np.uint8)
    residual, iterations = literal_support(p, b)
    rows = export_rows(p, b, residual)
    graph = graph_certificate(rows)
    oracle = exact_cost(rows, residual)
    full_rows = export_rows(p, b, list(range(len(p))))
    full_oracle = exact_cost(full_rows, list(range(len(p))))
    assert oracle['complete'] and full_oracle['complete']
    assert oracle['cost'] == full_oracle['cost'] == expected_cost, name
    assert graph['G'] == expected_g and graph['G_prime'] == expected_gp, name
    witness = witness_success(full_rows, list(range(len(p))), oracle['pins'])
    assert witness['checked'] and witness['success'], name
    return dict(name=name, cost=expected_cost, residual=residual,
                G=expected_g, G_prime=expected_gp, iterations=iterations,
                physical_duplicate_rows=sum(a == c for a, c in parents),
                optimum_candidates=oracle['candidates'],
                optimum_compatible_states=oracle['states'])


def run():
    exhaustive = self_check()
    fixtures = [
        fixture('constant_zero', [[0, 0]], [0], 0, True, True),
        fixture('duplicate_parent_XOR_is_physically_zero', [[0, 0]], [3],
                0, False, False),
        fixture('duplicate_parent_AND_is_physical_copy', [[0, 0]], [4],
                1, False, False),
        fixture('duplicate_parent_COPY1', [[0, 0]], [5], 1, False, False),
        fixture('nominally_coupled_independent_self_copies', [[0, 1], [0, 1]],
                [5, 6], 2, False, False),
        fixture('cooperative_AND_pair', [[0, 1], [0, 1]], [4, 4],
                1, False, False),
        fixture('three_cycle_with_zero_exterior', [[2, 3], [0, 3], [1, 3], [3, 3]],
                [5, 5, 5, 0], 1, True, True),
    ]
    rows = [[0, 0, 0, 1, 0, 1]]
    guards = {
        'size': exact_cost(rows, [0], max_size=0),
        'candidate': exact_cost(rows, [0], max_candidates=1),
        'compatible_state': exact_cost(rows, [0], max_compatible_states=2),
        'witness_size': witness_success(rows, [0], [], max_size=0),
    }
    assert not guards['size']['complete']
    assert guards['size']['reason'] == 'independent_size_guard'
    for name in ['candidate', 'compatible_state']:
        assert not guards[name]['complete']
        assert guards[name]['reason'] == 'independent_enumeration_guard'
        assert guards[name]['lower'] == 1
        assert guards[name]['candidates'] == 1 and guards[name]['states'] == 2
    assert not guards['witness_size']['checked']
    assert guards['witness_size']['remaining_dimension'] == 1
    empty = exact_cost([], [])
    assert empty['complete'] and empty['cost'] == 0
    return dict(status='PASS', exhaustive_physical_oracle_check=exhaustive,
                edge_fixtures=fixtures, bounded_exit_fixtures=guards,
                empty_residual=empty,
                scope='Independent oracle self-check and physical repeated-parent and search-limit fixtures.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    result['source_sha256'] = {p.name: sha(p) for p in [
        Path(__file__), Path(__file__).with_name('independent_network_oracle.py')]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'output': str(args.output),
                      'exhaustive_networks': result['exhaustive_physical_oracle_check']['networks'],
                      'edge_fixtures': len(result['edge_fixtures']),
                      'bounded_exit_fixtures': len(result['bounded_exit_fixtures'])}), flush=True)


if __name__ == '__main__':
    main()
