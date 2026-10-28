"""Repeat the eleven frozen contexts and all full-state target computations."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor'))
from scc_control import BooleanNetwork, support_residual, solve_network
from biology_core import (read_model, enumerate_fixed_points, discover_zero_orbit,
                          restricted_rows, successors, success, minimum_pinsets)

CAPS = {'max_component_size': 18, 'max_pins': 18, 'max_work': 100_000_000}
MARKERS = ('v_TBET', 'v_GATA3', 'v_Dec2', 'v_FOXP3', 'v_RORGT', 'v_proliferation',
           'v_STAT5_HIGH', 'v_IL4', 'v_IRF4', 'v_SMAD3', 'v_TGFB', 'v_TGFBR')


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze_target(model, target, assignment, target_id, discovery_target):
    baseline = tuple(model.names.index(n) for n in model.inputs)
    net = BooleanNetwork(model.parents, model.tables, target)
    support = support_residual(net, baseline, mode='physical')
    R = support['residual']; rpos = {v: j for j, v in enumerate(R)}
    if len(R) > 20:
        raise RuntimeError('Declared whole-residual witness dimension cap of 20 exceeded')
    rows = restricted_rows(model.parents, model.tables, R, target, target)
    reduced = BooleanNetwork([[rpos[u] for u in row['parents']] for row in rows],
                             [row['table'] for row in rows])
    solved = solve_network(reduced, support_mode='physical', **CAPS)
    assert solved['residual'] == list(range(len(R)))
    if solved['exact_cost'] is None: raise RuntimeError('Declared local search budget did not resolve target')
    components = []
    global_pins = []
    functional_parents = [()] * len(model.names)
    functional_tables = [(0,)] * len(model.names)
    for row in rows:
        functional_parents[row['node']] = tuple(row['parents'])
        functional_tables[row['node']] = tuple(row['table'])
    for cindex, comp in enumerate(solved['components']):
        C = [R[v] for v in comp['nodes']]
        localrows = restricted_rows(functional_parents, functional_tables, C, [0] * len(model.names))
        nxt = successors(localrows, C)
        minimum = minimum_pinsets(nxt, C, known_minimum=comp['exact_cost'])
        assert minimum['status'] == 'exact'
        canonical = minimum['minimum_pinsets'][0]
        global_pins.extend(canonical)
        record = {'component_id': cindex, 'nodes': C, 'names': [model.names[v] for v in C],
                  'cyclic': comp['cyclic'], 'restricted_rows': localrows,
                  'exact_cost': comp['exact_cost'], 'pins': canonical,
                  'pin_names': [model.names[v] for v in canonical],
                  'minimum_pinsets': minimum['minimum_pinsets'],
                  'minimum_pinsets_by_name': [[model.names[v] for v in K] for K in minimum['minimum_pinsets']],
                  'all_minimum_pinsets_enumerated': True,
                  'minimum_cardinality_enumeration_candidates': minimum['tested_candidates'],
                  'solver_record_reduced_indices': comp}
        components.append(record)
    global_pins.sort()
    nxt = successors(rows, R)
    ok, failed = success(nxt, sum(1 << rpos[v] for v in global_pins))
    assert ok and len(global_pins) == solved['exact_cost']
    marker = {name.removeprefix('v_'): target[model.names.index(name)] for name in MARKERS if name in model.names}
    eligibility = {'status': 'not_applicable', 'cost': None}
    if 'v_proliferation' in model.names:
        p = model.names.index('v_proliferation')
        if p not in rpos:
            eligibility = {'status': 'exact', 'cost': solved['exact_cost'],
                           'global_pins': global_pins,
                           'reason': 'forbidden_coordinate_outside_target_residual'}
        else:
            constrained_parts = []
            constrained_pins = []
            feasible = True
            for comp in components:
                C = comp['nodes']
                if p not in C:
                    constrained_parts.append({'component_id': comp['component_id'], 'status': 'exact',
                                              'cost': comp['exact_cost'], 'pins': comp['pins']})
                    constrained_pins.extend(comp['pins'])
                else:
                    con = minimum_pinsets(successors(comp['restricted_rows'], C), C,
                                          allowed=[v for v in C if v != p])
                    con['component_id'] = comp['component_id']
                    constrained_parts.append(con)
                    if con['status'] != 'exact': feasible = False
                    else: constrained_pins.extend(con['minimum_pinsets'][0])
            eligibility = {'status': 'exact' if feasible else 'infeasible',
                           'cost': len(constrained_pins) if feasible else None,
                           'global_pins': sorted(constrained_pins) if feasible else None,
                           'reason': 'exact_component_search_with_proliferation_ineligible',
                           'components': constrained_parts}
            if feasible:
                assert success(nxt, sum(1 << rpos[v] for v in constrained_pins))[0]
        eligibility.update(forbidden_node=p, forbidden_name='v_proliferation',
                           target_proliferation=target[p], target_STAT5_HIGH=marker['STAT5_HIGH'])
    return {'target_id': target_id, 'target': list(target), 'target_by_name': dict(zip(model.names, target)),
            'target_bitstring': ''.join(map(str, target)),
            'selected_zero_orbit': list(target) == discovery_target,
            'marker_signature': marker, 'residual': R, 'residual_names': [model.names[v] for v in R],
            'residual_size': len(R), 'support_record': support, 'functional_rows': rows,
            'components': components, 'exact_cost': solved['exact_cost'], 'global_pins': global_pins,
            'global_pin_names': [model.names[v] for v in global_pins],
            'global_witness_compatible_residual_states': 1 << (len(R) - len(global_pins)),
            'global_witness_exhaustively_verified': True,
            'functional_G_prime': solved['G_prime'],
            'proliferation_forbidden': eligibility}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results')
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    ppath = ROOT / 'larger_networks/protocol/contexts.json'
    protocol = json.loads(ppath.read_text())
    data = {'schema_version': 'BiologicalTargetCensus_v1', 'protocol_sha256': digest(ppath),
            'bbm_immutable_commit': protocol['bbm_immutable_commit'],
            'semantics': 'synchronous; full-state fixed target; permanent target-value pins; every compatible initial state; unit costs; baseline environmental clamps excluded from additional cost',
            'graph_semantics': 'exact functional dependence on target-residual cube; all inessential parents removed by complete Boolean-difference tests',
            'budgets': CAPS, 'fixed_point_remaining_dimension_cap': 24,
            'all_minimizer_component_size_cap': 12, 'global_witness_residual_dimension_cap': 20,
            'source_sha256': {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), ROOT/'src/biology_core.py', ROOT/'vendor/scc_control.py']},
            'models': {}, 'contexts': []}
    models = {}
    for key, info in protocol['models'].items():
        path = ROOT / 'larger_networks/data' / key / 'model.bnet'
        assert digest(path) == info['bnet_sha256']
        model = read_model(path); models[key] = model
        assert len(model.names) == info['total_nodes'] and list(model.inputs) == info['input_names']
        data['models'][key] = {'node_order': list(model.names), 'regulated_names': list(model.dynamic),
                               'input_names': list(model.inputs), 'bnet_sha256': digest(path),
                               'source_rule_count': len(model.dynamic),
                               'maximum_compiled_physical_arity': max(map(len, model.parents)),
                               'compiled_truth_table_entries': sum(map(len, model.tables))}
        if 'v_proliferation' in model.names:
            p = model.names.index('v_proliferation'); s = model.names.index('v_STAT5_HIGH')
            assert set(model.parents[p]) == {p, s}
            assert model.tables[p] == (0, 1, 1, 1)
            data['models'][key]['proliferation_rule_verified'] = 'proliferation_next = STAT5_HIGH OR proliferation'
    summary, targetrows = [], []
    for ctx in protocol['contexts']:
        model = models[ctx['model_id']]
        targets, forced, free, rounds = enumerate_fixed_points(model, ctx['input_assignment'])
        discovery = discover_zero_orbit(model, ctx['input_assignment'])
        archive = json.loads((ROOT / 'larger_networks/results' / ctx['id'] / 'certificate.json').read_text())
        if discovery['status'] == 'fixed':
            assert dict(zip(model.names, discovery['target'])) == archive['target_by_name']
        entry = {'id': ctx['id'], 'model_id': ctx['model_id'],
                 'input_assignment': ctx['input_assignment'], 'on_inputs': ctx['on_inputs'],
                 'baseline_count': len(model.inputs),
                 'forced_values': {model.names[i]: value for i, value in sorted(forced.items())},
                 'forced_value_rounds': rounds, 'free_nodes': free,
                 'free_names': [model.names[i] for i in free],
                 'enumerated_assignments': 1 << len(free), 'fixed_point_count': len(targets),
                 'zero_orbit_discovery': discovery, 'fixed_points': []}
        for j, target in enumerate(targets):
            tid = f"{ctx['id']}_FP{j+1:02d}"
            rec = analyze_target(model, target, ctx['input_assignment'], tid, discovery['target'])
            if rec['selected_zero_orbit']:
                assert rec['residual_size'] == archive['residual_size']
                assert rec['exact_cost'] == archive['additional_cost_bounds'][0] == archive['additional_cost_bounds'][1]
            entry['fixed_points'].append(rec)
            targetrows.append({'context_id': ctx['id'], 'target_id': tid,
                               'selected_zero_orbit': rec['selected_zero_orbit'],
                               'residual_size': rec['residual_size'], 'additional_cost': rec['exact_cost'],
                               'functional_component_sizes': '+'.join(str(len(c['nodes'])) for c in rec['components'] if c['cyclic']),
                               'cyclic_component_costs': '+'.join(str(c['exact_cost']) for c in rec['components'] if c['cyclic']),
                               'global_pin_names': ';'.join(rec['global_pin_names']),
                               'proliferation_forbidden_status': rec['proliferation_forbidden']['status'],
                               'proliferation_forbidden_cost': rec['proliferation_forbidden']['cost'],
                               **rec['marker_signature']})
        selected = [r for r in entry['fixed_points'] if r['selected_zero_orbit']]
        assert len(selected) == (1 if discovery['status'] == 'fixed' else 0)
        entry['selected_target_id'] = selected[0]['target_id'] if selected else None
        entry['cost_histogram'] = dict(sorted(Counter(r['exact_cost'] for r in entry['fixed_points']).items()))
        entry['proliferation_forbidden_status_counts'] = dict(Counter(r['proliferation_forbidden']['status'] for r in entry['fixed_points']))
        data['contexts'].append(entry)
        summary.append({'context_id': ctx['id'], 'model_id': ctx['model_id'], 'fixed_point_count': len(targets),
                        'enumerated_assignments': entry['enumerated_assignments'],
                        'free_dimension': len(free), 'cost_histogram': json.dumps(entry['cost_histogram'], sort_keys=True),
                        'selected_cost': selected[0]['exact_cost'] if selected else None,
                        'selected_residual_size': selected[0]['residual_size'] if selected else None,
                        'selected_proliferation_forbidden_status': selected[0]['proliferation_forbidden']['status'] if selected else None})
        print(ctx['id'], 'fixed_points', len(targets), 'costs', entry['cost_histogram'], flush=True)
    data['total_fixed_targets'] = sum(c['fixed_point_count'] for c in data['contexts'])
    data['status'] = 'completed'
    (args.output / 'all_targets.json').write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    for filename, rows in [('context_summary.csv', summary), ('all_target_summary.csv', targetrows)]:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (args.output / filename).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


if __name__ == '__main__': main()
