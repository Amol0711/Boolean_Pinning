"""Verify frozen larger-network certificates independently of their producer.

Example: python verification/verify_results.py --protocol protocol.json \
    --results results --output verification/independent_results.json
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

from independent_core import (read_bnet, step, support_trace, parse_expression,
                              variables, evaluate, residual_successors, controlled_graph)

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def named_trace(trace, order):
    return [{order[j] for j in support} for support in trace]


def verify_case(path, context, protocol, network):
    result = json.loads(path.read_text())
    rows, regulated, inputs = network
    order = result['node_order']
    assert len(order) == len(set(order)) == len(rows)
    assert set(order) == {row[0] for row in rows}
    assert set(regulated) == set(result['dynamic_node_names'])
    assert inputs == result['input_node_names']
    assert result['context'] == context
    assert result['budget'] == protocol['budgets']
    assert result['baseline_input_assignment'] == context['input_assignment']
    assert set(context['input_assignment']) == set(inputs)
    assert all(type(v) is int and v in (0, 1)
               for v in context['input_assignment'].values())
    assert set(context['on_inputs']) == {name for name,value in context['input_assignment'].items() if value == 1}
    source = ROOT / 'data' / result['model_id']
    assert sha(source/'model.bnet') == result['source']['model_bnet_sha256']
    assert sha(source/'model.bnet') == protocol['models'][result['model_id']]['bnet_sha256']
    assert sha(source/'metadata.json') == result['source']['metadata_sha256']
    metadata = json.loads((source/'metadata.json').read_text())
    assert set(metadata['input-names']) == set(inputs)
    assert metadata['variables'] == len(regulated)
    assert metadata['inputs'] == len(inputs)
    corners = sum(len(row[2]) for row in rows)
    assert corners == result['source']['local_truth_table_corners']
    rowmap = {row[0]: row for row in rows}
    for name, formula in context.get('model_rule_mapping', {}).items():
        tree = parse_expression(formula)
        _, scope, table = rowmap[name]
        assert variables(tree) == set(scope)
        assert all(evaluate(tree,dict(zip(scope,bits))) == value for bits,value in table.items())
    identity = []
    for name, scope, table in rows:
        if name in scope and all(value == bits[scope.index(name)]
                                 for bits, value in table.items()):
            identity.append(name)
    assert set(inputs) <= set(identity)
    assignment = context['input_assignment']
    state = {name: assignment.get(name, 0) for name in order}
    discovery = result['discovery']
    assert [state[name] for name in order] == discovery['initial_state']
    seen = {}; trajectory = []
    for update in range(protocol['budgets']['target_orbit_max_updates'] + 1):
        code = tuple(state[name] for name in order)
        if code in seen:
            first = seen[code]; cycle = trajectory[first:]
            break
        seen[code] = update; trajectory.append(code)
        state = step(rows, state)
    else:
        cycle = []
    record = {'context_id': context['id'], 'model_id': result['model_id'],
              'source_bnet_sha256': sha(source/'model.bnet'),
              'certificate_sha256': sha(path), 'nodes': len(rows),
              'baseline_inputs': len(inputs), 'local_truth_corners': corners,
              'status': discovery['status']}
    if discovery['target'] is None:
        if discovery['status'] == 'target_not_found_by_declared_orbit':
            assert len(cycle) > 1 and len(cycle) == discovery['cycle_length']
            assert [list(s) for s in cycle] == discovery['cycle']
            record['verified_seed_orbit_cycle_length'] = len(cycle)
        elif discovery['status'] == 'target_discovery_update_budget':
            assert not cycle
        else:
            assert discovery['status'] == 'target_discovery_wall_budget'
            record['wall_cutoff_scope'] = 'No assertion that a fixed point does not exist.'
        return record
    assert len(cycle) == 1
    assert list(cycle[0]) == discovery['target']
    assert first == discovery['updates']
    if 'candidate_orbit' in discovery:
        assert [list(s) for s in trajectory] == discovery['candidate_orbit']
    target = dict(zip(order, discovery['target']))
    assert target == result['target_by_name']
    assert step(rows, target) == target
    assert all(target[name] == assignment[name] for name in inputs)
    assert set(identity) == set(result['identity_rows'])
    assert {order[j] for j in result['baseline_pin_indices']} == set(inputs)
    trace = support_trace(rows, target, inputs)
    assert trace == named_trace(result['baseline_support_trace'], order)
    assert [len(s) for s in trace] == result['baseline_support_sizes']
    nodes = result['residual_node_names']
    assert nodes == [order[j] for j in result['residual_node_indices']]
    assert set(nodes) == trace[-1]
    assert len(nodes) == result['residual_size']
    search = result['search']
    lower = search['lower_bound']; upper = search['upper_bound']
    assert 0 <= lower <= upper <= len(nodes)
    assert [lower, upper] == result['additional_cost_bounds']
    assert result['baseline_cost'] == len(inputs)
    assert result['total_cost_bounds_counting_inputs'] == [len(inputs)+lower, len(inputs)+upper]
    checked = {}
    graph = None
    cycles_checked = 0
    if result['residual_successor_table'] is not None:
        graph = residual_successors(rows, target, nodes)
        assert graph == result['residual_successor_table']
        for candidate in search['candidate_results']:
            pins = tuple(candidate['pin_indices_local'])
            assert pins == tuple(sorted(set(pins)))
            assert pins not in checked
            assert candidate['cardinality'] == len(pins)
            assert candidate['pin_names'] == [nodes[j] for j in pins]
            assert candidate['pin_values'] == {nodes[j]: target[nodes[j]] for j in pins}
            independent = controlled_graph(graph, pins)
            checked[pins] = independent
            assert independent['compatible_states'] == candidate['compatible_states']
            assert independent['convergent_states'] == candidate['convergent_states']
            assert independent['success'] == candidate['globally_convergent']
            mask = sum(1 << j for j in pins)
            if independent['success']:
                assert independent['maximum_hitting_time'] == candidate['max_residual_hitting_time']
                route = candidate['longest_residual_trajectory']
                assert len(route)-1 == independent['maximum_hitting_time']
                assert route[-1] == 0
                assert all(not s & mask for s in route)
                assert all((graph[a] & ~mask) == b for a,b in zip(route, route[1:]))
            else:
                cycle = candidate['failure_cycle']
                assert cycle and 0 not in cycle and len(cycle) == len(set(cycle))
                assert all(not s & mask for s in cycle)
                for a, b in zip(cycle, cycle[1:]+cycle[:1]):
                    assert (graph[a] & ~mask) == b
                    # Lift the failed cycle to the original source network.
                    full = target.copy()
                    for j, name in enumerate(nodes):
                        full[name] ^= (a >> j) & 1
                    following = step(rows, full)
                    for name in inputs + candidate['pin_names']:
                        following[name] = target[name]
                    expected = target.copy()
                    for j, name in enumerate(nodes):
                        expected[name] ^= (b >> j) & 1
                    assert following == expected
                    cycles_checked += 1
        # Every declared complete level must contain each subset exactly once.
        for cost in search['completed_cardinality_levels']:
            expected = set(itertools.combinations(range(len(nodes)), cost))
            observed = {p for p in checked if len(p) == cost}
            assert observed == expected
        # This is the substantive lower-bound check, not just a witness check.
        for cost in range(lower):
            for pins in itertools.combinations(range(len(nodes)), cost):
                assert pins in checked and not checked[pins]['success']
        if search['all_minimum_sets_complete']:
            assert lower == upper
            optima = sorted(tuple(nodes[j] for j in p) for p, c in checked.items()
                            if len(p) == lower and c['success'])
            assert optima == sorted(tuple(x) for x in search['all_minimum_additional_sets'])
            assert lower in search['completed_cardinality_levels']
        declared_work = search['transition_evaluations_reserved']
        audited_work = len(graph) + sum(c['compatible_states'] for c in checked.values())
        assert audited_work <= declared_work <= protocol['budgets']['max_transition_evaluations']
    elif nodes:
        if lower:
            assert search['status'] == 'exact_copy_certificate'
            copies = {}
            for name, scope, table in rows:
                if name not in nodes:
                    continue
                compatible_corners = [(bits,value) for bits,value in table.items()
                                      if all(p in nodes or bit == target[p]
                                             for p,bit in zip(scope,bits))]
                candidates = [p for p in scope if p in nodes and all(
                    (value ^ target[name]) == (bits[scope.index(p)] ^ target[p])
                    for bits,value in compatible_corners)]
                assert len(candidates) == 1
                copies[name] = candidates[0]
            seen = set(); copy_cycles = []
            for name in nodes:
                route = []; position = {}; current = name
                while current not in seen and current not in position:
                    position[current] = len(route); route.append(current)
                    current = copies[current]
                if current in position:
                    copy_cycles.append(set(route[position[current]:]))
                seen.update(route)
            assert lower == upper == len(copy_cycles)
            assert all(cycle.intersection(search['feasible_additional_set'])
                       for cycle in copy_cycles)
    else:
        assert lower == upper == 0
        assert search['all_minimum_additional_sets'] == [[]]
    feasible = search['feasible_additional_set']
    assert len(feasible) == upper and len(set(feasible)) == len(feasible)
    assert set(feasible) <= set(nodes)
    assert search['feasible_additional_values'] == {name: target[name] for name in feasible}
    controlled_trace = support_trace(rows, target, inputs + feasible)
    if 'first_feasible_full_support_trace' in result:
        assert controlled_trace == named_trace(result['first_feasible_full_support_trace'], order)
    justified_bounds = []
    if not controlled_trace[-1]:
        justified_bounds.append(len(controlled_trace)-1)
    if graph is not None:
        pin_indices = tuple(nodes.index(name) for name in feasible)
        feasible_check = checked.get(pin_indices) or controlled_graph(graph, pin_indices)
        assert feasible_check['success']
        justified_bounds.append(len(trace)-1 + feasible_check['maximum_hitting_time'])
    assert justified_bounds, 'Feasible upper bound lacks an independently checked witness'
    if 'global_convergence_upper_bound' in result:
        assert result['global_convergence_upper_bound'] >= min(justified_bounds)
    resources = json.loads(path.with_name('resources.json').read_text())
    assert resources['analysis_wall_seconds'] > 0 and resources['peak_worker_rss_bytes'] > 0
    record.update(status=search['status'], target_verified=True,
                  input_completion_identity_verified=True,
                  support_sizes=[len(s) for s in trace], residual_size=len(nodes),
                  residual_states_independently_evaluated=0 if graph is None else len(graph),
                  candidates_independently_classified=len(checked),
                  compatible_candidate_states=sum(c['compatible_states'] for c in checked.values()),
                  full_network_failure_cycle_edges_checked=cycles_checked,
                  additional_cost_bounds=[lower,upper],
                  optimum_list_complete=search['all_minimum_sets_complete'],
                  feasible_pins=feasible,
                  independently_justified_convergence_bound=min(justified_bounds),
                  declared_convergence_bound=result.get('global_convergence_upper_bound'))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    protocol = json.loads(args.protocol.read_text())
    assert protocol['frozen'] is True
    contexts = protocol['contexts']
    assert len({c['id'] for c in contexts}) == len(contexts)
    manifest = json.loads((args.results/'run_manifest.json').read_text())
    assert manifest['protocol_sha256'] == sha(args.protocol)
    assert manifest['contexts_executed_in_frozen_order'] == [c['id'] for c in contexts]
    expected = {c['id'] for c in contexts}
    assert {p.parent.name for p in args.results.glob('*/certificate.json')} == expected
    networks = {model: read_bnet(ROOT/'data'/model/'model.bnet')
                for model in {c['model_id'] for c in contexts}}
    cases = []
    for context in contexts:
        key = context['id']+'/certificate.json'
        assert manifest['scientific_result_sha256'][key] == sha(args.results/key)
        case = verify_case(args.results/context['id']/'certificate.json',
                           context, protocol, networks[context['model_id']])
        cases.append(case)
        print(context['id'], 'verified', case.get('additional_cost_bounds', case['status']), flush=True)
    output = {'status': 'ALL_FROZEN_CASES_INDEPENDENTLY_VERIFIED',
              'method': 'Native BNET recursive-descent parser; complete local truth-table evaluation; exact coordinate-cube propagation; reverse reachability from zero in controlled residual graphs; original-network lifting of every reported failure-cycle edge. No production parser, compiled table, certificate engine, or search imports.',
              'protocol_sha256': sha(args.protocol), 'cases': cases,
              'totals': {
                  'contexts': len(cases),
                  'candidate_sets': sum(c.get('candidates_independently_classified',0) for c in cases),
                  'candidate_states': sum(c.get('compatible_candidate_states',0) for c in cases),
                  'residual_successors': sum(c.get('residual_states_independently_evaluated',0) for c in cases),
                  'lifted_failure_cycle_edges': sum(c.get('full_network_failure_cycle_edges_checked',0) for c in cases)},
              'elapsed_seconds': time.perf_counter()-started,
              'limits': ['Exactness is for the declared synchronous, identity-completed model and unrestricted model-coordinate pin eligibility.',
                         'This audit is not biological phenotype validation or a complexity guarantee for all networks.',
                         'Runtime and RSS measurements are descriptive single-worker measurements, not a cross-method performance comparison.']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    main()
