"""Independent fair-coordinate asynchronous and absorbing-memory checks.

Run from workflows/biological:
  python verification/async_memory_check.py --output results/async_memory_check.json

Standard library only. No production solver, archived compiled rule table,
external validation implementation, NumPy, or SciPy is imported. Source BNET expressions
are parsed by shunting yard and evaluated as postfix Boolean programs.
"""
from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
import hashlib
import itertools
import json
from pathlib import Path
import re


@dataclass
class Rule:
    name: str
    parents: tuple
    table: tuple


def postfix(expression):
    tokens = re.findall(r"[A-Za-z_][A-Za-z_0-9]*|[01!&|()]", expression)
    assert ''.join(tokens) == re.sub(r'\s+', '', expression)
    precedence = {'!': 3, '&': 2, '|': 1}
    output, operators = [], []
    for token in tokens:
        if token == '(':
            operators.append(token)
        elif token == ')':
            while operators[-1] != '(':
                output.append(operators.pop())
            operators.pop()
        elif token in precedence:
            while (operators and operators[-1] in precedence and token != '!'
                   and precedence[operators[-1]] >= precedence[token]):
                output.append(operators.pop())
            operators.append(token)
        else:
            output.append(token)
    while operators:
        token = operators.pop()
        assert token != '('
        output.append(token)
    return output


def evaluate_program(program, values):
    stack = []
    for token in program:
        if token == '!':
            stack.append(1 - stack.pop())
        elif token in ('&', '|'):
            right, left = stack.pop(), stack.pop()
            stack.append(left & right if token == '&' else left | right)
        else:
            stack.append(int(token) if token in ('0', '1') else values[token])
    assert len(stack) == 1
    return stack[0]


def read_model(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('targets,'):
            continue
        name, expression = line.split(',', 1)
        program = postfix(expression)
        parents = tuple(sorted({token for token in program
                                if re.fullmatch('[A-Za-z_][A-Za-z_0-9]*', token)}))
        table = tuple(evaluate_program(program, {p: (s >> j) & 1
                                                for j, p in enumerate(parents)})
                      for s in range(1 << len(parents)))
        rows.append(Rule(name.strip(), parents, table))
    regulated = {row.name for row in rows}
    assert len(regulated) == len(rows)
    inputs = sorted({p for row in rows for p in row.parents} - regulated)
    return rows, inputs


def row_value(row, values):
    return row.table[sum(values[p] << j for j, p in enumerate(row.parents))]


def support_trace(rows, target, baseline):
    alive = set(target) - set(baseline)
    trace = [sorted(alive)]
    while True:
        following = set()
        for row in rows:
            if row.name in baseline:
                continue
            for pattern, value in enumerate(row.table):
                if value == target[row.name]:
                    continue
                if all(p in alive or ((pattern >> j) & 1) == target[p]
                       for j, p in enumerate(row.parents)):
                    following.add(row.name)
                    break
        assert following <= alive
        if following == alive:
            return trace
        trace.append(sorted(following))
        alive = following


def restricted_rules(rowmap, target, nodes):
    node_set = set(nodes)
    result = {}
    for name in nodes:
        row = rowmap[name]
        scope = tuple(p for p in row.parents if p in node_set)
        table = []
        for state in range(1 << len(scope)):
            values = {p: target[p] for p in row.parents}
            for j, parent in enumerate(scope):
                values[parent] ^= (state >> j) & 1
            table.append(row_value(row, values) ^ target[name])
        assert table[0] == 0
        result[name] = (scope, tuple(table))
    return result


def dependency_components(rowmap, target, residual):
    restricted = restricted_rules(rowmap, target, residual)
    edges = []
    for name, (scope, table) in restricted.items():
        for j, parent in enumerate(scope):
            if any(table[s] != table[s ^ (1 << j)]
                   for s in range(len(table)) if not s & (1 << j)):
                edges.append((parent, name))
    adjacent = {name: [] for name in residual}
    for source, destination in edges:
        adjacent[source].append(destination)
    reach = {}
    for name in residual:
        seen, stack = {name}, [name]
        while stack:
            for following in adjacent[stack.pop()]:
                if following not in seen:
                    seen.add(following)
                    stack.append(following)
        reach[name] = seen
    unused, components = set(residual), []
    while unused:
        first = min(unused)
        component = sorted(v for v in unused if v in reach[first] and first in reach[v])
        unused -= set(component)
        components.append(component)
    return sorted(edges), components


def graph_data(rowmap, target, nodes, pins=(), keep_graph=False):
    free = [name for name in nodes if name not in pins]
    position = {name: j for j, name in enumerate(free)}
    restricted = restricted_rules(rowmap, target, free)
    compiled = [(tuple(position[p] for p in restricted[name][0]), restricted[name][1])
                for name in free]
    count = 1 << len(free)
    graph = [[] for _ in range(count)]
    indegree = [0] * count
    edge_hash = hashlib.sha256()
    edge_count = 0
    for state in range(count):
        for coordinate, (parents, table) in enumerate(compiled):
            pattern = sum(((state >> source) & 1) << j for j, source in enumerate(parents))
            if table[pattern] != ((state >> coordinate) & 1):
                following = state ^ (1 << coordinate)
                graph[state].append(following)
                indegree[following] += 1
                edge_count += 1
                edge_hash.update(f'{state},{coordinate},{following}\n'.encode())
    sinks = [state for state, edges in enumerate(graph) if not edges]
    assert 0 in sinks
    queue = deque(state for state, degree in enumerate(indegree) if degree == 0)
    topological = []
    while queue:
        state = queue.popleft()
        topological.append(state)
        for following in graph[state]:
            indegree[following] -= 1
            if indegree[following] == 0:
                queue.append(following)
    acyclic = len(topological) == count
    path_bound = None
    if acyclic:
        length = [0] * count
        for state in reversed(topological):
            length[state] = max((1 + length[nxt] for nxt in graph[state]), default=0)
        path_bound = max(length)
    nonzero_fixed = next((state for state in sinks if state), None)
    result = {'free_node_order': free, 'compatible_states': count,
              'state_changing_edges': edge_count, 'state_changing_graph_sha256': edge_hash.hexdigest(),
              'acyclic': acyclic, 'topologically_removed_states': len(topological),
              'sink_count': len(sinks), 'target_is_only_sink': sinks == [0],
              'maximum_state_changes_after_capture': path_bound,
              'fair_convergence_certificate': acyclic and sinks == [0]}
    if nonzero_fixed is not None:
        result['nonzero_fixed_point_witness'] = {
            'error_state_index': nonzero_fixed,
            'physical_values': {name: target[name] ^ ((nonzero_fixed >> j) & 1)
                                for j, name in enumerate(free)},
            'pinned_values': {name: target[name] for name in pins},
        }
    if keep_graph:
        result['_graph'] = graph
    return result


def local_async_cost(rowmap, target, nodes):
    rejected = []
    for size in range(len(nodes) + 1):
        for pins in itertools.combinations(nodes, size):
            checked = graph_data(rowmap, target, nodes, pins)
            if checked['fair_convergence_certificate']:
                assert all('nonzero_fixed_point_witness' in rec for rec in rejected
                           if len(rec['pins']) < size)
                return {'nodes': nodes, 'exact_fair_coordinate_cost': size,
                        'feasible_pins': list(pins), 'upper_certificate': checked,
                        'lower_fixed_point_witnesses': [rec for rec in rejected
                                                       if len(rec['pins']) < size],
                        'tested_candidates_before_first_success': len(rejected) + 1}
            # A cycle alone is not assumed coordinate-fair. This checker only
            # rejects a candidate after exhibiting a distinct fixed point.
            if 'nonzero_fixed_point_witness' not in checked:
                raise AssertionError(f'Candidate needs a separate fair-cycle certificate: {nodes}, {pins}')
            rejected.append({'pins': list(pins),
                             'nonzero_fixed_point_witness': checked['nonzero_fixed_point_witness']})
    raise AssertionError('All-pinned isolated map must succeed')


def fairness_counterexample():
    # x1'=x1*x2, x2'=x1*(1-x2); physical target 00.
    rows = {'x1': Rule('x1', ('x1', 'x2'), (0, 0, 0, 1)),
            'x2': Rule('x2', ('x1', 'x2'), (0, 1, 0, 0))}
    target = {'x1': 0, 'x2': 0}
    info = graph_data(rows, target, ['x1', 'x2'], keep_graph=True)
    graph = info.pop('_graph')
    reachable = []
    for start in range(4):
        seen, stack = {start}, [start]
        while stack:
            for nxt in graph[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt); stack.append(nxt)
        reachable.append(seen)
    unused, components = set(range(4)), []
    while unused:
        first = min(unused)
        component = sorted(s for s in unused if s in reachable[first] and first in reachable[s])
        unused -= set(component); components.append(component)
    terminals = [c for c in components if all(nxt in c for s in c for nxt in graph[s])]
    assert terminals == [[0]] and not info['acyclic']
    state, path = 1, [1]
    schedule = [1, 0, 1]  # zero-based coordinate indices, repeated forever
    for coordinate in schedule:
        x1, x2 = state & 1, (state >> 1) & 1
        updated = x1 & x2 if coordinate == 0 else x1 & (1 - x2)
        state = (state & ~(1 << coordinate)) | (updated << coordinate)
        path.append(state)
    assert state == 1 and 0 not in path and set(schedule) == {0, 1}
    return {'rules': ["x1'=x1 AND x2", "x2'=x1 AND NOT x2"],
            'target': [0, 0], 'terminal_sccs': terminals,
            'coordinate_order': ['x1', 'x2'], 'periodic_update_indices': schedule,
            'periodic_state_indices': path, 'both_coordinates_updated_each_period': True,
            'converges_to_target': False}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(base):
    larger = base / 'larger_networks'
    model_path = larger / 'data/cd4_188/model.bnet'
    protocol_path = larger / 'protocol/contexts.json'
    rows, inputs = read_model(model_path)
    rowmap = {row.name: row for row in rows}
    protocol = json.loads(protocol_path.read_text())
    contexts = {context['id']: context for context in protocol['contexts']}
    results = []
    for context_id in ('CD4_OFF', 'CD4_APC'):
        certificate_path = larger / 'results' / context_id / 'certificate.json'
        archived = json.loads(certificate_path.read_text())
        assert sha(model_path) == archived['source']['model_bnet_sha256']
        target = archived['target_by_name']
        baseline = contexts[context_id]['input_assignment']
        assert set(baseline) == set(inputs)
        assert baseline == archived['baseline_input_assignment']
        assert all(target[name] == value for name, value in baseline.items())
        assert set(target) == set(rowmap) | set(inputs)
        assert all(row_value(row, target) == target[row.name] for row in rows)
        trace = support_trace(rows, target, baseline)
        residual = trace[-1]
        assert set(residual) == set(archived['residual_node_names'])
        edges, components = dependency_components(rowmap, target, residual)
        local_results = [local_async_cost(rowmap, target, component) for component in components]
        positive_components = [c for c in local_results if c['exact_fair_coordinate_cost']]
        optimum = sum(c['exact_fair_coordinate_cost'] for c in local_results)
        expected = 3 if context_id == 'CD4_OFF' else 4
        assert optimum == expected and all(c['exact_fair_coordinate_cost'] == 1
                                           for c in positive_components)
        variants = [['v_Dec2', 'v_TBET', 'v_proliferation'],
                    ['v_GATA3', 'v_TBET', 'v_proliferation']]
        if context_id == 'CD4_APC':
            variants = [pins + ['v_FOXP3'] for pins in variants]
        global_certificates = []
        for pins in variants:
            pins = sorted(pins)
            checked = graph_data(rowmap, target, residual, pins)
            assert checked['fair_convergence_certificate'] and len(pins) == optimum
            global_certificates.append({'pins': pins, 'pin_values': {p: target[p] for p in pins},
                                        **checked})
        results.append({'context_id': context_id, 'archive_certificate_sha256': sha(certificate_path),
                        'baseline_assignment': baseline, 'baseline_count': len(baseline),
                        'target_by_name': target, 'residual_node_order': residual,
                        'support_trace_sizes': list(map(len, trace)),
                        'functional_residual_edges': [list(edge) for edge in edges],
                        'isolated_components': local_results,
                        'exact_minimum_for_all_coordinate_fair_schedules': optimum,
                        'lower_bound_method': 'Each disjoint positive-cost SCC has a nonzero isolated fixed point when unpinned; strict ancestor initialization preserves it under any outside pins.',
                        'full_residual_upper_certificates': global_certificates})

    proliferation = rowmap['v_proliferation']
    assert set(proliferation.parents) == {'v_STAT5_HIGH', 'v_proliferation'}
    for state, value in enumerate(proliferation.table):
        assignment = {p: (state >> j) & 1 for j, p in enumerate(proliferation.parents)}
        assert value == (assignment['v_STAT5_HIGH'] | assignment['v_proliferation'])
    memory = []
    for context_id in sorted(k for k in contexts if k.startswith('CD4_')):
        target = json.loads((larger / 'results' / context_id / 'certificate.json').read_text())['target_by_name']
        memory.append({'context_id': context_id, 'target_proliferation': target['v_proliferation'],
                       'target_STAT5_HIGH': target['v_STAT5_HIGH'],
                       'forbidden_proliferation_infeasible_by_absorbing_one': target['v_proliferation'] == 0,
                       'scope': 'Selected archived target only; false here does not establish restricted-eligibility feasibility.'})
    off = next(r for r in results if r['context_id'] == 'CD4_OFF')
    apc = next(r for r in results if r['context_id'] == 'CD4_APC')
    positive = lambda r: {tuple(c['nodes']) for c in r['isolated_components']
                          if c['exact_fair_coordinate_cost']}
    new_components = sorted(positive(apc) - positive(off))
    assert new_components == [('v_FOXP3', 'v_SMAD3', 'v_TGFB', 'v_TGFBR')]
    return {'schema': 'FairCoordinateMemory_v1', 'status': 'PASS',
            'model_bnet_sha256': sha(model_path), 'protocol_sha256': sha(protocol_path),
            'independence': 'Own postfix BNET parser, truth tables, physical support, functional-dependency SCCs and asynchronous graph checks; no production or external validation imports.',
            'update_semantics': 'One unpinned coordinate selected per step, possibly unchanged; every unpinned coordinate is selected infinitely often; all pins permanent at target from time zero.',
            'certificate_logic': 'Schedule-dependent fair-round support capture; acyclic state-changing residual graph and unique sink zero; local competing fixed points certify all lower costs.',
            'contexts': results, 'memory_rule': "proliferation'=STAT5_HIGH OR proliferation",
            'memory_rule_verified_on_every_local_assignment': True,
            'selected_target_memory_checks': memory,
            'additional_APC_positive_cost_components_relative_to_OFF': [list(c) for c in new_components],
            'unique_terminal_scc_fairness_counterexample': fairness_counterexample(),
            'limitations': ['No uniform bound on asynchronous update count is claimed under unbounded fair delays.',
                            'The computation concerns the selected biological contexts.',
                            'A target with proliferation=1 requires separate eligibility analysis; the absorbing-one test alone does not decide it.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = args.output or args.base / 'results/async_memory_check.json'
    result = run(args.base)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    summary = [{'context': r['context_id'], 'exact_fair_cost': r['exact_minimum_for_all_coordinate_fair_schedules'],
                'states_per_full_residual_check': [x['compatible_states'] for x in r['full_residual_upper_certificates']]}
               for r in result['contexts']]
    print(json.dumps({'status': result['status'], 'contexts': summary, 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
