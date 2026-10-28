"""Frozen-context exact Boolean benchmark; Python standard library only.

Coordinates are physical model variables. The only dynamical completion is an
explicit identity update for each free input omitted from the original BNET.
Each worker analyzes one context in a fresh process. No global 2^n enumeration
and no monotonicity assumption on success under additional pins are used.
"""

import argparse
import ast
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

from boolean_certificates import (BooleanModel, LocalRule, identity_rows,
                                  physical_residual_model, support_trace,
                                  pure_copy_certificate)

ROOT = Path(__file__).resolve().parents[1]
FROZEN_PROTOCOL_SHA256 = '2d8f64156c5643cd0559a149b5ccde1cc7bed26d06bf0f517963a7db188aa47d'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def parse_boolean(formula):
    tree = ast.parse(formula.replace('!', ' not ').replace('&', ' and ')
                    .replace('|', ' or ').strip(), mode='eval')
    whitelist = (ast.Expression, ast.Name, ast.Load, ast.BoolOp, ast.And, ast.Or,
                 ast.UnaryOp, ast.Not, ast.Constant)
    for node in ast.walk(tree):
        if not isinstance(node, whitelist):
            raise ValueError('Unsupported BNET Boolean syntax')
        if isinstance(node, ast.Constant) and node.value not in (0, 1):
            raise ValueError('Unsupported BNET constant')
    return tree


def evaluate(tree, values):
    if isinstance(tree, ast.Expression):
        return evaluate(tree.body, values)
    if isinstance(tree, ast.Name):
        return values[tree.id]
    if isinstance(tree, ast.Constant):
        return int(tree.value)
    if isinstance(tree, ast.UnaryOp):
        return 1 - evaluate(tree.operand, values)
    if isinstance(tree, ast.BoolOp) and isinstance(tree.op, ast.And):
        return int(all(evaluate(child, values) for child in tree.values))
    if isinstance(tree, ast.BoolOp) and isinstance(tree.op, ast.Or):
        return int(any(evaluate(child, values) for child in tree.values))
    raise ValueError('Unexpected BNET expression node')


def load_network(model_id):
    directory = ROOT / 'data' / model_id
    formulas = {}
    for line in (directory/'model.bnet').read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#') or line.replace(' ', '') == 'targets,factors':
            continue
        name, formula = line.split(',', 1)
        name, formula = name.strip(), formula.strip()
        if name in formulas:
            raise ValueError('Duplicate BNET row')
        formulas[name] = parse_boolean(formula)
    dynamic = tuple(sorted(formulas))
    referenced = {node.id for tree in formulas.values() for node in ast.walk(tree)
                  if isinstance(node, ast.Name)}
    inputs = tuple(sorted(referenced - set(dynamic)))
    names = dynamic + inputs
    location = {name: i for i, name in enumerate(names)}
    metadata = json.loads((directory/'metadata.json').read_text())
    if set(metadata['input-names']) != set(inputs):
        raise ValueError('Missing BNET input rules do not match input metadata')
    if metadata['variables'] != len(dynamic) or metadata['inputs'] != len(inputs):
        raise ValueError('BNET node counts do not match metadata')
    for name in inputs:
        formulas[name] = parse_boolean(name)
    rules = []
    for name in names:
        tree = formulas[name]
        parents = tuple(sorted({location[node.id] for node in ast.walk(tree)
                                if isinstance(node, ast.Name)}))
        if len(parents) > 20:
            raise ValueError('Local truth-table arity exceeds the declared compilation guard')
        table = tuple(evaluate(tree, {names[p]: (code >> j) & 1
                                     for j, p in enumerate(parents)})
                      for code in range(1 << len(parents)))
        rules.append(LocalRule(parents, table))
    source_info = {
        'model_bnet_sha256': digest(directory/'model.bnet'),
        'metadata_sha256': digest(directory/'metadata.json'),
        'source_publication': metadata['url-publication'],
        'source_model': metadata['url-model'],
        'dynamic_node_count_from_bnet': len(dynamic),
        'input_node_count_from_missing_rules': len(inputs),
        'metadata_variable_names_count': len(metadata['variable-names']),
        'variable_names_note': 'Node order is derived from BNET left-hand sides and referenced inputs. The metadata variable-names list may omit output coordinates and is not used as the complete node list.',
        'local_truth_table_corners': sum(len(rule.table) for rule in rules),
        'maximum_physical_arity': max(len(rule.parents) for rule in rules),
    }
    return tuple(rules), names, dynamic, inputs, formulas, source_info


def discover_target(rules, names, inputs, assignment, maximum_updates, deadline):
    if set(assignment) != set(inputs) or any(type(value) is not int or value not in (0, 1)
                                           for value in assignment.values()):
        raise ValueError('Context must explicitly assign every input exactly once')
    pins = {names.index(name): value for name, value in assignment.items()}
    initial = tuple(pins.get(i, 0) for i in range(len(names)))
    trajectory, seen = [initial], {initial: 0}
    state = initial
    for update in range(maximum_updates + 1):
        if time.perf_counter() > deadline:
            return {'status': 'target_discovery_wall_budget', 'initial_state': list(initial),
                    'updates': update, 'target': None}
        next_state = tuple(pins[i] if i in pins else rule.evaluate(state)
                           for i, rule in enumerate(rules))
        if next_state == state:
            return {'status': 'fixed_target_found', 'initial_state': list(initial),
                    'updates': update, 'target': list(state),
                    'candidate_orbit': [list(x) for x in trajectory]}
        if update == maximum_updates:
            break
        if next_state in seen:
            start = seen[next_state]
            cycle = trajectory[start:]
            return {'status': 'target_not_found_by_declared_orbit', 'initial_state': list(initial),
                    'updates': update + 1, 'target': None, 'cycle_start_update': start,
                    'cycle_length': len(cycle), 'cycle': [list(x) for x in cycle]}
        state = next_state
        trajectory.append(state)
        seen[state] = update + 1
    return {'status': 'target_discovery_update_budget', 'initial_state': list(initial),
            'updates': maximum_updates, 'target': None}


class BudgetExceeded(Exception):
    pass


class WorkBudget:
    def __init__(self, maximum, deadline):
        self.maximum, self.deadline = maximum, deadline
        self.used = 0

    def check(self):
        if time.perf_counter() > self.deadline:
            raise BudgetExceeded('wall_budget')

    def reserve(self, amount):
        self.check()
        if self.used + amount > self.maximum:
            raise BudgetExceeded('transition_budget')
        self.used += amount


def successor_table(model, budget):
    count = 1 << model.n
    budget.reserve(count)
    table = []
    for code in range(count):
        if code % 1024 == 0:
            budget.check()
        state = tuple((code >> i) & 1 for i in range(model.n))
        table.append(sum(rule.evaluate(state) << i for i, rule in enumerate(model.rules)))
    if table[0] != 0:
        raise AssertionError('Residual target is not fixed')
    return table


def classify_pins(successor, r, pins, budget):
    mask = sum(1 << p for p in pins)
    compatible_count = 1 << (r - len(pins))
    budget.reserve(compatible_count)
    # Classify the compatible functional graph with cached forward paths.
    distance = {0: 0}
    failure_cycle = None
    convergent_count = 1
    for start in range(1 << r):
        if start % 1024 == 0:
            budget.check()
        if start & mask or start in distance:
            continue
        path, position = [], {}
        state = start
        while state not in distance and state not in position:
            if len(path) % 1024 == 0:
                budget.check()
            position[state] = len(path)
            path.append(state)
            state = successor[state] & ~mask
        if state in position:
            cycle = path[position[state]:]
            if failure_cycle is None:
                failure_cycle = cycle
            for code in path:
                distance[code] = -1
        else:
            d = distance[state]
            for code in reversed(path):
                d = d + 1 if d >= 0 else -1
                distance[code] = d
                if d >= 0:
                    convergent_count += 1
    succeeds = failure_cycle is None
    result = {'cardinality': len(pins), 'pin_indices_local': list(pins),
              'compatible_states': compatible_count, 'convergent_states': convergent_count,
              'globally_convergent': succeeds}
    if succeeds:
        maximum = max(distance.values())
        witness = min(code for code, d in distance.items() if d == maximum)
        route = [witness]
        while route[-1] != 0:
            route.append(successor[route[-1]] & ~mask)
        result.update({'max_residual_hitting_time': maximum,
                       'longest_residual_trajectory': route})
    else:
        result['failure_cycle'] = failure_cycle
    return result


def exact_search(model, baseline, node_order, budgets, deadline):
    trace = support_trace(model, baseline)
    residual_model, residual_nodes = physical_residual_model(model, trace[-1])
    names = [node_order[i] for i in residual_nodes]
    r = len(residual_nodes)
    output = {'baseline_support_trace': [sorted(s) for s in trace],
              'baseline_support_sizes': [len(s) for s in trace],
              'residual_node_indices': list(residual_nodes), 'residual_node_names': names,
              'residual_size': r, 'residual_successor_table': None}
    search = {'lower_bound': 0, 'upper_bound': r, 'status': 'not_started',
              'completed_cardinality_levels': [], 'candidate_results': [],
              'all_minimum_additional_sets': [], 'all_minimum_sets_complete': False,
              'feasible_additional_set': names.copy(),
              'feasible_additional_values': {node_order[i]: model.target[i] for i in residual_nodes},
              'transition_evaluations_reserved': 0,
              'work_definition': 'One unit per original residual successor evaluated plus one per compatible state in each candidate pin-set classification. Reserved work includes interrupted candidate work conservatively.',
              'feasible_upper_bound_reason': 'Pinning every original residual coordinate is globally sufficient by physical residual equivalence.'}
    output['search'] = search
    if r == 0:
        search.update({'status': 'exact_empty_residual', 'all_minimum_sets_complete': True,
                       'all_minimum_additional_sets': [[]], 'completed_cardinality_levels': [0]})
        output['first_feasible_full_support_trace'] = output['baseline_support_trace']
        output['global_convergence_upper_bound'] = len(trace) - 1
        output['global_convergence_upper_bound_by_empty_full_support'] = len(trace) - 1
        return output
    if r > budgets['max_residual_nodes']:
        certificate = pure_copy_certificate(model, baseline)
        if certificate is not None:
            chosen = [node_order[i] for i in certificate['one_optimal_additional_set']]
            count = certificate['conditional_optimum']
            search.update({'status': 'exact_copy_certificate', 'lower_bound': count,
                           'upper_bound': count, 'feasible_additional_set': chosen,
                           'feasible_additional_values': {name: model.target[node_order.index(name)] for name in chosen},
                           'copy_certificate': certificate})
        else:
            search['status'] = 'residual_size_cap'
        return output
    budget = WorkBudget(budgets['max_transition_evaluations'], deadline)
    try:
        successor = successor_table(residual_model, budget)
        output['residual_successor_table'] = successor
        for cardinality in range(min(r, budgets['max_additional_cost']) + 1):
            successes = []
            for pins in itertools.combinations(range(r), cardinality):
                result = classify_pins(successor, r, pins, budget)
                result['pin_names'] = [names[i] for i in pins]
                result['pin_values'] = {names[i]: model.target[residual_nodes[i]] for i in pins}
                search['candidate_results'].append(result)
                if result['globally_convergent']:
                    successes.append(result)
                    search['upper_bound'] = cardinality
                    search['feasible_additional_set'] = result['pin_names']
                    search['feasible_additional_values'] = result['pin_values']
            search['completed_cardinality_levels'].append(cardinality)
            if successes:
                search.update({'status': 'exact_residual_enumeration',
                               'lower_bound': cardinality, 'upper_bound': cardinality,
                               'all_minimum_sets_complete': True,
                               'all_minimum_additional_sets': [x['pin_names'] for x in successes]})
                first = successes[0]
                search['feasible_additional_set'] = first['pin_names']
                search['feasible_additional_values'] = first['pin_values']
                chosen = set(baseline) | {node_order.index(name) for name in first['pin_names']}
                output['first_feasible_full_support_trace'] = [sorted(s) for s in support_trace(model, chosen)]
                output['global_convergence_upper_bound'] = len(trace) - 1 + first['max_residual_hitting_time']
                output['global_convergence_upper_bound_by_residual_reduction'] = output['global_convergence_upper_bound']
                break
            search['lower_bound'] = cardinality + 1
        else:
            search['status'] = 'additional_cost_cap'
    except BudgetExceeded as error:
        search['status'] = str(error)
        if search['upper_bound'] == search['lower_bound']:
            search['status'] += '_optimum_known_minimum_set_list_incomplete'
    search['transition_evaluations_reserved'] = budget.used
    return output


def run_context(context, protocol):
    started = time.perf_counter()
    cpu_started = time.process_time()
    budgets = protocol['budgets']
    deadline = started + budgets['wall_seconds_per_case']
    rules, names, dynamic, inputs, formulas, source_info = load_network(context['model_id'])
    expected_source = protocol['models'][context['model_id']]['bnet_sha256']
    if source_info['model_bnet_sha256'] != expected_source:
        raise ValueError('BNET source differs from the frozen protocol hash')
    assignment = context['input_assignment']
    discovery = discover_target(rules, names, inputs, assignment,
                                budgets['target_orbit_max_updates'], deadline)
    result = {'schema': 'LargerNetworkCertificate_v1', 'context_id': context['id'],
              'model_id': context['model_id'], 'context': context,
              'source': source_info, 'node_order': list(names),
              'dynamic_node_names': list(dynamic), 'input_node_names': list(inputs),
              'baseline_input_assignment': assignment, 'discovery': discovery,
              'eligibility': 'All regulated coordinates, including outputs and phenotypic nodes, are eligible for additional target-value pins. No gene-only or molecular-actuator restriction is imposed.',
              'dynamics': 'Synchronous Boolean update; omitted input rules explicitly completed as identities; static target pins applied from time zero.',
              'target_scope': 'Computational fixed point found by the declared orbit for the frozen input context; not asserted to be a validated biological phenotype.',
              'budget': budgets}
    if discovery['target'] is not None:
        target = tuple(discovery['target'])
        model = BooleanModel(rules, target)
        target_map = dict(zip(names, target))
        if not all(evaluate(formulas[name], target_map) == target_map[name] for name in names):
            raise AssertionError('Candidate target failed independent source-expression evaluation')
        baseline = tuple(names.index(name) for name in inputs)
        result.update({'target_by_name': target_map,
                       'target_fixed_checked': True,
                       'baseline_pin_indices': list(baseline),
                       'identity_rows': [names[i] for i in identity_rows(model)],
                       'baseline_pins_mandatory': True})
        if not set(baseline) <= set(identity_rows(model)):
            raise AssertionError('Completed baseline inputs are not identity rows')
        result.update(exact_search(model, baseline, names, budgets, deadline))
        # Always expose the chosen feasible control, including cap-limited runs.
        feasible = result['search']['feasible_additional_set']
        complete_pins = set(baseline) | {names.index(name) for name in feasible}
        full_trace = support_trace(model, complete_pins)
        result['first_feasible_full_support_trace'] = [sorted(s) for s in full_trace]
        result['first_feasible_full_support_sizes'] = [len(s) for s in full_trace]
        if not full_trace[-1]:
            cube_bound = len(full_trace) - 1
            result['global_convergence_upper_bound_by_empty_full_support'] = cube_bound
            result['global_convergence_upper_bound'] = min(
                result.get('global_convergence_upper_bound', cube_bound), cube_bound)
        result['transient_bound_interpretation'] = 'A proved uniform upper bound from an empty full support and/or baseline capture plus exact residual hitting time; not claimed to be the attained worst transient of the original full network.'
        result['additional_cost_bounds'] = [result['search']['lower_bound'], result['search']['upper_bound']]
        result['baseline_cost'] = len(inputs)
        result['total_cost_bounds_counting_inputs'] = [len(inputs) + x for x in result['additional_cost_bounds']]
        result['total_cost_interpretation'] = 'Identity inputs are individually mandatory for convergence from every compatible initial state before input clamping; total cost adds those mandatory input clamps to the conditional additional optimum/bounds.'
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = int(rss if sys.platform == 'darwin' else rss * 1024)
    metrics = {'analysis_wall_seconds': time.perf_counter() - started,
               'analysis_cpu_seconds': time.process_time() - cpu_started,
               'peak_worker_rss_bytes': rss_bytes,
               'wall_cpu_scope': 'From model parsing through target discovery, support reduction, residual search, and certificate object construction; excludes interpreter startup, imports, and file serialization.',
               'peak_memory_scope': 'Peak resident-set size of the isolated context worker process, including interpreter/import overhead. Linux ru_maxrss converted from KiB to bytes; macOS value is already bytes.',
               'python': sys.version.split()[0], 'platform': platform.platform(),
               'processor': platform.processor(), 'measurement_repeats': 1}
    return result, metrics


def summary_row(result, resources):
    search = result.get('search', {})
    discovery = result['discovery']
    return {'context_id': result['context_id'], 'model_id': result['model_id'],
            'classification': result['context'].get('classification', ''),
            'nodes': len(result['node_order']), 'dynamic_nodes': len(result['dynamic_node_names']),
            'baseline_inputs': len(result['input_node_names']),
            'target_status': discovery['status'], 'residual_size': result.get('residual_size', ''),
            'additional_lower': search.get('lower_bound', ''), 'additional_upper': search.get('upper_bound', ''),
            'total_lower': result.get('total_cost_bounds_counting_inputs', ['', ''])[0],
            'total_upper': result.get('total_cost_bounds_counting_inputs', ['', ''])[1],
            'status': search.get('status', discovery['status']),
            'all_optima_complete': search.get('all_minimum_sets_complete', False),
            'minimum_sets_count': (len(search.get('all_minimum_additional_sets', []))
                                   if search.get('all_minimum_sets_complete', False) else ''),
            'successful_candidates_found': sum(bool(item['globally_convergent'])
                                                for item in search.get('candidate_results', [])),
            'candidate_sets_checked': len(search.get('candidate_results', [])),
            'work_reserved': search.get('transition_evaluations_reserved', 0),
            'uniform_convergence_bound': result.get('global_convergence_upper_bound', ''),
            'analysis_wall_seconds': resources['analysis_wall_seconds'],
            'analysis_cpu_seconds': resources['analysis_cpu_seconds'],
            'peak_worker_rss_bytes': resources['peak_worker_rss_bytes']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--worker-context', help=argparse.SUPPRESS)
    args = parser.parse_args()
    protocol_path = Path(args.protocol).resolve()
    if digest(protocol_path) != FROZEN_PROTOCOL_SHA256:
        parser.error('Protocol bytes do not match the frozen protocol SHA256')
    protocol = json.loads(protocol_path.read_text())
    if not protocol.get('frozen', False):
        parser.error('Protocol must be explicitly frozen before outcome execution')
    contexts = protocol['contexts']
    ids = [item['id'] for item in contexts]
    if len(set(ids)) != len(ids):
        parser.error('Context IDs must be unique')
    output = Path(args.output).resolve()
    if output.exists():
        parser.error('Refusing to overwrite an existing output directory')
    output.mkdir(parents=True)
    if args.worker_context is not None:
        context = next(item for item in contexts if item['id'] == args.worker_context)
        result, metrics = run_context(context, protocol)
        write_json(output/'certificate.json', result)
        write_json(output/'resources.json', metrics)
        return
    summaries = []
    for context in contexts:
        if '/' in context['id'] or '\\' in context['id'] or context['id'] in ('.', '..'):
            parser.error('Unsafe context ID')
        destination = output/context['id']
        started = time.perf_counter()
        completed = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                    '--protocol', str(protocol_path), '--output', str(destination),
                                    '--worker-context', context['id']], capture_output=True, text=True)
        if completed.returncode:
            write_json(output/'execution_failure.json', {'context_id': context['id'],
                       'returncode': completed.returncode, 'stderr': completed.stderr})
            raise RuntimeError('Context worker failed; see execution_failure.json')
        resources = json.loads((destination/'resources.json').read_text())
        resources['isolated_process_wall_seconds_including_startup_and_export'] = time.perf_counter() - started
        write_json(destination/'resources.json', resources)
        result = json.loads((destination/'certificate.json').read_text())
        summaries.append(summary_row(result, resources))
        print(context['id'], summaries[-1]['status'],
              summaries[-1]['additional_lower'], summaries[-1]['additional_upper'], flush=True)
    with (output/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    source_files = [p for parent in (ROOT/'src', ROOT/'data') for p in parent.rglob('*')
                    if p.is_file() and '__pycache__' not in p.parts]
    write_json(output/'run_manifest.json', {
        'schema': 'LargerNetworkBenchmarkRun_v1', 'protocol_sha256': digest(protocol_path),
        'protocol_file_name': protocol_path.name,
        'contexts_executed_in_frozen_order': ids,
        'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in sorted(source_files)},
        'scientific_result_sha256': {str(p.relative_to(output)): digest(p)
                                    for p in sorted(output.rglob('certificate.json'))},
        'timing_note': 'Resource measurements and summary timing columns vary between reruns. Scientific certificate files are deterministic unless a wall-time cutoff changes the reached search point.'})


if __name__ == '__main__':
    main()
