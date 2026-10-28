"""Independent BNET/cube/functional-graph oracle (standard library only).

This module neither imports the production engine nor reads compiled tables.
It parses the source Boolean grammar directly, using recursive descent with
precedence ! > & > |. Controlled graph convergence is checked by reverse
reachability from zero, rather than by the engine's trajectory algorithm.
"""
from collections import deque
from itertools import combinations, product
from pathlib import Path
import re


def parse_expression(text):
    tokens = re.findall(r"[A-Za-z_][A-Za-z_0-9]*|[01!&|()]", text)
    if ''.join(tokens) != re.sub(r'\s+', '', text):
        raise ValueError('Unsupported BNET token')
    position = 0

    def atom():
        nonlocal position
        if position == len(tokens):
            raise ValueError('Truncated expression')
        token = tokens[position]
        position += 1
        if token == '!':
            return ('not', atom())
        if token == '(':
            value = disjunction()
            if position == len(tokens) or tokens[position] != ')':
                raise ValueError('Unbalanced expression')
            position += 1
            return value
        if token in ('0', '1'):
            return ('constant', int(token))
        if re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', token):
            return ('variable', token)
        raise ValueError('Unexpected token')

    def conjunction():
        nonlocal position
        value = atom()
        while position < len(tokens) and tokens[position] == '&':
            position += 1
            value = ('and', value, atom())
        return value

    def disjunction():
        nonlocal position
        value = conjunction()
        while position < len(tokens) and tokens[position] == '|':
            position += 1
            value = ('or', value, conjunction())
        return value

    result = disjunction()
    if position != len(tokens):
        raise ValueError('Trailing tokens')
    return result


def variables(tree):
    if tree[0] == 'variable':
        return {tree[1]}
    if tree[0] == 'constant':
        return set()
    return set().union(*(variables(child) for child in tree[1:]))


def evaluate(tree, values):
    operator = tree[0]
    if operator == 'variable':
        return values[tree[1]]
    if operator == 'constant':
        return tree[1]
    if operator == 'not':
        return 1 - evaluate(tree[1], values)
    left = evaluate(tree[1], values)
    right = evaluate(tree[2], values)
    return left * right if operator == 'and' else max(left, right)


def read_bnet(path):
    rows = []
    names = set()
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('targets,'):
            continue
        name, text = line.split(',', 1)
        name = name.strip()
        if name in names:
            raise ValueError('Repeated target')
        names.add(name)
        tree = parse_expression(text)
        scope = sorted(variables(tree))
        truth = {bits: evaluate(tree, dict(zip(scope, bits)))
                 for bits in product((0, 1), repeat=len(scope))}
        rows.append((name, scope, truth))
    external = sorted(set().union(*(set(row[1]) for row in rows)) - names)
    regulated = [row[0] for row in rows]
    rows.extend((name, [name], {(0,): 0, (1,): 1}) for name in external)
    return rows, regulated, external


def step(rows, state):
    return {name: table[tuple(state[parent] for parent in scope)]
            for name, scope, table in rows}


def support_trace(rows, target, pins):
    pinned = set(pins)
    support = set(target) - pinned
    trace = [support]
    while True:
        following = set()
        for name, scope, table in rows:
            if name in pinned:
                continue
            for bits, value in table.items():
                if value != target[name] and all(
                    parent in support or bit == target[parent]
                    for parent, bit in zip(scope, bits)
                ):
                    following.add(name)
                    break
        if not following <= support:
            raise AssertionError('Support sequence is not decreasing')
        if following == support:
            return trace
        trace.append(following)
        support = following


def residual_successors(rows, target, nodes):
    order = {name: j for j, name in enumerate(nodes)}
    rowmap = {row[0]: row for row in rows}
    successors = []
    for state in range(1 << len(nodes)):
        following = 0
        for j, name in enumerate(nodes):
            _, scope, table = rowmap[name]
            corner = tuple(target[parent] ^ ((state >> order[parent]) & 1)
                           if parent in order else target[parent]
                           for parent in scope)
            following |= (table[corner] ^ target[name]) << j
        successors.append(following)
    assert successors[0] == 0
    return successors


def controlled_graph(successors, pin_indices):
    """Return exact graph statistics and one nonzero cycle for failure."""
    mask = sum(1 << j for j in pin_indices)
    compatible = [s for s in range(len(successors)) if not s & mask]
    reverse = [[] for _ in successors]
    for state in compatible:
        reverse[successors[state] & ~mask].append(state)
    distance = {0: 0}
    queue = deque([0])
    while queue:
        state = queue.popleft()
        for predecessor in reverse[state]:
            if predecessor not in distance:
                distance[predecessor] = distance[state] + 1
                queue.append(predecessor)
    result = {'compatible_states': len(compatible),
              'convergent_states': len(distance),
              'success': len(compatible) == len(distance),
              'maximum_hitting_time': max(distance.values())}
    if not result['success']:
        state = next(s for s in compatible if s not in distance)
        result['failure_initial_state'] = state
        trail = []; seen = {}
        while state not in seen:
            seen[state] = len(trail); trail.append(state)
            state = successors[state] & ~mask
        result['failure_cycle'] = trail[seen[state]:]
        assert 0 not in result['failure_cycle']
    return result


def exact_optimum(successors, nodes, eligible=None, max_cardinality=None):
    eligible = list(range(len(nodes))) if eligible is None else [
        j for j, name in enumerate(nodes) if name in eligible]
    max_cardinality = len(eligible) if max_cardinality is None else max_cardinality
    records = []; by_cardinality = []
    for size in range(max_cardinality + 1):
        successful = []
        count = 0
        for subset in combinations(eligible, size):
            count += 1
            checked = controlled_graph(successors, subset)
            checked.update(pin_indices=list(subset),
                           pins=[nodes[j] for j in subset], cardinality=size)
            records.append(checked)
            if checked['success']:
                successful.append(checked)
        by_cardinality.append({'cardinality': size, 'tested': count,
                               'successful': len(successful)})
        if successful:
            return {'optimum': size, 'all_optimal_sets': successful,
                    'by_cardinality': by_cardinality, 'candidates': records}
    return {'optimum': None, 'by_cardinality': by_cardinality,
            'candidates': records}
