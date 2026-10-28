"""Source-rule and exact finite-state operations for the biological study."""
from __future__ import annotations
import ast
from dataclasses import dataclass
import itertools
import numpy as np


def parse_formula(text):
    tree = ast.parse(text.replace('!', ' not ').replace('&', ' and ').replace('|', ' or ').strip(), mode='eval')
    allowed = (ast.Expression, ast.Name, ast.Load, ast.BoolOp, ast.And, ast.Or,
               ast.UnaryOp, ast.Not, ast.Constant)
    if any(not isinstance(node, allowed) for node in ast.walk(tree)):
        raise ValueError('Unsupported Boolean syntax')
    if any(isinstance(node, ast.Constant) and node.value not in (0, 1) for node in ast.walk(tree)):
        raise ValueError('Unsupported Boolean constant')
    return tree


def expression_value(node, values):
    if isinstance(node, ast.Expression): return expression_value(node.body, values)
    if isinstance(node, ast.Name): return values[node.id]
    if isinstance(node, ast.Constant): return int(node.value)
    if isinstance(node, ast.UnaryOp): return 1 - expression_value(node.operand, values)
    if isinstance(node, ast.BoolOp):
        vals = [expression_value(child, values) for child in node.values]
        return int(all(vals) if isinstance(node.op, ast.And) else any(vals))
    raise ValueError('Unexpected AST node')


@dataclass
class Model:
    names: tuple
    dynamic: tuple
    inputs: tuple
    parents: tuple
    tables: tuple
    expressions: dict

    def step(self, state):
        return tuple(table[sum(state[u] << j for j, u in enumerate(row))]
                     for row, table in zip(self.parents, self.tables))


def read_model(path):
    expressions = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#') or line.replace(' ', '') == 'targets,factors':
            continue
        name, formula = line.split(',', 1)
        name = name.strip()
        if name in expressions: raise ValueError('Duplicate BNET row')
        expressions[name] = parse_formula(formula)
    dynamic = tuple(sorted(expressions))
    refs = {node.id for tree in expressions.values() for node in ast.walk(tree) if isinstance(node, ast.Name)}
    inputs = tuple(sorted(refs - set(dynamic)))
    names = dynamic + inputs
    pos = {v: i for i, v in enumerate(names)}
    for name in inputs: expressions[name] = parse_formula(name)
    parents, tables = [], []
    for name in names:
        scope = tuple(sorted({pos[node.id] for node in ast.walk(expressions[name]) if isinstance(node, ast.Name)}))
        if len(scope) > 20: raise ValueError('Declared source-arity guard exceeded')
        table = tuple(expression_value(expressions[name], {names[u]: (code >> j) & 1 for j, u in enumerate(scope)})
                      for code in range(1 << len(scope)))
        parents.append(scope); tables.append(table)
    return Model(names, dynamic, inputs, tuple(parents), tuple(tables), expressions)


def forced_values(model, assignment):
    if set(assignment) != set(model.inputs): raise ValueError('Context must assign every external input')
    if any(type(v) is not int or v not in (0, 1) for v in assignment.values()):
        raise ValueError('Context values must be Boolean integers')
    forced = {model.names.index(n): v for n, v in assignment.items()}
    rounds = []
    while True:
        additions = {}
        for i, (row, table) in enumerate(zip(model.parents, model.tables)):
            if i in forced: continue
            vals = {value for code, value in enumerate(table)
                    if all(u not in forced or ((code >> j) & 1) == forced[u] for j, u in enumerate(row))}
            if len(vals) == 1: additions[i] = vals.pop()
        if not additions: break
        forced.update(additions)
        rounds.append({model.names[i]: val for i, val in sorted(additions.items())})
    return forced, rounds


def enumerate_fixed_points(model, assignment, chunk_size=8192):
    forced, rounds = forced_values(model, assignment)
    free = [i for i in range(len(model.names)) if i not in forced]
    if len(free) > 24: raise ValueError('Declared remaining-fixed-point dimension guard exceeded')
    pos = {v: i for i, v in enumerate(free)}
    targets = []
    for start in range(0, 1 << len(free), chunk_size):
        codes = np.arange(start, min(1 << len(free), start + chunk_size), dtype=np.uint64)
        alive = np.ones(len(codes), dtype=bool)
        for v in free:
            pattern = np.zeros(len(codes), dtype=np.uint64)
            for j, u in enumerate(model.parents[v]):
                bits = ((codes >> pos[u]) & 1) if u in pos else forced[u]
                pattern |= np.asarray(bits, dtype=np.uint64) << j
            expected = ((codes >> pos[v]) & 1)
            alive &= np.asarray(model.tables[v], dtype=np.uint8)[pattern] == expected
        for code in codes[alive]:
            value = int(code)
            x = tuple(forced[i] if i in forced else (value >> pos[i]) & 1 for i in range(len(model.names)))
            assert model.step(x) == x
            # Also evaluate every original parsed expression, independently of
            # the source-table index, for each retained full target.
            env = dict(zip(model.names, x))
            assert all(expression_value(model.expressions[name], env) == x[i] for i, name in enumerate(model.names))
            targets.append(x)
    return sorted(targets), forced, free, rounds


def discover_zero_orbit(model, assignment, max_updates=2000):
    inputs = {model.names.index(n): v for n, v in assignment.items()}
    x = tuple(inputs.get(i, 0) for i in range(len(model.names)))
    seen, trajectory = {}, []
    for t in range(max_updates + 1):
        if x in seen:
            period = t - seen[x]
            return {'status': 'fixed' if period == 1 else 'cycle', 'period': period,
                    'cycle_start': seen[x], 'target': list(x) if period == 1 else None,
                    'cycle': [list(y) for y in trajectory[seen[x]:]]}
        seen[x] = t; trajectory.append(x)
        x = model.step(x)
    return {'status': 'update_cap', 'target': None}


def restricted_rows(parents, tables, selected, outside, output_xor=None):
    """Exact restriction with essential-variable testing on the whole cube.

    ``outside`` fixes physical variables outside selected. With output_xor,
    selected variables denote deviations from that same physical target.
    Returned parents are original coordinate indices, with all inessential
    restricted inputs removed; truth-table index bit order is parents order.
    """
    selected = tuple(sorted(selected)); selected_set = set(selected)
    answer = []
    for v in selected:
        free = tuple(u for u in parents[v] if u in selected_set)
        position = {u: k for k, u in enumerate(free)}
        table = []
        for code in range(1 << len(free)):
            pattern = 0
            for j, u in enumerate(parents[v]):
                if u in position:
                    val = (code >> position[u]) & 1
                    if output_xor is not None: val ^= outside[u]
                else: val = outside[u]
                pattern |= val << j
            val = tables[v][pattern]
            table.append(val ^ output_xor[v] if output_xor is not None else val)
        essential = [j for j in range(len(free)) if any(table[k] != table[k ^ (1 << j)] for k in range(len(table)))]
        compressed = [table[sum(((code >> k) & 1) << j for k, j in enumerate(essential))]
                      for code in range(1 << len(essential))]
        answer.append({'node': v, 'parents': [free[j] for j in essential], 'table': compressed})
    return answer


def successors(rows, nodes):
    """Synchronous error-state successor array; all exterior errors zero."""
    nodes = tuple(nodes); pos = {v: j for j, v in enumerate(nodes)}
    codes = np.arange(1 << len(nodes), dtype=np.uint64)
    nxt = np.zeros(len(codes), dtype=np.uint64)
    by_node = {row['node']: row for row in rows}
    for v in nodes:
        row = by_node[v]
        pattern = np.zeros(len(codes), dtype=np.uint64)
        for j, u in enumerate(row['parents']):
            if u in pos: pattern |= ((codes >> pos[u]) & 1) << j
        nxt |= np.asarray(row['table'], dtype=np.uint64)[pattern] << pos[v]
    assert nxt[0] == 0
    return nxt


def success(nxt, pin_mask):
    """Every compatible state reaches zero; return counterexample when false."""
    nstates = len(nxt)
    good = bytearray(nstates); good[0] = 1
    keep = (nstates - 1) ^ pin_mask
    for start in range(nstates):
        if start & pin_mask or good[start]: continue
        path, active, x = [], set(), start
        while not good[x]:
            if x in active: return False, start
            active.add(x); path.append(x)
            x = int(nxt[x]) & keep
        for x in path: good[x] = 1
    return True, None


def minimum_pinsets(nxt, nodes, allowed=None, known_minimum=None):
    nodes = tuple(nodes)
    if len(nodes) > 12:
        raise RuntimeError('Declared all-minimizer component-size cap of 12 exceeded')
    eligible = nodes if allowed is None else tuple(v for v in nodes if v in set(allowed))
    pos = {v: j for j, v in enumerate(nodes)}
    tested = 0
    lower, upper = (known_minimum, known_minimum) if known_minimum is not None else (0, len(eligible))
    for cost in range(lower, upper + 1):
        found = []
        for K in itertools.combinations(eligible, cost):
            tested += 1
            ok, _ = success(nxt, sum(1 << pos[v] for v in K))
            if ok: found.append(list(K))
        if found: return {'status': 'exact', 'cost': cost, 'minimum_pinsets': found,
                          'all_minimum_pinsets_enumerated': True, 'tested_candidates': tested}
    return {'status': 'infeasible', 'cost': None, 'minimum_pinsets': [],
            'all_eligible_subsets_rejected': known_minimum is None,
            'tested_candidates': tested}
