"""Independent JointRules reference routines; no production modules are imported.

The prescribed NumPy stream is reconstructed directly from the fixed protocol.
Support uses synchronous set iteration; SCCs use Tarjan; physical success uses
indegree peeling, so its reference algorithms differ from production routines.
"""
from __future__ import annotations
from collections import deque
import hashlib
import itertools
import sys

import numpy as np

sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))


def reconstruct(protocol, family, n, replicate):
    spec = protocol['families'][family]
    counts = []
    for row in spec['types']:
        numerator, denominator = row['fraction']
        assert (n*numerator) % denominator == 0
        counts.append(n*numerator//denominator)
    assert sum(counts) == n
    degrees = np.concatenate([np.full(c, d, dtype=np.int64)
                              for c, d in zip(counts, [t['degree'] for t in spec['types']])])
    sources = np.repeat(np.arange(n, dtype=np.int64), degrees)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([
        protocol['master_seed'], spec['id'], n, replicate])))
    parents = rng.permutation(sources).reshape(n, 2)
    codes = np.concatenate([rng.choice(8,size=count,
        p=np.asarray(row['weights'],dtype=np.float64)/row['probability_denominator'])
        for count,row in zip(counts,spec['types'])])
    corners = np.array([(codes >> j) & 1 for j in range(3)],dtype=np.uint8).T.copy()
    digest = hashlib.sha256(parents.astype('<i8').tobytes() + corners.tobytes()).hexdigest()
    return parents, corners, digest


def literal_support(parents, corners):
    """Direct synchronous application of the stated nominal-slot support map."""
    active = np.ones(len(parents), dtype=bool)
    b1, b2, b12 = (corners[:, j].astype(bool) for j in range(3))
    iterations = 0
    while True:
        a, b = active[parents[:, 0]], active[parents[:, 1]]
        updated = (a & b1) | (b & b2) | (a & b & b12)
        iterations += 1
        if np.array_equal(active, updated):
            return np.flatnonzero(active).tolist(), iterations
        assert not np.any(updated & ~active)
        active = updated


def export_rows(parents, corners, residual):
    return [[int(v), int(parents[v, 0]), int(parents[v, 1]),
             *map(int, corners[v])] for v in residual]


def graph_certificate(rows):
    table = {row[0]: row for row in rows}
    assert len(table) == len(rows)
    vertices = sorted(table)
    adjacency = {v: sorted(set(u for u in table[v][1:3] if u in table))
                 for v in vertices}
    # Reversing every edge preserves SCCs. Here receiver -> parent is convenient.
    index, low, stack, on_stack, groups = {}, {}, [], set(), []

    def visit(v):
        index[v] = low[v] = len(index)
        stack.append(v)
        on_stack.add(v)
        for u in adjacency[v]:
            if u not in index:
                visit(u)
                low[v] = min(low[v], low[u])
            elif u in on_stack:
                low[v] = min(low[v], index[u])
        if low[v] == index[v]:
            group = []
            while True:
                u = stack.pop()
                on_stack.remove(u)
                group.append(u)
                if u == v:
                    break
            groups.append(sorted(group))

    for v in vertices:
        if v not in index:
            visit(v)
    components = []
    for nodes in sorted(groups, key=lambda c: c[0]):
        node_set = set(nodes)
        cyclic = len(nodes) > 1 or nodes[0] in adjacency[nodes[0]]
        marked = cyclic
        slots = []
        for v in nodes:
            row = table[v]
            internal = [j for j, u in enumerate(row[1:3]) if u in node_set]
            slots.extend((row[1+j], v, j) for j in internal)
            if len(internal) != 1 or row[3+internal[0]] != 1:
                marked = False
        components.append({'nodes': nodes, 'cyclic': cyclic,
                           'marked_slot_cycle': marked, 'internal_slots': slots})
    double_rows = sum(sum(u in table for u in row[1:3]) == 2 for row in rows)
    G = all(sum(u in table for u in row[1:3]) == 1 for row in rows)
    G_prime = all(not c['cyclic'] or c['marked_slot_cycle'] for c in components)
    assert not G or G_prime
    return {'G': G, 'G_prime': G_prime, 'double_rows': double_rows,
            'components': components,
            'cyclic_component_count': sum(c['cyclic'] for c in components),
            'largest_complex': max((len(c['nodes']) for c in components
                                    if c['cyclic'] and not c['marked_slot_cycle']), default=0)}


def physical_successors(rows, nodes):
    """Full physical map on listed coordinates, exterior coordinates fixed zero."""
    lookup = {row[0]: row for row in rows}
    position = {v: j for j, v in enumerate(nodes)}
    states = np.arange(1 << len(nodes), dtype=np.uint64)
    successor = np.zeros(len(states), dtype=np.uint64)
    for j, v in enumerate(nodes):
        row = lookup[v]
        a = (states >> position[row[1]]) & 1 if row[1] in position else 0
        b = (states >> position[row[2]]) & 1 if row[2] in position else 0
        truth = np.array([0, *row[3:]], dtype=np.uint64)
        successor |= truth[a + 2*b] << j
    assert successor[0] == 0
    return successor


def succeeds(successor, pin_mask):
    """A compatible functional graph converges to zero iff no other cycle remains."""
    count = len(successor)
    keep = (count - 1) ^ pin_mask
    compatible = np.flatnonzero((np.arange(count, dtype=np.uint64) & pin_mask) == 0)
    image = (successor & keep).astype(np.int64)
    # A nonzero compatible fixed point is already a complete failure witness.
    # This vector test is independent of production trajectory traversal.
    if np.any(image[compatible[1:]] == compatible[1:]):
        return False
    indegree = np.bincount(image[compatible], minlength=count)
    queue = deque(int(x) for x in compatible if indegree[x] == 0)
    removed = 0
    while queue:
        v = queue.popleft()
        removed += 1
        u = int(image[v])
        indegree[u] -= 1
        if indegree[u] == 0:
            queue.append(u)
    return removed == len(compatible) - 1


def exact_cost(rows, nodes, max_size=18, max_candidates=100000,
               max_compatible_states=20_000_000):
    if len(nodes) > max_size:
        return {'complete': False, 'reason': 'independent_size_guard', 'candidates': 0,
                'states': 0}
    successor = physical_successors(rows, nodes)
    checked = states = 0
    for cost in range(len(nodes) + 1):
        for pins in itertools.combinations(range(len(nodes)), cost):
            candidate_states = 1 << (len(nodes)-cost)
            if checked >= max_candidates or states+candidate_states > max_compatible_states:
                return {'complete': False, 'reason': 'independent_enumeration_guard',
                        'lower': cost, 'candidates': checked, 'states': states}
            mask = sum(1 << j for j in pins)
            checked += 1
            states += candidate_states
            if succeeds(successor, mask):
                return {'complete': True, 'cost': cost,
                        'pins': [nodes[j] for j in pins],
                        'candidates': checked, 'states': states}
    raise AssertionError('all coordinates pinned must succeed')


def witness_success(rows, nodes, pins, max_size=18):
    # Physical support capture after pins can itself certify success and avoids
    # exponential enumeration of transient coordinates in a witness check.
    alive = set(nodes) - set(pins)
    lookup = {row[0]: row for row in rows}
    while True:
        new = set()
        for v in alive:
            _, a, b, b1, b2, b12 = lookup[v]
            if a == b:
                survives = a in alive and b12
            else:
                survives = (a in alive and b1) or (b in alive and b2) or (
                    a in alive and b in alive and b12)
            if survives:
                new.add(v)
        if new == alive:
            break
        alive = new
    if not alive:
        return {'checked': True, 'method': 'physical_support_capture',
                'states': 0, 'success': True}
    nodes = sorted(alive)
    if len(nodes) > max_size:
        return {'checked': False, 'reason': 'independent_size_guard',
                'remaining_dimension': len(nodes)}
    successor = physical_successors(rows, nodes)
    return {'checked': True, 'method': 'physical_state_enumeration',
            'states': 1 << len(nodes), 'success': succeeds(successor, 0)}


def normalized_key(rows, nodes):
    lookup = {row[0]: row for row in rows}
    position = {v: j for j, v in enumerate(nodes)}
    return tuple((position.get(lookup[v][1], -1), position.get(lookup[v][2], -1),
                  *lookup[v][3:]) for v in nodes)


def baseline_refined_certificate(rows):
    """Independently recover the former double-row deletion/cycle interval."""
    vertices = {row[0] for row in rows}
    doubles = {row[0] for row in rows if row[1] in vertices and row[2] in vertices}
    alive = vertices - doubles
    while True:
        updated = {v for v, a, b, b1, b2, b12 in rows if v in alive and (
            (a in alive and b1) or (b in alive and b2) or
            (a in alive and b in alive and b12))}
        if updated == alive:
            break
        alive = updated
    surviving = [row for row in rows if row[0] in alive]
    graph = graph_certificate(surviving)
    assert graph['G']
    cycles = [c for c in graph['components'] if c['cyclic']]
    # The lower witnesses consist only of original unary residual coordinates.
    for c in cycles:
        assert not (set(c['nodes']) & doubles)
    pins = sorted(doubles | {min(c['nodes']) for c in cycles})
    return {'lower': len(cycles), 'upper': len(cycles)+len(doubles),
            'pins': pins, 'double_vertices': sorted(doubles)}


def self_check():
    """Exhaustive tiny oracle check against a separate literal trajectory evaluator."""
    networks = pin_sets = compatible_states = 0
    for parents in itertools.product(range(2), repeat=4):
        for truth in itertools.product(range(2), repeat=6):
            rows = [[i, parents[2*i], parents[2*i+1], *truth[3*i:3*i+3]] for i in range(2)]
            successor = physical_successors(rows, [0, 1])
            for mask in range(4):
                expected = True
                for initial in range(4):
                    if initial & mask:
                        continue
                    compatible_states += 1
                    seen, state = set(), initial
                    while state and state not in seen:
                        seen.add(state)
                        result = 0
                        for i, a, b, b1, b2, b12 in rows:
                            x, y = (state >> a) & 1, (state >> b) & 1
                            value = 0 if not (x or y) else b1 if x and not y else b2 if y and not x else b12
                            if not (mask >> i) & 1:
                                result |= value << i
                        state = result
                    expected &= state == 0
                assert succeeds(successor, mask) == expected
                pin_sets += 1
            networks += 1
    return {'status': 'PASS', 'networks': networks, 'pin_sets': pin_sets,
            'compatible_states': compatible_states,
            'scope': 'All 1024 two-coordinate ordered two-slot zero-target networks; independent literal trajectory evaluator checks all physical pin sets.'}


if __name__ == '__main__':
    import json
    print(json.dumps(self_check(), indent=2, sort_keys=True))
