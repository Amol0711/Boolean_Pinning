"""Independent finite-state oracle: no support reduction or cycle certificate.

Enumerates every compatible PHYSICAL state, builds its synchronous successor,
and classifies the full functional graph. Intended only for small networks.
"""


def all_subsets(nodes):
    nodes = tuple(sorted(nodes))
    return [frozenset(nodes[j] for j in range(len(nodes)) if (mask >> j) & 1)
            for mask in range(1 << len(nodes))]


def reachable_support_series(rules, target, pins, steps, state_limit=1 << 18):
    """Exact successive images of the full compatible physical state set.

    Counts the UNION of coordinates that deviate in at least one reachable
    state; this is not the largest Hamming distance in an individual state.
    """
    n = len(rules)
    pins = frozenset(pins)
    free = tuple(i for i in range(n) if i not in pins)
    if 1 << len(free) > state_limit:
        raise ValueError("Reachable-image oracle limit exceeded")
    target_integer = sum(x << i for i, x in enumerate(target))
    states = set()
    for code in range(1 << len(free)):
        state = target_integer
        for j, i in enumerate(free):
            state = (state & ~(1 << i)) | (((code >> j) & 1) << i)
        states.add(state)
    records = []
    for time in range(steps + 1):
        union = 0
        for state in states:
            union |= state ^ target_integer
        records.append({"update": time, "distinct_reachable_states": len(states),
                        "exact_reachable_union_support_size": union.bit_count(),
                        "exact_reachable_union_support_zero_based":
                        [i for i in range(n) if (union >> i) & 1]})
        if time < steps:
            next_states = set()
            for state in states:
                successor = target_integer
                for i in free:
                    corner = sum(((state >> p) & 1) << j
                                 for j, p in enumerate(rules[i].parents))
                    successor = ((successor & ~(1 << i))
                                 | (rules[i].table[corner] << i))
                next_states.add(successor)
            states = next_states
    return records


def classify(rules, target, pins=(), state_limit=1 << 18):
    n = len(rules)
    pins = frozenset(pins)
    free = tuple(i for i in range(n) if i not in pins)
    count = 1 << len(free)
    if count > state_limit:
        raise ValueError("Full-state oracle limit exceeded")
    target_code = sum(target[i] << j for j, i in enumerate(free))
    successors = []
    for code in range(count):
        physical = list(target)
        for j, i in enumerate(free):
            physical[i] = (code >> j) & 1
        next_code = 0
        for j, i in enumerate(free):
            parents, table = rules[i].parents, rules[i].table
            # Written independently of LocalRule.evaluate and model.step.
            multiplier, corner = 1, 0
            for parent in parents:
                corner += multiplier * physical[parent]
                multiplier *= 2
            next_code += (1 << j) * table[corner]
        successors.append(next_code)
    if successors[target_code] != target_code:
        raise AssertionError("Target was not fixed in the physical oracle")
    # Each entry is the number of updates to the target, or -1 for failure.
    distances = {target_code: 0}
    failure_cycle = None
    for start in range(count):
        if start in distances:
            continue
        path, local = [], {}
        state = start
        while state not in distances and state not in local:
            local[state] = len(path)
            path.append(state)
            state = successors[state]
        if state in local:
            cycle = path[local[state]:]
            if failure_cycle is None:
                failure_cycle = cycle
            for code in path:
                distances[code] = -1
        else:
            distance = distances[state]
            for code in reversed(path):
                distance = distance + 1 if distance >= 0 else -1
                distances[code] = distance
    maximum = max(distances.values())
    witness = min(code for code, distance in distances.items() if distance == maximum)

    def decode(code):
        result = list(target)
        for j, i in enumerate(free):
            result[i] = (code >> j) & 1
        return result

    trajectory = [witness]
    while trajectory[-1] != target_code:
        trajectory.append(successors[trajectory[-1]])
    return {
        "globally_convergent": failure_cycle is None,
        "compatible_states_enumerated": count,
        "maximum_convergent_transient": maximum,
        "all_state_maximum_transient": maximum if failure_cycle is None else None,
        "longest_convergent_trajectory": [decode(code) for code in trajectory],
        "failure_cycle": None if failure_cycle is None else [decode(c) for c in failure_cycle],
        "target": list(target),
        "pins": sorted(pins),
    }
