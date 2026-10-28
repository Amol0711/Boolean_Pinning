#!/usr/bin/env python3
"""Independent finite-state checks of SCC localization and the slot certificate.

Uses no production solver, graph library, or model parser. Coverage comprises
finite exhaustive cases and seeded physical-network cases.
Run from any directory. The JSON report is deterministic and contains no timings.
"""
import json
from itertools import product
from pathlib import Path
from random import Random


def components(n, edges):
    """SCCs by transitive closure, independently of a Tarjan/Kosaraju solver."""
    reach = [[i == j or (i, j) in edges for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    remaining = set(range(n))
    result = []
    while remaining:
        i = min(remaining)
        block = {j for j in remaining if reach[i][j] and reach[j][i]}
        result.append(block)
        remaining -= block
    return result


def dependencies(n, functions):
    return {(j, i) for i in range(n) for j in range(n)
            if any(functions[i][s] != functions[i][s ^ (1 << j)]
                   for s in range(1 << n))}


def succeeds(n, functions, pins):
    keep = ((1 << n) - 1) ^ pins
    for start in range(1 << n):
        if start & pins:
            continue
        state, seen = start, set()
        while state and state not in seen:
            seen.add(state)
            state = sum(f[state] << i for i, f in enumerate(functions)) & keep
        if state:
            return False
    return True


def restrict(functions, indices):
    return tuple(tuple(functions[i][sum(((x >> k) & 1) << j
                                       for k, j in enumerate(indices))]
                       for x in range(1 << len(indices))) for i in indices)


def check_localization(n, functions, edges):
    local = []
    for block in components(n, edges):
        indices = sorted(block)
        f_local = restrict(functions, indices)
        local.append((indices, [succeeds(len(indices), f_local, k)
                                for k in range(1 << len(indices))]))
    global_success = []
    for pins in range(1 << n):
        global_ok = succeeds(n, functions, pins)
        local_ok = all(success[sum(((pins >> j) & 1) << k
                                  for k, j in enumerate(indices))]
                       for indices, success in local)
        assert global_ok == local_ok, (functions, edges, pins, global_ok, local_ok)
        global_success.append(global_ok)
    for eligible in range(1 << n):
        global_min = min([pins.bit_count()
                          for pins, ok in enumerate(global_success)
                          if ok and pins & ~eligible == 0] or [float("inf")])
        local_sum = sum(min([pins.bit_count() for pins, ok in enumerate(success)
                             if ok and all(not (pins >> k & 1) or (eligible >> j & 1)
                                           for k, j in enumerate(indices))]
                            or [float("inf")]) for indices, success in local)
        assert global_min == local_sum, ("eligibility", functions, edges, eligible)
    return 1 << n


def check_slot_realization(n, parents, tables):
    functions = tuple(tuple(tables[i][((s >> parents[i][0]) & 1)
                                       + 2 * ((s >> parents[i][1]) & 1)]
                            for s in range(1 << n)) for i in range(n))

    def support(coordinates):
        return {i for i, (a, b) in enumerate(parents)
                if (tables[i][1] and a in coordinates)
                or (tables[i][2] and b in coordinates)
                or (tables[i][3] and a in coordinates and b in coordinates)}

    residual = set(range(n))
    while support(residual) != residual:
        residual = support(residual)
    edges = {(p, i) for i in residual for p in parents[i] if p in residual}
    blocks = [block for block in components(n, edges) if block <= residual]
    cyclic = [block for block in blocks
              if len(block) > 1 or any((i, i) in edges for i in block)]
    good_prime = all(all(sum(p in block for p in parents[i]) == 1
                        and all(p not in block or tables[i][1 << slot]
                                for slot, p in enumerate(parents[i]))
                        for i in block) for block in cyclic)

    # Count simple directed cycles with slot multiplicity; smallest vertex fixes
    # the rotation. No production marked-cycle or graph routines are imported.
    marked_count = 0
    for start in range(n):
        def visit(path):
            nonlocal marked_count
            u = path[-1]
            for v in range(n):
                for slot, parent in enumerate(parents[v]):
                    if parent != u or not tables[v][1 << slot]:
                        continue
                    if v == start:
                        marked_count += 1
                    elif v > start and v not in path:
                        visit(path + [v])
        visit([start])
    cost = min(pins.bit_count() for pins in range(1 << n)
               if succeeds(n, functions, pins))
    if good_prime:
        assert cost == len(cyclic) == marked_count, (parents, tables, residual)
    indices = sorted(residual)
    f_residual = restrict(functions, indices)
    cost_residual = min(pins.bit_count() for pins in range(1 << len(indices))
                        if succeeds(len(indices), f_residual, pins))
    assert cost == cost_residual, ("residual", parents, tables)
    return bool(good_prime)


def main():
    report = {"status": "PASS", "seed": 20261015,
              "semantics": "synchronous; fixed zero target; permanent zero pins; all compatible initial states",
              "scope": "finite exhaustive and seeded physical-state, reduction, component and eligibility checks"}
    pairs = pin_checks = 0
    for values in product(range(8), repeat=2):
        functions = tuple(tuple([0] + [(a >> i) & 1 for i in range(3)])
                          for a in values)
        required = dependencies(2, functions)
        optional = [e for e in product(range(2), repeat=2) if e not in required]
        for mask in range(1 << len(optional)):
            edges = required | {e for k, e in enumerate(optional) if mask >> k & 1}
            pin_checks += check_localization(2, functions, edges)
            pairs += 1
    report["exhaustive_two_variable_maps"] = {
        "maps": 64, "map_dependency_supergraph_pairs": pairs,
        "pin_sets": pin_checks, "eligibility_sets": pin_checks, "failures": 0}
    rng = Random(report["seed"])
    pin_checks = 0
    for iteration in range(1200):
        n = 3 if iteration < 1000 else 4
        edges = {(j, i) for i in range(n) for j in range(n) if rng.random() < .35}
        functions = []
        for i in range(n):
            indices = [j for j in range(n) if (j, i) in edges]
            table = [0] + [rng.randrange(2) for _ in range((1 << len(indices)) - 1)]
            functions.append(tuple(table[sum(((s >> j) & 1) << k
                                             for k, j in enumerate(indices))]
                                   for s in range(1 << n)))
        if iteration % 3 == 0:
            edges |= {(j, i) for i in range(n) for j in range(n) if rng.random() < .1}
        pin_checks += check_localization(n, functions, edges)
    report["seeded_three_four_variable_maps"] = {
        "three_variable_maps": 1000, "four_variable_maps": 200,
        "pin_sets": pin_checks, "eligibility_sets": pin_checks, "failures": 0}
    count = good = 0
    for flat in product(range(2), repeat=4):
        parents = [flat[:2], flat[2:]]
        for values in product(range(8), repeat=2):
            tables = [(0, *(a >> i & 1 for i in range(3))) for a in values]
            good += check_slot_realization(2, parents, tables)
            count += 1
    report["exhaustive_two_variable_two_slot_models"] = {
        "realizations": count, "G_prime_passes": good, "failures": 0,
        "checks": ["conservative residual cost equivalence", "G-prime implies exact marked-cycle cost"]}
    count = good = 0
    for _ in range(1000):
        n = 3
        parents = [(rng.randrange(n), rng.randrange(n)) for _ in range(n)]
        tables = [(0, rng.randrange(2), rng.randrange(2), rng.randrange(2)) for _ in range(n)]
        good += check_slot_realization(n, parents, tables)
        count += 1
    report["seeded_three_variable_two_slot_models"] = {
        "realizations": count, "G_prime_passes": good, "failures": 0}
    output = Path(__file__).resolve().parents[1] / "results" / "reduction_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
