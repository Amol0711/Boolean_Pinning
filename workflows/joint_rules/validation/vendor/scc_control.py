"""Deterministic SCC localization for permanent target-value Boolean pins.

No external dependencies.  JSON CLI: python scc_control.py input.json -o result.json
All returned costs exclude baseline clamps.  Exact-search caps are per component.
"""
from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
import itertools
import json
from numbers import Integral
from pathlib import Path
from typing import Sequence


def _bit(value, name):
    if not isinstance(value, Integral) or value not in (0, 1):
        raise ValueError(f"{name} must be an integer bit")
    return int(value)


@dataclass(frozen=True, init=False)
class BooleanNetwork:
    """Ordered input slots and truth tables; slot j is index bit j.

    Repeated slots refer to the same physical variable during dynamics.  The
    target must be a fixed point of the supplied uncontrolled physical map.
    """

    parents: tuple[tuple[int, ...], ...]
    tables: tuple[tuple[int, ...], ...]
    target: tuple[int, ...]
    g_tables: tuple[tuple[int, ...], ...]
    n: int

    def __init__(self, parents: Sequence[Sequence[int]],
                 tables: Sequence[Sequence[int]], target=None):
        p = tuple(tuple(row) for row in parents)
        n = len(p)
        if len(tables) != n:
            raise ValueError("one truth table is required per coordinate")
        for row in p:
            for u in row:
                if not isinstance(u, Integral) or isinstance(u, bool) or not 0 <= u < n:
                    raise ValueError("parent indices must be integers in [0,n)")
        p = tuple(tuple(int(u) for u in row) for row in p)
        t = tuple(tuple(_bit(v, "truth-table entry") for v in row) for row in tables)
        for row, tab in zip(p, t):
            if len(tab) != 1 << len(row):
                raise ValueError("a k-slot truth table must contain exactly 2**k entries")
        x = (0,) * n if target is None else tuple(_bit(v, "target entry") for v in target)
        if len(x) != n:
            raise ValueError("target length differs from network size")
        g = []
        for i, (row, tab) in enumerate(zip(p, t)):
            mask = sum(x[u] << j for j, u in enumerate(row))
            if tab[mask] != x[i]:
                raise ValueError(f"target is not a fixed point at coordinate {i}")
            g.append(tuple(tab[e ^ mask] ^ x[i] for e in range(len(tab))))
        for name, value in (("parents", p), ("tables", t), ("target", x),
                            ("g_tables", tuple(g)), ("n", n)):
            object.__setattr__(self, name, value)


def from_two_input_arrays(parents, rules_nonzero_corners, target=None):
    """Adapt two-input simulation (n,2) parents and (n,3) error-corner rule arrays.

    The three entries mean g(1,0), g(0,1), g(1,1), not physical F corners.
    For a nonzero target the supplied error rules are conjugated back to F.
    This function copies arrays and never modifies inputs or engine records.
    """
    parents = tuple(tuple(row) for row in parents)
    rules = tuple(tuple(row) for row in rules_nonzero_corners)
    if len(parents) != len(rules) or any(len(row) != 2 for row in parents):
        raise ValueError("adapter requires (n,2) parents and (n,3) corner rules")
    if any(len(row) != 3 for row in rules):
        raise ValueError("adapter requires three nonzero error-corner entries")
    if target is None:
        return BooleanNetwork(parents, [(0,) + row for row in rules])
    # Validate physical fixed point and indices first using zero error rules.
    x = tuple(_bit(v, "target entry") for v in target)
    if len(x) != len(parents):
        raise ValueError("target length differs from network size")
    validation = BooleanNetwork(parents, [(x[i],) * 4 for i in range(len(parents))], x)
    tables = []
    for i, row in enumerate(validation.parents):
        g = (0,) + tuple(_bit(v, "error-corner entry") for v in rules[i])
        mask = sum(x[u] << j for j, u in enumerate(row))
        tables.append(tuple(g[state ^ mask] ^ x[i] for state in range(4)))
    return BooleanNetwork(parents, tables, x)


def _indices(values, n, name="indices"):
    values = tuple(values)
    if any(not isinstance(v, Integral) or isinstance(v, bool) or not 0 <= v < n for v in values):
        raise ValueError(f"{name} must be integers in [0,n)")
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must not contain duplicates")
    return tuple(sorted(int(v) for v in values))


def evaluate_state(network: BooleanNetwork, state: int, pins=()) -> int:
    """One physical synchronous update; clamp input and output of pinned rows.

    State bit i is coordinate i.  Pins always have the physical target value.
    Compatible input states are sufficient for global-convergence tests.
    """
    if not isinstance(state, Integral) or not 0 <= state < (1 << network.n):
        raise ValueError("state is outside the network state space")
    state = int(state)
    pinned = set(_indices(pins, network.n, "pins"))
    for i in pinned:
        state = (state & ~(1 << i)) | (network.target[i] << i)
    result = 0
    for i, (row, tab) in enumerate(zip(network.parents, network.tables)):
        val = network.target[i] if i in pinned else tab[
            sum(((state >> u) & 1) << j for j, u in enumerate(row))]
        result |= val << i
    return result


def _subset_or(table):
    result = list(table)
    bit = 1
    while bit < len(result):
        for mask in range(len(result)):
            if mask & bit:
                result[mask] |= result[mask ^ bit]
        bit <<= 1
    return tuple(result)


def support_residual(network: BooleanNetwork, baseline=(), mode="slot"):
    """Greatest support using a deletion worklist.

    ``slot`` allows nominal slots to vary independently, hence is conservative
    with repeated parents. ``physical`` first identifies repeated variables.
    Evaluations count worklist row-support queries, not truth-table preprocessing.
    """
    if mode not in ("slot", "physical"):
        raise ValueError("support mode must be 'slot' or 'physical'")
    baseline = _indices(baseline, network.n, "baseline")
    pin = set(baseline)
    variables, closures = [], []
    for row, g in zip(network.parents, network.g_tables):
        if mode == "slot":
            variables.append(row)
            closures.append(_subset_or(g))
        else:
            unique = tuple(dict.fromkeys(row))
            pos = {u: j for j, u in enumerate(unique)}
            table = []
            for state in range(1 << len(unique)):
                idx = sum(((state >> pos[u]) & 1) << j for j, u in enumerate(row))
                table.append(g[idx])
            variables.append(unique)
            closures.append(_subset_or(table))
    alive = [i not in pin for i in range(network.n)]
    dependents = [set() for _ in range(network.n)]
    for i, row in enumerate(variables):
        for u in row:
            dependents[u].add(i)
    queue = deque(i for i in range(network.n) if alive[i])
    queued = set(queue)
    evaluations = 0
    while queue:
        i = queue.popleft()
        queued.remove(i)
        if not alive[i]:
            continue
        mask = sum(int(alive[u]) << j for j, u in enumerate(variables[i]))
        evaluations += 1
        if closures[i][mask]:
            continue
        alive[i] = False
        for v in sorted(dependents[i]):
            if alive[v] and v not in queued:
                queue.append(v)
                queued.add(v)
    residual = [i for i in range(network.n) if alive[i]]
    return {"residual": residual, "mode": mode, "evaluations": evaluations,
            "deletion_count": network.n - len(baseline) - len(residual)}


def _sccs(network, residual):
    """Kosaraju on all nominal residual dependencies; iterative to avoid stack caps."""
    rset = set(residual)
    adj = {v: set() for v in residual}
    rev = {v: set() for v in residual}
    for v in residual:
        for u in network.parents[v]:
            if u in rset:
                adj[u].add(v)
                rev[v].add(u)
    adj = {v: tuple(sorted(row)) for v, row in adj.items()}
    rev = {v: tuple(sorted(row)) for v, row in rev.items()}
    visited, order = set(), []
    for root in residual:
        if root in visited:
            continue
        visited.add(root)
        stack = [(root, 0)]
        while stack:
            v, index = stack[-1]
            if index == len(adj[v]):
                stack.pop()
                order.append(v)
            else:
                u = adj[v][index]
                stack[-1] = (v, index + 1)
                if u not in visited:
                    visited.add(u)
                    stack.append((u, 0))
    assigned, components = set(), []
    for root in reversed(order):
        if root in assigned:
            continue
        component, stack = [], [root]
        assigned.add(root)
        while stack:
            v = stack.pop()
            component.append(v)
            for u in rev[v]:
                if u not in assigned:
                    assigned.add(u)
                    stack.append(u)
        components.append(sorted(component))
    return sorted(components, key=lambda c: c[0])


def _cycle_info(network, nodes):
    cset = set(nodes)
    edges = [(u, v, j) for v in nodes for j, u in enumerate(network.parents[v]) if u in cset]
    cyclic = len(nodes) > 1 or any(u == v for u, v, j in edges)
    good = cyclic and len(edges) == len(nodes)
    if good:
        incoming_counts = {v: 0 for v in nodes}
        for u, v, j in edges:
            incoming_counts[v] += 1
            if network.g_tables[v][1 << j] != 1:
                good = False
                break
        good = good and all(incoming_counts[v] == 1 for v in nodes)
    return {"cyclic": cyclic, "marked_slot_cycle": good,
            "internal_slots": [{"source": u, "receiver": v, "slot": j} for u, v, j in edges]}


class _BudgetExceeded(Exception):
    pass


def _cap(value, name):
    if value is not None and (not isinstance(value, Integral) or isinstance(value, bool) or value < 0):
        raise ValueError(f"{name} must be a nonnegative integer or null")
    return None if value is None else int(value)


def solve_component(network: BooleanNetwork, nodes, max_component_size=18,
                    max_pins=8, max_work=2_000_000):
    """Exact isolated subsystem search with exterior errors held at zero.

    Searches each candidate directly. Work units are one local row evaluation
    in successor construction or one controlled successor-edge traversal.
    All component coordinates are eligible; pinning them all is always the
    fallback feasible witness, including when construction exhausts its budget.
    """
    nodes = _indices(nodes, network.n, "component nodes")
    max_component_size = _cap(max_component_size, "max_component_size")
    max_pins = _cap(max_pins, "max_pins")
    max_work = _cap(max_work, "max_work")
    s = len(nodes)
    result = {"nodes": list(nodes), "exact_cost": None, "lower": 0, "upper": s,
              "pins": list(nodes), "reason": None, "work": 0,
              "row_evaluations": 0, "successor_traversals": 0,
              "completed_cardinality": -1, "next_cardinality": 0,
              "candidate_counts": [], "rejected_candidates": 0,
              "failed_candidate_last": None, "last_failure_state": None,
              "limits": {"max_component_size": max_component_size,
                         "max_pins": max_pins, "max_work": max_work}}

    def finish(reason, exact=None):
        result["reason"] = reason
        if exact is not None:
            result.update(exact_cost=exact, lower=exact, upper=exact, next_cardinality=None)
        return result

    if not s:
        return finish("empty_component", 0)
    if max_component_size is not None and s > max_component_size:
        return finish("component_size_cap")

    def charge(kind):
        if max_work is not None and result["work"] >= max_work:
            raise _BudgetExceeded
        result["work"] += 1
        result[kind] += 1

    pos = {v: k for k, v in enumerate(nodes)}
    local_slots = [tuple((j, pos[u]) for j, u in enumerate(network.parents[v]) if u in pos)
                   for v in nodes]
    count = 1 << s
    result["successor_construction_work_required"] = count * s
    if max_work is not None and count * s > max_work:
        # Preflight prevents an allocation whose complete table cannot fit
        # the declared work budget. No row work is consumed in this branch.
        return finish("work_cap_successor_construction")
    nxt = [0] * count
    try:
        for state in range(count):
            successor = 0
            for k, v in enumerate(nodes):
                charge("row_evaluations")
                pattern = sum(((state >> bit) & 1) << j for j, bit in local_slots[k])
                successor |= network.g_tables[v][pattern] << k
            nxt[state] = successor
    except _BudgetExceeded:
        return finish("work_cap_successor_construction")
    assert nxt[0] == 0

    def succeeds(mask):
        # Each candidate uses fresh visitation state. No monotonicity of
        # success under adding or deleting pins is assumed.
        good = bytearray(count)
        good[0] = 1
        keep = (count - 1) ^ mask
        for start in range(count):
            if start & mask or good[start]:
                continue
            path, index = [], {}
            state = start
            while not good[state]:
                if state in index:
                    return False, start
                index[state] = len(path)
                path.append(state)
                charge("successor_traversals")
                state = nxt[state] & keep
            for state in path:
                good[state] = 1
        return True, None

    limit = s if max_pins is None else min(s, max_pins)
    for cardinality in range(limit + 1):
        result["next_cardinality"] = cardinality
        stats = {"cardinality": cardinality, "tested": 0, "rejected": 0,
                 "all_rejected": False}
        result["candidate_counts"].append(stats)
        # Full clamping is a symbolic witness; no graph traversal needed.
        if cardinality == s:
            result["pins"] = list(nodes)
            return finish("all_pinned_witness", s)
        for candidate in itertools.combinations(range(s), cardinality):
            mask = sum(1 << k for k in candidate)
            try:
                ok, failure = succeeds(mask)
            except _BudgetExceeded:
                result["lower"] = cardinality
                return finish("work_cap_candidate_search")
            stats["tested"] += 1
            if ok:
                result["pins"] = [nodes[k] for k in candidate]
                return finish("exact_search", cardinality)
            stats["rejected"] += 1
            result["rejected_candidates"] += 1
            result["failed_candidate_last"] = [nodes[k] for k in candidate]
            # Witness is an error-state index in the declared nodes order.
            result["last_failure_state"] = failure
        stats["all_rejected"] = True
        result["completed_cardinality"] = cardinality
        result["lower"] = cardinality + 1
    if result["lower"] == s:
        return finish("exact_from_bounds", s)
    result["next_cardinality"] = limit + 1
    return finish("pin_count_cap")


def solve_network(network: BooleanNetwork, baseline=(), support_mode="slot",
                  max_component_size=18, max_pins=8, max_work=2_000_000):
    """Localize to SCCs and return exact cost or complete certified bounds.

    G and G_prime are distinct nominal-slot certificates. Repeated slots are
    counted separately. G_prime ignores acyclic singleton SCCs and requires
    every cyclic SCC to be a marked single-slot cycle after exterior target
    substitution. All costs and witnesses are ADDITIONAL to baseline clamps.
    """
    baseline = _indices(baseline, network.n, "baseline")
    for value, name in ((max_component_size, "max_component_size"),
                        (max_pins, "max_pins"), (max_work, "max_work")):
        _cap(value, name)
    support = support_residual(network, baseline, support_mode)
    residual = support["residual"]
    rset = set(residual)
    G = all(sum(u in rset for u in network.parents[v]) == 1 for v in residual)
    components = []
    for nodes in _sccs(network, residual):
        info = _cycle_info(network, nodes)
        if not info["cyclic"]:
            record = {"nodes": nodes, "exact_cost": 0, "lower": 0, "upper": 0,
                      "pins": [], "reason": "acyclic_singleton", "work": 0,
                      "completed_cardinality": None, "next_cardinality": None}
        elif info["marked_slot_cycle"]:
            record = {"nodes": nodes, "exact_cost": 1, "lower": 1, "upper": 1,
                      "pins": [min(nodes)], "reason": "marked_slot_cycle", "work": 0,
                      "completed_cardinality": None, "next_cardinality": None}
        else:
            record = solve_component(network, nodes, max_component_size, max_pins, max_work)
        record.update(info)
        components.append(record)
    lower = sum(c["lower"] for c in components)
    upper = sum(c["upper"] for c in components)
    Gp = all(not c["cyclic"] or c["marked_slot_cycle"] for c in components)
    cycles = sum(c["cyclic"] for c in components)
    pins = sorted(v for c in components for v in c["pins"])
    assert len(pins) == upper
    assert not G or Gp
    return {"schema_version": "1.0", "n": network.n,
            "target": list(network.target), "baseline": list(baseline),
            "support_mode": support_mode, "support": support,
            "residual": residual, "residual_size": len(residual),
            "G": G, "G_prime": Gp, "cyclic_component_count": cycles,
            "certificate_semantics": "nominal_slot_residual",
            "original_two_slot_unpinned_scope": (
                support_mode == "slot" and not baseline
                and all(len(row) == 2 for row in network.parents)),
            "certificate_cycle_count": cycles if Gp else None,
            "exact_cost": lower if lower == upper else None,
            "lower": lower, "upper": upper, "feasible_pins": pins,
            "total_feasible_pins": sorted(set(baseline) | set(pins)),
            "unresolved_components": sum(c["exact_cost"] is None for c in components),
            "total_work": sum(c["work"] for c in components),
            "work_budget_scope": "per_complex_component",
            "components": components}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args(argv)
    try:
        data = json.loads(args.input.read_text())
        allowed = {"parents", "tables", "target", "baseline", "options"}
        if not isinstance(data, dict) or set(data) - allowed:
            raise ValueError("input must use only parents, tables, target, baseline, options")
        network = BooleanNetwork(data["parents"], data["tables"], data.get("target"))
        result = solve_network(network, baseline=data.get("baseline", ()), **data.get("options", {}))
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(str(exc))
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
