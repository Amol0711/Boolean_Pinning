"""Independent finite validation of the separate post-run control bounds.

The base engine validator and its results are verified without modification.
This script runs no main-size experiment and reads no main-run outcome table.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

import engine
from refine_control import control_bounds
from validate_engine import (all_pinset_oracle, deterministic_examples,
                             literal_support, mask_of)


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR_SHA256 = "31f3779053ddda1a76016f14d3b2de3b679b8ce7eb17445e0d6ddb3e00ff3d05"
VALIDATION_CONTENT_SHA256 = "3ea7e56de8fcfaaca9e26f3297f8efdffc06b099e30d654293ef7480625ef3f2"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_inputs():
    assert digest(ROOT/"simulation/validate_engine.py") == VALIDATOR_SHA256, "Base validator changed"
    report = json.loads((ROOT/"results/engine_validation.json").read_text())
    content = json.dumps(report, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(content).hexdigest() == VALIDATION_CONTENT_SHA256, "Base validation content changed"


def independent_reference(parents, rules):
    residual, _ = literal_support(parents, rules)
    doubles = {i for i in residual
               if int(parents[i][0]) in residual and int(parents[i][1]) in residual}
    pinned_rules = rules.copy()
    for i in doubles:
        pinned_rules[i] = [0, 0, 0]
    refined, _ = literal_support(parents, pinned_rules)
    assert refined <= residual
    assert refined.isdisjoint(doubles)
    parent_map = {}
    for i in refined:
        surviving = [int(p) for p in parents[i] if int(p) in refined]
        assert len(surviving) == 1
        parent_map[i] = surviving[0]
    cycles = set()
    for start in refined:
        position, path, vertex = {}, [], start
        while vertex not in position:
            position[vertex] = len(path)
            path.append(vertex)
            vertex = parent_map[vertex]
        cycles.add(frozenset(path[position[vertex]:]))
    return residual, doubles, refined, cycles


def validate_one(parents, rules, label):
    p = np.asarray(parents, dtype=np.int64)
    b = np.asarray(rules, dtype=np.uint8)
    degree = np.bincount(p.ravel(), minlength=len(p)).astype(np.int64)
    before = tuple(a.copy() for a in (p, b, degree))
    result = control_bounds(p, b, degree)
    assert all(np.array_equal(a, saved) for a, saved in zip((p, b, degree), before)), (label, "input mutation")
    residual, doubles, refined, cycles = independent_reference(p, b)
    oracle = all_pinset_oracle(p, b)
    c = len(cycles)
    lower, upper = int(result["lower"]), int(result["upper"])
    pins = list(map(int, result["pins"]))
    assert lower == c and upper == c+len(doubles), label
    assert lower <= oracle["cost"] <= upper <= len(residual), (label, "invalid interval")
    assert bool(result["original_certified"]) == (not doubles), label
    assert bool(result["exact"]) == (lower == upper), label
    assert int(result["cycle_cost"]) == c, label
    assert int(result["double_rows"]) == len(doubles), label
    assert set(map(int, result["double_vertices"])) == doubles, label
    assert int(result["residual_size"]) == len(residual), label
    assert int(result["refined_residual_size"]) == len(refined), label
    assert len(pins) == len(set(pins)) == upper, label
    assert doubles <= set(pins) <= residual, label
    assert mask_of(pins) in oracle["successful"], (label, "upper candidate physically fails")
    for cycle in cycles:
        assert len(cycle.intersection(pins)) == 1, (label, "cycle representative")
    if not doubles:
        assert lower == oracle["cost"] == upper, (label, "certified exactness")
    return {"cost": oracle["cost"], "lower": lower, "upper": upper,
            "residual_size": len(residual), "double_rows": len(doubles),
            "refined_residual_size": len(refined), "certified": not doubles}


def main():
    check_inputs()
    counts = Counter()
    boundaries = {}
    for name, (p, b, expected) in deterministic_examples().items():
        row = validate_one(p, b, name)
        assert row["cost"] == expected
        boundaries[name] = row
        counts["deterministic_examples"] += 1
    for sources in itertools.product(range(2), repeat=4):
        p = np.array(sources, dtype=np.int64).reshape(2, 2)
        for bits in itertools.product((0, 1), repeat=6):
            b = np.array(bits, dtype=np.uint8).reshape(2, 3)
            row = validate_one(p, b, f"exhaustive:{sources}:{bits}")
            counts["exhaustive_n2_networks"] += 1
            counts["exhaustive_n2_certified"] += row["certified"]
            counts["exhaustive_n2_uncertified"] += not row["certified"]
            counts["exhaustive_n2_strictly_narrower_than_zero_to_R"] += (
                row["lower"] > 0 or row["upper"] < row["residual_size"])
    random_rows = []
    for m_index, mixture in enumerate(("A", "B", "C", "D", "E")):
        for n in (2, 4, 6):
            for rep in range(16):
                seed = 970_001 + 10_000*m_index + 100*n + rep
                p, b, _, _ = engine.generate_network(n, mixture, np.random.default_rng(seed))
                row = validate_one(p, b, f"random:{mixture}:{n}:{rep}")
                counts["random_networks"] += 1
                counts["random_certified"] += row["certified"]
                counts["random_uncertified"] += not row["certified"]
                random_rows.append({"mixture": mixture, "n": n, "rep": rep,
                                    "seed": seed, **row})
    check_inputs()
    report = {
        "status": "passed",
        "scope": "Independent finite validation of a separately identified post-run deterministic refinement; no main data regenerated or changed.",
        "checked_property": "For every tested network, the returned interval contains the independently enumerated minimum cost and its upper-endpoint pin set succeeds from every compatible initial state.",
        "independence": "Literal synchronous support, independent parent-map cycle counting, and all physical pin subsets and compatible states of the original network.",
        "counts": dict(counts),
        "checks": ["original D reconstruction", "modified support contained in original R",
                   "D excluded after zeroing its rules", "modified support satisfies G",
                   "independent copy-cycle count", "interval contains original physical optimum",
                   "upper-bound candidate succeeds for every physical initial state",
                   "upper candidate cardinality equals bound", "one selected vertex per retained cycle",
                   "exact on original certified event", "no source-array mutation",
                   "base validator and canonical validation content unchanged"],
        "boundary_examples": boundaries,
        "random_checks": random_rows,
        "checked_inputs": {"validator_sha256": VALIDATOR_SHA256, "canonical_validation_sha256": VALIDATION_CONTENT_SHA256},
        "source_sha256": {name: digest(ROOT/"simulation"/name) for name in (
            "validate_refined_bounds.py", "refine_control.py", "engine.py", "engine_core.cpp")},
    }
    output = ROOT/"results"/"refined_bounds_validation.json"
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": "passed", "counts": dict(counts), "output": str(output)}))


if __name__ == "__main__":
    main()
