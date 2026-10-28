"""Independent, finite physical-state validation of the Baseline simulation engine.

The reference oracle evaluates the original physical Boolean network and all
permanent pin sets. It does not call the engine's support or transition code.
These checks compare implementations on complete finite state spaces.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

import engine


def literal_support(parents, rules):
    """Synchronous set iteration with the three terms written literally."""
    n = len(parents)
    current = set(range(n))
    history = [current.copy()]
    while True:
        following = set()
        for i in range(n):
            a, b = (int(j) for j in parents[i])
            b10, b01, b11 = (int(j) for j in rules[i])
            if ((b10 and a in current) or (b01 and b in current)
                    or (b11 and a in current and b in current)):
                following.add(i)
        assert following <= current
        history.append(following.copy())
        if following == current:
            return current, history
        current = following


def physical_successor(parents, rules, state, pinmask=0):
    """Physical inputs sharing a source always receive the identical bit."""
    output = 0
    for i in range(len(parents)):
        if (pinmask >> i) & 1:
            continue
        first = (state >> int(parents[i][0])) & 1
        second = (state >> int(parents[i][1])) & 1
        if first == 0 and second == 0:
            value = 0
        elif first == 1 and second == 0:
            value = int(rules[i][0])
        elif first == 0 and second == 1:
            value = int(rules[i][1])
        else:
            value = int(rules[i][2])
        output |= value << i
    return output


def all_pinset_oracle(parents, rules):
    """All physical states and all subsets of ALL vertices, including outside R."""
    n = len(parents)
    full = (1 << n) - 1
    unpinned_next = [physical_successor(parents, rules, s) for s in range(1 << n)]
    successful = set()
    for pinmask in range(1 << n):
        good = True
        zero_basin = {0}
        for initial in range(1 << n):
            if initial & pinmask:
                continue
            state, visited = initial, set()
            while state not in zero_basin and state not in visited:
                visited.add(state)
                state = unpinned_next[state] & (full ^ pinmask)
            if state not in zero_basin:
                good = False
                break
            zero_basin.update(visited)
        if good:
            successful.add(pinmask)
    assert full in successful
    return {
        "successful": successful,
        "cost": min(k.bit_count() for k in successful),
        "fixed_points": [s for s, nxt in enumerate(unpinned_next) if s == nxt],
        "unpinned_next": unpinned_next,
    }


def mask_of(vertices):
    mask = 0
    for vertex in vertices:
        mask |= 1 << int(vertex)
    return mask


def validate_instance(parents, rules, *, label):
    parents = np.asarray(parents, dtype=np.int64)
    rules = np.asarray(rules, dtype=np.uint8)
    n = len(parents)
    degrees = np.bincount(parents.ravel(), minlength=n).astype(np.int64)
    reference_R, history = literal_support(parents, rules)
    ref_mask = mask_of(reference_R)
    expected_double = sum(sum(int(p) in reference_R for p in parents[i]) == 2
                          for i in reference_R)
    expected_good = all(sum(int(p) in reference_R for p in parents[i]) == 1
                        for i in reference_R)
    oracle = all_pinset_oracle(parents, rules)
    result = engine.certify(parents, rules, degrees)
    python_result = engine.certify(parents, rules, degrees, backend="python")
    for key in result:
        assert np.array_equal(result[key], python_result[key]), (label, key, "backend disagreement")
    assert set(np.flatnonzero(result["R"])) == reference_R, label
    assert int(result["residual_size"]) == len(reference_R), label
    assert int(result["double_rows"]) == expected_double, label
    assert bool(result["certified"]) == expected_good, label
    assert ref_mask in oracle["successful"], (label, "pin-all-R bound")

    # Direct check of uniform capture for every pin set and initial state.
    full = (1 << n) - 1
    support_masks = [mask_of(s) for s in history]
    for pinmask in range(1 << n):
        for initial in range(1 << n):
            if initial & pinmask:
                continue
            state = initial
            for support_mask in support_masks:
                assert state & ~support_mask == 0, (label, "uniform capture")
                assert state & pinmask == 0, (label, "persistent pin")
                state = oracle["unpinned_next"][state] & (full ^ pinmask)

    if expected_good:
        pins = list(map(int, result["pins"]))
        assert len(set(pins)) == len(pins), label
        assert mask_of(pins) in oracle["successful"], (label, "certificate success")
        assert int(result["cost"]) == len(pins) == oracle["cost"], (label, "optimality")
        assert len(oracle["fixed_points"]) == 2 ** oracle["cost"], (label, "fixed points")
        for i in reference_R:
            unique = [int(p) for p in parents[i] if int(p) in reference_R]
            assert len(unique) == 1
            assert int(result["residual_parent"][i]) == unique[0], label
    else:
        assert int(result["cost"]) == -1, (label, "no false certification")

    fallback = engine.exact_residual_cost(parents, rules, result["R"],
                                          max_residual=12, work_cap=2_000_000)
    assert bool(fallback["exact"]), (label, fallback)
    assert int(fallback["cost"]) == oracle["cost"], (label, "fallback minimum")
    fallback_pins = list(map(int, fallback["pins"]))
    assert len(fallback_pins) == oracle["cost"], label
    assert mask_of(fallback_pins) in oracle["successful"], (label, "fallback success")
    return {"certified": expected_good, "cost": oracle["cost"],
            "residual_size": len(reference_R), "double_rows": expected_double,
            "fixed_point_count": len(oracle["fixed_points"])}


def deterministic_examples():
    # Copies use [1,0,1] and ignore the second input. Slot order is explicit.
    return {
        "cooperative_and": ([[1, 2], [0, 2], [0, 1]], [[0, 0, 1]] * 3, 1),
        "pruning_not_minimum": ([[1, 2], [0, 2], [0, 1]],
                                [[1, 1, 1], [1, 0, 1], [1, 0, 1]], 1),
        "nonmonotone_success": ([[0, 1], [0, 1]], [[1, 0, 0], [1, 0, 1]], 0),
        "copy_cycle_with_zero_second_input": ([[2, 3], [0, 3], [1, 3], [0, 1]],
                                               [[1, 0, 1]] * 3 + [[0, 0, 0]], 1),
        "repeated_source_physically_zero": ([[0, 0]], [[1, 0, 0]], 0),
        "two_marked_slots_one_physical_self_loop": ([[0, 0]], [[1, 1, 1]], 1),
        "repeated_input_and": ([[0, 0]], [[0, 0, 1]], 1),
        "empty_support": ([[0, 1], [0, 1]], [[0, 0, 0], [0, 0, 0]], 0),
    }


def validate_generator():
    mixtures = {
        "A": ((1, .5), (3, .1)),
        "B": ((1, .1), (3, .5)),
        "C": ((1, .8), (3, .2)),
        "D": ((1, .35), (3, .35)),
        "E": ((0, .8), (4, .3)),
    }
    rows = []
    counts = Counter()
    for m_index, (mixture, types) in enumerate(mixtures.items()):
        for n in (2, 4, 6):
            for repetition in range(16):
                seed = 970_001 + 10_000*m_index + 100*n + repetition
                generated = engine.generate_network(n, mixture, np.random.default_rng(seed))
                parents, rules, degree, q = generated
                replay = engine.generate_network(n, mixture, np.random.default_rng(seed))
                assert all(np.array_equal(a, b) for a, b in zip(generated, replay))
                assert parents.shape == (n, 2) and rules.shape == (n, 3)
                assert degree.shape == (n,) and q.shape == (n,)
                assert np.all((0 <= parents) & (parents < n))
                assert np.all((rules == 0) | (rules == 1))
                assert np.array_equal(np.bincount(parents.ravel(), minlength=n), degree)
                assert sum(map(int, degree)) == 2*n
                for d_value, q_value in types:
                    assert int(np.sum((degree == d_value) & (q == q_value))) == n//2
                row = validate_instance(parents, rules,
                                        label=f"mixture={mixture},n={n},rep={repetition}")
                counts["random_networks"] += 1
                counts["random_certified"] += row["certified"]
                counts["random_uncertified"] += not row["certified"]
                rows.append({"mixture": mixture, "n": n, "rep": repetition,
                             "seed": seed, **row})
    return counts, rows


class PrescribedDraws:
    """Expose each labeled-stub permutation and each corner draw explicitly."""
    def __init__(self, permutation, uniforms):
        self.order = np.asarray(permutation, dtype=np.int64)
        self.uniforms = np.asarray(uniforms, dtype=float)
        self.calls = []

    def permutation(self, values):
        self.calls.append("permutation")
        assert len(values) == len(self.order)
        return np.asarray(values)[self.order]

    def random(self, shape):
        self.calls.append("random")
        assert tuple(shape) == self.uniforms.shape
        return self.uniforms.copy()


def validate_generator_construction():
    """Exact projection multiplicities, without a fragile random p-value test."""
    projection_results = []
    for label, mixture, degrees, probabilities, expected_patterns, multiplicity in (
        ("degrees_1_3", "A", [1, 3], [.5, .1], 4, 6),
        ("degrees_0_4", "E", [0, 4], [.8, .3], 1, 24),
        ("degrees_2_2", ((2, .3, 1.0),), [2, 2], [.3, .3], 6, 4),
    ):
        projected = Counter()
        expected_bits = np.array([[1, 0, 1], [0, 1, 0]], dtype=np.uint8)
        q = np.array(probabilities)[:, None]
        uniforms = np.where(expected_bits, q/2, (1+q)/2)
        for permutation in itertools.permutations(range(4)):
            rng = PrescribedDraws(permutation, uniforms)
            parents, bits, actual_degrees, actual_q = engine.generate_network(2, mixture, rng)
            assert rng.calls == ["permutation", "random"]
            assert np.array_equal(bits, expected_bits)
            assert np.array_equal(actual_degrees, degrees)
            assert np.array_equal(actual_q, probabilities)
            projected[tuple(map(int, parents.ravel()))] += 1
        assert len(projected) == expected_patterns
        assert set(projected.values()) == {multiplicity}
        projection_results.append({"case": label, "labeled_matchings": 24,
                                   "physical_parent_patterns": expected_patterns,
                                   "preimages_per_pattern": multiplicity})

    # Every one of the 64 two-node corner patterns must survive the independent
    # thresholding of the full (n,3) uniform draw, in the declared column order.
    for bits in itertools.product((0, 1), repeat=6):
        expected = np.array(bits, dtype=np.uint8).reshape(2, 3)
        q = np.array([.5, .1])[:, None]
        rng = PrescribedDraws((0, 1, 2, 3), np.where(expected, q/2, (1+q)/2))
        generated = engine.generate_network(2, "A", rng)
        assert np.array_equal(generated[1], expected)
        assert rng.calls == ["permutation", "random"]

    # Actual production SeedSequence streams are replayable independently of
    # replicate processing order; this is separate from default_rng seed checks.
    first = {}
    for rep in (0, 1, 2, 3):
        first[rep] = engine.generate_network(6, "C", engine.make_rng(8301, engine.MIXTURE_IDS["C"], 6, rep))
    for rep in (3, 1, 0, 2):
        replay = engine.generate_network(6, "C", engine.make_rng(8301, engine.MIXTURE_IDS["C"], 6, rep))
        assert all(np.array_equal(a, b) for a, b in zip(first[rep], replay))
    return {"matching_projection": projection_results,
            "corner_patterns_checked": 64, "production_seed_replays": 4,
            "argument": "Each physical source pattern has product_i(d_i!) labeled-stub preimages, so the source-array permutation is exactly the projection of a uniform labeled matching. A separate (n,3) uniform array is thresholded coordinatewise by receiver q, without corner sharing or conditioning on the matching."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1]/"results"/"engine_validation.json")
    args = parser.parse_args()
    totals = Counter()
    boundary_results = {}
    for name, (parents, rules, expected) in deterministic_examples().items():
        result = validate_instance(parents, rules, label=name)
        assert result["cost"] == expected, name
        boundary_results[name] = result
        totals["deterministic_examples"] += 1

    # Exhaust ALL two-node physical wirings, including repeated parents, and all
    # planted truth tables. This newly checks the new queue/certificate/fallback
    # implementation against an independent full-network physical oracle.
    for sources in itertools.product(range(2), repeat=4):
        parents = np.array(sources, dtype=np.int64).reshape(2, 2)
        for bits in itertools.product((0, 1), repeat=6):
            rules = np.array(bits, dtype=np.uint8).reshape(2, 3)
            result = validate_instance(parents, rules, label=f"n2:{sources}:{bits}")
            totals["exhaustive_n2_networks"] += 1
            totals["exhaustive_n2_certified"] += result["certified"]
            totals["exhaustive_n2_uncertified"] += not result["certified"]

    random_counts, random_rows = validate_generator()
    totals.update(random_counts)
    generator_construction = validate_generator_construction()

    # Deterministic exhaustion limits must return unresolved, never a guessed
    # optimum. Empty residuals are allowed to resolve with zero work.
    parents = np.array([[1, 2], [0, 2], [0, 1]], dtype=np.int64)
    rules = np.array([[0, 0, 1]]*3, dtype=np.uint8)
    R = np.ones(3, dtype=bool)
    capped = engine.exact_residual_cost(parents, rules, R, max_residual=12, work_cap=0)
    assert not capped["exact"] and int(capped["cost"]) == -1
    assert int(capped["work"]) <= 0
    oversized = engine.exact_residual_cost(parents, rules, R, max_residual=2,
                                           work_cap=2_000_000)
    assert not oversized["exact"] and int(oversized["cost"]) == -1
    uncapped = engine.exact_residual_cost(parents, rules, R, max_residual=12,
                                          work_cap=2_000_000)
    needed_work = int(uncapped["work"])
    assert needed_work > 0
    one_short = engine.exact_residual_cost(parents, rules, R, max_residual=12,
                                           work_cap=needed_work - 1)
    exactly_enough = engine.exact_residual_cost(parents, rules, R, max_residual=12,
                                                work_cap=needed_work)
    assert not one_short["exact"] and int(one_short["cost"]) == -1
    assert int(one_short["work"]) <= needed_work - 1
    assert exactly_enough["exact"] and int(exactly_enough["cost"]) == 1
    assert int(exactly_enough["work"]) == needed_work
    totals["work_cap_checks"] = 4

    report = {
        "status": "passed",
        "scope": "Independent finite implementation validation, not an asymptotic experiment",
        "independence": "Literal synchronous set support; physical truth-table oracle; all pin subsets of all physical vertices and all compatible states. No engine transitions used in reference calculations.",
        "counts": dict(totals),
        "backend": engine.backend_info(),
        "checks": ["prescribed output degrees and type probabilities",
                   "seed replay", "support fixed set", "slot-count certificate event",
                   "all-pin/all-state uniform capture", "certificate physical success and global minimum",
                   "physical fixed-point count on certified event",
                   "exact residual fallback agrees with full physical optimum",
                   "repeated-source semantics", "self-loops", "zero output degrees",
                   "empty residual", "deterministic work/residual caps"],
        "generator_law": generator_construction,
        "boundary_examples": boundary_results,
        "random_checks": random_rows,
        "source_sha256": {name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                          for name, path in {"engine.py": engine.__file__,
                                             "engine_core.cpp": Path(engine.__file__).with_name("engine_core.cpp"),
                                             "validate_engine.py": __file__}.items()},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": "passed", "counts": dict(totals), "output": str(args.output)}))


if __name__ == "__main__":
    main()
