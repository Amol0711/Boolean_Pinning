"""Targeted cap/mode checks against a separate physical-state trajectory oracle.

From scc_solver run:
  python verification/targeted_cap_check.py --output results/targeted_cap_check.json

The oracle uses original physical tables and never calls production evaluation,
support, SCC, candidate-search, or convergence helpers. The independent check
covers 720 cap/mode cases; it does not replay the production simulation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scc_control import BooleanNetwork, solve_network


def succeeds(network, pins):
    """Check each compatible physical trajectory directly for a repeated state."""
    pin = set(pins)
    target = sum(bit << i for i, bit in enumerate(network.target))
    pin_mask = sum(1 << i for i in pin)
    for start in range(1 << network.n):
        if (start ^ target) & pin_mask:
            continue
        current, seen = start, set()
        while current != target:
            if current in seen:
                return False
            seen.add(current)
            following = 0
            for i, (parents, table) in enumerate(zip(network.parents, network.tables)):
                pattern = sum(((current >> source) & 1) << slot
                              for slot, source in enumerate(parents))
                bit = network.target[i] if i in pin else table[pattern]
                following |= bit << i
            current = following
    return True


def optimum(network, baseline):
    eligible = [i for i in range(network.n) if i not in baseline]
    for cardinality in range(len(eligible) + 1):
        for pins in itertools.combinations(eligible, cardinality):
            if succeeds(network, tuple(baseline) + pins):
                return cardinality
    raise AssertionError("Full clamping must reach the prescribed target")


def fixtures():
    return [
        ("repeated_xor", BooleanNetwork([[0, 0]], [[0, 1, 1, 0]]), ()),
        ("repeated_or", BooleanNetwork([[0, 0]], [[0, 1, 1, 1]]), ()),
        ("or_complex", BooleanNetwork([[0, 1], [0, 1]],
                                      [[0, 1, 1, 1], [0, 1, 1, 1]]), ()),
        ("pinning_nonmonotone", BooleanNetwork([[0, 1], [0, 1]],
                                               [[0, 1, 0, 0], [0, 1, 0, 0]]), ()),
        ("mixed_target_baseline", BooleanNetwork([[], [2], [2]],
                                                 [[1], [1, 0], [0, 1]],
                                                 [1, 0, 1]), (2,)),
    ]


def run():
    checked = 0
    records = []
    all_reasons = Counter()
    for name, network, baseline in fixtures():
        exact = optimum(network, baseline)
        local_count = 0
        reasons = Counter()
        stopped_after_rejection = 0
        for mode, size, pins, work in itertools.product(
                ["slot", "physical"], [0, 1, None], [0, 1, None],
                [0, 1, 7, 8, 9, 12, 100, None]):
            result = solve_network(network, baseline, mode, size, pins, work)
            context = (name, mode, size, pins, work)
            assert result["lower"] <= exact <= result["upper"], (context, result)
            assert succeeds(network, result["total_feasible_pins"]), (context, result)
            assert result["exact_cost"] is None or result["exact_cost"] == exact, (
                context, result)
            if work is not None:
                assert all(component["work"] <= work for component in result["components"])
            for component in result["components"]:
                reasons[component["reason"]] += 1
                if component["reason"] in ("work_cap_candidate_search", "pin_count_cap"):
                    if component.get("completed_cardinality", -1) >= 0:
                        stopped_after_rejection += 1
            local_count += 1
            checked += 1
        all_reasons.update(reasons)
        records.append({"name": name, "n": network.n, "baseline": list(baseline),
                        "oracle_exact_additional_cost": exact,
                        "solver_calls": local_count,
                        "component_reason_counts": dict(sorted(reasons.items())),
                        "stopped_after_completed_rejected_cardinality": stopped_after_rejection})
    assert checked == 720
    return {
        "status": "PASS", "targeted_cap_mode_cases": checked,
        "support_modes": ["slot", "physical"],
        "max_component_size_values": [0, 1, None],
        "max_pins_values": [0, 1, None],
        "max_work_values": [0, 1, 7, 8, 9, 12, 100, None],
        "oracle": "Direct original-physical-table trajectories from every compatible state",
        "checks": ["lower <= exact optimum <= upper", "every returned full pin witness succeeds",
                   "every exact solver answer equals the oracle optimum",
                   "each component respects its finite charged-work cap"],
        "fixtures": records,
        "component_reason_counts": dict(sorted(all_reasons.items())),
        "scope": "Targeted finite checks; no replay of the original 250000-network experiment",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "results" / "targeted_cap_check.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {result['targeted_cap_mode_cases']} targeted cap/mode calls; {args.output}")


if __name__ == "__main__":
    main()
