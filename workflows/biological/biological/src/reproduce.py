"""Reproduce exact model certificates and validation without third-party packages.

Run `python src/reproduce.py --output /path/to/a/new/directory` from any working
directory. Existing output directories are refused. Original learning scripts
are never imported or executed and the 29-node global state space is not scanned.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tests"))
from boolean_certificates import (all_minimum_copy_sets, cube_row_table, identity_rows,
                                  load_json_model, physical_residual_model,
                                  pure_copy_certificate, support_map, support_trace)
from independent_oracle import classify, reachable_support_series
from validate_engine import validate_all


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def trace_record(model, pins):
    trace = support_trace(model, pins)
    return {"pins_zero_based": sorted(pins), "pins_one_based": [i + 1 for i in sorted(pins)],
            "support_trace_zero_based": [sorted(s) for s in trace],
            "support_sizes": [len(s) for s in trace],
            "residual_zero_based": sorted(trace[-1]),
            "capture_bound_updates": len(trace) - 1}


def residual_record(model, baseline):
    residual = support_trace(model, baseline)[-1]
    zero = model.relabeled()
    rows = []
    for i in sorted(residual):
        active, table = cube_row_table(zero, i, residual)
        rows.append({"row_zero_based": i, "row_one_based": i + 1,
                     "parents_zero_based": list(active), "deviation_truth_table": list(table)})
    return {"coordinate_convention": "y_i=x_i XOR target_i; coordinates outside residual are zero.",
            "rows": rows}


def trajectory(model, initial, pins, steps):
    states = [model.clamp(initial, pins)]
    for _ in range(steps):
        states.append(model.step(states[-1], pins))
    return [list(s) for s in states]


def verify29(path):
    document = json.loads(path.read_text())
    model = load_json_model(document)
    mandatory = identity_rows(model)
    if mandatory != (0, 8):
        raise AssertionError("Unexpected identity-row lower bound")
    certificate = pure_copy_certificate(model, mandatory)
    if certificate is None or certificate["residual"] != [14, 15, 16, 17, 18]:
        raise AssertionError("Unexpected 29-node residual")
    if certificate["cycles"] != [[14, 15]]:
        raise AssertionError("Unexpected residual cycle")
    minimum_sets = [tuple(sorted(set(mandatory) | set(k)))
                    for k in all_minimum_copy_sets(certificate)]
    if minimum_sets != [(0, 8, 14), (0, 8, 15)]:
        raise AssertionError("Unexpected optimal pin sets")
    competitors = list(model.target)
    for i in certificate["residual"]:
        competitors[i] ^= 1
    assert tuple(competitors) != model.target
    assert model.step(competitors, mandatory) == tuple(competitors)
    controlled_records = []
    for pins in minimum_sets:
        record = trace_record(model, pins)
        if record["residual_zero_based"] or record["capture_bound_updates"] != 4:
            raise AssertionError("Expected a four-update global support bound")
        initial = list(model.target)
        initial[9] ^= 1
        path_states = trajectory(model, initial, pins, 4)
        if path_states[-1] != list(model.target) or path_states[-2] == list(model.target):
            raise AssertionError("Four-update lower-bound witness failed")
        record["four_update_witness_physical_trajectory"] = path_states
        record["all_state_maximum_transient"] = 4
        controlled_records.append(record)
    # Independent local/global checks on the 5-node physical residual (32 states),
    # never the original 2^29 state space.
    reduced, order = physical_residual_model(model, certificate["residual"])
    reduced_unpinned = classify(reduced.rules, reduced.target)
    assert not reduced_unpinned["globally_convergent"]
    reduced_pinned = []
    for node in (14, 15):
        check = classify(reduced.rules, reduced.target, {order.index(node)})
        assert check["globally_convergent"]
        reduced_pinned.append({"original_pin_zero_based": node, "oracle": check})
    return {"status": "PASS", "model_id": document["model_id"],
            "model_json_sha256": digest(path), "target_is_fixed": True,
            "claim_scope": "Exact optimum and global convergence for the explicit synchronous coded Boolean implementation; no claim of biological reconstruction fidelity or novelty of these control sets.",
            "proof_method": "Mandatory physical identity rows, exact support-cube descent, exact pure-copy residual, one pin per residual cycle. Full global state enumeration was not used.",
            "mandatory_identity_nodes_zero_based": list(mandatory),
            "mandatory_identity_nodes_one_based": [i + 1 for i in mandatory],
            "mandatory_lower_bound_explanation": "An unpinned identity coordinate initially different from its target never changes, regardless of any pins on other coordinates.",
            "baseline": trace_record(model, mandatory),
            "copy_certificate": certificate,
            "residual_map": residual_record(model, mandatory),
            "baseline_competing_fixed_point_physical": competitors,
            "baseline_competing_fixed_point_checked": True,
            "additional_optimum_given_free_baseline": 1,
            "total_optimum_counting_all_pins": 3,
            "all_minimum_total_sets_zero_based": [list(s) for s in minimum_sets],
            "all_minimum_total_sets_one_based": [[i + 1 for i in s] for s in minimum_sets],
            "controlled_support_certificates": controlled_records,
            "residual_oracle_unpinned": reduced_unpinned,
            "residual_oracle_one_pin": reduced_pinned}


def verify17(path):
    document = json.loads(path.read_text())
    model = load_json_model(document)
    mandatory = identity_rows(model)
    assert mandatory == (0, 1, 2)
    baseline = trace_record(model, mandatory)
    assert baseline["residual_zero_based"] == [3, 4, 9, 10, 12, 15, 16]
    assert pure_copy_certificate(model, mandatory) is None
    reduced, order = physical_residual_model(model, baseline["residual_zero_based"])
    reduced_oracle = classify(reduced.rules, reduced.target)
    physical_oracle = classify(model.rules, model.target, mandatory)
    assert reduced_oracle["globally_convergent"] and physical_oracle["globally_convergent"]
    assert reduced_oracle["compatible_states_enumerated"] == 128
    assert physical_oracle["compatible_states_enumerated"] == 16384
    assert reduced_oracle["all_state_maximum_transient"] == 8
    assert physical_oracle["all_state_maximum_transient"] == 11
    reachable = reachable_support_series(model.rules, model.target, mandatory, 11)
    assert reachable[-1]["distinct_reachable_states"] == 1
    assert reachable[-1]["exact_reachable_union_support_size"] == 0
    assert reachable[-2]["exact_reachable_union_support_size"] > 0
    for states, pinned_model, pins in ((physical_oracle["longest_convergent_trajectory"], model, mandatory),
                                       (reduced_oracle["longest_convergent_trajectory"], reduced, ())):
        assert all(pinned_model.step(a, pins) == tuple(b) for a, b in zip(states, states[1:]))
    return {"status": "PASS", "model_id": document["model_id"],
            "model_json_sha256": digest(path), "target_is_fixed": True,
            "claim_scope": "Exact optimum for the coded synchronous Boolean map before automatic freezing, at its supplied fixed-point target. This is not an assertion that the coded map exactly reproduces a published biological model.",
            "proof_method": "Three mandatory identity coordinates give the cost lower bound. Exact physical-state enumeration proves that pinning just those three coordinates globally converges. Exact residual enumeration independently agrees.",
            "mandatory_identity_nodes_zero_based": list(mandatory),
            "mandatory_identity_nodes_one_based": [i + 1 for i in mandatory],
            "baseline": baseline,
            "pure_copy_certificate_available": False,
            "residual_map": residual_record(model, mandatory),
            "additional_optimum_given_free_baseline": 0,
            "total_optimum_counting_all_pins": 3,
            "unique_minimum_total_set_zero_based": list(mandatory),
            "unique_minimum_total_set_one_based": [i + 1 for i in mandatory],
            "residual_all_state_oracle": reduced_oracle,
            "original_all_compatible_state_oracle": physical_oracle,
            "exact_reachable_set_images": reachable,
            "support_capture_plus_residual_transient_upper_bound": baseline["capture_bound_updates"] + reduced_oracle["all_state_maximum_transient"],
            "interpretation": "A nonempty invariant cube support need not imply recurrent nonzero dynamics. Here the seven-node nonlinear residual is globally transient, so adding a pin merely to force empty cube support would overestimate the optimum."}


def write_json(path, document):
    path.write_text(json.dumps(document, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="New output directory; existing paths are refused")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    if output.exists():
        parser.error("Refusing to overwrite an existing output directory")
    output.mkdir(parents=True, exist_ok=False)
    validation = validate_all(ROOT / "models")
    result29 = verify29(ROOT / "models" / "tlgl_coded_29.json")
    result17 = verify17(ROOT / "models" / "drosophila_coded_17.json")
    write_json(output / "validation.json", validation)
    write_json(output / "tlgl_coded_29_certificate.json", result29)
    write_json(output / "drosophila_coded_17_certificate.json", result17)
    with (output / "convergence_curves.csv").open("w", newline="") as stream:
        fields = ["update", "cube_support_size", "exact_reachable_union_support_size",
                  "distinct_reachable_states"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        sizes = result17["baseline"]["support_sizes"]
        for record in result17["exact_reachable_set_images"]:
            writer.writerow({"update": record["update"],
                             "cube_support_size": sizes[min(record["update"], len(sizes) - 1)],
                             "exact_reachable_union_support_size": record["exact_reachable_union_support_size"],
                             "distinct_reachable_states": record["distinct_reachable_states"]})
    source_files = sorted(p for directory in (ROOT / "src", ROOT / "tests", ROOT / "models")
                          for p in directory.iterdir() if p.is_file() and p.suffix in (".py", ".json", ".bnet"))
    write_json(output / "run_manifest.json", {
        "schema": "ExactBooleanCertificateRun_v1",
        "random_data_generated": False,
        "engine_tests_seed": 2026101401,
        "dependencies": "Python 3.10+ standard library only",
        "source_sha256": {str(p.relative_to(ROOT)): digest(p) for p in source_files},
        "result_sha256": {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()},
    })
    print("PASS: exact 29-node and 17-node certificates and independent engine checks reproduced.")


if __name__ == "__main__":
    main()
