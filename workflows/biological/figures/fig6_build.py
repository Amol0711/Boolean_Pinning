#!/usr/bin/env python3
"""Export/check the exact biological figure data and compile its editable vector PDF.

Run from any directory. The complete computational tree is needed to regenerate
data. The clean LaTeX archive can rebuild the figure using --compile-only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

FIGURES = Path(__file__).resolve().parent
ROOT = FIGURES.parent
NAME = "fig6_cd4_target_mechanisms"
CONTEXTS = ["CD4_OFF", "CD4_ON", "CD4_APC", "CD4_ECM", "CD4_IL2", "CD4_GALPHAQ"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


def export_and_check():
    source = ROOT / "results/all_targets.json"
    production = json.loads(source.read_text())
    contexts = {c["id"]: c for c in production["contexts"]}
    independent_path = ROOT / "results/independent_biological_validation.json"
    independent = json.loads(independent_path.read_text())
    assert independent["status"] == "PASS"
    independent_contexts = {c["id"]: c for c in independent["expected"]["contexts"]}
    data = FIGURES / "data"
    summaries = []
    labels = []
    for j, context_id in enumerate(CONTEXTS):
        c = contexts[context_id]
        counts = {str(k): 0 for k in range(2, 6)}
        for target in c["fixed_points"]:
            cost = target["exact_cost"]
            assert cost in range(2, 6)
            assert cost == sum(component["exact_cost"] for component in target["components"])
            counts[str(cost)] += 1
        assert {k: v for k, v in counts.items() if v} == c["cost_histogram"]
        assert sum(counts.values()) == c["fixed_point_count"]
        assert independent_contexts[context_id]["cost_histogram"] == c["cost_histogram"]
        assert independent_contexts[context_id]["fixed_point_count"] == c["fixed_point_count"]
        summaries.append({"context": context_id, "row": j,
                          **{"cost" + k: v for k, v in counts.items()},
                          "total": c["fixed_point_count"],
                          "target_costs": [{"target_id": t["target_id"],
                                            "exact_cost": t["exact_cost"]}
                                           for t in c["fixed_points"]]})
        start = 0
        for k, count in counts.items():
            if count:
                color = "white" if int(k) >= 4 else "black!85"
                labels.append(r"\node[text=%s] at (axis cs:%s,%s) {%s};" %
                              (color, start + count / 2, j, count))
            start += count
        labels.append(r"\node[anchor=west,text=black!75,xshift=2pt] at (axis cs:%s,%s) {%s};" %
                      (start, j, start))
    with (data / "fig6_cd4_cost_distribution.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["context", "row", "cost2", "cost3", "cost4", "cost5", "total"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)
    (data / "fig6_bar_labels.tex").write_text("\n".join(labels) + "\n")

    selected = [t for t in contexts["CD4_APC"]["fixed_points"] if t["selected_zero_orbit"]]
    assert len(selected) == 1
    target = selected[0]
    assert target["target_id"] == contexts["CD4_APC"]["selected_target_id"]
    assert target["residual_size"] == 17 and target["exact_cost"] == 4
    independent_target = next(t for t in independent_contexts["CD4_APC"]["fixed_points"]
                              if t["selected_zero_orbit"])
    assert independent_target["target_by_name"] == target["target_by_name"]
    cyclic = {tuple(c["names"]): c for c in target["components"] if c["cyclic"]}
    expected_costs = {
        ("v_Dec2", "v_GATA3"): 1,
        ("v_FOXP3", "v_SMAD3", "v_TGFB", "v_TGFBR"): 1,
        ("v_RORGT",): 0, ("v_TBET",): 1, ("v_proliferation",): 1,
    }
    assert set(cyclic) == set(expected_costs)
    assert {tuple(c["names"]): c["cost"] for c in independent_target["components"] if c["cyclic"]} == expected_costs
    acyclic = [c for c in target["components"] if not c["cyclic"]]
    assert len(acyclic) == 8
    assert all(len(c["names"]) == 1 and c["exact_cost"] == 0 for c in acyclic)

    # These exact truth tables correspond to the displayed symbolic rules.
    expected_rows = {
        "v_Dec2": (["v_GATA3"], [0, 1]),
        "v_GATA3": (["v_Dec2"], [0, 1]),
        "v_FOXP3": (["v_FOXP3", "v_SMAD3"], [0, 1, 1, 1]),
        "v_SMAD3": (["v_TGFBR"], [0, 1]),
        "v_TGFB": ([], [0]),
        "v_TGFBR": (["v_TGFB"], [0, 1]),
        "v_RORGT": ([], [0]),
        "v_TBET": (["v_TBET"], [0, 1]),
        "v_proliferation": (["v_proliferation"], [0, 1]),
    }
    index_names = {i: n for i, n in zip(target["residual"], target["residual_names"])}
    modules = []
    for names, c in cyclic.items():
        assert c["exact_cost"] == expected_costs[names]
        module_rows = []
        for row in c["restricted_rows"]:
            name = index_names[row["node"]]
            parents = [index_names[i] for i in row["parents"]]
            assert (parents, row["table"]) == expected_rows[name]
            module_rows.append({"name": name, "parents": parents, "table": row["table"]})
        modules.append({"names": list(names), "cost": c["exact_cost"],
                        "isolated_target_error_rows": module_rows})

    # Independent source grammar evaluation, no production model evaluator.
    checker_path = ROOT / "larger_networks/verification/independent_core.py"
    spec = importlib.util.spec_from_file_location("fig6_independent_core", checker_path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    bnet_path = ROOT / "larger_networks/data/cd4_188/model.bnet"
    source_rows, _, _ = checker.read_bnet(bnet_path)
    source_rows = {row[0]: row for row in source_rows}
    truth_checks = 0
    for module in modules:
        names = module["names"]
        for state in range(1 << len(names)):
            values = dict(target["target_by_name"])
            for j, name in enumerate(names):
                values[name] ^= (state >> j) & 1
            for row in module["isolated_target_error_rows"]:
                name = row["name"]
                _, scope, truth = source_rows[name]
                actual = truth[tuple(values[parent] for parent in scope)] ^ target["target_by_name"][name]
                corner = sum((values[parent] ^ target["target_by_name"][parent]) << j
                             for j, parent in enumerate(row["parents"]))
                assert actual == row["table"][corner]
                truth_checks += 1
    assert all(target["target_by_name"][name] == 0 for name in expected_rows)
    record = {
        "schema": "CD4TargetMechanismsFigure_v1",
        "source_json": "results/all_targets.json",
        "source_json_sha256": digest(source),
        "independent_validation_sha256": digest(independent_path),
        "independent_validation_status": "PASS",
        "model_bnet_sha256": digest(bnet_path),
        "independent_parser_sha256": digest(checker_path),
        "context_order": CONTEXTS,
        "panel_a": summaries,
        "panel_b": {"target_id": target["target_id"],
                    "target_by_name": target["target_by_name"],
                    "residual_names": target["residual_names"],
                    "cyclic_components": modules,
                    "acyclic_singleton_names": [c["names"][0] for c in acyclic],
                    "exact_additional_cost": 4,
                    "independently_checked_source_rule_values": truth_checks},
        "interpretation": {
            "scc_partition": "Functional dependencies on the baseline residual coordinate cube.",
            "displayed_rules": "Component isolation fixes all coordinates outside that component at the selected target.",
            "cost": "Additional permanent target-valued pins; environmental input clamps are excluded.",
            "admissible_interventions": "All regulated coordinates, including the modeled proliferation output.",
            "scope": "Exact synchronous Boolean-model statements; no intervention-feasibility or clinical claim.",
        },
    }
    write_json(data / "fig6_source_records.json", record)
    with (data / "fig6_apc_isolated_rules.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["component", "component_cost", "node", "parents_little_endian", "table"])
        for j, module in enumerate(modules):
            for row in module["isolated_target_error_rows"]:
                writer.writerow([j, module["cost"], row["name"], ";".join(row["parents"]), "".join(map(str, row["table"]))])
    return record


def compile_figure():
    env = dict(os.environ, SOURCE_DATE_EPOCH="1792022400", FORCE_SOURCE_DATE="1")
    with tempfile.TemporaryDirectory(prefix="cd4_target_figure_") as directory:
        temporary = Path(directory)
        shutil.copy2(FIGURES / (NAME + ".tex"), temporary)
        shutil.copytree(FIGURES / "data", temporary / "data")
        for _ in range(2):
            compiled = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", NAME + ".tex"],
                                      cwd=temporary, env=env, text=True, capture_output=True)
            if compiled.returncode:
                raise RuntimeError(compiled.stdout[-4000:] + compiled.stderr)
        log = (temporary / (NAME + ".log")).read_text()
        assert "Overfull" not in log and "Undefined" not in log
        shutil.copy2(temporary / (NAME + ".pdf"), FIGURES / (NAME + ".pdf"))
        qa = ROOT / "qa/fig6"
        qa.mkdir(parents=True, exist_ok=True)
        (qa / "build.log").write_text(log)
        subprocess.run(["pdftoppm", "-png", "-singlefile", "-r", "180", str(FIGURES / (NAME + ".pdf")), str(qa / "fig6_preview")], check=True)
    metadata = {
        "source_date_epoch": 1792022400,
        "latex_passes": 2,
        "hashes": {str(p.relative_to(FIGURES)): digest(p) for p in
                   [FIGURES / (NAME + ".tex"), FIGURES / (NAME + ".pdf"),
                    *sorted(p for p in (FIGURES / "data").glob("fig6_*")
                            if p.name != "fig6_build_metadata.json")]},
    }
    write_json(FIGURES / "data/fig6_build_metadata.json", metadata)
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compile-only", action="store_true")
    args = parser.parse_args()
    if not args.compile_only:
        record = export_and_check()
        print("Exported 46 exact CD4 target costs and verified", record["panel_b"]["independently_checked_source_rule_values"], "isolated rule values against the BNET source.")
    metadata = compile_figure()
    print("Compiled", NAME + ".pdf", metadata["hashes"][NAME + ".pdf"])
