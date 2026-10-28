"""Reproduce the complete simulation study in a new empty destination.

The delivered data are never overwritten.  Use --prepare-only to inspect the
staged sources, or --validate-only for small-network tests without the main run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


SOURCE_NAMES = (
    "engine.py", "engine_core.cpp", "build_core.py", "validate_engine.py",
    "validate_refined_bounds.py", "run_experiments.py", "refine_control.py",
    "refine_results.py", "analyze_results.py", "build_figures.py",
    "pilot_timing.py", "reproduce.py",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(destination, script, *arguments):
    command = [sys.executable, str(destination/"simulation"/script)]
    command.extend(map(str, arguments))
    print("Running " + script, flush=True)
    subprocess.run(command, cwd=destination, check=True)


def verify_raw(destination, original, refined):
    checks = []
    for cell in original["cells"]:
        name = cell["file"]
        actual = digest(destination/name)
        checks.append({"file": name, "expected_sha256": cell["sha256"],
                       "actual_sha256": actual, "match": actual == cell["sha256"]})
    for cell in refined["files"]:
        name = "results/refined/raw/" + cell["file"]
        actual = digest(destination/name)
        checks.append({"file": name, "expected_sha256": cell["refined_sha256"],
                       "actual_sha256": actual,
                       "match": actual == cell["refined_sha256"]})
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, required=True,
                        help="New directory, or existing empty directory; never overwritten")
    parser.add_argument("--workers", type=int, default=4,
                        help="Main-run worker processes (default 4)")
    parser.add_argument("--figures", action="store_true",
                        help="Also rebuild the three TikZ/PGFPlots figures")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only", action="store_true",
                      help="Only copy portable sources and reference manifests")
    mode.add_argument("--validate-only", action="store_true",
                      help="Stage and run small-network validation; no main-size draws")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    source = Path(__file__).resolve().parents[1]
    destination = args.destination.resolve()
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        parser.error("destination must be a new or empty directory")
    # Locate reference manifests before creating the destination.  A staged
    # copy remains usable as a source for another new reproduction destination.
    original_path = source/"results/run_manifest.json"
    refined_path = source/"results/refined/refinement_manifest.json"
    if not original_path.exists():
        original_path = source/"reference/run_manifest.json"
        refined_path = source/"reference/refinement_manifest.json"
    original = json.loads(original_path.read_text())
    refined = json.loads(refined_path.read_text())
    if original["status"] != "completed" or refined["status"] != "completed":
        raise RuntimeError("Reference manifests must describe complete runs")
    if len(original["cells"]) != 25 or len(refined["files"]) != 25:
        raise RuntimeError("Expected 25 original and 25 refined reference cells")
    source_checks = {}
    for name in ("engine.py", "engine_core.cpp", "run_experiments.py"):
        actual = digest(source/"simulation"/name)
        if actual != original["sources_sha256"][name]:
            raise RuntimeError("Frozen source does not match run manifest: " + name)
        source_checks[name] = actual
    for name in ("refine_control.py", "refine_results.py"):
        actual = digest(source/"simulation"/name)
        if actual != refined["source_sha256"][name]:
            raise RuntimeError("Refinement source does not match manifest: " + name)
        source_checks[name] = actual
    if digest(source/"Protocol_Baseline.json") != original["protocol_sha256"]:
        raise RuntimeError("Protocol hash does not match original run")
    if digest(source/"Bound_Refinement_Protocol.json") != refined["analysis_addendum_sha256"]:
        raise RuntimeError("Analysis addendum hash does not match refinement")
    destination.mkdir(parents=True, exist_ok=True)
    (destination/"simulation").mkdir()
    (destination/"results").mkdir()
    (destination/"reference").mkdir()
    for name in SOURCE_NAMES:
        shutil.copy2(source/"simulation"/name, destination/"simulation"/name)
    for name in ("Protocol_Baseline.json", "Bound_Refinement_Protocol.json", "README.txt", "requirements.txt"):
        shutil.copy2(source/name, destination/name)
    shutil.copy2(original_path, destination/"reference/run_manifest.json")
    shutil.copy2(refined_path, destination/"reference/refinement_manifest.json")
    report = {"status": "prepared", "source_verification": source_checks,
              "main_draws_requested": 0 if (args.prepare_only or args.validate_only) else 250000,
              "worker_processes": args.workers, "figures_requested": args.figures,
              "new_destination": str(destination)}
    report_path = destination/"reproduction_report.json"
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    if args.prepare_only:
        print(json.dumps(report, indent=2))
        return
    try:
        run(destination, "validate_engine.py")
        run(destination, "validate_refined_bounds.py")
        run(destination, "analyze_results.py", "--self-test")
        if args.validate_only:
            report["status"] = "small-network validation passed"
        else:
            run(destination, "run_experiments.py", "--root", destination,
                "--workers", args.workers)
            run(destination, "analyze_results.py", "--root", destination)
            run(destination, "refine_results.py")
            run(destination, "analyze_results.py", "--root", destination,
                "--raw-dir", destination/"results/refined/raw",
                "--out-dir", destination/"results/refined")
            checks = verify_raw(destination, original, refined)
            report["raw_sha256_checks"] = checks
            report["raw_files_matching"] = sum(x["match"] for x in checks)
            if not all(x["match"] for x in checks):
                raise RuntimeError("Raw-file SHA256 mismatch; inspect reproduction_report.json")
            if args.figures:
                run(destination, "build_figures.py", "--analysis-dir", destination/"results/refined")
            report["status"] = "complete; all 50 raw files match delivered manifests"
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        report_path.write_text(json.dumps(report, indent=2)+"\n")
        raise
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
