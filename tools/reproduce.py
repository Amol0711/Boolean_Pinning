#!/usr/bin/env python3
"""Run declared workflows in an isolated copy after checking file integrity.

The manifest checks distributed file sizes and SHA256 digests. The CLI has no
option to disable verification or select a different repository root.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from typing import Any

WORKFLOWS = frozenset({"scc", "biological", "bernoulli", "joint_rules", "joint_replication"})


class RunnerError(ValueError):
    """A declared workflow or protected path is invalid."""


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _no_symlinks(path: Path) -> None:
    absolute = Path(os.path.abspath(path))
    for item in [*reversed(absolute.parents), absolute]:
        if item.is_symlink():
            raise RunnerError(f"Symlinks are not permitted in this path: {item}")


def _relative_path(value: str, *, allow_dot: bool = False) -> PurePosixPath:
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise RunnerError(f"Invalid relative path: {value!r}")
    if allow_dot and value == ".":
        return PurePosixPath(".")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in ("", ".", "..") for part in value.split("/")):
        raise RunnerError(f"Unsafe relative path: {value!r}")
    if str(path) != value or re.match(r"^[A-Za-z]:", value):
        raise RunnerError(f"Noncanonical relative path: {value!r}")
    return path


def _inside(root: Path, relative: str, *, allow_dot: bool = False) -> Path:
    path = root.joinpath(*_relative_path(relative, allow_dot=allow_dot).parts)
    _no_symlinks(path)
    if not path.resolve().is_relative_to(root.resolve()):
        raise RunnerError(f"Path escapes its root: {relative!r}")
    return path


def _load_manifest(root: Path) -> tuple[dict[str, dict[str, Any]], str]:
    root = Path(root).absolute()
    _no_symlinks(root)
    manifest_path = _inside(root, "MANIFEST.json")
    if not manifest_path.is_file():
        raise RunnerError("MANIFEST.json is missing or is not a regular file")
    raw = manifest_path.read_bytes()
    try:
        manifest = json.loads(raw)
    except (ValueError, UnicodeError) as exc:
        raise RunnerError(f"Cannot parse MANIFEST.json: {exc}") from exc
    if not isinstance(manifest, dict) or not isinstance(manifest.get("files"), dict):
        raise RunnerError("MANIFEST.json must contain a files mapping")
    files = manifest["files"]
    if not files:
        raise RunnerError("MANIFEST.json has no files")
    for relative, record in files.items():
        _relative_path(relative)
        if relative == "MANIFEST.json":
            raise RunnerError("The manifest must not list itself")
        if not isinstance(record, dict):
            raise RunnerError(f"Invalid manifest entry: {relative}")
        size = record.get("bytes")
        digest = record.get("sha256")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise RunnerError(f"Invalid byte count: {relative}")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise RunnerError(f"Invalid SHA256: {relative}")
    return files, hashlib.sha256(raw).hexdigest()


def verify_repository(root: Path, *, expected_manifest_sha256: str | None = None) -> dict[str, Any]:
    """Check every listed regular file; do not execute any repository code."""
    report: dict[str, Any] = {
        "ok": False, "manifest_sha256": None, "checked_files": 0,
        "verified_bytes": 0, "errors": [],
        "scope": "File sizes and SHA256 digests against the distributed manifest.",
    }
    try:
        files, manifest_hash = _load_manifest(Path(root))
        report["manifest_sha256"] = manifest_hash
        if expected_manifest_sha256 is not None and manifest_hash != expected_manifest_sha256:
            report["errors"].append("MANIFEST.json changed since the baseline check")
        for relative, record in sorted(files.items()):
            try:
                path = _inside(Path(root), relative)
                if not path.is_file():
                    raise RunnerError("missing or not a regular file")
                actual_size = path.stat().st_size
                actual_hash = _sha256(path)
                report["checked_files"] += 1
                if actual_size != record["bytes"]:
                    raise RunnerError(f"byte count differs ({actual_size} != {record['bytes']})")
                if actual_hash != record["sha256"]:
                    raise RunnerError("SHA256 differs")
                report["verified_bytes"] += actual_size
            except (OSError, RunnerError) as exc:
                report["errors"].append(f"{relative}: {exc}")
    except (OSError, RunnerError) as exc:
        report["errors"].append(str(exc))
    report["ok"] = not report["errors"]
    return report


def load_tasks(root: Path) -> dict[str, dict[str, Any]]:
    """Read and validate the declarative registry; no imports or evaluation."""
    registry_path = _inside(Path(root), "config/tasks.json")
    try:
        registry = json.loads(registry_path.read_bytes())
    except (OSError, ValueError, UnicodeError) as exc:
        raise RunnerError(f"Cannot read config/tasks.json: {exc}") from exc
    if not isinstance(registry, dict) or registry.get("schema") != 1:
        raise RunnerError("Task registry schema must be 1")
    tasks = registry.get("tasks")
    if not isinstance(tasks, dict):
        raise RunnerError("Task registry must contain a tasks mapping")
    for name, task in tasks.items():
        if not isinstance(name, str) or not name or not isinstance(task, dict):
            raise RunnerError("Invalid task name or record")
        if task.get("workflow") not in WORKFLOWS:
            raise RunnerError(f"Unsupported workflow for {name}")
        if not isinstance(task.get("description"), str) or not isinstance(task.get("tier"), str):
            raise RunnerError(f"Task {name} needs string description and tier")
        requirements = task.get("requirements", [])
        if not isinstance(requirements, list) or any(not isinstance(x, str) or not x for x in requirements):
            raise RunnerError(f"Invalid requirements for {name}")
        commands = task.get("commands")
        if not isinstance(commands, list) or not commands:
            raise RunnerError(f"Task {name} needs at least one command")
        for command in commands:
            if not isinstance(command, dict):
                raise RunnerError(f"Invalid command for {name}")
            argv = command.get("argv")
            if not isinstance(argv, list) or not argv or any(not isinstance(x, str) or "\x00" in x for x in argv):
                raise RunnerError(f"Invalid argv for {name}")
            if not argv[0]:
                raise RunnerError(f"Empty executable for {name}")
            for value in argv:
                if any(token in value and value != token for token in ("{python}", "{workers}")):
                    raise RunnerError(f"Placeholders must occupy whole argv tokens: {value!r}")
            _relative_path(command.get("cwd", "."), allow_dot=True)
            timeout = command.get("timeout_seconds")
            if timeout is not None and (
                not isinstance(timeout, (int, float)) or isinstance(timeout, bool)
                or not math.isfinite(timeout) or timeout <= 0
            ):
                raise RunnerError(f"Invalid timeout for {name}")
    return tasks


def _binary_info(name: str) -> dict[str, Any]:
    executable = shutil.which(name)
    result: dict[str, Any] = {"available": executable is not None, "path": executable, "version": None}
    if executable:
        try:
            completed = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=10, check=False)
            lines = (completed.stdout or completed.stderr).splitlines()
            result["version"] = lines[0] if lines else None
            result["version_exit_code"] = completed.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            result["version_error"] = str(exc)
    return result


def _package_info(name: str) -> dict[str, Any]:
    try:
        present = importlib.util.find_spec(name) is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        present = False
    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"available": present, "version": version}


def assertion_policy() -> dict[str, Any]:
    """Keep assertion-based validators active in this process and its children."""
    value = os.environ.get("PYTHONOPTIMIZE")
    # Only unambiguous zero-valued settings are accepted. Do not edit user env.
    environment_safe = value is None or value == "" or value == "0"
    enabled = sys.flags.optimize == 0 and environment_safe
    return {
        "enabled": enabled, "current_optimization_level": sys.flags.optimize,
        "PYTHONOPTIMIZE": value,
        "required": "Run without -O/-OO; PYTHONOPTIMIZE must be unset, empty, or 0.",
    }


def python_policy() -> dict[str, Any]:
    """Require the Python version used by the preserved scientific workflows."""
    return {"supported": sys.version_info >= (3, 11), "minimum_version": "3.11"}


def environment_report(requirements: list[str] | None = None) -> dict[str, Any]:
    report: dict[str, Any] = {
        "python": {"available": True, "executable": sys.executable, "version": sys.version, **python_policy()},
        "platform": platform.platform(), "machine": platform.machine(),
        "assertions": assertion_policy(),
        "numpy": _package_info("numpy"), "scipy": _package_info("scipy"),
        "g++": _binary_info("g++"), "pdflatex": _binary_info("pdflatex"),
    }
    for requirement in requirements or []:
        if requirement in report:
            continue
        if re.fullmatch(r"[A-Za-z_]\w*", requirement) and requirement not in {"pdftoppm", "latexmk", "xelatex", "lualatex", "make"}:
            report[requirement] = _package_info(requirement)
        else:
            report[requirement] = _binary_info(requirement)
    return report


def _output_destination(root: Path, output: Path) -> Path:
    raw = Path(output).expanduser().absolute()
    _no_symlinks(raw)
    path = raw.resolve()
    root = Path(root).resolve()
    if path.exists():
        raise RunnerError(f"Output already exists: {path}")
    if path.is_relative_to(root):
        runs = root / "runs"
        if path == runs or not path.is_relative_to(runs):
            raise RunnerError("Output inside the repository must be a new directory below runs/")
    return path


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _stop_process(process: subprocess.Popen[bytes]) -> None:
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.wait()
    finally:
        # A child may ignore SIGTERM after its parent has exited. The session
        # belongs to this command; also remove any remaining group members.
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def _execute_command(command: dict[str, Any], workspace: Path, logs: Path, index: int, workers: int) -> dict[str, Any]:
    substitutions = {"{python}": sys.executable, "{workers}": str(workers)}
    argv = [substitutions.get(value, value) for value in command["argv"]]
    cwd = _inside(workspace, command.get("cwd", "."), allow_dot=True)
    if not cwd.is_dir():
        raise RunnerError(f"Command working directory is absent: {cwd}")
    stdout_path, stderr_path = logs / f"{index:02d}.stdout.log", logs / f"{index:02d}.stderr.log"
    result: dict[str, Any] = {
        "argv": argv, "cwd": str(cwd), "started_at": _utc_now(),
        "timeout_seconds": command.get("timeout_seconds"), "timed_out": False,
        "stdout_log": str(stdout_path.relative_to(logs.parent)),
        "stderr_log": str(stderr_path.relative_to(logs.parent)),
    }
    start = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        try:
            process = subprocess.Popen(argv, cwd=cwd, stdout=stdout, stderr=stderr, start_new_session=(os.name == "posix"))
        except OSError as exc:
            stderr.write((str(exc) + "\n").encode())
            result.update(exit_code=127, error=str(exc))
        else:
            try:
                result["exit_code"] = process.wait(timeout=command.get("timeout_seconds"))
            except subprocess.TimeoutExpired:
                result["timed_out"] = True
                _stop_process(process)
                result["exit_code"] = process.returncode
            except KeyboardInterrupt:
                _stop_process(process)
                result.update(exit_code=130, interrupted=True, error="Command interrupted")
                stderr.write(b"\nRunner interrupted the command and stopped its process group.\n")
            except BaseException:
                _stop_process(process)
                raise
    result.update(finished_at=_utc_now(), runtime_seconds=time.monotonic() - start)
    result["status"] = "passed" if result["exit_code"] == 0 and not result["timed_out"] else "failed"
    return result


def run_task(root: Path, name: str, output: Path, *, workers: int = 1) -> dict[str, Any]:
    """Internal API for the CLI and isolated fixture tests; verification is mandatory."""
    if not isinstance(workers, int) or isinstance(workers, bool) or workers < 1:
        raise RunnerError("Workers must be a positive integer")
    if not python_policy()["supported"]:
        raise RunnerError("Python 3.11 or later is required for the preserved workflows")
    if not assertion_policy()["enabled"]:
        raise RunnerError("Assertions must remain enabled: run without -O/-OO and leave PYTHONOPTIMIZE unset, empty, or 0")
    root = Path(root).absolute()
    _no_symlinks(root)
    root = root.resolve()
    destination = _output_destination(root, Path(output))
    baseline_before = verify_repository(root)
    if not baseline_before["ok"]:
        raise RunnerError("Baseline integrity check failed: " + "; ".join(baseline_before["errors"]))
    files, manifest_hash = _load_manifest(root)
    if manifest_hash != baseline_before["manifest_sha256"]:
        raise RunnerError("MANIFEST.json changed while preparing the run")
    if "config/tasks.json" not in files:
        raise RunnerError("config/tasks.json must be covered by MANIFEST.json")
    tasks = load_tasks(root)
    if name not in tasks:
        raise RunnerError(f"Unknown task: {name}")
    task = tasks[name]
    prefix = f"workflows/{task['workflow']}/"
    selected = {path: record for path, record in files.items() if path.startswith(prefix)}
    if not selected:
        raise RunnerError(f"No manifest-covered files for workflow {task['workflow']}")
    destination.mkdir(parents=True, exist_ok=False)
    workspace, logs = destination / "workspace", destination / "logs"
    workspace.mkdir()
    logs.mkdir()
    summary: dict[str, Any] = {
        "schema": 1, "task": name, "workflow": task["workflow"], "tier": task["tier"],
        "result_scope": "Declared command completion and baseline integrity; scientific conclusions remain in workflow reports.",
        "description": task["description"], "status": "running", "exit_code": None,
        "started_at": _utc_now(), "repository": str(root), "output": str(destination),
        "workspace": str(workspace), "workers": workers,
        "baseline_before": baseline_before, "baseline_after": None, "baseline_intact": None,
        "workflow_files": len(selected), "workspace_verified": False, "commands": [], "errors": [],
        "workflow_manifest_sha256": hashlib.sha256(json.dumps(selected, sort_keys=True).encode()).hexdigest(),
    }
    start = time.monotonic()
    summary_path = destination / "run_summary.json"
    _atomic_json(summary_path, summary)
    try:
        environment = environment_report(task.get("requirements", []))
        summary["environment"] = environment
        missing = [item for item in task.get("requirements", []) if not environment[item].get("available", False)]
        if missing:
            raise RunnerError("Missing task requirements: " + ", ".join(missing))
        for relative, record in sorted(selected.items()):
            source = _inside(root, relative)
            target = workspace / relative[len(prefix):]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            if target.stat().st_size != record["bytes"] or _sha256(target) != record["sha256"]:
                raise RunnerError(f"Copied workflow file failed integrity check: {relative}")
        summary["workspace_verified"] = True
        _atomic_json(summary_path, summary)
        for index, command in enumerate(task["commands"], 1):
            result = _execute_command(command, workspace, logs, index, workers)
            summary["commands"].append(result)
            _atomic_json(summary_path, summary)
            if result["status"] != "passed":
                summary["status"] = "failed"
                summary["exit_code"] = 124 if result["timed_out"] else (result["exit_code"] if 0 < result["exit_code"] < 256 else 1)
                summary["errors"].append(f"Command {index} failed; see its stdout/stderr logs")
                break
        else:
            summary.update(status="passed", exit_code=0)
    except BaseException as exc:
        summary.update(status="failed", exit_code=130 if isinstance(exc, KeyboardInterrupt) else 1)
        summary["errors"].append(f"{type(exc).__name__}: {exc}")
        (logs / "runner.stderr.log").write_text(traceback.format_exc(), encoding="utf-8")
    finally:
        summary["baseline_after"] = verify_repository(root, expected_manifest_sha256=manifest_hash)
        summary["baseline_intact"] = summary["baseline_after"]["ok"]
        if not summary["baseline_intact"]:
            summary.update(status="failed", exit_code=1)
            summary["errors"].append("Repository baseline changed during the run")
        summary.update(finished_at=_utc_now(), runtime_seconds=time.monotonic() - start)
        _atomic_json(summary_path, summary)
    return summary


def _positive_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("list", help="List declared tasks, workflows and tiers")
    doctor = subparsers.add_parser("doctor", help="Report environment availability; optionally check a task")
    doctor.add_argument("task", nargs="?")
    subparsers.add_parser("verify", help="Check all manifest-listed files without running a workflow")
    run = subparsers.add_parser("run", help="Run a verified workflow in a new isolated directory")
    run.add_argument("task")
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--workers", type=_positive_integer, default=1)
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    try:
        if args.action == "verify":
            report = verify_repository(root)
            print(json.dumps(report, indent=2))
            return 0 if report["ok"] else 1
        if args.action == "list":
            for name, task in sorted(load_tasks(root).items()):
                print(f"{name}\t{task['workflow']}\t{task['tier']}\t{task['description']}")
            return 0
        if args.action == "doctor":
            task = None
            if args.task:
                tasks = load_tasks(root)
                if args.task not in tasks:
                    raise RunnerError(f"Unknown task: {args.task}")
                task = tasks[args.task]
            requirements = task.get("requirements", []) if task else []
            report = environment_report(requirements)
            missing = [item for item in requirements if not report[item].get("available", False)]
            print(json.dumps({"environment": report, "task": args.task, "requirements": requirements, "missing": missing}, indent=2))
            return 1 if missing or not report["assertions"]["enabled"] or not report["python"]["supported"] else 0
        summary = run_task(root, args.task, args.output, workers=args.workers)
        print(json.dumps({"status": summary["status"], "exit_code": summary["exit_code"], "baseline_intact": summary["baseline_intact"], "run_summary": str(Path(summary["output"]) / "run_summary.json")}, indent=2))
        return summary["exit_code"]
    except (OSError, RunnerError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
