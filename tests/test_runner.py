"""Synthetic integration tests for runner protection and isolation.

These fixtures contain no scientific models and perform no study reruns.
Run with: python -m unittest discover -s tests -p 'test_runner.py' -v
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

RUNNER_PATH = Path(__file__).resolve().parents[1] / "tools" / "reproduce.py"
SPEC = importlib.util.spec_from_file_location("repository_reproduction_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)
BUILDER_PATH = RUNNER_PATH.with_name("build_manifest.py")
BUILDER_SPEC = importlib.util.spec_from_file_location("repository_manifest_builder", BUILDER_PATH)
assert BUILDER_SPEC is not None and BUILDER_SPEC.loader is not None
builder = importlib.util.module_from_spec(BUILDER_SPEC)
BUILDER_SPEC.loader.exec_module(builder)


class RunnerIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="reproduction-runner-test-")
        self.addCleanup(self.temporary.cleanup)
        self.parent = Path(self.temporary.name)
        self.root = self.parent / "repository"
        (self.root / "config").mkdir(parents=True)
        (self.root / "workflows" / "scc").mkdir(parents=True)
        self.workflow = self.root / "workflows" / "scc"
        self.output = self.parent / "output"

    def fixture(self, source: str, *, commands=None, requirements=None):
        (self.workflow / "task.py").write_text(source, encoding="utf-8")
        (self.workflow / "input.txt").write_text("baseline\n", encoding="utf-8")
        task = {
            "workflow": "scc", "description": "Synthetic fixture only", "tier": "test",
            "requirements": requirements or [],
            "commands": commands or [{"argv": ["{python}", "task.py"], "cwd": "."}],
        }
        (self.root / "config" / "tasks.json").write_text(json.dumps({"schema": 1, "tasks": {"fixture": task}}), encoding="utf-8")
        self.write_manifest()

    def write_manifest(self):
        files = {}
        for path in self.root.rglob("*"):
            if path.is_file() and path.name != "MANIFEST.json":
                raw = path.read_bytes()
                files[path.relative_to(self.root).as_posix()] = {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
        (self.root / "MANIFEST.json").write_text(json.dumps({"files": files}, sort_keys=True), encoding="utf-8")

    def test_success_is_isolated_and_argv_is_literal(self):
        self.fixture(
            "from pathlib import Path\nimport json, sys\n"
            "Path('input.txt').write_text('workspace changed')\n"
            "Path('result.json').write_text(json.dumps(sys.argv[1:]))\n"
            "print('complete stdout')\nprint('complete stderr', file=sys.stderr)\n",
            commands=[{"argv": ["{python}", "task.py", "{workers}", "$(touch shell_should_not_run)", "semi;colon"], "cwd": "."}],
        )
        (self.workflow / "not_in_manifest.txt").write_text("do not copy")
        result = runner.run_task(self.root, "fixture", self.output, workers=3)
        self.assertEqual(result["exit_code"], 0)
        self.assertEqual(result["status"], "passed")
        self.assertTrue(result["baseline_intact"])
        self.assertTrue(result["workspace_verified"])
        self.assertEqual((self.workflow / "input.txt").read_text(), "baseline\n")
        workspace = self.output / "workspace"
        self.assertEqual((workspace / "input.txt").read_text(), "workspace changed")
        self.assertFalse((workspace / "not_in_manifest.txt").exists())
        self.assertFalse((workspace / "shell_should_not_run").exists())
        self.assertEqual(json.loads((workspace / "result.json").read_text()), ["3", "$(touch shell_should_not_run)", "semi;colon"])
        self.assertEqual(result["commands"][0]["argv"][0], sys.executable)
        self.assertIn("complete stdout", (self.output / "logs/01.stdout.log").read_text())
        self.assertIn("complete stderr", (self.output / "logs/01.stderr.log").read_text())
        self.assertEqual(json.loads((self.output / "run_summary.json").read_text()), result)
        self.assertFalse(list(self.output.glob(".run_summary.json.*.tmp")))

    def test_protected_and_existing_outputs_are_refused(self):
        self.fixture("print('fixture')\n")
        for path in [self.root / "new_results", self.workflow / "new_results", self.root / "runs"]:
            with self.subTest(path=path), self.assertRaises(runner.RunnerError):
                runner.run_task(self.root, "fixture", path)
            self.assertFalse(path.exists())
        self.output.mkdir()
        (self.output / "keep.txt").write_text("untouched")
        with self.assertRaises(runner.RunnerError):
            runner.run_task(self.root, "fixture", self.output)
        self.assertEqual((self.output / "keep.txt").read_text(), "untouched")
        permitted = self.root / "runs" / "new_run"
        result = runner.run_task(self.root, "fixture", permitted)
        self.assertTrue(result["baseline_intact"])
        self.assertEqual(result["exit_code"], 0)

    def test_symlink_output_parent_is_refused(self):
        self.fixture("print('fixture')\n")
        alias = self.parent / "alias"
        try:
            alias.symlink_to(self.workflow, target_is_directory=True)
        except (OSError, NotImplementedError):
            self.skipTest("Symlinks unavailable on this platform")
        with self.assertRaises(runner.RunnerError):
            runner.run_task(self.root, "fixture", alias / "new_result")
        self.assertFalse((self.workflow / "new_result").exists())

    def test_tamper_prevents_execution(self):
        self.fixture("from pathlib import Path\nPath('ran.txt').write_text('bad')\n")
        (self.workflow / "input.txt").write_text("tampered\n")
        verified = runner.verify_repository(self.root)
        self.assertFalse(verified["ok"])
        self.assertTrue(any("input.txt" in item for item in verified["errors"]))
        with self.assertRaises(runner.RunnerError):
            runner.run_task(self.root, "fixture", self.output)
        self.assertFalse(self.output.exists())

    def test_unsafe_manifest_paths_and_symlinks_fail_verification(self):
        self.fixture("print('fixture')\n")
        manifest_path = self.root / "MANIFEST.json"
        original = json.loads(manifest_path.read_text())
        record = {"bytes": 0, "sha256": hashlib.sha256(b"").hexdigest()}
        for unsafe in ["../outside.txt", "/absolute.txt", "workflows/../escape.txt", "C:/escape.txt"]:
            with self.subTest(path=unsafe):
                edited = {"files": {unsafe: record}}
                manifest_path.write_text(json.dumps(edited))
                self.assertFalse(runner.verify_repository(self.root)["ok"])
        manifest_path.write_text(json.dumps(original))
        target = self.workflow / "input.txt"
        target.unlink()
        external = self.parent / "external.txt"
        external.write_text("baseline\n")
        try:
            target.symlink_to(external)
        except (OSError, NotImplementedError):
            self.skipTest("Symlinks unavailable on this platform")
        self.assertFalse(runner.verify_repository(self.root)["ok"])
        with self.assertRaises(runner.RunnerError):
            runner.run_task(self.root, "fixture", self.output)

    def test_task_failure_stops_following_commands_and_keeps_logs(self):
        self.fixture(
            "import sys\nprint('before failure', flush=True)\n"
            "print('reason for failure', file=sys.stderr, flush=True)\nsys.exit(7)\n",
            commands=[{"argv": ["{python}", "task.py"], "cwd": "."}, {"argv": ["{python}", "should_never_run.py"], "cwd": "."}],
        )
        result = runner.run_task(self.root, "fixture", self.output)
        self.assertEqual(result["exit_code"], 7)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(len(result["commands"]), 1)
        self.assertTrue(result["baseline_intact"])
        self.assertIn("before failure", (self.output / "logs/01.stdout.log").read_text())
        self.assertIn("reason for failure", (self.output / "logs/01.stderr.log").read_text())
        self.assertEqual(json.loads((self.output / "run_summary.json").read_text())["exit_code"], 7)

    def test_timeout_is_reported_as_failure_with_partial_output(self):
        self.fixture("import time\nprint('started', flush=True)\ntime.sleep(30)\n", commands=[{"argv": ["{python}", "task.py"], "cwd": ".", "timeout_seconds": 0.15}])
        result = runner.run_task(self.root, "fixture", self.output)
        self.assertEqual(result["exit_code"], 124)
        self.assertEqual(result["status"], "failed")
        self.assertTrue(result["commands"][0]["timed_out"])
        self.assertTrue(result["baseline_intact"])
        self.assertIn("started", (self.output / "logs/01.stdout.log").read_text())

    def test_manifest_change_during_run_fails_postcheck(self):
        self.fixture("from pathlib import Path\np=Path(" + repr(str(self.root / "MANIFEST.json")) + ")\np.write_text(p.read_text()+'\\n')\n")
        result = runner.run_task(self.root, "fixture", self.output)
        self.assertEqual(result["commands"][0]["exit_code"], 0)
        self.assertNotEqual(result["exit_code"], 0)
        self.assertFalse(result["baseline_intact"])
        self.assertTrue(any("MANIFEST.json changed" in item for item in result["baseline_after"]["errors"]))

    def test_registry_must_be_manifest_covered(self):
        self.fixture("print('fixture')\n")
        path = self.root / "MANIFEST.json"
        manifest = json.loads(path.read_text())
        del manifest["files"]["config/tasks.json"]
        path.write_text(json.dumps(manifest))
        with self.assertRaises(runner.RunnerError):
            runner.run_task(self.root, "fixture", self.output)
        self.assertFalse(self.output.exists())

    def test_optimized_interpreters_and_child_environment_are_refused(self):
        self.fixture("from pathlib import Path\nPath('ran.txt').write_text('bad')\n")
        for value in ["1", "2", "enabled"]:
            with self.subTest(environment=value), mock.patch.dict(os.environ, {"PYTHONOPTIMIZE": value}):
                self.assertFalse(runner.assertion_policy()["enabled"])
                with self.assertRaises(runner.RunnerError):
                    runner.run_task(self.root, "fixture", self.output)
                self.assertEqual(os.environ["PYTHONOPTIMIZE"], value)
                self.assertFalse(self.output.exists())
        (self.root / "tools").mkdir()
        copied_runner = self.root / "tools/reproduce.py"
        copied_runner.write_bytes(RUNNER_PATH.read_bytes())
        self.write_manifest()
        child_env = os.environ.copy()
        child_env.pop("PYTHONOPTIMIZE", None)
        completed = subprocess.run(
            [sys.executable, "-O", str(copied_runner), "run", "fixture", "--output", str(self.output)],
            capture_output=True, text=True, env=child_env, check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Assertions must remain enabled", completed.stderr)
        self.assertFalse(self.output.exists())
        doctor = subprocess.run(
            [sys.executable, "-O", str(copied_runner), "doctor", "fixture"],
            capture_output=True, text=True, env=child_env, check=False,
        )
        self.assertNotEqual(doctor.returncode, 0)
        self.assertFalse(json.loads(doctor.stdout)["environment"]["assertions"]["enabled"])

    def test_missing_requirement_and_invalid_workers_never_run(self):
        self.fixture("print('must not run')\n", requirements=["definitely_missing_runner_fixture_module_20261017"])
        for workers in [0, -1, True]:
            with self.subTest(workers=workers), self.assertRaises(runner.RunnerError):
                runner.run_task(self.root, "fixture", self.output, workers=workers)
        self.assertFalse(self.output.exists())
        result = runner.run_task(self.root, "fixture", self.output)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["commands"], [])
        self.assertFalse(result["workspace_verified"])
        self.assertTrue(result["baseline_intact"])
        self.assertIn("Missing task requirements", (self.output / "logs/runner.stderr.log").read_text())

    def test_python_version_guard_prevents_execution_and_is_reported(self):
        self.fixture("from pathlib import Path\nPath('ran.txt').write_text('bad')\n")
        with mock.patch.object(runner.sys, "version_info", (3, 10, 99)):
            self.assertFalse(runner.python_policy()["supported"])
            self.assertFalse(runner.environment_report()["python"]["supported"])
            with self.assertRaisesRegex(runner.RunnerError, "Python 3.11 or later"):
                runner.run_task(self.root, "fixture", self.output)
            self.assertFalse(self.output.exists())
        with mock.patch.object(runner.sys, "version_info", (3, 11, 0)):
            self.assertTrue(runner.python_policy()["supported"])


class ManifestBuilderTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="manifest-builder-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "repository"
        (self.root / "provenance").mkdir(parents=True)
        (self.root / "workflows/scc").mkdir(parents=True)
        (self.root / "docs").mkdir()
        self.original = self.root / "workflows/scc/model.txt"
        self.document = self.root / "docs/preserved.txt"
        self.original.write_text("original model\n")
        self.document.write_text("original qualification\n")
        files = {}
        for path in (self.original, self.document):
            raw = path.read_bytes()
            files[path.relative_to(self.root).as_posix()] = {
                "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
            }
        for inventory in builder.PROVENANCE:
            (self.root / inventory).write_text(json.dumps({"schema": 1, "files": files}))
        self.manifest = self.root / "MANIFEST.json"

    def test_preserved_original_and_document_tamper_or_loss_refuses_rebuild(self):
        builder.build_manifest(self.root)
        original_manifest = self.manifest.read_bytes()
        for path in (self.original, self.document):
            original_bytes = path.read_bytes()
            for operation in ("modify", "remove"):
                with self.subTest(path=path.name, operation=operation):
                    if operation == "modify":
                        path.write_bytes(original_bytes + b"alteration")
                    else:
                        path.unlink()
                    with self.assertRaises(builder.ManifestBuildError):
                        builder.build_manifest(self.root)
                    self.assertEqual(self.manifest.read_bytes(), original_manifest)
                    path.write_bytes(original_bytes)

    def test_check_detects_unexpected_file_and_leaves_manifest_unchanged(self):
        builder.build_manifest(self.root)
        original_manifest = self.manifest.read_bytes()
        extra = self.root / "new_document.txt"
        extra.write_text("new deliverable\n")
        with self.assertRaisesRegex(builder.ManifestBuildError, "unlisted=.*new_document.txt"):
            builder.build_manifest(self.root, check=True)
        self.assertEqual(self.manifest.read_bytes(), original_manifest)
        builder.build_manifest(self.root)
        self.assertTrue(builder.build_manifest(self.root, check=True)["ok"])
        self.assertIn("new_document.txt", json.loads(self.manifest.read_bytes())["files"])
        extra.unlink()
        with self.assertRaisesRegex(builder.ManifestBuildError, "absent=.*new_document.txt"):
            builder.build_manifest(self.root, check=True)

    def test_deterministic_inventory_excludes_generated_files_and_retains_nested_manifest(self):
        nested = self.root / "workflows/scc/MANIFEST.json"
        nested.write_text('{"files": {}}\n')
        (self.root / "shared.so.notes.txt").write_text("documentation, not a shared object\n")
        builder.build_manifest(self.root)
        original_manifest = self.manifest.read_bytes()
        exclusions = [".git/index", ".venv/installed.py", "runs/output/data.txt", "tools/__pycache__/cached.pyc",
                      "orphan.pyc", "workflows/scc/kernel.so", "workflows/scc/kernel.so.1", "native.pyd", "native.dll", "native.dylib"]
        for relative in exclusions:
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("generated\n")
        self.assertTrue(builder.build_manifest(self.root, check=True)["ok"])
        builder.build_manifest(self.root)
        self.assertEqual(self.manifest.read_bytes(), original_manifest)
        data = json.loads(original_manifest)
        self.assertIn("workflows/scc/MANIFEST.json", data["files"])
        self.assertIn("shared.so.notes.txt", data["files"])
        self.assertEqual(list(data["files"]), sorted(data["files"]))
        self.assertEqual(set(data), {"schema", "scope", "files"})
        self.assertFalse(list(self.root.glob(".MANIFEST.json.*.tmp")))

    def test_included_symlink_and_unsafe_provenance_are_rejected(self):
        link = self.root / "docs/alias.txt"
        try:
            link.symlink_to(self.document)
        except (OSError, NotImplementedError):
            self.skipTest("Symlinks unavailable on this platform")
        with self.assertRaisesRegex(builder.ManifestBuildError, "Symlinks"):
            builder.build_manifest(self.root)
        self.assertFalse(self.manifest.exists())
        link.unlink()
        inventory = self.root / builder.PROVENANCE[0]
        original = json.loads(inventory.read_text())
        record = next(iter(original["files"].values()))
        original["files"] = {"../outside.txt": record}
        inventory.write_text(json.dumps(original))
        with self.assertRaisesRegex(builder.ManifestBuildError, "Unsafe relative path"):
            builder.build_manifest(self.root)
        self.assertFalse(self.manifest.exists())


if __name__ == "__main__":
    unittest.main()
