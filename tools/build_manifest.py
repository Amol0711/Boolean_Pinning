#!/usr/bin/env python3
"""Maintain the repository integrity manifest against computational file records.

Run after repository edits have stopped. The computational inventory is checked
before inventorying deliverable files. --check performs no writes and requires
the existing manifest to describe the exact current deliverable inventory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import tempfile
from typing import Any

PROVENANCE = ("provenance/COMPUTATIONAL_INVENTORY.json",)
SCOPE = "Current distributed file sizes and SHA256 digests, with computational payloads checked against their separate inventory."
EXCLUDED_ROOT_DIRECTORIES = frozenset({".git", ".venv", "runs"})


class ManifestBuildError(ValueError):
    """The preserved payload or current inventory is unsafe or inconsistent."""


def _no_symlinks(path: Path) -> None:
    absolute = Path(os.path.abspath(path))
    for item in [*reversed(absolute.parents), absolute]:
        if item.is_symlink():
            raise ManifestBuildError(f"Symlinks are not permitted in included paths: {item}")


def _safe_path(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
        raise ManifestBuildError(f"Invalid relative path: {relative!r}")
    parsed = PurePosixPath(relative)
    if parsed.is_absolute() or any(part in ("", ".", "..") for part in relative.split("/")):
        raise ManifestBuildError(f"Unsafe relative path: {relative!r}")
    if str(parsed) != relative or re.match(r"^[A-Za-z]:", relative):
        raise ManifestBuildError(f"Noncanonical relative path: {relative!r}")
    path = root.joinpath(*parsed.parts)
    _no_symlinks(path)
    return path


def _file_record(path: Path) -> dict[str, Any]:
    _no_symlinks(path)
    if not stat.S_ISREG(path.stat().st_mode):
        raise ManifestBuildError(f"Not a regular file: {path}")
    digest = hashlib.sha256()
    count = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            count += len(block)
            digest.update(block)
    return {"bytes": count, "sha256": digest.hexdigest()}


def _read_json(path: Path) -> Any:
    _no_symlinks(path)
    if not stat.S_ISREG(path.stat().st_mode):
        raise ManifestBuildError(f"Not a regular file: {path}")
    try:
        return json.loads(path.read_bytes())
    except (ValueError, UnicodeError) as exc:
        raise ManifestBuildError(f"Cannot parse {path}: {exc}") from exc


def verify_preserved_files(root: Path) -> dict[str, dict[str, Any]]:
    """Verify the current computational inventory before any manifest write."""
    expected: dict[str, dict[str, Any]] = {}
    for provenance in PROVENANCE:
        data = _read_json(_safe_path(root, provenance))
        if not isinstance(data, dict) or data.get("schema") != 1 or not isinstance(data.get("files"), dict) or not data["files"]:
            raise ManifestBuildError(f"Invalid or empty provenance inventory: {provenance}")
        for relative, record in sorted(data["files"].items()):
            path = _safe_path(root, relative)
            if not isinstance(record, dict):
                raise ManifestBuildError(f"Invalid provenance record: {relative}")
            size, digest = record.get("bytes"), record.get("sha256")
            if not isinstance(size, int) or isinstance(size, bool) or size < 0:
                raise ManifestBuildError(f"Invalid provenance byte count: {relative}")
            if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
                raise ManifestBuildError(f"Invalid provenance SHA256: {relative}")
            canonical = {"bytes": size, "sha256": digest}
            if relative in expected and expected[relative] != canonical:
                raise ManifestBuildError(f"Conflicting provenance records: {relative}")
            try:
                actual = _file_record(path)
            except OSError as exc:
                raise ManifestBuildError(f"Preserved file missing or unreadable: {relative}: {exc}") from exc
            if actual != canonical:
                raise ManifestBuildError(f"Preserved file differs from {provenance}: {relative}")
            expected[relative] = canonical
    return expected


def _excluded_file(name: str) -> bool:
    # Native build products are rebuilt in isolated workspaces. Versioned .so
    # names and Python/Windows/macOS shared-library suffixes are also excluded.
    return name.endswith((".pyc", ".so", ".pyd", ".dll", ".dylib")) or re.search(r"\.so(?:\.\d+)+$", name) is not None


def collect_manifest(root: Path) -> dict[str, Any]:
    """Read-only inventory API; useful before a maintainer freezes the tree."""
    root = Path(root).absolute()
    _no_symlinks(root)
    if not root.is_dir():
        raise ManifestBuildError(f"Repository is not a directory: {root}")
    expected = verify_preserved_files(root)
    files: dict[str, dict[str, Any]] = {}

    def walk_error(error: OSError) -> None:
        raise error

    for directory, directories, filenames in os.walk(root, topdown=True, followlinks=False, onerror=walk_error):
        current = Path(directory)
        included_directories = []
        for name in sorted(directories):
            if name == "__pycache__" or (current == root and name in EXCLUDED_ROOT_DIRECTORIES):
                continue
            _no_symlinks(current / name)
            included_directories.append(name)
        directories[:] = included_directories
        for name in sorted(filenames):
            if name == "__pycache__" or (current == root and (name == "MANIFEST.json" or name in EXCLUDED_ROOT_DIRECTORIES)) or _excluded_file(name):
                continue
            path = current / name
            relative = path.relative_to(root).as_posix()
            files[relative] = _file_record(_safe_path(root, relative))
    # Catch a preserved file changing between its first check and inventorying.
    # A preserved payload must also be covered by the deliverable manifest.
    for relative, record in expected.items():
        if files.get(relative) != record:
            raise ManifestBuildError(f"Preserved file changed during inventory or was excluded: {relative}")
    return {"schema": 1, "scope": SCOPE, "files": dict(sorted(files.items()))}


def _atomic_write(path: Path, manifest: dict[str, Any]) -> None:
    _no_symlinks(path)
    descriptor, temporary = tempfile.mkstemp(prefix=".MANIFEST.json.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def build_manifest(root: Path, *, check: bool = False) -> dict[str, Any]:
    """Build or check a manifest; root is an internal API, not a CLI override."""
    root = Path(root).absolute()
    manifest = collect_manifest(root)
    path = _safe_path(root, "MANIFEST.json")
    if check:
        existing = _read_json(path)
        if existing != manifest:
            previous = existing.get("files", {}) if isinstance(existing, dict) else {}
            if not isinstance(previous, dict):
                previous = {}
            current = manifest["files"]
            missing = sorted(set(current) - set(previous))
            stale = sorted(set(previous) - set(current))
            changed = sorted(key for key in set(previous) & set(current) if previous[key] != current[key])
            raise ManifestBuildError(
                "Manifest is not the exact current inventory; "
                f"unlisted={missing}, absent={stale}, changed={changed}; schema/scope must also match"
            )
    else:
        _atomic_write(path, manifest)
    return {"ok": True, "action": "checked" if check else "written", "files": len(manifest["files"]),
            "bytes": sum(record["bytes"] for record in manifest["files"].values()), "scope": SCOPE}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify exact inventory without writing")
    args = parser.parse_args(argv)
    try:
        print(json.dumps(build_manifest(Path(__file__).resolve().parents[1], check=args.check), indent=2))
        return 0
    except (OSError, ManifestBuildError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
