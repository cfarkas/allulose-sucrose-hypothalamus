#!/usr/bin/env python3
"""Snapshot or verify immutable raw-data and script inputs in the Paper tree.

The all-figure launcher renders in temporary stages and installs only validated
generated outputs. This utility makes that contract fail closed: every file in
every directory named ``raw`` or ``raw_data`` is inventoried by relative path,
kind, byte size and nanosecond mtime, while active scripts are SHA-256 hashed.
The post-install inventory must exactly match the pre-run snapshot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable


ACTIVE_FIGURES = (
    "Fig1", "Fig2", "Fig3", "Fig4", "Fig5",
    "FigS1", "FigS2", "FigS3", "FigS4", "FigS5", "FigS6",
)
SCRIPT_SUFFIXES = {".py", ".sh", ".R"}
PRUNED_TOP_LEVEL = {".claude"}
PRUNED_SCRIPT_DIRS = {"raw", "raw_data", "__pycache__"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def raw_roots(root: Path) -> list[Path]:
    found: list[Path] = []
    for current, dirs, _files in os.walk(root, topdown=True, followlinks=False):
        here = Path(current)
        if here == root:
            dirs[:] = sorted(
                name for name in dirs
                if name not in PRUNED_TOP_LEVEL
                and not name.startswith("apotome_rebuild_")
            )
            continue
        if here.name in {"raw", "raw_data"}:
            found.append(here)
            dirs[:] = []
            continue
        dirs.sort()
    return sorted(found, key=lambda path: relative(path, root))


def inventory_raw(root: Path, roots: Iterable[Path]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for raw_root in roots:
        for current, dirs, files in os.walk(raw_root, topdown=True, followlinks=False):
            here = Path(current)
            dirs.sort()
            files.sort()
            for name in files:
                path = here / name
                stat = path.lstat()
                record: dict[str, Any] = {
                    "path": relative(path, root),
                    "kind": "symlink" if path.is_symlink() else "file",
                    "bytes": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
                if path.is_symlink():
                    record["target"] = os.readlink(path)
                entries.append(record)
    return sorted(entries, key=lambda item: item["path"])


def script_candidates(root: Path) -> list[Path]:
    starts = [root / "scripts", *(root / figure for figure in ACTIVE_FIGURES)]
    scripts: list[Path] = []
    for start in starts:
        if not start.is_dir():
            continue
        for current, dirs, files in os.walk(start, topdown=True, followlinks=False):
            dirs[:] = sorted(name for name in dirs if name not in PRUNED_SCRIPT_DIRS)
            here = Path(current)
            for name in sorted(files):
                path = here / name
                if path.suffix in SCRIPT_SUFFIXES and not path.is_symlink():
                    scripts.append(path)
    unique = {relative(path, root): path for path in scripts}
    return [unique[name] for name in sorted(unique)]


def inventory_scripts(root: Path) -> list[dict[str, Any]]:
    return [
        {
            "path": relative(path, root),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in script_candidates(root)
    ]


def package_requirements(root: Path) -> dict[str, Any]:
    figures: list[dict[str, Any]] = []
    for figure in ACTIVE_FIGURES:
        figure_root = root / figure
        numbered = sorted(
            path.name for path in figure_root.iterdir()
            if path.is_file()
            and len(path.name) >= 3
            and path.name[:2].isdigit()
            and path.name[2] in {"_", "a", "b"}
            and path.suffix in SCRIPT_SUFFIXES
        )
        raw_dirs = sorted(
            relative(path, root) for path in figure_root.iterdir()
            if path.is_dir() and path.name in {"raw", "raw_data"}
        )
        if not numbered:
            raise SystemExit(f"No active numbered script at {figure_root}")
        if not raw_dirs:
            raise SystemExit(f"No raw/raw_data directory at {figure_root}")
        figures.append({
            "figure": figure,
            "numbered_scripts": numbered,
            "raw_directories": raw_dirs,
        })
    return {"active_figures": figures}


def make_snapshot(root: Path) -> dict[str, Any]:
    roots = raw_roots(root)
    raw = inventory_raw(root, roots)
    scripts = inventory_scripts(root)
    payload: dict[str, Any] = {
        "schema": "paper_canonical_immutable_inputs_v1",
        "root_basename": root.name,
        "raw_roots": [relative(path, root) for path in roots],
        "raw_entries": raw,
        "raw_file_count": len(raw),
        "raw_total_bytes": sum(int(item["bytes"]) for item in raw),
        "scripts": scripts,
        "script_count": len(scripts),
    }
    payload.update(package_requirements(root))
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--snapshot", type=Path)
    group.add_argument("--verify", type=Path)
    args = parser.parse_args()

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Paper root not found: {root}")
    current = make_snapshot(root)

    if args.snapshot is not None:
        destination = args.snapshot.expanduser().resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(current, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            f"[PASS] canonical input snapshot: {current['raw_file_count']} raw files, "
            f"{current['raw_total_bytes']} raw bytes, {current['script_count']} scripts"
        )
        return 0

    baseline_path = args.verify.expanduser().resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if current != baseline:
        raise SystemExit(
            "Canonical raw-data/script audit changed during the rebuild. "
            f"Baseline: {baseline_path}"
        )
    print(
        f"[PASS] canonical inputs unchanged: {current['raw_file_count']} raw files, "
        f"{current['raw_total_bytes']} raw bytes, {current['script_count']} scripts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

