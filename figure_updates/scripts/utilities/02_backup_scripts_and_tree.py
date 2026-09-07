#!/usr/bin/env python3
"""Back up every Paper script plus a full directory listing into one ZIP.

The archive holds source files only — no raw data, no figures, no binaries —
plus DIRECTORY_TREE.txt and MANIFEST.csv. It never traverses outside the Paper
tree, and never touches Minion_Data or storage_disk2.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import zipfile
from datetime import datetime, timezone
from pathlib import Path

PAPER = next(
    path for path in Path(__file__).resolve().parents
    if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                                 and (path / "README.txt").is_file()))
)
# Written beside the tree, never inside it. A scripts archive kept within the
# tree it was taken from goes stale against that tree and the two then disagree,
# which is why the old Paper-root copy was dropped on 2026-08-27.
DEFAULT_ZIP = PAPER.parent / f"{PAPER.name}_scripts_and_tree.zip"
SCRIPT_SUFFIXES = {".py", ".sh", ".R", ".r", ".ipynb", ".yml", ".yaml", ".patch"}
SKIP_DIRECTORIES = {"__pycache__", ".git", ".ipynb_checkpoints"}
MAX_SCRIPT_BYTES = 8 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def walk(root: Path):
    for path in sorted(root.rglob("*")):
        if any(part in SKIP_DIRECTORIES for part in path.parts):
            continue
        if path.is_symlink():
            continue
        yield path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_ZIP)
    args = parser.parse_args()
    output = args.output.resolve()
    if PAPER == output.parent or PAPER in output.parents:
        raise RuntimeError(
            "Refusing to write the scripts backup inside the Paper tree it "
            f"describes; choose a path outside {PAPER}"
        )

    tree_lines: list[str] = []
    rows: list[dict[str, object]] = []
    scripts: list[Path] = []
    directories = files = 0
    for path in walk(PAPER):
        if path == output:
            continue
        relative = path.relative_to(PAPER).as_posix()
        if path.is_dir():
            directories += 1
            tree_lines.append(f"{relative}/")
            continue
        files += 1
        size = path.stat().st_size
        tree_lines.append(f"{relative}\t{size}")
        if path.suffix in SCRIPT_SUFFIXES and size <= MAX_SCRIPT_BYTES:
            scripts.append(path)
            rows.append({"path": relative, "bytes": size, "sha256": sha256(path)})

    stamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    header = (
        f"Paper script and directory backup\n"
        f"created_utc\t{stamp}\n"
        f"root\t{PAPER}\n"
        f"directories\t{directories}\n"
        f"files\t{files}\n"
        f"scripts_archived\t{len(scripts)}\n"
        f"note\tsource files only; no raw data, figures, Minion_Data or storage_disk2\n\n"
    )
    manifest = io.StringIO()
    writer = csv.DictWriter(manifest, fieldnames=["path", "bytes", "sha256"])
    writer.writeheader()
    writer.writerows(rows)

    temporary = output.with_suffix(".zip.partial")
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        archive.writestr("DIRECTORY_TREE.txt", header + "\n".join(tree_lines) + "\n")
        archive.writestr("MANIFEST.csv", manifest.getvalue())
        for path in scripts:
            archive.write(path, f"Paper/{path.relative_to(PAPER).as_posix()}")
    temporary.replace(output)
    print(f"[OK] {output} — {len(scripts)} scripts, {files} files and {directories} directories listed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
