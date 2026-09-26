#!/usr/bin/env python3
"""Refresh an existing GitHub clone's figure_updates/ payload from this Paper tree.

The public repository carries a `figure_updates/` mirror of the figure
directories and a hashed `figure-updates.json` manifest beside it.
`export_figure_updates.py` builds that payload from scratch and refuses to run
when one already exists. This script does the other half: it brings an existing
payload up to date with the current published figures, file by file, and
rewrites the manifest.

What it touches is decided by the manifest, not by a hard-coded list:

  * every tracked path is re-copied from this Paper tree, and recorded as
    `unchanged` or `updated` by comparing SHA-256 before and after;
  * a tracked path that no longer exists in the Paper tree is deleted from the
    payload and recorded as `removed`;
  * a file that is new in the Paper tree is added only when its directory is
    already tracked and its suffix already appears in that directory, so the
    payload grows with the figures it already carries and never with a
    directory the repository deliberately excludes.

Nothing is committed or pushed; the result is a dirty working tree for review.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import date
from pathlib import Path
import shutil

PAPER = Path(__file__).resolve().parents[1]
SCHEMA = "apotome-figure-updates-v1"
PAYLOAD = "figure_updates"
MANIFEST = "figure-updates.json"
# GitHub rejects blobs above 100 MB; stay clear of the limit.
MAXIMUM_BYTES = 90 * 1024 * 1024
SKIP_PARTS = {"__pycache__", "legacy", ".ipynb_checkpoints"}
# Files at the Paper root are curated by hand: README.txt and the launcher are
# published, the author's machine-specific notes and incident logs are not. A
# new root file is therefore never added automatically.
ADD_ROOT_FILES = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_of(path: Path) -> str:
    return "0755" if path.stat().st_mode & 0o111 else "0644"


def candidate_additions(tracked: set[str]) -> list[str]:
    """New Paper files in directories the payload already carries."""
    directories: dict[str, set[str]] = {}
    for relative in tracked:
        parent = str(Path(relative).parent)
        directories.setdefault(parent, set()).add(Path(relative).suffix.lower())
    found = []
    for parent, suffixes in directories.items():
        folder = PAPER / parent if parent != "." else PAPER
        if not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            if not path.is_file() or path.is_symlink():
                continue
            if any(part in SKIP_PARTS for part in path.relative_to(PAPER).parts):
                continue
            relative = path.relative_to(PAPER).as_posix()
            if parent == "." and not ADD_ROOT_FILES:
                continue
            if relative not in tracked and path.suffix.lower() in suffixes:
                found.append(relative)
    return sorted(found)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repository", type=Path, required=True, help="Path to the GitHub clone")
    parser.add_argument("--description", required=True, help="One paragraph for the manifest's description field")
    parser.add_argument("--zenodo-status", default=("Existing records still describe the archived release; new versions "
                                                    "have not been published for this update."))
    parser.add_argument("--no-additions", action="store_true", help="Refresh tracked files only; add nothing new")
    args = parser.parse_args()

    repository = args.repository.expanduser().resolve()
    payload, manifest_path = repository / PAYLOAD, repository / MANIFEST
    if not payload.is_dir() or not manifest_path.is_file():
        raise SystemExit(f"Not an exported repository: {repository}")
    manifest = json.loads(manifest_path.read_text())
    tracked = {row["path"] for row in manifest["files"]}
    previous = {row["path"]: row["sha256"] for row in manifest["files"]}

    additions = [] if args.no_additions else candidate_additions(tracked)
    rows, counts = [], {"unchanged": 0, "updated": 0, "added": 0, "removed": 0}
    for relative in sorted(tracked | set(additions)):
        source, target = PAPER / relative, payload / relative
        if not source.is_file():
            if target.exists():
                target.unlink()
            counts["removed"] += 1
            rows.append({"path": relative, "sha256": None, "size": 0, "mode": "0644", "change": "removed"})
            continue
        if source.stat().st_size > MAXIMUM_BYTES:
            raise SystemExit(f"Refusing a file above {MAXIMUM_BYTES // 2**20} MiB: {relative}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        digest = sha256(target)
        target.chmod(0o755 if source.stat().st_mode & 0o111 else 0o644)
        change = ("added" if relative not in previous else
                  "unchanged" if previous[relative] == digest else "updated")
        counts[change] += 1
        rows.append({"path": relative, "sha256": digest, "size": target.stat().st_size,
                     "mode": mode_of(target), "change": change})

    kept = [row for row in rows if row["change"] != "removed"]
    manifest.update(schema=SCHEMA, date=date.today().isoformat(), description=args.description,
                    zenodo_status=args.zenodo_status,
                    updated_figure_count=len({row["path"].split("/")[0] for row in kept
                                              if row["change"] in {"updated", "added"}
                                              and row["path"].startswith("Fig")}),
                    files=kept)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    summary = dict(repository=str(repository), paper=str(PAPER), counts=counts,
                   tracked_after=len(kept),
                   payload_mib=round(sum(row["size"] for row in kept) / 2**20, 1),
                   manifest=str(manifest_path), manifest_sha256=sha256(manifest_path))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
