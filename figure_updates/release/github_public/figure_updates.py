#!/usr/bin/env python3
"""Verify or install the small GitHub figure update over the archived Paper tree."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile

SCHEMA = "apotome-figure-updates-v1"
MANIFEST = "figure-updates.json"
PAYLOAD = "figure_updates"
RECEIPT = ".figure-updates.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_path(root, relative):
    pure = PurePosixPath(relative)
    if (not relative or pure.is_absolute() or pure.as_posix() != relative
            or any(p in {".", ".."} for p in pure.parts) or "\\" in relative):
        raise ValueError(f"Unsafe update path: {relative!r}")
    path = root
    for part in pure.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f"Linked update path: {path}")
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Update path escapes its root: {relative}")
    return path


def load_updates(root):
    root = Path(root).resolve()
    path = safe_path(root, MANIFEST)
    update = json.loads(path.read_text())
    if update.get("schema") != SCHEMA or not update.get("files"):
        raise ValueError("Invalid figure update manifest")
    if update["base_manifest_sha256"] != digest(root / "manifest.json"):
        raise ValueError("Figure updates are bound to a different archived release")
    payload = safe_path(root, PAYLOAD)
    seen = set()
    for row in update["files"]:
        relative = row["path"]
        if relative in seen:
            raise ValueError(f"Duplicate update path: {relative}")
        seen.add(relative)
        source = safe_path(payload, relative)
        if (not source.is_file() or source.stat().st_size != row["size"]
                or digest(source) != row["sha256"] or row["mode"] not in {"0644", "0755"}):
            raise ValueError(f"Update file is absent or changed: {relative}")
        if bool(source.stat().st_mode & 0o111) != (row["mode"] == "0755"):
            raise ValueError(f"Update executable mode differs: {relative}")
    return update


def installed(root, paper):
    receipt = safe_path(paper, RECEIPT)
    if not receipt.exists():
        return False
    data = json.loads(receipt.read_text())
    if data.get("manifest_sha256") != digest(root / MANIFEST):
        raise ValueError("Installed figure update differs from this checkout; use its matching Git revision")
    return True


def effective_manifest(root, base):
    """Keep archived rows, substituting only authenticated, installed updates."""
    root = Path(root).resolve()
    paper = root / "Paper"
    if not (paper / RECEIPT).exists() and not (paper / RECEIPT).is_symlink():
        return base
    update = load_updates(root)
    installed(root, paper)
    rows = {row["path"]: row for row in base["tree"]["files"]}
    for row in update["files"]:
        rows[row["path"]] = {**rows.get(row["path"], {}), **row}
        # Frozen benchmark inputs must also be authenticated before replay.
        if row["path"].startswith("FigS6/raw_data/"):
            path = safe_path(paper, row["path"])
            if not path.is_file() or digest(path) != row["sha256"]:
                raise ValueError(f"Installed S6 input is absent or changed: {row['path']}")
    return {**base, "tree": {**base["tree"], "files": list(rows.values())}}


def atomic_copy(source, target, mode):
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".figure-update-", delete=False) as stream:
        temp = Path(stream.name)
        try:
            with source.open("rb") as incoming:
                shutil.copyfileobj(incoming, stream)
            stream.flush()
            os.fchmod(stream.fileno(), mode)
            os.replace(temp, target)
        finally:
            temp.unlink(missing_ok=True)


def apply_updates(root, paper=None):
    root = Path(root).resolve()
    paper = Path(paper) if paper is not None else root / "Paper"
    if paper.is_symlink() or not paper.is_dir():
        raise ValueError("Reconstruct an existing, non-linked Paper directory first")
    paper = paper.resolve()
    payload = root / PAYLOAD
    if paper == payload or paper.is_relative_to(payload) or payload.is_relative_to(paper):
        raise ValueError("The update destination overlaps its source")
    update = load_updates(root)
    if installed(root, paper):
        return {"status": "already_installed", "files": len(update["files"])}
    base = json.loads((root / "manifest.json").read_text())
    originals = {row["path"]: row for row in base["tree"]["files"]}
    planned = []
    # Validate the complete destination before creating or replacing any file.
    for row in update["files"]:
        target = safe_path(paper, row["path"])
        if target.exists():
            if not target.is_file():
                raise ValueError(f"Update destination is not a file: {target}")
            actual = digest(target)
            if actual == row["sha256"] and stat.S_IMODE(target.stat().st_mode) == int(row["mode"], 8):
                continue
            allowed = {row["sha256"], originals.get(row["path"], {}).get("sha256")}
            if actual not in allowed:
                raise ValueError(f"Preserving locally changed file; update stopped: {target}")
        elif row["path"] in originals:
            raise ValueError(f"Archived file is missing; complete reconstruction first: {target}")
        planned.append((row, target))
    backup_parent = safe_path(root, ".figure-update-backups")
    backup_parent.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-")
    backup = Path(tempfile.mkdtemp(prefix=stamp, dir=backup_parent))
    for row, target in planned:
        if target.exists():
            saved = backup / row["path"]
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(target, saved)
        atomic_copy(payload / row["path"], target, int(row["mode"], 8))
        if digest(target) != row["sha256"]:
            raise ValueError(f"Installed update hash differs: {target}")
    receipt = {"schema": SCHEMA, "manifest_sha256": digest(root / MANIFEST),
               "base_manifest_sha256": update["base_manifest_sha256"],
               "backup": str(backup), "installed_files": len(planned)}
    with tempfile.NamedTemporaryFile(mode="w", dir=paper, prefix=".figure-update-receipt-", delete=False) as stream:
        temp = Path(stream.name)
        try:
            json.dump(receipt, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.replace(temp, paper / RECEIPT)
        finally:
            temp.unlink(missing_ok=True)
    return {"status": "installed", **receipt}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--apply", action="store_true", help="Install after verifying all files; preserve replaced archived files in a backup.")
    parser.add_argument("--paper", type=Path, help="Existing reconstructed Paper directory (default: ROOT/Paper).")
    args = parser.parse_args()
    try:
        if args.apply:
            print(json.dumps(apply_updates(args.root, args.paper), indent=2))
        else:
            update = load_updates(args.root)
            print(f"PASS: {len(update['files'])} figure update files match their manifest.")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"STOPPED: {exc}\n")


if __name__ == "__main__":
    main()
