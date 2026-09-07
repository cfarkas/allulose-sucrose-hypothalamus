#!/usr/bin/env python3
"""Export reviewed figure code, small inputs and reference renders to a Git clone.

The original Zenodo manifests and their DOIs remain unchanged. The result is a
separately hashed source update that can be installed with figure_updates.py.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil

SCHEMA = "apotome-figure-updates-v1"


def selected_files(paper):
    selected = set()
    for name in ("Fig1", "FigS6"):
        for path in (paper / name).rglob("*"):
            if path.is_file() and not any(part in {"legacy", "__pycache__"} for part in path.relative_to(paper).parts):
                selected.add(path)
    for name in ("Fig5", "scripts/Fig1", "scripts/Fig5", "scripts/FigS6", "scripts/shared", "scripts/utilities", "scripts/setup"):
        for path in (paper / name).glob("*"):
            if path.is_file() and path.suffix in {".py", ".sh", ".txt", ".yml"}:
                selected.add(path)
    selected.update((paper / "scripts/shared/tests").glob("test_microglia_*.py"))
    selected.update(p for p in (paper / "Fig5/microglial_review_data").glob("*") if p.is_file())
    for name in ("README.txt", "scripts/README.txt", "reproduce_all_figures.sh", "release/build_archives.py",
                 "release/reconstruct_paper.py", "release/README.md", "release/export_figure_updates.py",
                 "release/build_microglial_review_data.py",
                 "release/tests/test_build_archives.py", "release/github_public/figure_updates.py",
                 "release/github_public/build_archives.py"):
        selected.add(paper / name)
    return sorted(selected)


def export(paper, repository):
    paper, repository = paper.resolve(), repository.resolve()
    destination = repository / "figure_updates"
    manifest_path = repository / "figure-updates.json"
    if destination.exists() or manifest_path.exists():
        raise ValueError("Export requires absent figure_updates/ and figure-updates.json; preserve an existing update first")
    base_path = repository / "manifest.json"
    base = json.loads(base_path.read_text())
    old = {row["path"]: row for row in base["tree"]["files"]}
    rows = []
    files = selected_files(paper)
    for source in files:
        if source.is_symlink() or not source.is_file() or source.stat().st_size > 10 * 1024 * 1024:
            raise ValueError(f"Absent, linked or unexpectedly large source: {source}")
    for source in files:
        relative = source.relative_to(paper).as_posix()
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        mode = 0o755 if source.stat().st_mode & 0o111 else 0o644
        target.chmod(mode)
        sha = hashlib.sha256(target.read_bytes()).hexdigest()
        rows.append({"path": relative, "sha256": sha, "size": target.stat().st_size,
                     "mode": f"{mode:04o}", "change": "added" if relative not in old else
                     "unchanged" if sha == old[relative]["sha256"] else "updated"})
    manifest = {"schema": SCHEMA, "date": "2026-09-07",
                "description": "Figure 1 Holm terminology, Figure 5 reviewed segmentation/classifier code, and supplementary Figure S6",
                "base_manifest_sha256": hashlib.sha256(base_path.read_bytes()).hexdigest(),
                "archived_figure_count": 10, "updated_figure_count": 11,
                "zenodo_status": "Existing records still describe the archived release; new versions have not been published for this update.",
                "files": rows}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    helper = paper / "release/github_public/figure_updates.py"
    shutil.copyfile(helper, repository / "figure_updates.py")
    (repository / "figure_updates.py").chmod(0o644)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--repository", type=Path, required=True)
    args = parser.parse_args()
    manifest = export(args.paper_root, args.repository)
    print(f"Exported {len(manifest['files'])} files ({sum(r['size'] for r in manifest['files']) / 2**20:.1f} MiB)")


if __name__ == "__main__":
    main()
