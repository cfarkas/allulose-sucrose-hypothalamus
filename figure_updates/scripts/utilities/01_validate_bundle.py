#!/usr/bin/env python3
"""Validate the self-contained Paper/Zenodo bundle contract.

This checks the shape of what actually gets uploaded. The directories in
LOCAL_ONLY are working state that README2.txt says to leave out of the archive,
so what is inside them is not the bundle's problem and is skipped: an earlier
version walked the whole tree and failed on two relative symlinks inside
recovery_archive, which meant this utility could never pass.

It reads and reports; it writes nothing and executes no figure step.
scripts/utilities/03_check_reproducibility.py is the complementary check on the
numbered scripts themselves.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from collections import defaultdict
from pathlib import Path, PurePosixPath


HERE = Path(__file__).resolve()
PAPER = next(
    (
        parent
        for parent in HERE.parents
        if (
            parent.name == "Paper"
            or (
                (parent / "scripts" / "setup").is_dir()
                and (parent / "README.txt").is_file()
            )
        )
    ),
    None,
)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

# Working state and material that the public release contract excludes.
LOCAL_ONLY = {
    "analyses",
    "legacy",
    "quarantine",
    "recovery_archive",
    "recovery_support",
    "docs",
    ".claude",
    ".mypy_cache",
    "__pycache__",
}
LOCAL_ONLY_ROOT_PREFIXES = ("apotome_rebuild_",)
LOCAL_ONLY_FILES = {"README2.txt"}
LOCAL_ONLY_FILE_SUFFIXES = (".orig", ".rej")
# The Allen atlas experiment and former Figure S3 cross-slide organ-alignment
# experiment were not used for the intended paper analyses. Exclude both exact
# prefixes while retaining microscope tile/channel stitching receipts required
# to reconstruct acquisitions for Figures 2–4.
PUBLIC_EXCLUDED_RELATIVE_PREFIXES = (
    "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot",
    "FigS3/registration",
    "literature_webscrap_30_08_2026",
)
# These downloaded article copies are not required for figure reproduction.
# An index-provided open-access URL is not, by itself, sufficient evidence that
# every third-party PDF can be redistributed under this package's licenses.
# Normal bibliographic references remain in the manuscript; the working scrape
# and its adjacent third-party metadata are excluded as one exact root prefix.
REDISTRIBUTION_EXCLUDED_BASENAMES = {"full_text_open_access.pdf"}
# scripts/shared and scripts/setup hold named helpers that the numbered scripts
# import or call; only the per-figure mirrors carry the NN_ naming rule.
NUMBERED_SCRIPT_ROOTS = (
    "Fig1",
    "Fig2",
    "Fig3",
    "Fig4",
    "Fig5",
    "FigS1",
    "FigS2",
    "FigS3",
    "FigS4",
    "FigS5",
    "FigS6",
)
FIGURES = NUMBERED_SCRIPT_ROOTS
SOURCE_ONLY_FIGURES: dict[str, str] = {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_public_relative_path(relative: str, *, is_directory: bool) -> bool:
    """Return whether one normalized Paper-relative path is public."""
    pure = PurePosixPath(relative)
    if (
        not relative
        or "\\" in relative
        or pure.is_absolute()
        or pure.as_posix() != relative
        or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        return False
    parts = pure.parts
    directory_parts = parts if is_directory else parts[:-1]
    if any(part in LOCAL_ONLY for part in directory_parts):
        return False
    if (
        is_directory
        and len(parts) == 1
        and any(parts[0].startswith(prefix) for prefix in LOCAL_ONLY_ROOT_PREFIXES)
    ):
        return False
    if not is_directory:
        if len(parts) == 1 and parts[0] in LOCAL_ONLY_FILES:
            return False
        if parts[-1].endswith(LOCAL_ONLY_FILE_SUFFIXES):
            return False
        if parts[-1] in REDISTRIBUTION_EXCLUDED_BASENAMES:
            return False
    return not any(
        relative == prefix or relative.startswith(prefix + "/")
        for prefix in PUBLIC_EXCLUDED_RELATIVE_PREFIXES
    )


def bundle_entries(root: Path):
    """Every path that would be packed into the archive."""
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        dirnames[:] = [
            name
            for name in dirnames
            if is_public_relative_path(
                (base / name).relative_to(root).as_posix(), is_directory=True
            )
        ]
        for name in dirnames:
            yield base / name
        for name in filenames:
            path = base / name
            relative = path.relative_to(root).as_posix()
            if not is_public_relative_path(relative, is_directory=False):
                continue
            yield path


def check_symlinks(entries) -> list[str]:
    links = [path for path in entries if path.is_symlink()]
    return [f"symbolic link in the bundle: {path}" for path in links]


def check_hard_links(entries, *, verify_hashes: bool) -> tuple[list[str], list[str]]:
    """Hard links are fine only if the bytes also exist independently in the bundle.

    Fig3/raw is a second view of the Figure 3 acquisition in which most files are
    hard links to originals elsewhere on the authoring machine. That is recorded
    in Fig3/raw_data/RAW_LINK_MANIFEST.csv, and Fig3/raw_data holds the same
    bytes as independent files, so packing the archive loses nothing. A hard link
    with no independent counterpart would be a real hole in self-containment.
    """
    problems: list[str] = []
    notes: list[str] = []
    linked = [
        path
        for path in entries
        if path.is_file() and not path.is_symlink() and path.stat().st_nlink > 1
    ]
    if not linked:
        return problems, notes

    independent: dict[tuple[str, int], list[Path]] = defaultdict(list)
    for path in entries:
        if path.is_file() and not path.is_symlink() and path.stat().st_nlink == 1:
            independent[(path.name, path.stat().st_size)].append(path)

    per_root: dict[str, int] = defaultdict(int)
    for path in linked:
        relative = path.relative_to(PAPER)
        per_root[relative.parts[0]] += 1
        candidates = independent.get((path.name, path.stat().st_size), [])
        if not candidates:
            problems.append(
                f"hard-linked file with no independent copy in the bundle: {relative}"
            )
            continue
        if verify_hashes and sha256(path) not in {
            sha256(other) for other in candidates
        }:
            problems.append(
                f"hard-linked file whose independent copies differ in content: {relative}"
            )
    for root, count in sorted(per_root.items()):
        notes.append(
            f"{count} hard-linked files under {root}/, each with an independent "
            f"copy elsewhere in the bundle; an archive materialises them as bytes"
        )
    return problems, notes


def check_script_naming(entries) -> tuple[list[str], list[str]]:
    """Only things you run need NN_ prefixes; imported helpers do not.

    Each of FigS1, FigS2 and FigS4 keeps a copy of spanish_panel_text.py beside
    its renderer, which imports it by name. Those copies are byte-identical to
    scripts/shared/spanish_panel_text.py, so they are a helper sitting where its
    importer can see it, not an unnumbered entry point.
    """
    problems, notes = [], []
    shared = PAPER / "scripts" / "shared"
    shared_helpers = (
        {
            path.name: sha256(path)
            for path in shared.glob("*.py")
            if not path.name[:2].isdigit()
        }
        if shared.is_dir()
        else {}
    )
    figure_local_helpers = {
        ("FigS3", "figs3_common.py"),
    }
    for path in entries:
        if not path.is_file() or path.suffix.lower() not in {".py", ".sh"}:
            continue
        relative = path.relative_to(PAPER)
        parts = relative.parts
        numbered_zone = (parts[0] in NUMBERED_SCRIPT_ROOTS and len(parts) == 2) or (
            parts[0] == "scripts"
            and len(parts) > 1
            and parts[1] in NUMBERED_SCRIPT_ROOTS
            and len(parts) == 3
        )
        if not numbered_zone or path.name[:2].isdigit():
            continue
        if len(parts) == 2 and (parts[0], parts[1]) in figure_local_helpers:
            notes.append(f"{relative} is an imported figure-local helper")
            continue

        expected = shared_helpers.get(path.name)
        if expected is not None and sha256(path) == expected:
            notes.append(
                f"{relative} is a byte-exact copy of scripts/shared/{path.name}, "
                f"placed where its importer can find it"
            )
            continue
        problems.append(
            f"unnumbered script at the front of a figure, and not a copy of a "
            f"scripts/shared helper: {relative}"
        )
    return problems, notes


def check_figure_shape(entries_set: set[Path]) -> tuple[list[str], list[str]]:
    problems, notes = [], []
    for figure in FIGURES:
        directory = PAPER / figure
        if not directory.is_dir():
            problems.append(f"missing figure directory: {figure}")
            continue
        if not (directory / "README.txt").is_file():
            problems.append(f"{figure} has no README.txt")
        numbered = sorted(
            path.name
            for path in directory.iterdir()
            if path.is_file() and path.name[:2].isdigit()
        )
        if not numbered:
            problems.append(f"{figure} has no numbered scripts at its front")
        graphics = sorted(
            path.name
            for path in directory.iterdir()
            if path.is_file()
            and path.name.startswith("Figure_")
            and path.suffix.lower() in {".pdf", ".png"}
        )
        if not graphics:
            if figure in SOURCE_ONLY_FIGURES:
                notes.append(
                    f"{figure} is the source dataset for {SOURCE_ONLY_FIGURES[figure]}"
                )
            else:
                notes.append(f"{figure} has no promoted Figure_* graphic")
    return problems, notes


def validate_manifest(manifest: Path, verify_hashes: bool) -> tuple[int, int]:
    rows = manifest.read_text(encoding="utf-8").splitlines()
    if not rows or rows[0] != (
        "source_relative_to_Apotome\ttarget_relative_to_Paper\tsize_bytes\tsha256"
    ):
        raise RuntimeError(f"Unexpected Figure 2 manifest header: {manifest}")
    count = 0
    total = 0
    for line_no, line in enumerate(rows[1:], start=2):
        source, target, size_text, expected_hash = line.split("\t")
        del source
        path = PAPER / target
        expected_size = int(size_text)
        if not path.is_file():
            raise FileNotFoundError(f"Manifest line {line_no} is missing: {path}")
        observed_size = path.stat().st_size
        if observed_size != expected_size:
            raise RuntimeError(
                f"Size mismatch at line {line_no}: {path}: "
                f"{observed_size} != {expected_size}"
            )
        if verify_hashes and sha256(path) != expected_hash:
            raise RuntimeError(f"SHA-256 mismatch at line {line_no}: {path}")
        count += 1
        total += observed_size
    return count, total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hash-large-inputs",
        action="store_true",
        help="Also hash every large Figure 2 manifest input, and every hard-link "
        "counterpart. Slow.",
    )
    parser.add_argument(
        "--strict-hard-links",
        action="store_true",
        help="Treat any hard-linked file as a problem, even one with an "
        "independent copy in the bundle.",
    )
    args = parser.parse_args()

    entries = list(bundle_entries(PAPER))
    files = [path for path in entries if path.is_file() and not path.is_symlink()]
    total_bytes = sum(path.stat().st_size for path in files)

    problems: list[str] = []
    notes: list[str] = []

    problems += check_symlinks(entries)
    link_problems, link_notes = check_hard_links(
        entries, verify_hashes=args.hash_large_inputs
    )
    if args.strict_hard_links:
        problems += link_problems + [
            f"hard links present: {note}" for note in link_notes
        ]
    else:
        problems += link_problems
        notes += link_notes
    naming_problems, naming_notes = check_script_naming(entries)
    problems += naming_problems
    notes += naming_notes
    shape_problems, shape_notes = check_figure_shape(set(entries))
    problems += shape_problems
    notes += shape_notes

    manifest = PAPER / "analyses/Fig2/provenance/minimal_required_files.tsv"
    if manifest.is_file():
        count, total = validate_manifest(manifest, args.hash_large_inputs)
        notes.append(f"Figure 2 minimal-input manifest: {count} files, {total} bytes")
    else:
        notes.append(
            "Figure 2 minimal-input manifest is absent; it lived under analyses/, "
            "which is working state and is not part of the bundle"
        )

    print(f"Bundle root: {PAPER}")
    print(f"Bundle contents: {len(files)} files, {total_bytes / 1024**3:.1f} GB")
    excluded = sorted(LOCAL_ONLY | LOCAL_ONLY_FILES)
    excluded.extend(f"{prefix}*" for prefix in LOCAL_ONLY_ROOT_PREFIXES)
    excluded.extend(f"*{suffix}" for suffix in LOCAL_ONLY_FILE_SUFFIXES)
    excluded.extend(sorted(REDISTRIBUTION_EXCLUDED_BASENAMES))
    excluded.extend(f"{prefix}/**" for prefix in PUBLIC_EXCLUDED_RELATIVE_PREFIXES)
    print(f"Excluded from the public archive: {', '.join(excluded)}")
    print()
    for note in notes:
        print(f"  note     {note}")
    for problem in problems:
        print(f"  PROBLEM  {problem}")
    print()
    print(
        "OVERALL: PASS" if not problems else f"OVERALL: FAIL ({len(problems)} problems)"
    )
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
