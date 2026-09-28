#!/usr/bin/env python3
"""Validate that every figure-required Cellpose mask accompanies its raw image."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve()
DEFAULT_PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if DEFAULT_PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")


@dataclass(frozen=True)
class Coverage:
    figure: str
    root: Path
    expected: int
    images: tuple[Path, ...]


def figure_2_coverage(paper: Path) -> Coverage:
    root = paper / "analyses" / "Fig2" / "raw" / "datasets"
    images = tuple(
        sorted(
            path
            for path in root.rglob("*")
            if path.is_file()
            and path.suffix.casefold() in {".tif", ".tiff"}
            and path.stem.casefold().endswith(("dapi", "fos"))
        )
    )
    return Coverage("Figure 2", root, 50, images)


def figure_4_coverage(paper: Path) -> Coverage:
    root = paper / "analyses" / "Fig4" / "raw" / "input"
    plane_pattern = re.compile(r"^S\d+_(?:left|right)(?:_T\d+)?$", re.IGNORECASE)
    marker_pattern = re.compile(r"_(?:DAPI|NPY|POMC|cFOS)$", re.IGNORECASE)
    images: list[Path] = []
    if root.is_dir():
        for sample in sorted(path for path in root.iterdir() if path.is_dir()):
            for plane in sorted(path for path in sample.iterdir() if path.is_dir()):
                if not plane_pattern.fullmatch(plane.name):
                    continue
                images.extend(
                    path
                    for path in sorted(plane.iterdir())
                    if path.is_file()
                    and path.suffix.casefold() in {".tif", ".tiff"}
                    and marker_pattern.search(path.stem)
                )
    return Coverage("Figure 4", root, 491, tuple(images))


def figure_5_coverage(paper: Path) -> Coverage:
    root = paper / "analyses" / "Fig5" / "raw" / "input"
    section_pattern = re.compile(r"^S\d+$", re.IGNORECASE)
    marker_pattern = re.compile(r"_(?:DAPI|Iba1)$", re.IGNORECASE)
    images: list[Path] = []
    if root.is_dir():
        for sample in sorted(path for path in root.iterdir() if path.is_dir()):
            for section in sorted(path for path in sample.iterdir() if path.is_dir()):
                if not section_pattern.fullmatch(section.name):
                    continue
                images.extend(
                    path
                    for path in sorted(section.iterdir())
                    if path.is_file()
                    and path.suffix.casefold() in {".tif", ".tiff"}
                    and marker_pattern.search(path.stem)
                )
    return Coverage("Figure 5", root, 122, tuple(images))


BUILDERS = {
    "fig2": figure_2_coverage,
    "fig4": figure_4_coverage,
    "fig5": figure_5_coverage,
}


def validate(coverage: Coverage) -> list[str]:
    errors: list[str] = []
    if not coverage.root.is_dir():
        return [f"raw root is missing: {coverage.root}"]
    if len(coverage.images) != coverage.expected:
        errors.append(
            f"required raw-image inventory is {len(coverage.images)}, expected {coverage.expected}"
        )
    seen: set[Path] = set()
    for image in coverage.images:
        mask = image.with_name(f"{image.stem}_seg.npy")
        if mask in seen:
            errors.append(f"duplicate mask resolution: {mask}")
            continue
        seen.add(mask)
        if not mask.is_file():
            errors.append(f"missing mask beside raw image: {mask}")
        elif mask.stat().st_size <= 0:
            errors.append(f"empty mask beside raw image: {mask}")
        elif mask.is_symlink():
            errors.append(f"mask must be a physical raw-data file, not a symlink: {mask}")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figures", nargs="+", choices=(*BUILDERS, "all"))
    parser.add_argument("--paper-root", type=Path, default=DEFAULT_PAPER_ROOT)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paper = args.paper_root.expanduser().resolve()
    requested = list(BUILDERS) if "all" in args.figures else list(dict.fromkeys(args.figures))
    failed = False
    for name in requested:
        coverage = BUILDERS[name](paper)
        errors = validate(coverage)
        total_masks = sum(1 for path in coverage.root.rglob("*_seg.npy") if path.is_file())
        if errors:
            failed = True
            print(f"[FAIL] {coverage.figure}: {len(errors)} raw-local Cellpose coverage error(s)", file=sys.stderr)
            for error in errors[:20]:
                print(f"  - {error}", file=sys.stderr)
            if len(errors) > 20:
                print(f"  - ... and {len(errors) - 20} more", file=sys.stderr)
        else:
            print(
                f"[OK] {coverage.figure}: {len(coverage.images)}/{coverage.expected} "
                f"required masks accompany raw images ({total_masks} total *_seg.npy files under raw)"
            )
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
