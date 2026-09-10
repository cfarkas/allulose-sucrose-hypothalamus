#!/usr/bin/env python3
"""Fail-closed, no-delete staging driver for recovered Figure 4 and Figure 5.

This driver deliberately separates raw/HIL checking, interactive review,
preparation, analysis, and rendering.  It never passes a force/replacement flag
to an analyzer or renderer and never removes an existing path.  Generated
outputs must be placed below an explicit staging root outside ``Paper``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


# Repinned 2026-08-25 for the Make Figure output-location fix, the Figure 5
# default paths and the name-independent paper-root discovery. Superseded
# values, kept so an older staged run can still be identified:
#   scripts/Fig4/02_analyze_pomc_cfos.py
#     d6071b74d6bc0ce333ea3f13af65798e58b6974ae4ba9ba89b24d5d3f1bc3c20
#   scripts/Fig4/03_make_figure_4_pomc_cfos.py
#     b7a4029aa1e3a6101b0e4db08600988a8e12b3105ad0eda22f3864bf73c75cc1
#   scripts/Fig5/01_analyze_gfap_iba1_microglia.py
#     e0da7f2a2b62328328a5377a6355b81fb9f0d1472b3063fb96840703bff780a8
#   scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py
#     085e8538e5e487fd33dc33f4fc860f5e6c81fa3163155bb5e55199a102459607
#   scripts/shared/01_annotate_regions.py
#     ece6aef4890ea8aa7d50476a558673e3547a75f535a49e8947aba8f7be9c742a
#   scripts/shared/02_validate_raw_cellpose_masks.py
#     b806220de779e433ebfacead65da12954dd272fbf2bae67080e10d96c5f6346b
# Repinned 2026-08-25 (second pass) for the Figure 5 render defaults: its
# --output-dir no longer points at the live Fig5 directory and Panel B's
# representative moved from FR8-3 to FR8-1. Superseded values:
#   scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py
#     85f281738bf7509485edb6a2be533b2f3e3749b95798be160c461dd4e002bd08
# Repinned 2026-08-25 (third pass) for the microglia nomenclature change,
# the second Figure 5 reviewer, and the Figure 4 NCC rounding. Superseded:
#   scripts/Fig4/02_analyze_pomc_cfos.py
#     8fc5c28ae74fceb1578763f8a5e22fd004554c0c319f5b4af300cfc09a0f871e
#   scripts/Fig5/01_analyze_gfap_iba1_microglia.py
#     bd2e18ea1af260dcb6c922833f24443c59de8119a81996268f0972e646e7c81d
#   scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py
#     e79f6c4e70b95794b8fe165cda3f5035e86317ac61c68b0430f6fb1942aa2ac8
# Repinned 2026-08-27 for the Figure 4 E/F reconstruction: E is the exact
# A+B+C microscopy merge with its cell magnification, and F adds the registered
# NPY-GFP intensity channel plus the enlarged distribution inset. Superseded:
#   scripts/Fig4/03_make_figure_4_pomc_cfos.py
#     eeabd498ea3fe89697d49e6cdb7fa2d02401f98dfb5c33d86893348f7470c191
# Repinned 2026-08-27 for the ROI-painted Figure 4F, human-HIL-ARC upper
# ventricle correction in D, FR722 S02 condition labels in A-C, and wider
# occurrence heatmaps in I/J. Superseded values:
#   scripts/Fig4/03_make_figure_4_pomc_cfos.py
#     aea0595ab54e842bd92159e490ccdb882f75d19902b06afe28988bbec0e18d55
#   scripts/Fig4/05_analyze_spatial_distributions.py
#     2aceef80c581af516e14298adc23337ad76f3af782115ac805d02c059ccd997b
# Repinned 2026-08-27 for real POMC intensity strictly masked to accepted POMC
# ROIs in Figure 4F and dashed white tissue/ventricle guides in its ROI inset.
# Superseded value:
#   scripts/Fig4/03_make_figure_4_pomc_cfos.py
#     743bee282de68469f9bfcbe89659ba39ea3a0b7e4f1d8084f6e8f4dc3449db64
# Repinned 2026-08-27 after the former Figure 5 workflow was renumbered to
# Figure 5. Active scripts and shared helpers now use Fig5 paths/names; the
# analyzer alone translates immutable pre-renumbering Fig5 paths recorded in
# historical HIL receipts after verifying the relocated files.
# Repinned 2026-08-31 for calibrated scale bars in all microscopy-derived
# Figure 4/5 views and Spanish-only collision repairs in Figure 5.
# Repinned 2026-08-31 after presentation terminology was standardized to HIL
# in English and Spanish figure labels and legends.
# Repinned 2026-09-01 after active workflow and QC presentation strings were normalized to HIL; stable CLI and receipt identifiers remain unchanged.
EXPECTED_HASHES = {
    "scripts/Fig4/02_analyze_pomc_cfos.py": "4f847cdd0521d9ef0f3c037b167ae06abf2f54f7927869867ab9e798c356b05c",
    "scripts/Fig4/03_make_figure_4_pomc_cfos.py": "471f8b974ad35f4251f7c27deb6b895b3e8b2e7c4a65ab09b4565109755bbe2f",
    "scripts/Fig4/05_analyze_spatial_distributions.py": "224e60de4fbb1049e070c717f9597916336bdce32440264ecadb987e90648051",
    "scripts/Fig5/01_analyze_gfap_iba1_microglia.py": "6e889ca1d822a408ddefc303bf18989628186e828565a8cc431762990614a698",
    "scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py": "c826581833f411c74f2cbaa27e264a9e3103850f31460f86b5bc7806f82d493c",
    "scripts/shared/01_annotate_regions.py": "189fbba65700a4b20db030f0ec0c018fd3e877fc6d5c594a43eb8e28c349c03d",
    "scripts/shared/02_validate_raw_cellpose_masks.py": "020b64b1d11638dd5601702696bc482d59f0dc139e02be7457b4e89871c10dab",
}

FIG4_ANALYZER = "scripts/Fig4/02_analyze_pomc_cfos.py"
FIG4_RENDERER = "scripts/Fig4/03_make_figure_4_pomc_cfos.py"
FIG4_SPATIAL = "scripts/Fig4/05_analyze_spatial_distributions.py"
FIG5_ANALYZER = "scripts/Fig5/01_analyze_gfap_iba1_microglia.py"
FIG5_RENDERER = "scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py"
REVIEWER = "scripts/shared/01_annotate_regions.py"
MASK_VALIDATOR = "scripts/shared/02_validate_raw_cellpose_masks.py"


class ContractError(RuntimeError):
    """A prerequisite or safety contract was not met."""


@dataclass(frozen=True)
class RuntimePaths:
    figure: str
    paper: Path
    stage: Path
    raw: Path
    analysis: Path
    hil: Path
    publication: Path
    python: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def lexical_absolute(path: Path) -> Path:
    expanded = path.expanduser()
    return expanded if expanded.is_absolute() else Path.cwd() / expanded


def reject_symlink_components(path: Path) -> None:
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        if current.is_symlink():
            raise ContractError(f"symbolic-link path component is forbidden: {current}")


def validate_stage_root(stage_arg: Path, paper: Path, *, create: bool) -> Path:
    lexical = lexical_absolute(stage_arg)
    reject_symlink_components(lexical)
    stage = lexical.resolve(strict=False)
    filesystem_root = Path(stage.anchor).resolve()
    if stage == filesystem_root:
        raise ContractError(f"filesystem root cannot be a staging root: {stage}")
    if stage == paper or is_within(stage, paper) or is_within(paper, stage):
        raise ContractError(
            f"staging root must be outside and must not contain Paper: stage={stage}, Paper={paper}"
        )
    if stage.exists() and not stage.is_dir():
        raise ContractError(f"staging root exists but is not a directory: {stage}")
    if create:
        stage.mkdir(parents=True, exist_ok=True)
    return stage


def require_stage_output(path_arg: Path, stage: Path, purpose: str, *, must_be_absent: bool) -> Path:
    lexical = lexical_absolute(path_arg)
    reject_symlink_components(lexical)
    path = lexical.resolve(strict=False)
    if path == stage or not is_within(path, stage):
        raise ContractError(f"{purpose} must be a strict descendant of staging root {stage}: {path}")
    if must_be_absent and path.exists():
        raise ContractError(
            f"{purpose} must not exist; this wrapper never replaces or cleans outputs: {path}"
        )
    return path


def require_regular_independent_file(path: Path, purpose: str) -> None:
    if path.is_symlink():
        raise ContractError(f"{purpose} is a symbolic link: {path}")
    if not path.is_file():
        raise ContractError(f"{purpose} is missing or not a regular file: {path}")
    stat = path.stat()
    if stat.st_size <= 0:
        raise ContractError(f"{purpose} is empty: {path}")
    if stat.st_nlink != 1:
        raise ContractError(
            f"{purpose} is not an independent copy (link count {stat.st_nlink}): {path}"
        )


def verify_recovered_scripts(paper: Path, figure: str) -> dict[str, str]:
    required = [REVIEWER, MASK_VALIDATOR]
    required.extend(
        [FIG4_ANALYZER, FIG4_SPATIAL, FIG4_RENDERER]
        if figure == "fig4"
        else [FIG5_ANALYZER, FIG5_RENDERER]
    )
    verified: dict[str, str] = {}
    for relative in required:
        path = paper / relative
        if path.is_symlink() or not path.is_file():
            raise ContractError(f"required recovered script is missing or linked: {path}")
        actual = sha256_file(path)
        expected = EXPECTED_HASHES[relative]
        if actual != expected:
            raise ContractError(
                f"recovered-script hash mismatch for {relative}: expected {expected}, got {actual}"
            )
        verified[relative] = actual
    return verified


def marker_from_stem(stem: str, allowed: Sequence[str]) -> str | None:
    for marker in allowed:
        if re.search(rf"_{re.escape(marker)}$", stem, flags=re.IGNORECASE):
            return marker
    return None


def validate_fig4_raw(raw: Path) -> dict[str, int]:
    if raw.is_symlink() or not raw.is_dir():
        raise ContractError(f"Figure 4 raw input root is missing, linked, or not a directory: {raw}")
    plane_re = re.compile(r"^S(?P<section>\d+)_(?P<side>left|right)(?:_T(?P<tile>\d+))?$", re.I)
    images: list[Path] = []
    identities: dict[tuple[str, int, str], set[str]] = {}
    for sample in sorted(path for path in raw.iterdir() if path.is_dir() and not path.is_symlink()):
        for plane in sorted(path for path in sample.iterdir() if path.is_dir() and not path.is_symlink()):
            match = plane_re.fullmatch(plane.name)
            if not match:
                continue
            identity = (sample.name, int(match.group("section")), match.group("tile") or "")
            side = match.group("side").lower()
            if side in identities.setdefault(identity, set()):
                raise ContractError(f"duplicate Figure 4 plane side for {identity}: {side}")
            identities[identity].add(side)
            for path in sorted(plane.iterdir()):
                if not path.is_file() or path.suffix.casefold() not in {".tif", ".tiff"}:
                    continue
                if marker_from_stem(path.stem, ("DAPI", "NPY", "POMC", "cFOS")):
                    images.append(path)
    if len(identities) != 76:
        raise ContractError(
            f"Figure 4 reconstructed-section identity inventory is {len(identities)}, expected 76"
        )
    # Repinned 2026-08-25 from 491 to 481. 491 was measured against the author's
    # private analyses/Fig4/raw/input; the acquisition that ships as Fig4/raw_data
    # holds 481, which its own RAW_DATA_MANIFEST.csv confirms sample by sample,
    # and which is what the 76-section human-final analysis reads.
    if len(images) != 481:
        raise ContractError(f"Figure 4 marker-TIFF inventory is {len(images)}, expected 481")
    masks: set[Path] = set()
    for image in images:
        require_regular_independent_file(image, "Figure 4 raw microscopy image")
        mask = image.with_name(f"{image.stem}_seg.npy")
        require_regular_independent_file(mask, "Figure 4 raw-local Cellpose mask")
        if mask in masks:
            raise ContractError(f"duplicate Figure 4 mask resolution: {mask}")
        masks.add(mask)
    return {"sections": len(identities), "marker_tiffs": len(images), "required_masks": len(masks)}


def validate_fig5_raw(raw: Path) -> dict[str, int]:
    if raw.is_symlink() or not raw.is_dir():
        raise ContractError(f"Figure 5 raw input root is missing, linked, or not a directory: {raw}")
    section_re = re.compile(r"^S\d+$", re.I)
    sections = 0
    images: list[Path] = []
    masks: set[Path] = set()
    for sample in sorted(path for path in raw.iterdir() if path.is_dir() and not path.is_symlink()):
        for section in sorted(path for path in sample.iterdir() if path.is_dir() and not path.is_symlink()):
            if not section_re.fullmatch(section.name):
                continue
            by_marker: dict[str, list[Path]] = {"DAPI": [], "GFAP": [], "Iba1": []}
            for path in sorted(section.iterdir()):
                if not path.is_file() or path.suffix.casefold() not in {".tif", ".tiff"}:
                    continue
                marker = marker_from_stem(path.stem, tuple(by_marker))
                if marker:
                    by_marker[marker].append(path)
            for marker, matches in by_marker.items():
                if len(matches) != 1:
                    raise ContractError(
                        f"Figure 5 {sample.name}/{section.name} has {len(matches)} {marker} TIFFs, expected 1"
                    )
                image = matches[0]
                require_regular_independent_file(image, f"Figure 5 {marker} raw microscopy image")
                images.append(image)
                if marker in {"DAPI", "Iba1"}:
                    mask = image.with_name(f"{image.stem}_seg.npy")
                    require_regular_independent_file(mask, f"Figure 5 {marker} raw-local Cellpose mask")
                    if mask in masks:
                        raise ContractError(f"duplicate Figure 5 mask resolution: {mask}")
                    masks.add(mask)
            sections += 1
    if sections != 61:
        raise ContractError(f"Figure 5 section inventory is {sections}, expected 61")
    if len(images) != 183:
        raise ContractError(f"Figure 5 raw channel-TIFF inventory is {len(images)}, expected 183")
    if len(masks) != 122:
        raise ContractError(f"Figure 5 required Cellpose-mask inventory is {len(masks)}, expected 122")
    return {"sections": sections, "channel_tiffs": len(images), "required_masks": len(masks)}


def validate_raw(paths: RuntimePaths) -> dict[str, int]:
    return validate_fig4_raw(paths.raw) if paths.figure == "fig4" else validate_fig5_raw(paths.raw)


def run_command(command: Sequence[str], *, capture: bool = False, stdin_null: bool = False) -> subprocess.CompletedProcess[str]:
    print("[COMMAND]", shlex.join([str(item) for item in command]), flush=True)
    completed = subprocess.run(
        [str(item) for item in command],
        check=False,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        stdin=subprocess.DEVNULL if stdin_null else None,
    )
    if completed.returncode:
        detail = ""
        if capture:
            detail = f"\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        raise ContractError(f"command failed with status {completed.returncode}: {shlex.join(command)}{detail}")
    return completed


def reviewer_command(paths: RuntimePaths, *, validate_only: bool) -> list[str]:
    input_root = (
        paths.analysis / "reconstructed_sections"
        if paths.figure == "fig4"
        else paths.raw
    )
    command = [
        str(paths.python), str(paths.paper / REVIEWER),
        "--dataset", "pomc" if paths.figure == "fig4" else "gfap",
        "--input-root", str(input_root),
        "--analysis-dir", str(paths.analysis),
        "--output-dir", str(paths.hil),
        "--host", "127.0.0.1", "--port", "0",
    ]
    if validate_only:
        command.append("--validate-only")
        return command
    # An interactive review session exists to segment sections, including ones
    # that were accepted earlier, so it always runs in replacement mode. This is
    # the HIL annotation writer, not an analyzer or renderer: it rewrites its
    # own JSON/TIFF receipt pair atomically and never deletes an output.
    command.append("--force")
    command.extend([
        "--figure-command", shlex.join(figure_button_command(paths)),
        "--figure-analysis-dir", str(paths.analysis),
        "--figure-output-root", str(reviewer_render_root(paths)),
    ])
    return command


def reviewer_render_root(paths: RuntimePaths) -> Path:
    """Staging parent for figures rendered from the reviewer's Make Figure button."""
    fig_name = "Fig4" if paths.figure == "fig4" else "Fig5"
    return paths.stage / fig_name / "reviewer_renders"


def figure_button_command(paths: RuntimePaths) -> list[str]:
    """Renderer the reviewer's Make Figure 4/5 button runs, into {output_dir}."""
    placeholder = Path("{output_dir}")
    staged = RuntimePaths(
        figure=paths.figure, paper=paths.paper, stage=paths.stage, raw=paths.raw,
        analysis=paths.analysis, hil=paths.hil, publication=placeholder,
        python=paths.python,
    )
    return fig4_render_command(staged) if paths.figure == "fig4" else fig5_render_command(staged)


def parse_reviewer_json(text: str) -> dict[str, object]:
    try:
        record = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ContractError(f"reviewer validation did not emit one JSON document: {exc}\n{text}") from exc
    if not isinstance(record, dict):
        raise ContractError("reviewer validation JSON is not an object")
    return record


def validate_hil(paths: RuntimePaths, *, require_complete: bool) -> dict[str, object]:
    completed = run_command(reviewer_command(paths, validate_only=True), capture=True, stdin_null=True)
    record = parse_reviewer_json(completed.stdout)
    expected = 76 if paths.figure == "fig4" else 61
    sections = int(record.get("sections", -1))
    accepted = int(record.get("accepted", -1))
    invalid = int(record.get("invalid_existing_pairs", -1))
    if sections != expected:
        raise ContractError(f"{paths.figure} HIL source inventory is {sections}, expected {expected}")
    if invalid != 0 or record.get("ok") is not True:
        raise ContractError(f"{paths.figure} HIL has invalid receipt pairs: {invalid}")
    if require_complete and accepted != expected:
        raise ContractError(f"{paths.figure} accepted HIL coverage is {accepted}/{expected}, expected complete coverage")
    if paths.figure == "fig4" and int(record.get("native_automated_mask_available", -1)) != expected:
        raise ContractError("Figure 4 HIL sources do not have 76/76 native automated masks")
    return record


def common_fig4_analysis_args(paths: RuntimePaths) -> list[str]:
    return [
        "--input", str(paths.raw), "--raw-input", str(paths.raw),
        "--output", str(paths.analysis),
        "--no-segment-missing", "--no-auto-conda",
        "--region-mode", "arc-me-walls", "--stitch-align", "dapi-overlap",
        "--stitch-align-max-shift-px", "180",
        "--stitch-align-medial-fraction", "0.55",
        "--stitch-align-smooth-px", "10",
        "--stitch-align-ventral-quantile", "0.985",
        "--stitch-align-anchor-max-corr-loss", "0.035",
        "--stitch-overlap-min-px", "24",
        "--stitch-overlap-max-fraction", "0.40",
        "--stitch-sift-ratio", "0.72",
        "--stitch-sift-ransac-px", "2.5",
        "--stitch-sift-min-inliers", "6",
        "--stitch-sift-min-inlier-fraction", "0.45",
        "--stitch-sift-max-residual-p95-px", "3",
        "--stitch-sift-min-y-span-px", "25",
        "--stitch-sift-min-scale", "0.98",
        "--stitch-sift-max-scale", "1.02",
        "--stitch-sift-max-rotation-deg", "1.5",
        "--stitch-overlap-min-raw-ncc", "0.85",
        "--stitch-overlap-min-gradient-ncc", "0.70",
        "--stitch-ncc-template-width-px", "24",
        "--stitch-ncc-min-score", "0.72",
        "--stitch-ncc-min-uniqueness-margin", "0.08",
        "--stitch-seam-search-fraction", "0.60",
        "--stitch-seam-max-touching-unmatched", "0",
        "--scalebar-um", "100", "--plot-suite", "final", "--verbose",
    ]


def fig4_prepare_command(paths: RuntimePaths) -> list[str]:
    return [
        str(paths.python), str(paths.paper / FIG4_ANALYZER),
        "--prepare-human-review", *common_fig4_analysis_args(paths),
    ]


def fig4_registration_check_command(paths: RuntimePaths) -> list[str]:
    return [
        str(paths.python), str(paths.paper / FIG4_ANALYZER),
        "--validate-stitch-registration-only", *common_fig4_analysis_args(paths),
    ]


def fig4_render_command(paths: RuntimePaths) -> list[str]:
    return [
        str(paths.python), str(paths.paper / FIG4_RENDERER),
        "--analysis-dir", str(paths.analysis),
        "--manual-region-dir", str(paths.hil), "--require-manual-regions",
        "--outdir", str(paths.publication), "--scalebar-um", "100",
        "--inset-field-um", "50", "--inset-scalebar-um", "20",
        "--um-per-px", "0.755", "--dpi", "600",
        "--spatial-dir", str(paths.analysis / "spatial"),
    ]


def fig4_spatial_command(paths: RuntimePaths) -> list[str]:
    return [
        str(paths.python), str(paths.paper / FIG4_SPATIAL),
        "--analysis-dir", str(paths.analysis),
        "--output-dir", str(paths.analysis / "spatial"),
        "--dpi", "600",
    ]


def fig5_analysis_command(paths: RuntimePaths, args: argparse.Namespace) -> list[str]:
    command = [
        str(paths.python), str(paths.paper / FIG5_ANALYZER),
        "--full-run", "--input", str(paths.raw), "--raw-input", str(paths.raw),
        "--output", str(paths.analysis), "--no-segment-missing", "--no-auto-conda",
        "--um-per-px", "1.51", "--human_in_the_loop", "--require-manual-regions",
        "--manual-region-dir", str(paths.hil),
        "--reviewed-dapi-root", str(paths.paper / "Fig5" / "hil_review" / "dapi_final_20260830_v1"),
        "--reviewed-iba1-root", str(paths.paper / "Fig5" / "hil_review" / "iba1_final_20260830_v1"),
        "--microglia-classifier-mode", str(args.classifier_mode),
        "--cnn-device", str(args.cnn_device), "--cnn-epochs", str(args.cnn_epochs),
        "--cnn-min-training-animals", "3", "--scalebar-um", "100",
        "--plot-suite", "final", "--verbose",
    ]
    if args.labels_csv is not None:
        command.extend(["--microglia-classifier-mode", "supervised", "--microglia-labels-csv", str(args.labels_csv)])
    return command


def fig5_render_command(paths: RuntimePaths) -> list[str]:
    return [
        str(paths.python), str(paths.paper / FIG5_RENDERER),
        "--analysis-dir", str(paths.analysis), "--input-root", str(paths.raw),
        "--manual-region-dir", str(paths.hil), "--output-dir", str(paths.publication),
        "--um-per-px", "1.51", "--dpi", "600", "--formats", "pdf,png",
    ]


def assert_no_replacement_flag(command: Sequence[str]) -> None:
    """Guard analyzer and renderer commands only; see reviewer_command."""
    forbidden = "--" + "force"
    if forbidden in command:
        raise ContractError(f"internal safety error: replacement flag appeared in command: {shlex.join(command)}")


def require_outputs(paths: RuntimePaths) -> dict[str, int]:
    if paths.figure == "fig4":
        required_analysis = [
            "per_reconstructed_section_summary.csv",
            "per_cell_reconstructed_measurements.csv",
            "definitive_endpoint_values_long.csv",
            "definitive_statistics.csv",
        ]
        panel_suffixes = {
            "A": "dapi", "B": "cfos", "C": "pomc",
            "D": "cellpose_cartoon", "E": "microscopy",
            "F": "microscopy_npy_gfp",
            "G": "cfos_dapi", "H": "pomc_activation",
            "I": "spatial_cfos", "J": "spatial_cfos_pomc",
        }
        required_publication = [
            "Figure_ARC_ME_multipanel.pdf", "Figure_ARC_ME_multipanel.png",
            *[
                f"panels/Figure_ARC_ME_panel_{letter}_{suffix}{locale}.{extension}"
                for letter, suffix in panel_suffixes.items()
                for locale in ("", "_spanish")
                for extension in ("pdf", "png")
            ],
            "provenance/all_section_stitch_registration_audit.csv",
            "provenance/multipanel_source_manifest.csv",
            "provenance/multipanel_spanish_translation_receipt.json",
        ]
        geometry_glob = "reconstructed_sections/**/stitch_registration.json"
        expected_geometry = 76
    else:
        required_analysis = [
            "human_review_coverage.csv", "discovered_sections.csv",
            "per_section_summary.csv", "per_section_region_summary.csv",
            "per_microglia_cell_measurements.csv", "per_animal_region_summary.csv",
            "gfap_iba1_endpoint_values_long.csv", "gfap_iba1_definitive_statistics.csv",
            "microglia_classifier_metadata.json",
        ]
        stem = "Figure_GFAP_Iba1_microglia_ARC_ME"
        required_publication = [
            f"{stem}.pdf", f"{stem}.png",
            *[f"panels/{stem}_Panel_{letter}.pdf" for letter in "ABCDEFGHI"],
            f"provenance/{stem}_manifest.csv",
            f"source_data/{stem}_all_HIL_final_mask_audit.csv",
        ]
        geometry_glob = "processed_sections/**/region_geometry.json"
        expected_geometry = 61
    for relative in required_analysis:
        path = paths.analysis / relative
        if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
            raise ContractError(f"required staged analysis artifact is absent/empty/linked: {path}")
    for relative in required_publication:
        path = paths.publication / relative
        if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
            raise ContractError(f"required staged publication artifact is absent/empty/linked: {path}")
    geometry = [path for path in paths.analysis.glob(geometry_glob) if path.is_file() and not path.is_symlink()]
    if len(geometry) != expected_geometry:
        raise ContractError(
            f"staged geometry/registration receipt count is {len(geometry)}, expected {expected_geometry}"
        )
    for root in (paths.analysis, paths.publication):
        linked = next((path for path in root.rglob("*") if path.is_symlink()), None)
        if linked is not None:
            raise ContractError(f"staged output contains a symbolic link: {linked}")
    return {
        "analysis_required_files": len(required_analysis),
        "publication_required_files": len(required_publication),
        "geometry_receipts": len(geometry),
    }


def write_receipt(paths: RuntimePaths, mode: str, payload: dict[str, object]) -> Path:
    receipt_dir = paths.stage / "receipts"
    receipt_dir.mkdir(parents=True, exist_ok=True)
    target = receipt_dir / f"{paths.figure}_{mode}_receipt.json"
    if target.exists():
        raise ContractError(f"receipt already exists; refusing to replace it: {target}")
    record = {
        "schema_version": "fig45_safe_stage_receipt_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "figure": paths.figure,
        "mode": mode,
        "paper_root": str(paths.paper),
        "stage_root": str(paths.stage),
        "raw_root": str(paths.raw),
        "analysis_dir": str(paths.analysis),
        "hil_dir": str(paths.hil),
        "publication_dir": str(paths.publication),
        "no_delete_policy": True,
        "analyzer_replacement_flag_used": False,
        **payload,
    }
    temporary = target.with_name(target.name + f".tmp.{os.getpid()}")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    return target


def resolve_paths(args: argparse.Namespace) -> RuntimePaths:
    paper = lexical_absolute(args.paper_root).resolve()
    if paper.is_symlink() or not paper.is_dir():
        raise ContractError(f"Paper root is missing, linked, or not a directory: {paper}")
    stage = validate_stage_root(args.stage_root, paper, create=args.mode != "check")
    fig_name = "Fig4" if args.figure == "fig4" else "Fig5"
    # analyses/FigN/raw/input is the author's private working copy of the
    # acquisition and is absent in a freshly unpacked tree; the acquisition
    # itself ships as FigN/raw_data, so the default has to fall through to it.
    raw_default = next(
        (candidate for candidate in (
            paper / "analyses" / fig_name / "raw" / "input",
            paper / fig_name / "raw_data",
        ) if candidate.is_dir()),
        paper / fig_name / "raw_data",
    )
    analysis_default = stage / fig_name / "analysis"
    hil_default = stage / fig_name / "hil" / f"accepted_annotations_rebuilt_{date.today():%Y%m%d}"
    publication_default = stage / fig_name / "publication"
    python = lexical_absolute(args.python).resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise ContractError(f"Python interpreter is missing or not executable: {python}")
    return RuntimePaths(
        figure=args.figure,
        paper=paper,
        stage=stage,
        raw=lexical_absolute(args.raw_root or raw_default).resolve(strict=False),
        analysis=lexical_absolute(args.analysis_dir or analysis_default).resolve(strict=False),
        hil=lexical_absolute(args.hil_dir or hil_default).resolve(strict=False),
        publication=lexical_absolute(args.figure_dir or publication_default).resolve(strict=False),
        python=python,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Fail-closed Figure 4/5 staging driver. It verifies exact recovered scripts, "
            "raw microscopy/Cellpose coverage, and complete HIL receipts; it never deletes "
            "or replaces an existing output."
        ),
    )
    parser.add_argument("--figure", choices=("fig4", "fig5"), required=True)
    parser.add_argument("--mode", choices=("check", "prepare-hil", "review", "build", "render"), default="check")
    parser.add_argument("--paper-root", type=Path, required=True)
    parser.add_argument("--stage-root", type=Path, required=True, help="External, non-Paper staging root")
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument("--hil-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--classifier-mode", choices=("hybrid", "cnn", "morphometry", "supervised"), default="hybrid")
    parser.add_argument("--cnn-device", choices=("cpu", "cuda", "auto"), default="cpu")
    parser.add_argument("--cnn-epochs", type=int, default=40)
    parser.add_argument("--labels-csv", type=Path)
    return parser


def execute(args: argparse.Namespace) -> int:
    paths = resolve_paths(args)
    scripts = verify_recovered_scripts(paths.paper, paths.figure)

    if args.mode in {"prepare-hil", "review", "build"}:
        mutable: list[tuple[Path, str, bool]] = []
        if args.mode == "prepare-hil":
            mutable.append((paths.analysis, "analysis preparation output", True))
        elif args.mode == "review":
            mutable.append((paths.hil, "versioned HIL output", False))
        elif args.mode == "build":
            mutable.extend([
                (paths.analysis, "analysis output", True),
                (paths.publication, "publication output", True),
            ])
        for path, purpose, absent in mutable:
            require_stage_output(path, paths.stage, purpose, must_be_absent=absent)
        if args.mode == "review" and not re.fullmatch(r"accepted_annotations_rebuilt_\d{8}", paths.hil.name):
            raise ContractError(
                "new HIL output must be versioned as accepted_annotations_rebuilt_YYYYMMDD"
            )
    elif args.mode == "render":
        require_stage_output(paths.publication, paths.stage, "publication output", must_be_absent=True)

    raw = validate_raw(paths)

    if args.mode == "prepare-hil":
        if paths.figure != "fig4":
            raise ContractError(
                "Figure 5 reviews native DAPI directly; use --mode review. Its exact analyzer "
                "does not expose a safe proposal-only phase before HIL."
            )
        paths.analysis.parent.mkdir(parents=True, exist_ok=True)
        command = fig4_prepare_command(paths)
        assert_no_replacement_flag(command)
        run_command(command, stdin_null=True)
        registration = fig4_registration_check_command(paths)
        assert_no_replacement_flag(registration)
        run_command(registration, stdin_null=True)
        hil_status = validate_hil(paths, require_complete=False)
        receipt = write_receipt(paths, args.mode, {"script_sha256": scripts, "raw": raw, "hil_status": hil_status})
        print(f"[OK] Figure 4 HIL preparation is staged; receipt: {receipt}")
        return 0

    if args.mode == "review":
        if paths.figure == "fig4" and not paths.analysis.is_dir():
            raise ContractError("Figure 4 review requires a completed staged prepare-hil analysis directory")
        paths.hil.parent.mkdir(parents=True, exist_ok=True)
        reviewer_render_root(paths).mkdir(parents=True, exist_ok=True)
        # assert_no_replacement_flag deliberately does not apply here: the guard
        # exists to keep a replacement flag out of the analyzers and renderers,
        # and the reviewer needs it to let a human segment a section again.
        command = reviewer_command(paths, validate_only=False)
        run_command(command)
        hil_status = validate_hil(paths, require_complete=True)
        receipt = write_receipt(paths, args.mode, {"script_sha256": scripts, "raw": raw, "hil_status": hil_status})
        print(f"[OK] Complete HIL validated; receipt: {receipt}")
        return 0

    hil_status = validate_hil(paths, require_complete=True)
    if args.mode == "check":
        print(json.dumps({
            "ok": True, "figure": paths.figure, "script_sha256": scripts,
            "raw": raw, "hil": hil_status,
            "analysis_dir": str(paths.analysis), "publication_dir": str(paths.publication),
        }, indent=2, sort_keys=True))
        return 0

    if args.mode == "build":
        if paths.figure == "fig4":
            raise ContractError(
                "Figure 4 fresh strict-HIL analysis cannot be run by this no-replacement wrapper: "
                "the exact analyzer validates receipts against reconstructed files in its already-"
                "nonempty output and then requires replacement mode. Use prepare-hil/check now; "
                "an analyzer transactional-output change is required before safe final inference."
            )
        paths.analysis.parent.mkdir(parents=True, exist_ok=True)
        paths.publication.parent.mkdir(parents=True, exist_ok=True)
        analysis_command = fig5_analysis_command(paths, args)
        render_command = fig5_render_command(paths)
        for command in (analysis_command, render_command):
            assert_no_replacement_flag(command)
        run_command(analysis_command, stdin_null=True)
        run_command(render_command, stdin_null=True)
    elif args.mode == "render":
        if not paths.analysis.is_dir():
            raise ContractError(f"analysis input is missing: {paths.analysis}")
        paths.publication.parent.mkdir(parents=True, exist_ok=True)
        if paths.figure == "fig4":
            spatial_dir = paths.analysis / "spatial"
            if not spatial_dir.exists():
                spatial_command = fig4_spatial_command(paths)
                assert_no_replacement_flag(spatial_command)
                run_command(spatial_command, stdin_null=True)
            command = fig4_render_command(paths)
        else:
            command = fig5_render_command(paths)
        assert_no_replacement_flag(command)
        run_command(command, stdin_null=True)

    outputs = require_outputs(paths)
    receipt = write_receipt(
        paths, args.mode,
        {"script_sha256": scripts, "raw": raw, "hil_status": hil_status, "outputs": outputs},
    )
    print(f"[OK] {paths.figure} staged {args.mode} completed; receipt: {receipt}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cnn_epochs <= 0:
        parser.error("--cnn-epochs must be positive")
    try:
        return execute(args)
    except ContractError as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
