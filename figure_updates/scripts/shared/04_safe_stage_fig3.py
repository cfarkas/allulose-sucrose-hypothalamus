#!/usr/bin/env python3
"""Plan or execute the recovered Figure 3 workflow in a new external stage.

Planning is the default and is read-only. Execution requires ``--execute``.
The launcher never removes, replaces, or reuses a path. A run receives one
absent direct child of an external staging parent, and a failed run is retained
in place for inspection.
"""

from __future__ import annotations

import argparse
import ast
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
from typing import Iterable, Mapping, Sequence


FIGURE_NAME = "Figure_3_cFos_NPY"
KNOWN_LIVE_PAPER = Path("/media/server/STORAGE/Apotome/Paper")

# Pinned identities of the recovered Figure 3 scripts.
#
# 2026-08-24: two entries were re-pinned. 02_make_cfos_npy_cartoons.py and
# 03_make_figure_3_cfos_npy.py resolved the Paper root by fixed parent depth,
# which breaks as soon as the tree is relocated or a copy of the script is run
# from the figure directory. Both now search their parents for the directory
# named "Paper". No analysis, statistic, drawing or output path changed.
# Superseded pins, for the record:
#   02_make_cfos_npy_cartoons.py
#     4b2542e09e7262ade16374917884e0ab8e86d6344a593696df13f4ead217069c
#   03_make_figure_3_cfos_npy.py
#     12edc1a991af769f4439ec2c7fa83865cf6b489edf3afe34832a566590cfe157
#
# 2026-08-27: the two display renderers were re-pinned after a user-requested
# Figure 3 correction. The raw-DAPI lumen source masks remain separately
# hash-bound, while the final Water and Allulose ventricular polygons now come
# from source-bound HIL acceptance receipts. Panel A NPY ROI delimiters changed
# from black to purple (#A020F0) without changing ROI membership or counts.
# Superseded pins:
#   02_make_cfos_npy_cartoons.py
#     5e7d4c6180b1365d26aa82e5d5878008d96ee6522e34c5a3511ddd7789ebae1a
#   03_make_figure_3_cfos_npy.py
#     cdc023eb944b847a60305d2ab92dc046418553a7725fb35d728695e0b2270434
# The pre-HIL safety snapshot retains those exact bytes and the older pin
# history. The 00/01/05 entries below were also reconciled to their already
# active front/mirror copies after the plan-only gate exposed stale pins.
#
# 2026-08-28: the compositor was re-pinned after adding a protected D/E tick-
# label gutter and a render-time check that prevents the F/G row from covering
# those labels. Plot values, statistics, and panel source images are unchanged.
# Superseded pin:
#   03_make_figure_3_cfos_npy.py
#     96a193a51df5b3db4bd1b1afd4ba61d20352e4c3e002e215224df9fed7d09208
#
# 2026-08-28: the compositor was re-pinned after the user requested that the
# Panel A ROI strokes be visibly violet rather than black. The renderer now uses
# CSS violet (#EE82EE); ROI geometry, membership, counts, and source intensities
# are unchanged. The English reproduction copy was synchronized at the same time.
# Superseded pin:
#   03_make_figure_3_cfos_npy.py
#     a9fb31c96504e5edda5d487e1ff65afce348caa6e41a0d2ff365e611dcb7f942
#
# 2026-08-28: the Panel A ROI stroke was darkened to #C05ACD and is now
# clipped to the exact microscopy-image extent before black axes letterboxing.
# ROI geometry, membership, counts, and source intensities remain unchanged.
# Superseded pin:
#   03_make_figure_3_cfos_npy.py
#     50243b8be13739316e3aaf4f4425901abf6cc0724213a2ffdd8f9694bd0d9ca0
# The finalized clip disables scaling per ROI line, leaving the image-set limits
# intact without changing axis-wide autoscaling or emitting aspect warnings.
# Superseded implementation pins:
#   03_make_figure_3_cfos_npy.py
#     bec2066b75a6d5ba4d8015b493c69b91f2188bfdd2435739de20136045452c7d
#     92f3176b7a829b246076413b9b76a57bf43ffbe6f995695af87f30b15a00847e
#
# 2026-08-28: the compositor was re-pinned after widening each F/G spatial
# heatmap from a 7.42-in allocation to 7.69 in. The wider native-aspect boxes
# reclaim unused side and bottom margins while preserving the D/E tick-label
# clearance guard, all source pixels, values, and statistics.
# Superseded pin:
#   03_make_figure_3_cfos_npy.py
#     b92cfb4f7b97e1bef1f37cd130205c1149c0c5104c3c37e5b6385205e965ae7c
#
# 2026-08-28: Panel A ROI contours returned to one black (#000000) 1.56-pt
# stroke while retaining the accepted masks, assignments, counts, complete-frame
# contour geometry, and microscopy-extent clipping. The spatial analyzer now
# places condition means on the left of a widened 3.934-in radial-profile axis
# on a 6.40 x 3.12-in source canvas. The compositor preserves native aspect in
# 7.721-in F/G master boxes and retains the D/E tick-clearance guard.
# Superseded pins:
#   03_make_figure_3_cfos_npy.py
#     635c02d6360bd0ed88bbdce478f524622e5404842cc3d9402a1d9ef64fce043c
#   05_analyze_spatial_distributions.py
#     267c6dee448046caf34d41aefdaa95eebc946481c7f764c110d10206464d73c1
#
# 2026-08-28: F/G were restored to their accepted semantic order: the
# inferential radial-profile plot is on the left and the descriptive
# condition-means heatmap is on the right. A moderate 2.55:1.45 width
# allocation yields measured axes of 3.152 and 1.541 inches, respectively,
# while preserving the 6.40 x 3.12-in source canvas, black Panel A ROI strokes,
# native master aspect, all values/statistics, and the D/E tick-clearance guard.
# Superseded pins:
#   03_make_figure_3_cfos_npy.py
#     584cdae66eb014652ea6146d3cfa1866d6ce79e1f53edaf9ac60523295bcd4d1
#   05_analyze_spatial_distributions.py
#     13e30b1eed6a431833794767a1608118bc5219788eb61f766f8a2422b9da0e96
RECOVERED_SCRIPT_SHA256: Mapping[str, str] = {
    "scripts/Fig3/00_extract_water3_neun_display.py":
        "bf2e7ced97766627ca0eee66f0a63f3eca9254ef6cf64c7d370b311f35310b9d",
    "scripts/Fig3/01_analyze_cfos_npy.py":
        "8d6f328de6857ec39e233761c73e34af63afa00ee5dcb27852b5c3d4f91fcb91",
    "scripts/Fig3/02_make_cfos_npy_cartoons.py":
        "6f8557521291fecdfef14a1372d41d90d2347ba762dffc106e4553a3dd0f1967",
    "scripts/Fig3/03_make_figure_3_cfos_npy.py":
        "d855867632cdeb02440a7a83935531486522806760e8f8ffc09ddad00bed57b3",
    "scripts/Fig3/05_analyze_spatial_distributions.py":
        "ad101c2991f979d60a00ab7e932b458efdd10103fbb0451786a10aa6ac887570",
}

ANALYSIS_FILENAMES = (
    "per_animal_cfos_npy.csv",
    "anova_mwu_statistics.csv",
    "statistical_test_receipts.csv",
    "segmentation_qc.csv",
    "source_manifest.csv",
    "README.txt",
)

CARTOON_FILENAMES = (
    "Figure3_Panel_B_Water_marker_positive_nuclei.png",
    "Figure3_Panel_B_Water_marker_positive_nuclei.pdf",
    "Figure3_Panel_C_Allulose_marker_positive_nuclei.png",
    "Figure3_Panel_C_Allulose_marker_positive_nuclei.pdf",
    "Figure3_cFOS_NPY_cartoon_provenance.json",
    "Figure3_cFOS_NPY_cartoon_provenance.csv",
)

PANEL_A_FILENAMES = (
    "Figure3_Panel_A_raw_microscopy.png",
    "Figure3_Panel_A_raw_microscopy.pdf",
    "Figure3_Panel_A_raw_microscopy_provenance.json",
)

SPATIAL_FILENAMES = (
    "animal_spatial_counts.csv",
    "animal_radial_occurrence_features.csv",
    "exact_spatial_permanova.csv",
    "exact_spatial_permutation_receipts.csv",
    "exact_spatial_dispersion.csv",
    "exact_spatial_dispersion_permutation_receipts.csv",
    "segmentation_model_setting_audit.csv",
    "Figure3_Spatial_cFOS_occurrence.png",
    "Figure3_Spatial_cFOS_occurrence.pdf",
    "Figure3_Spatial_cFOS_NPY_occurrence.png",
    "Figure3_Spatial_cFOS_NPY_occurrence.pdf",
    "spatial_distribution_provenance.json",
    "output_manifest.csv",
)

SOURCE_DATA_FILENAMES = (
    f"{FIGURE_NAME}_animal_plot_values.csv",
    f"{FIGURE_NAME}_anova_mwu_statistics_used.csv",
    f"{FIGURE_NAME}_supplemental_NPY_DAPI_animal_values.csv",
    f"{FIGURE_NAME}_supplemental_NPY_DAPI_anova_mwu_statistics.csv",
    f"{FIGURE_NAME}_supplemental_NPY_DAPI_README.txt",
    f"{FIGURE_NAME}_spatial_exact_permanova.csv",
    f"{FIGURE_NAME}_spatial_exact_permutation_receipts.csv",
    f"{FIGURE_NAME}_spatial_exact_dispersion.csv",
    f"{FIGURE_NAME}_spatial_exact_dispersion_permutation_receipts.csv",
    f"{FIGURE_NAME}_spatial_animal_radial_occurrence_features.csv",
    f"{FIGURE_NAME}_spatial_animal_counts.csv",
    f"{FIGURE_NAME}_spatial_segmentation_model_setting_audit.csv",
)

RAW_METADATA_RELPATHS = (
    "apotome_2025_07_28_npy_samplesheet.csv",
    "water/apotome_2025_07_28_npy_water_samplesheet.csv",
    "water/WATER_NPY3_S01_extraction_provenance.json",
    "water/WATER_NPY3_S01_NeuN_display_only_provenance.json",
)

WATER3_BINDINGS = {
    "water3_dapi_tiff": "water/WATER_NPY3_S01_dapi.tif",
    "water3_cfos_tiff": "water/WATER_NPY3_S01_fos.tif",
    "water3_npy_tiff": "water/WATER_NPY3_S01_GFP.tif",
    "water3_dapi_mask": "water/WATER_NPY3_S01_dapi_seg.npy",
    "water3_cfos_mask": "water/WATER_NPY3_S01_fos_seg.npy",
    "water3_npy_mask": "water/WATER_NPY3_S01_GFP_seg.npy",
    "water3_extraction_provenance": "water/WATER_NPY3_S01_extraction_provenance.json",
}

VENTRICLE_HIL_SCHEMA = "fig3_cartoon_ventricle_hil_v1"
VENTRICLE_HIL_FILENAMES = (
    "water_ventricle_mask.png",
    "water_ventricle_receipt.json",
    "allulose_ventricle_mask.png",
    "allulose_ventricle_receipt.json",
)


class SafetyError(RuntimeError):
    """Raised before an unsafe or incomplete staging operation."""


@dataclass(frozen=True)
class Layout:
    stage_root: Path
    shadow_paper: Path
    shadow_scripts: Path
    analysis: Path
    panels: Path
    spatial: Path
    publication: Path
    logs: Path
    runtime: Path
    provenance: Path


@dataclass(frozen=True)
class Step:
    name: str
    argv: tuple[str, ...]


@dataclass(frozen=True)
class Plan:
    paper_root: Path
    stage_parent: Path
    raw_project_root: Path
    raw_root: Path
    ventricle_hil_dir: Path
    analysis_python: Path
    layout: Layout
    source_scripts: Mapping[str, Path]
    steps: tuple[Step, ...]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def lexical_absolute(value: Path | str) -> Path:
    return Path(os.path.abspath(os.path.expanduser(str(value))))


def is_within_or_equal(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def paths_overlap(first: Path, second: Path) -> bool:
    return is_within_or_equal(first, second) or is_within_or_equal(second, first)


def require_no_symlink_components(path: Path, *, description: str) -> None:
    absolute = lexical_absolute(path)
    current = Path(absolute.anchor)
    for component in absolute.parts[1:]:
        current = current / component
        if current.is_symlink():
            raise SafetyError(f"{description} contains a symbolic-link component: {current}")
        if not current.exists():
            break


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def exclusive_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def exclusive_json(path: Path, payload: object) -> None:
    exclusive_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def validate_recovered_scripts(paper_root: Path) -> dict[str, Path]:
    require_no_symlink_components(paper_root, description="Paper root")
    if not paper_root.is_dir():
        raise SafetyError(f"Paper root is not a directory: {paper_root}")
    observed: dict[str, Path] = {}
    for relative, expected_hash in RECOVERED_SCRIPT_SHA256.items():
        source = paper_root / relative
        if source.is_symlink() or not source.is_file():
            raise SafetyError(f"Recovered Figure 3 script is missing or linked: {source}")
        if source.resolve().parent != (paper_root / relative).parent.resolve():
            raise SafetyError(f"Recovered script escapes its expected directory: {source}")
        actual_hash = sha256_file(source)
        if actual_hash != expected_hash:
            raise SafetyError(
                f"Recovered script hash mismatch for {relative}: "
                f"expected {expected_hash}, observed {actual_hash}"
            )
        observed[relative] = source.resolve()
    return observed


def raw_relpaths_from_analysis_script(script: Path) -> tuple[str, ...]:
    """Extract only literal SampleSpec paths from the hash-pinned source AST."""
    tree = ast.parse(script.read_text(encoding="utf-8"), filename=str(script))
    relpaths: list[str] = []
    sample_calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Name) or node.func.id != "SampleSpec":
            continue
        if len(node.args) < 9:
            continue
        sample_calls.append(node)
        for item in node.args[3:9]:
            value = ast.literal_eval(item)
            if not isinstance(value, str):
                raise SafetyError("A recovered SampleSpec raw path is not a literal string")
            relpaths.append(value)
        for keyword in node.keywords:
            if keyword.arg in {"neun_seg", "neun_tif"}:
                value = ast.literal_eval(keyword.value)
                if value is not None:
                    if not isinstance(value, str):
                        raise SafetyError("A recovered NeuN path is not a literal string")
                    relpaths.append(value)
    if len(sample_calls) != 9:
        raise SafetyError(f"Expected 9 literal SampleSpec mappings; found {len(sample_calls)}")
    relpaths.extend(RAW_METADATA_RELPATHS)
    unique = tuple(sorted(set(relpaths)))
    if len(unique) != 75:
        raise SafetyError(f"Expected 75 mapped raw/provenance companions; found {len(unique)}")
    for relative in unique:
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise SafetyError(f"Unsafe mapped raw relative path: {relative}")
    return unique


def build_layout(stage_root: Path) -> Layout:
    shadow_paper = stage_root / "Paper"
    analysis = shadow_paper / "analyses/Fig3/results/final_run"
    return Layout(
        stage_root=stage_root,
        shadow_paper=shadow_paper,
        shadow_scripts=shadow_paper / "scripts/Fig3",
        analysis=analysis,
        panels=analysis / "panels",
        spatial=shadow_paper / "analyses/Fig3/results/spatial_distribution_run",
        publication=shadow_paper / "Fig3",
        logs=stage_root / "logs",
        runtime=stage_root / "runtime",
        provenance=stage_root / "provenance",
    )


def validate_water3_bindings(args: argparse.Namespace, raw_root: Path) -> None:
    for argument, relative in WATER3_BINDINGS.items():
        supplied = lexical_absolute(getattr(args, argument))
        expected = lexical_absolute(raw_root / relative)
        if supplied != expected:
            option = "--" + argument.replace("_", "-")
            raise SafetyError(f"{option} must equal {expected}; received {supplied}")


def validate_plan_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path, Path, Path]:
    paper_root = lexical_absolute(args.paper_root)
    stage_parent = lexical_absolute(args.stage_parent)
    stage_root = lexical_absolute(args.stage_root)
    raw_project_root = lexical_absolute(args.raw_project_root)
    raw_root = lexical_absolute(args.raw_root)
    ventricle_hil_dir = lexical_absolute(args.ventricle_hil_dir)

    require_no_symlink_components(stage_parent, description="Stage parent")
    if not stage_parent.is_dir():
        raise SafetyError(f"Stage parent must already exist as a normal directory: {stage_parent}")
    if stage_root.parent != stage_parent:
        raise SafetyError(
            "Stage root must be an absent direct child (and therefore strict descendant) "
            f"of stage parent {stage_parent}: {stage_root}"
        )
    if stage_root.exists() or stage_root.is_symlink():
        raise SafetyError(f"Stage root must be absent and is never reused: {stage_root}")

    protected_roots = {paper_root}
    if KNOWN_LIVE_PAPER.exists():
        protected_roots.add(KNOWN_LIVE_PAPER.resolve())
    for protected in protected_roots:
        if paths_overlap(stage_root, protected):
            raise SafetyError(f"Stage root must be disjoint from live Paper: {stage_root} vs {protected}")
        if paths_overlap(stage_parent, protected):
            raise SafetyError(f"Stage parent must be external to live Paper: {stage_parent} vs {protected}")
        if paths_overlap(raw_root, protected) or paths_overlap(raw_project_root, protected):
            raise SafetyError(f"Raw inputs must be an independent tree outside live Paper: {protected}")

    expected_raw = raw_project_root / "Fig3/raw/legacy_experiment_2025_07_28"
    if raw_root != expected_raw:
        raise SafetyError(f"Raw root must equal {expected_raw}; received {raw_root}")
    if paths_overlap(stage_root, raw_root) or paths_overlap(stage_parent, raw_project_root):
        raise SafetyError("Stage and independent raw trees must be disjoint")
    require_no_symlink_components(ventricle_hil_dir, description="Ventricle HIL directory")
    if ventricle_hil_dir.is_symlink() or not ventricle_hil_dir.is_dir():
        raise SafetyError(
            f"Accepted ventricle HIL directory is missing or linked: {ventricle_hil_dir}"
        )
    if paths_overlap(stage_root, ventricle_hil_dir):
        raise SafetyError("Ventricle HIL source and new stage must be disjoint")
    validate_water3_bindings(args, raw_root)
    return paper_root, stage_parent, stage_root, raw_project_root, raw_root, ventricle_hil_dir


def build_steps(layout: Layout, raw_project_root: Path, raw_root: Path, python: Path,
                dpi: int, panel_a_dpi: int) -> tuple[Step, ...]:
    script = lambda name: str(layout.shadow_scripts / name)
    steps = (
        Step("01_analysis", (
            str(python), script("01_analyze_cfos_npy.py"),
            "--raw-root", str(raw_root),
            "--output-dir", str(layout.analysis),
        )),
        Step("02_cartoons", (
            str(python), script("02_make_cfos_npy_cartoons.py"),
            "--root", str(raw_project_root),
            "--output-dir", str(layout.panels),
            "--dpi", str(dpi),
            "--ventricle-hil-dir", str(layout.provenance / "accepted_ventricle_hil"),
        )),
        Step("05_spatial", (
            str(python), script("05_analyze_spatial_distributions.py"),
            "--raw-root", str(raw_root),
            "--output-dir", str(layout.spatial),
            "--dpi", str(dpi),
        )),
        Step("03_compositor", (
            str(python), script("03_make_figure_3_cfos_npy.py"),
            "--analysis-dir", str(layout.analysis),
            "--output-dir", str(layout.publication),
            "--raw-root", str(raw_root),
            "--spatial-dir", str(layout.spatial),
            "--per-animal", str(layout.analysis / "per_animal_cfos_npy.csv"),
            "--statistics", str(layout.analysis / "anova_mwu_statistics.csv"),
            "--panel-b", str(layout.panels / CARTOON_FILENAMES[0]),
            "--panel-c", str(layout.panels / CARTOON_FILENAMES[2]),
            "--panel-f", str(layout.spatial / "Figure3_Spatial_cFOS_occurrence.png"),
            "--panel-g", str(layout.spatial / "Figure3_Spatial_cFOS_NPY_occurrence.png"),
            "--figure-name", FIGURE_NAME,
            "--dpi", str(dpi),
            "--panel-a-dpi", str(panel_a_dpi),
            "--seed", "31",
        )),
    )
    forbidden_option = "--" + "force"
    if any(forbidden_option in step.argv for step in steps):
        raise SafetyError("A generated command requested destructive replacement behavior")
    stage = layout.stage_root
    for step in steps:
        executable_script = lexical_absolute(step.argv[1])
        if not is_within_or_equal(executable_script, layout.shadow_scripts):
            raise SafetyError(f"Command does not use the staged script copy: {step.name}")
        for index, token in enumerate(step.argv[:-1]):
            if token in {"--output-dir", "--analysis-dir", "--spatial-dir"}:
                candidate = lexical_absolute(step.argv[index + 1])
                if not is_within_or_equal(candidate, stage):
                    raise SafetyError(f"Command path escapes stage: {step.name}: {candidate}")
    return steps


def audit_ventricle_hil_inputs(hil_dir: Path) -> dict[str, object]:
    """Validate the four accepted HIL files before any stage is created."""
    require_no_symlink_components(hil_dir, description="Ventricle HIL directory")
    if hil_dir.is_symlink() or not hil_dir.is_dir():
        raise SafetyError(f"Ventricle HIL directory is missing or linked: {hil_dir}")
    records: list[dict[str, object]] = []
    accepted: dict[str, dict[str, object]] = {}
    for condition in ("Water", "Allulose"):
        slug = condition.lower()
        mask_path = hil_dir / f"{slug}_ventricle_mask.png"
        receipt_path = hil_dir / f"{slug}_ventricle_receipt.json"
        for role, path in (("mask", mask_path), ("receipt", receipt_path)):
            require_no_symlink_components(path, description=f"{condition} HIL {role}")
            if path.is_symlink() or not path.is_file():
                raise SafetyError(f"Missing or linked accepted {condition} HIL {role}: {path}")
            info = path.stat()
            if not stat.S_ISREG(info.st_mode) or info.st_size <= 0:
                raise SafetyError(f"Empty or non-regular accepted HIL file: {path}")
            records.append({
                "condition": condition,
                "role": role,
                "filename": path.name,
                "bytes": int(info.st_size),
                "sha256": sha256_file(path),
            })
        with mask_path.open("rb") as handle:
            if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                raise SafetyError(f"Accepted {condition} HIL mask is not a PNG: {mask_path}")
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise SafetyError(f"Cannot parse accepted {condition} HIL receipt: {receipt_path}") from exc
        expected = {
            "schema": VENTRICLE_HIL_SCHEMA,
            "status": "accepted",
            "condition": condition,
            "coordinate_frame_yx": [866, 1374],
            "mask_filename": mask_path.name,
        }
        observed = {key: receipt.get(key) for key in expected}
        if observed != expected:
            raise SafetyError(f"Accepted {condition} HIL identity mismatch: {observed}")
        if receipt.get("mask_file_sha256") != sha256_file(mask_path):
            raise SafetyError(f"Accepted {condition} HIL mask hash receipt mismatch")
        reviewer = str(receipt.get("reviewer", "")).strip()
        session = str(receipt.get("session", "")).strip()
        accepted_at = str(receipt.get("accepted_at_utc", "")).strip()
        try:
            accepted_time = datetime.fromisoformat(accepted_at.replace("Z", "+00:00"))
        except ValueError as exc:
            raise SafetyError(f"Accepted {condition} HIL timestamp is invalid") from exc
        if not reviewer or not session or accepted_time.tzinfo is None:
            raise SafetyError(f"Accepted {condition} HIL lacks reviewer/session/timezone")
        accepted[condition] = {
            "sample_id": receipt.get("sample_id"),
            "animal_id": receipt.get("animal_id"),
            "reviewer": reviewer,
            "session": session,
            "accepted_at_utc": accepted_at,
            "mask_binary_sha256": receipt.get("mask_geometry", {}).get("binary_sha256"),
        }
    return {
        "schema_version": 1,
        "audited_at_utc": utc_now(),
        "source_directory": str(hil_dir.resolve()),
        "accepted_conditions": accepted,
        "required_file_count": len(records),
        "records": records,
    }


def prepare_plan(args: argparse.Namespace) -> Plan:
    if args.dpi < 300 or args.dpi > 1200:
        raise SafetyError("--dpi must be between 300 and 1200")
    if args.panel_a_dpi < 300 or args.panel_a_dpi > 1200:
        raise SafetyError("--panel-a-dpi must be between 300 and 1200")
    (
        paper_root, stage_parent, stage_root, raw_project_root, raw_root,
        ventricle_hil_dir,
    ) = validate_plan_paths(args)
    source_scripts = validate_recovered_scripts(paper_root)
    raw_relpaths_from_analysis_script(source_scripts["scripts/Fig3/01_analyze_cfos_npy.py"])
    audit_ventricle_hil_inputs(ventricle_hil_dir)
    analysis_python = lexical_absolute(args.analysis_python)
    layout = build_layout(stage_root)
    steps = build_steps(
        layout, raw_project_root, raw_root, analysis_python, args.dpi, args.panel_a_dpi,
    )
    return Plan(
        paper_root=paper_root,
        stage_parent=stage_parent,
        raw_project_root=raw_project_root,
        raw_root=raw_root,
        ventricle_hil_dir=ventricle_hil_dir,
        analysis_python=analysis_python,
        layout=layout,
        source_scripts=source_scripts,
        steps=steps,
    )


def audit_independent_raw_inputs(plan: Plan, args: argparse.Namespace) -> dict[str, object]:
    require_no_symlink_components(plan.raw_project_root, description="Raw project root")
    require_no_symlink_components(plan.raw_root, description="Raw root")
    project = plan.raw_project_root.resolve(strict=True)
    raw_root = plan.raw_root.resolve(strict=True)
    if not project.is_dir() or not raw_root.is_dir():
        raise SafetyError("Independent raw project/root must be normal directories")
    if raw_root != project / "Fig3/raw/legacy_experiment_2025_07_28":
        raise SafetyError("Resolved raw root no longer matches its project-root binding")
    if paths_overlap(raw_root, plan.paper_root.resolve()):
        raise SafetyError("Resolved raw inputs overlap live Paper")
    if paths_overlap(raw_root, plan.layout.stage_root):
        raise SafetyError("Resolved raw inputs overlap the stage")
    validate_water3_bindings(args, raw_root)

    relpaths = raw_relpaths_from_analysis_script(
        plan.source_scripts["scripts/Fig3/01_analyze_cfos_npy.py"]
    )
    records: list[dict[str, object]] = []
    seen_inodes: dict[tuple[int, int], Path] = {}
    for relative in relpaths:
        path = raw_root / relative
        require_no_symlink_components(path, description=f"Raw companion {relative}")
        if path.is_symlink() or not path.is_file():
            raise SafetyError(f"Missing or linked mapped raw companion: {path}")
        resolved = path.resolve(strict=True)
        if not is_within_or_equal(resolved, raw_root):
            raise SafetyError(f"Raw companion escapes independent root: {path}")
        info = path.stat()
        if not stat.S_ISREG(info.st_mode) or info.st_size <= 0:
            raise SafetyError(f"Raw companion is not a nonempty regular file: {path}")
        if info.st_nlink != 1:
            raise SafetyError(f"Raw companion is not an independent copy (link count {info.st_nlink}): {path}")
        inode = (int(info.st_dev), int(info.st_ino))
        if inode in seen_inodes:
            raise SafetyError(f"Two mapped inputs share one inode: {seen_inodes[inode]} and {path}")
        seen_inodes[inode] = path
        records.append({
            "relative_path": relative,
            "bytes": int(info.st_size),
            "sha256": sha256_file(path),
            "device": int(info.st_dev),
            "inode": int(info.st_ino),
            "link_count": int(info.st_nlink),
        })
    return {
        "schema_version": 1,
        "audited_at_utc": utc_now(),
        "raw_root": str(raw_root),
        "mapped_file_count": len(records),
        "recursive_traversal_used": False,
        "regular_nonempty_unlinked_files": True,
        "records": records,
    }


def expected_outputs(layout: Layout) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    for name in ANALYSIS_FILENAMES:
        outputs[f"analysis/{name}"] = layout.analysis / name
    for name in CARTOON_FILENAMES:
        outputs[f"analysis/panels/{name}"] = layout.panels / name
    for name in PANEL_A_FILENAMES:
        outputs[f"analysis/panels/{name}"] = layout.panels / name
    for name in SPATIAL_FILENAMES:
        outputs[f"spatial/{name}"] = layout.spatial / name
    outputs[f"publication/{FIGURE_NAME}.pdf"] = layout.publication / f"{FIGURE_NAME}.pdf"
    outputs[f"publication/{FIGURE_NAME}.png"] = layout.publication / f"{FIGURE_NAME}.png"
    outputs[f"publication/{FIGURE_NAME}_caption.txt"] = layout.publication / f"{FIGURE_NAME}_caption.txt"
    for letter in "ABCDEFG":
        for suffix in ("pdf", "png"):
            name = f"{FIGURE_NAME}_Panel_{letter}.{suffix}"
            outputs[f"publication/panels/{name}"] = layout.publication / "panels" / name
    outputs[f"publication/legends/{FIGURE_NAME}_LEGEND.txt"] = (
        layout.publication / "legends" / f"{FIGURE_NAME}_LEGEND.txt"
    )
    for letter in "ABCDEFG":
        name = f"{FIGURE_NAME}_Panel_{letter}_LEGEND.txt"
        outputs[f"publication/legends/{name}"] = layout.publication / "legends" / name
    for name in SOURCE_DATA_FILENAMES:
        outputs[f"publication/source_data/{name}"] = layout.publication / "source_data" / name
    for name in (f"{FIGURE_NAME}_provenance.json", f"{FIGURE_NAME}_manifest.csv"):
        outputs[f"publication/provenance/{name}"] = layout.publication / "provenance" / name
    if len(outputs) != 67:
        raise SafetyError(f"Internal Figure 3 postflight inventory is not 67 files: {len(outputs)}")
    return outputs


def _validate_signature(path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix == ".json":
        with path.open("r", encoding="utf-8") as handle:
            json.load(handle)
    elif suffix == ".png":
        with path.open("rb") as handle:
            if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                raise SafetyError(f"PNG signature is invalid: {path}")
    elif suffix == ".pdf":
        with path.open("rb") as handle:
            if handle.read(5) != b"%PDF-":
                raise SafetyError(f"PDF signature is invalid: {path}")
    elif suffix == ".csv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            if next(csv.reader(handle), None) is None:
                raise SafetyError(f"CSV has no header row: {path}")


def validate_postflight(layout: Layout) -> dict[str, object]:
    expected = expected_outputs(layout)
    records: list[dict[str, object]] = []
    for role, path in expected.items():
        require_no_symlink_components(path, description=f"Staged output {role}")
        if path.is_symlink() or not path.is_file():
            raise SafetyError(f"Required staged Figure 3 output is missing or linked: {path}")
        resolved = path.resolve(strict=True)
        if not is_within_or_equal(resolved, layout.stage_root):
            raise SafetyError(f"Staged output escapes stage: {path}")
        info = path.stat()
        if not stat.S_ISREG(info.st_mode) or info.st_size <= 0:
            raise SafetyError(f"Required staged output is empty or non-regular: {path}")
        _validate_signature(path)
        records.append({
            "role": role,
            "path": str(resolved),
            "bytes": int(info.st_size),
            "sha256": sha256_file(path),
        })

    expected_paths = {path.resolve() for path in expected.values()}
    observed_paths: set[Path] = set()
    for root in (layout.analysis, layout.spatial, layout.publication):
        if not root.is_dir() or root.is_symlink():
            raise SafetyError(f"Required staged output directory is missing or linked: {root}")
        for candidate in root.rglob("*"):
            if candidate.is_symlink():
                raise SafetyError(f"Linked path found in staged output tree: {candidate}")
            if candidate.is_file():
                observed_paths.add(candidate.resolve())
    if observed_paths != expected_paths:
        missing = sorted(str(path) for path in expected_paths - observed_paths)
        unexpected = sorted(str(path) for path in observed_paths - expected_paths)
        raise SafetyError(f"Staged output inventory mismatch; missing={missing}, unexpected={unexpected}")

    panel_letters = sorted(
        path.stem.rsplit("_", 1)[-1]
        for path in (layout.publication / "panels").glob(f"{FIGURE_NAME}_Panel_*.pdf")
    )
    if panel_letters != list("ABCDEFG"):
        raise SafetyError(f"Panel PDF registry is not exactly A-G: {panel_letters}")
    return {
        "schema_version": 1,
        "validated_at_utc": utc_now(),
        "stage_root": str(layout.stage_root),
        "required_output_count": len(records),
        "panel_registry": panel_letters,
        "scientific_and_provenance_inventory_complete": True,
        "records": records,
    }


def plan_payload(plan: Plan, *, execute: bool) -> dict[str, object]:
    return {
        "schema_version": 1,
        "created_at_utc": utc_now(),
        "mode": "execute" if execute else "plan-only",
        "paper_root": str(plan.paper_root),
        "stage_parent": str(plan.stage_parent),
        "stage_root": str(plan.layout.stage_root),
        "raw_project_root": str(plan.raw_project_root),
        "raw_root": str(plan.raw_root),
        "raw_audit": "required before execution; deferred in plan-only mode",
        "ventricle_hil_source": str(plan.ventricle_hil_dir),
        "ventricle_hil_staged_copy": str(
            plan.layout.provenance / "accepted_ventricle_hil"
        ),
        "analysis_python": str(plan.analysis_python),
        "expected_postflight_output_count": len(expected_outputs(plan.layout)),
        "safety_contract": {
            "stage_root_absent_at_plan_time": True,
            "stage_root_direct_child_of_external_parent": True,
            "live_paper_writes": False,
            "existing_path_reuse": False,
            "launcher_cleanup_or_replacement": False,
            "failed_stage_preserved": True,
            "recursive_raw_traversal": False,
            "accepted_hil_source_read_only": True,
            "cartoon_renderer_uses_staged_hil_copy": True,
        },
        "recovered_script_sha256": dict(RECOVERED_SCRIPT_SHA256),
        "commands": [{"name": step.name, "argv": list(step.argv)} for step in plan.steps],
    }


def copy_recovered_scripts(plan: Plan) -> None:
    plan.layout.shadow_scripts.mkdir(parents=True, exist_ok=False)
    for relative, source in plan.source_scripts.items():
        target = plan.layout.shadow_paper / relative
        if target.exists() or target.is_symlink():
            raise SafetyError(f"Staged script target unexpectedly exists: {target}")
        shutil.copyfile(source, target)
        if sha256_file(target) != RECOVERED_SCRIPT_SHA256[relative]:
            raise SafetyError(f"Staged script copy failed its exact hash gate: {target}")
        target.chmod(0o444)


def copy_ventricle_hil_inputs(plan: Plan, audit: Mapping[str, object]) -> Path:
    """Copy only the four validated accepted inputs into the fresh stage."""
    target_dir = plan.layout.provenance / "accepted_ventricle_hil"
    target_dir.mkdir(mode=0o750, parents=False, exist_ok=False)
    expected_hashes = {
        str(record["filename"]): str(record["sha256"])
        for record in audit["records"]  # type: ignore[index]
    }
    for filename in VENTRICLE_HIL_FILENAMES:
        source = plan.ventricle_hil_dir / filename
        target = target_dir / filename
        if target.exists() or target.is_symlink():
            raise SafetyError(f"Staged HIL target unexpectedly exists: {target}")
        shutil.copyfile(source, target)
        if sha256_file(target) != expected_hashes[filename]:
            raise SafetyError(f"Staged HIL copy failed its exact hash gate: {target}")
        target.chmod(0o444)
    return target_dir


def assert_absent(path: Path, description: str) -> None:
    if path.exists() or path.is_symlink():
        raise SafetyError(f"{description} must be absent and will not be replaced: {path}")


def run_step(step: Step, plan: Plan, env: Mapping[str, str]) -> None:
    log_path = plan.layout.logs / f"{step.name}.log"
    with log_path.open("xb") as log_handle:
        subprocess.run(
            step.argv,
            cwd=plan.layout.shadow_paper,
            env=dict(env),
            stdin=subprocess.DEVNULL,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def execute_plan(plan: Plan, args: argparse.Namespace) -> dict[str, object]:
    if not plan.analysis_python.is_file() or not os.access(plan.analysis_python, os.X_OK):
        raise SafetyError(f"Analysis Python is not an executable file: {plan.analysis_python}")
    raw_receipt = audit_independent_raw_inputs(plan, args)
    hil_receipt = audit_ventricle_hil_inputs(plan.ventricle_hil_dir)
    if plan.layout.stage_root.exists() or plan.layout.stage_root.is_symlink():
        raise SafetyError(f"Stage root appeared after planning: {plan.layout.stage_root}")
    if not os.access(plan.stage_parent, os.W_OK | os.X_OK):
        raise SafetyError(f"Stage parent is not writable/searchable: {plan.stage_parent}")

    plan.layout.stage_root.mkdir(mode=0o750, parents=False, exist_ok=False)
    try:
        plan.layout.logs.mkdir(exist_ok=False)
        plan.layout.runtime.mkdir(exist_ok=False)
        plan.layout.provenance.mkdir(exist_ok=False)
        for directory in (
            plan.layout.runtime / "tmp",
            plan.layout.runtime / "matplotlib",
            plan.layout.runtime / "xdg_cache",
            plan.layout.runtime / "numba_cache",
        ):
            directory.mkdir(exist_ok=False)
        exclusive_json(plan.layout.provenance / "STAGE_PLAN.json", plan_payload(plan, execute=True))
        exclusive_json(plan.layout.provenance / "RAW_INPUT_AUDIT.json", raw_receipt)
        exclusive_json(plan.layout.provenance / "VENTRICLE_HIL_INPUT_AUDIT.json", hil_receipt)
        copy_ventricle_hil_inputs(plan, hil_receipt)
        copy_recovered_scripts(plan)

        env = os.environ.copy()
        env.update({
            "PYTHONDONTWRITEBYTECODE": "1",
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": str(plan.layout.runtime / "matplotlib"),
            "XDG_CACHE_HOME": str(plan.layout.runtime / "xdg_cache"),
            "NUMBA_CACHE_DIR": str(plan.layout.runtime / "numba_cache"),
            "TMPDIR": str(plan.layout.runtime / "tmp"),
        })

        assert_absent(plan.layout.analysis, "Analysis output")
        run_step(plan.steps[0], plan, env)
        for filename in ANALYSIS_FILENAMES:
            if not (plan.layout.analysis / filename).is_file():
                raise SafetyError(f"Analysis did not produce {filename}")

        assert_absent(plan.layout.panels, "Cartoon panel output")
        run_step(plan.steps[1], plan, env)
        for filename in CARTOON_FILENAMES:
            if not (plan.layout.panels / filename).is_file():
                raise SafetyError(f"Cartoon renderer did not produce {filename}")

        assert_absent(plan.layout.spatial, "Spatial output")
        run_step(plan.steps[2], plan, env)
        for filename in SPATIAL_FILENAMES:
            if not (plan.layout.spatial / filename).is_file():
                raise SafetyError(f"Spatial analysis did not produce {filename}")

        assert_absent(plan.layout.publication, "Publication output")
        assert_absent(
            plan.layout.panels / "Figure3_Panel_A_reference_microscopy.png",
            "Legacy Panel A path",
        )
        run_step(plan.steps[3], plan, env)
        postflight = validate_postflight(plan.layout)
        exclusive_json(plan.layout.provenance / "STAGE_POSTFLIGHT.json", postflight)
        exclusive_text(
            plan.layout.provenance / "STAGE_COMPLETE.txt",
            "Figure 3 staged workflow completed and passed the 67-file postflight.\n",
        )
        return postflight
    except BaseException as exc:
        failure_path = plan.layout.provenance / "STAGE_FAILED.json"
        if plan.layout.provenance.is_dir() and not failure_path.exists():
            exclusive_json(failure_path, {
                "schema_version": 1,
                "failed_at_utc": utc_now(),
                "exception_type": type(exc).__name__,
                "message": str(exc),
                "stage_preserved": True,
            })
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, required=True,
                        help="Live recovered Paper root used only as the hash-pinned script source.")
    parser.add_argument("--stage-parent", type=Path, required=True,
                        help="Existing external directory that will contain one new stage.")
    parser.add_argument("--stage-root", type=Path, required=True,
                        help="Absent direct child of --stage-parent; never reused or cleaned.")
    parser.add_argument("--raw-project-root", type=Path, required=True,
                        help="Independent project tree containing Fig3/raw/legacy_experiment_2025_07_28.")
    parser.add_argument("--raw-root", type=Path, required=True,
                        help="Independent combined Figure 3 raw root at the required project-relative path.")
    parser.add_argument(
        "--ventricle-hil-dir", type=Path, required=True,
        help="Directory containing both accepted, source-bound Water/Allulose ventricle pairs.",
    )
    parser.add_argument("--water3-dapi-tiff", type=Path, required=True)
    parser.add_argument("--water3-cfos-tiff", type=Path, required=True)
    parser.add_argument("--water3-npy-tiff", type=Path, required=True)
    parser.add_argument("--water3-dapi-mask", type=Path, required=True)
    parser.add_argument("--water3-cfos-mask", type=Path, required=True)
    parser.add_argument("--water3-npy-mask", type=Path, required=True)
    parser.add_argument("--water3-extraction-provenance", type=Path, required=True)
    parser.add_argument("--analysis-python", type=Path, default=Path(sys.executable),
                        help="Python executable with the recovered Figure 3 environment.")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--panel-a-dpi", type=int, default=600)
    parser.add_argument("--execute", action="store_true",
                        help="Opt in to a real staged run; omission prints a read-only JSON plan.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        plan = prepare_plan(args)
        if not args.execute:
            print(json.dumps(plan_payload(plan, execute=False), indent=2, sort_keys=True))
            return 0
        result = execute_plan(plan, args)
        print(json.dumps({
            "status": "PASS",
            "stage_root": str(plan.layout.stage_root),
            "required_output_count": result["required_output_count"],
            "postflight_receipt": str(plan.layout.provenance / "STAGE_POSTFLIGHT.json"),
        }, indent=2, sort_keys=True))
        return 0
    except (SafetyError, FileNotFoundError, NotADirectoryError) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
