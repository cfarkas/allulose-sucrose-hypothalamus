#!/usr/bin/env python3
"""Quantify whole-field c-FOS, NPY, and NeuN phenotypes for Figure 3.

The nine biological animals are mapped explicitly to their raw-local Cellpose
``*_seg.npy`` files. NeuN uses the eight animals with valid curated masks;
Water3 remains missing because its display-only C2 channel failed. The
statistical unit is the animal; no pixel, ROI, or nucleus is an independent
replicate.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import shutil
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, mannwhitneyu, rankdata
from skimage.measure import regionprops
import tifffile


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

RAW_ROOT = PAPER_ROOT / "Fig3" / "raw" / "legacy_experiment_2025_07_28"
DEFAULT_OUTPUT = PAPER_ROOT / "analyses" / "Fig3" / "results" / "final_run"

CONDITIONS = ("Water", "Sucrose", "Allulose")
ENDPOINTS: tuple[tuple[str, str], ...] = (
    ("cfos_over_dapi", "c-FOS-positive nuclei / DAPI nuclei"),
    ("npy_over_dapi", "NPY-positive nuclei / DAPI nuclei"),
    ("cfos_npy_over_npy", "c-FOS and NPY double-positive nuclei / NPY-positive nuclei"),
    ("double_over_dapi", "c-FOS and NPY double-positive nuclei / DAPI nuclei"),
    ("cfos_neun_over_neun", "c-FOS and NeuN double-positive nuclei / NeuN-positive nuclei"),
)
PAIR_ORDER = (
    ("Water", "Sucrose"),
    ("Water", "Allulose"),
    ("Sucrose", "Allulose"),
)

CFOS_MIN_ROI_AREA_PX = 5
CFOS_DAPI_OVERLAP_FRACTION = 0.40
NPY_MIN_OVERLAP_PX = 1
NPY_DAPI_OVERLAP_FRACTION = 0.05
NEUN_MIN_OVERLAP_PX = 1
NEUN_DAPI_OVERLAP_FRACTION = 0.05
IQR_MULTIPLIER = 1.5
OUTPUT_NAMES = (
    "per_animal_cfos_npy.csv",
    "anova_mwu_statistics.csv",
    "statistical_test_receipts.csv",
    "segmentation_qc.csv",
    "source_manifest.csv",
    "README.txt",
)
OBSOLETE_OUTPUT_NAMES = (
    "exact_statistics.csv",
    "exact_test_receipts.csv",
)
PRESERVED_OUTPUT_NAMES = ("panels",)


@dataclass(frozen=True)
class SampleSpec:
    animal_id: str
    sample_id: str
    condition: str
    dapi_seg: str
    cfos_seg: str
    npy_seg: str
    dapi_tif: str
    cfos_tif: str
    npy_tif: str
    acquisition_batch: str
    animal_id_status: str = "confirmed"
    neun_seg: str | None = None
    neun_tif: str | None = None


SAMPLES: tuple[SampleSpec, ...] = (
    SampleSpec(
        "WATER_NPY1", "water-npy1-0001", "Water",
        "water/WATER_NPY1-0001_dapi_seg.npy",
        "water/WATER_NPY1-0001_fos_seg.npy",
        "water/WATER_NPY1-0001_GFP_seg.npy",
        "water/WATER_NPY1-0001_dapi.tif",
        "water/WATER_NPY1-0001_fos.tif",
        "water/WATER_NPY1-0001_GFP.tif",
        "legacy_water_2025",
        neun_seg="water/WATER_NPY1-0001_NeuN_seg.npy",
        neun_tif="water/WATER_NPY1-0001_NeuN.tif",
    ),
    SampleSpec(
        "WATER_NPY2", "water-npy2-0001", "Water",
        "water/WATER_NPY2-0001_dapi_seg.npy",
        "water/WATER_NPY2-0001_fos_seg.npy",
        "water/WATER_NPY2-0001_GFP_seg.npy",
        "water/WATER_NPY2-0001_dapi.tif",
        "water/WATER_NPY2-0001_fos.tif",
        "water/WATER_NPY2-0001_GFP.tif",
        "legacy_water_2025",
        neun_seg="water/WATER_NPY2-0001_NeuN_seg.npy",
        neun_tif="water/WATER_NPY2-0001_NeuN.tif",
    ),
    SampleSpec(
        "WATER_NPY3", "water-npy3-s01", "Water",
        "water/WATER_NPY3_S01_dapi_seg.npy",
        "water/WATER_NPY3_S01_fos_seg.npy",
        "water/WATER_NPY3_S01_GFP_seg.npy",
        "water/WATER_NPY3_S01_dapi.tif",
        "water/WATER_NPY3_S01_fos.tif",
        "water/WATER_NPY3_S01_GFP.tif",
        "leica_sp8_40x_2026_21_plane_zmax_rgb16",
        "provisional; replace with true animal ID before final inference",
        neun_tif="water/WATER_NPY3_S01_NeuN_display_only.tif",
    ),
    SampleSpec(
        "E7_FR7-1", "e7-fr7-1", "Sucrose",
        "E7-FR7-1_dapi_seg.npy", "E7-FR7-1_fos_seg.npy", "E7-FR7-1_gfp_seg.npy",
        "E7-FR7-1_dapi.tif", "E7-FR7-1_fos.tif", "E7-FR7-1_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E7-FR7-1_NeuN_seg.npy", neun_tif="E7-FR7-1_NeuN.tif",
    ),
    SampleSpec(
        "E7_FR7-5", "e7-fr7-5", "Sucrose",
        "E7-FR7-5_dapi_seg.npy", "E7-FR7-5_fos_seg.npy", "E7-FR7-5_gfp_seg.npy",
        "E7-FR7-5_dapi.tif", "E7-FR7-5_fos.tif", "E7-FR7-5_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E7-FR7-5_NeuN_seg.npy", neun_tif="E7-FR7-5_NeuN.tif",
    ),
    SampleSpec(
        "E9_FR8-1", "e9-fr8-1", "Sucrose",
        "E9-FR8-1_dapi_seg.npy", "E9-FR8-1_fos_seg.npy", "E9-FR8-1_gfp_seg.npy",
        "E9-FR8-1_dapi.tif", "E9-FR8-1_fos.tif", "E9-FR8-1_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E9-FR8-1_NeuN_seg.npy", neun_tif="E9-FR8-1_NeuN.tif",
    ),
    SampleSpec(
        "E8_FR6-1", "e8-fr6-1-alulosa-25x", "Allulose",
        "E8 FR6-1_ALULOSA_25X_dapi_seg.npy",
        "E8 FR6-1_ALULOSA_25X_fos_seg.npy",
        "E8 FR6-1_ALULOSA_25X_gfp_seg.npy",
        "E8 FR6-1_ALULOSA_25X_dapi.tif",
        "E8 FR6-1_ALULOSA_25X_fos.tif",
        "E8 FR6-1_ALULOSA_25X_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E8 FR6-1_ALULOSA_25X_NeuN_seg.npy",
        neun_tif="E8 FR6-1_ALULOSA_25X_NeuN.tif",
    ),
    SampleSpec(
        "E8_FR6-3", "e8-fr6-3-alulosa", "Allulose",
        "E8-FR6-3_ALULOSA_dapi_seg.npy",
        "E8-FR6-3_ALULOSA_fos_seg.npy",
        "E8-FR6-3_ALULOSA_gfp_seg.npy",
        "E8-FR6-3_ALULOSA_dapi.tif",
        "E8-FR6-3_ALULOSA_fos.tif",
        "E8-FR6-3_ALULOSA_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E8-FR6-3_ALULOSA_NeuN_seg.npy",
        neun_tif="E8-FR6-3_ALULOSA_NeuN.tif",
    ),
    SampleSpec(
        "E8_FR6-4", "e8-fr6-4-alulosa", "Allulose",
        "E8-FR6-4_ALULOSA_dapi_seg.npy",
        "E8-FR6-4_ALULOSA_fos_seg.npy",
        "E8-FR6-4_ALULOSA_gfp_seg.npy",
        "E8-FR6-4_ALULOSA_dapi.tif",
        "E8-FR6-4_ALULOSA_fos.tif",
        "E8-FR6-4_ALULOSA_gfp.tif",
        "legacy_apotome_25x_2025_rgb8",
        neun_seg="E8-FR6-4_ALULOSA_NeuN_seg.npy",
        neun_tif="E8-FR6-4_ALULOSA_NeuN.tif",
    ),
)


# Audited from the retired v10.1 result on 2026-08-21.  These constants are a
# regression oracle only; no retired/trash file is opened by this script.
ARCHIVED_FIRST_EIGHT_COUNTS: Mapping[str, Mapping[str, int]] = {
    "WATER_NPY1": {"dapi_nuclei": 7172, "cfos_nuclei": 615, "npy_nuclei": 272, "cfos_npy_nuclei": 92, "neun_nuclei": 1576, "cfos_neun_nuclei": 256},
    "WATER_NPY2": {"dapi_nuclei": 11608, "cfos_nuclei": 773, "npy_nuclei": 208, "cfos_npy_nuclei": 57, "neun_nuclei": 5876, "cfos_neun_nuclei": 657},
    "E7_FR7-1": {"dapi_nuclei": 3061, "cfos_nuclei": 43, "npy_nuclei": 403, "cfos_npy_nuclei": 12, "neun_nuclei": 1077, "cfos_neun_nuclei": 18},
    "E7_FR7-5": {"dapi_nuclei": 3241, "cfos_nuclei": 81, "npy_nuclei": 302, "cfos_npy_nuclei": 19, "neun_nuclei": 1464, "cfos_neun_nuclei": 56},
    "E9_FR8-1": {"dapi_nuclei": 2316, "cfos_nuclei": 33, "npy_nuclei": 301, "cfos_npy_nuclei": 5, "neun_nuclei": 1036, "cfos_neun_nuclei": 25},
    "E8_FR6-1": {"dapi_nuclei": 1976, "cfos_nuclei": 112, "npy_nuclei": 268, "cfos_npy_nuclei": 69, "neun_nuclei": 626, "cfos_neun_nuclei": 44},
    "E8_FR6-3": {"dapi_nuclei": 2976, "cfos_nuclei": 165, "npy_nuclei": 369, "cfos_npy_nuclei": 88, "neun_nuclei": 1062, "cfos_neun_nuclei": 79},
    "E8_FR6-4": {"dapi_nuclei": 3372, "cfos_nuclei": 305, "npy_nuclei": 494, "cfos_npy_nuclei": 131, "neun_nuclei": 1196, "cfos_neun_nuclei": 142},
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--force", action="store_true",
        help="Replace only this script's known files and remove its obsolete exact_* outputs.",
    )
    return parser


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_to_paper(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(PAPER_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def validate_inventory() -> None:
    if len(SAMPLES) != 9:
        raise RuntimeError(f"Figure 3 requires exactly 9 animals; mapped {len(SAMPLES)}")
    animals = [sample.animal_id for sample in SAMPLES]
    if len(set(animals)) != len(animals):
        raise RuntimeError("Duplicate animal IDs in explicit Figure 3 mapping")
    counts = {condition: sum(sample.condition == condition for sample in SAMPLES) for condition in CONDITIONS}
    if counts != {condition: 3 for condition in CONDITIONS}:
        raise RuntimeError(f"Expected 3 animals per condition, observed {counts}")


def prepare_output_dir(output: Path, force: bool) -> Path:
    requested = output.expanduser()
    if requested.exists() and requested.is_symlink():
        raise RuntimeError(f"Refusing output through symbolic link: {requested}")
    resolved = requested.resolve()
    allowed_parent = (PAPER_ROOT / "analyses" / "Fig3" / "results").resolve()
    if resolved != allowed_parent and allowed_parent not in resolved.parents:
        raise RuntimeError(f"Output must remain under {allowed_parent}: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    existing = [resolved / name for name in OUTPUT_NAMES if (resolved / name).exists()]
    obsolete = [resolved / name for name in OBSOLETE_OUTPUT_NAMES if (resolved / name).exists()]
    allowed_names = set(OUTPUT_NAMES) | set(PRESERVED_OUTPUT_NAMES) | set(OBSOLETE_OUTPUT_NAMES)
    unknown = [path for path in resolved.iterdir() if path.name not in allowed_names]
    if unknown:
        raise RuntimeError(
            "Refusing mixed Figure 3 output because unknown files are present: "
            + ", ".join(path.name for path in unknown)
        )
    # Known analysis files are atomically replaced on every run so the normal
    # publication build is deterministic.  The renderer-owned panels/ tree is
    # an explicit companion and is never deleted, including under --force.
    if force:
        for path in existing + obsolete:
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()
    elif obsolete:
        raise RuntimeError(
            "Obsolete exact_* statistics outputs are present; rerun with --force to remove them: "
            + ", ".join(path.name for path in obsolete)
        )
    return resolved


def load_cellpose_mask(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing or empty Cellpose file: {path}")
    payload = np.load(path, allow_pickle=True)
    try:
        record = payload.item()
    except Exception as exc:  # pragma: no cover - malformed external input
        raise ValueError(f"Cellpose file is not a scalar dictionary: {path}") from exc
    if not isinstance(record, dict) or "masks" not in record:
        raise ValueError(f"Cellpose dictionary lacks 'masks': {path}")
    masks = np.asarray(record["masks"])
    if masks.ndim != 2:
        raise ValueError(f"Cellpose masks must be 2D, got {masks.shape}: {path}")
    if not np.issubdtype(masks.dtype, np.integer):
        raise ValueError(f"Cellpose masks must use an integer dtype, got {masks.dtype}: {path}")
    if np.any(masks < 0):
        raise ValueError(f"Cellpose masks contain negative labels: {path}")
    positive_labels = np.unique(masks[masks > 0])
    contiguous = bool(
        positive_labels.size == 0
        or np.array_equal(positive_labels, np.arange(1, int(positive_labels[-1]) + 1))
    )
    ismanual = record.get("ismanual")
    manual_array = np.asarray(ismanual) if ismanual is not None else np.asarray([], dtype=bool)
    manual_changes = record.get("manual_changes")
    change_count = len(manual_changes) if hasattr(manual_changes, "__len__") else 0
    metadata = {
        "mask_shape_y": int(masks.shape[0]),
        "mask_shape_x": int(masks.shape[1]),
        "mask_dtype": str(masks.dtype),
        "object_count": int(positive_labels.size),
        "max_label": int(positive_labels[-1]) if positive_labels.size else 0,
        "labels_contiguous": contiguous,
        "manual_label_count": int(np.count_nonzero(manual_array)) if manual_array.size else 0,
        "manual_flag_count": int(manual_array.size),
        "manual_change_record_count": int(change_count),
        "model_path": str(record.get("model_path", "")),
        "flow_threshold": record.get("flow_threshold", np.nan),
        "cellprob_threshold": record.get("cellprob_threshold", np.nan),
        "diameter": record.get("diameter", np.nan),
        "cellpose_source_filename": str(record.get("filename", "")),
        "payload_keys": ";".join(sorted(map(str, record))),
    }
    return masks, metadata


def tiff_metadata(path: Path, expected_yx: tuple[int, int]) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing or empty TIFF companion: {path}")
    with tifffile.TiffFile(path) as handle:
        series = handle.series[0]
        shape = tuple(map(int, series.shape))
        dtype = str(series.dtype)
    shape_matches = shape == expected_yx or (len(shape) >= 3 and shape[:2] == expected_yx)
    if not shape_matches:
        raise ValueError(f"TIFF/mask shape mismatch for {path}: TIFF={shape}, mask={expected_yx}")
    return {"tiff_shape": "x".join(map(str, shape)), "tiff_dtype": dtype}


def cfos_positive_nuclei(dapi: np.ndarray, cfos: np.ndarray) -> tuple[set[int], int]:
    accepted_mappings: list[int] = []
    max_nucleus = int(dapi.max())
    if max_nucleus <= 0:
        raise ValueError("DAPI segmentation contains no nuclei")
    for prop in regionprops(cfos):
        if int(prop.area) < CFOS_MIN_ROI_AREA_PX:
            continue
        values = dapi[prop.coords[:, 0], prop.coords[:, 1]]
        counts = np.bincount(values.astype(np.int64), minlength=max_nucleus + 1)
        counts[0] = 0
        best = int(np.argmax(counts)) if counts.sum() else 0
        fraction = float(counts[best]) / float(prop.area) if best else 0.0
        if best > 0 and fraction >= CFOS_DAPI_OVERLAP_FRACTION:
            accepted_mappings.append(best)
    return set(accepted_mappings), len(accepted_mappings)


def marker_positive_nuclei(
    dapi: np.ndarray,
    marker_labels: np.ndarray,
    *,
    min_overlap_px: int,
    min_fraction: float,
) -> set[int]:
    marker = np.asarray(marker_labels > 0, dtype=bool)
    positives: set[int] = set()
    for prop in regionprops(dapi):
        coords = prop.coords
        overlap = int(marker[coords[:, 0], coords[:, 1]].sum())
        if (
            overlap >= int(min_overlap_px)
            and overlap / float(max(int(prop.area), 1)) >= float(min_fraction)
        ):
            positives.add(int(prop.label))
    return positives


def safe_ratio(numerator: int, denominator: int, label: str) -> float:
    if denominator <= 0:
        raise ValueError(f"Required denominator is zero for {label}")
    return float(numerator / denominator)


def quantify_sample(raw_root: Path, sample: SampleSpec) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    channel_specs: list[tuple[str, str, str]] = [
        ("DAPI", sample.dapi_seg, sample.dapi_tif),
        ("c-FOS", sample.cfos_seg, sample.cfos_tif),
        ("NPY", sample.npy_seg, sample.npy_tif),
    ]
    if sample.neun_seg is not None and sample.neun_tif is not None:
        channel_specs.append(("NeuN", sample.neun_seg, sample.neun_tif))
    labels: dict[str, np.ndarray] = {}
    qc_rows: list[dict[str, Any]] = []
    for channel, seg_rel, tif_rel in channel_specs:
        seg_path = raw_root / seg_rel
        tif_path = raw_root / tif_rel
        mask, metadata = load_cellpose_mask(seg_path)
        labels[channel] = mask
        tiff_info = tiff_metadata(tif_path, tuple(map(int, mask.shape)))
        qc_rows.append({
            "animal_id": sample.animal_id,
            "sample_id": sample.sample_id,
            "condition": sample.condition,
            "channel": channel,
            "segmentation_path": relative_to_paper(seg_path),
            "segmentation_sha256": sha256_file(seg_path),
            "segmentation_size_bytes": int(seg_path.stat().st_size),
            "tiff_path": relative_to_paper(tif_path),
            "tiff_sha256": sha256_file(tif_path),
            "tiff_size_bytes": int(tif_path.stat().st_size),
            "within_animal_shape_match": True,
            "segmentation_available": True,
            "numeric_analysis_eligible": True,
            "failure_reason": "",
            **metadata,
            **tiff_info,
        })

    shapes = {channel: tuple(mask.shape) for channel, mask in labels.items()}
    if len(set(shapes.values())) != 1:
        raise ValueError(f"Segmentation shape mismatch for {sample.animal_id}: {shapes}")
    for row in qc_rows:
        row["within_animal_shape_match"] = True

    dapi = labels["DAPI"]
    cfos = labels["c-FOS"]
    npy = labels["NPY"]
    dapi_labels = np.unique(dapi[dapi > 0])
    cfos_nuclei, accepted_cfos_rois = cfos_positive_nuclei(dapi, cfos)
    npy_nuclei = marker_positive_nuclei(
        dapi, npy,
        min_overlap_px=NPY_MIN_OVERLAP_PX,
        min_fraction=NPY_DAPI_OVERLAP_FRACTION,
    )
    double_nuclei = cfos_nuclei & npy_nuclei
    neun_available = "NeuN" in labels
    if neun_available:
        neun_nuclei: set[int] | None = marker_positive_nuclei(
            dapi, labels["NeuN"],
            min_overlap_px=NEUN_MIN_OVERLAP_PX,
            min_fraction=NEUN_DAPI_OVERLAP_FRACTION,
        )
        cfos_neun_nuclei: set[int] | None = cfos_nuclei & neun_nuclei
    else:
        neun_nuclei = None
        cfos_neun_nuclei = None

    if sample.neun_tif is None:
        raise RuntimeError(f"{sample.animal_id} has no NeuN TIFF provenance mapping")
    neun_tif_path = raw_root / sample.neun_tif
    if neun_available:
        neun_display_sha256 = ""
    else:
        tiff_metadata(neun_tif_path, tuple(map(int, dapi.shape)))
        neun_display_sha256 = sha256_file(neun_tif_path)

    counts = {
        "dapi_nuclei": int(dapi_labels.size),
        "cfos_nuclei": int(len(cfos_nuclei)),
        "npy_nuclei": int(len(npy_nuclei)),
        "cfos_npy_nuclei": int(len(double_nuclei)),
        "neun_nuclei": int(len(neun_nuclei)) if neun_nuclei is not None else np.nan,
        "cfos_neun_nuclei": int(len(cfos_neun_nuclei)) if cfos_neun_nuclei is not None else np.nan,
    }
    row: dict[str, Any] = {
        "animal_id": sample.animal_id,
        "sample_id": sample.sample_id,
        "condition": sample.condition,
        "condition_order": CONDITIONS.index(sample.condition) + 1,
        "experimental_unit": "biological animal",
        "field_count": 1,
        "image_height_px": int(dapi.shape[0]),
        "image_width_px": int(dapi.shape[1]),
        "dapi_mask_objects": int(qc_rows[0]["object_count"]),
        "cfos_mask_objects": int(qc_rows[1]["object_count"]),
        "npy_mask_objects": int(qc_rows[2]["object_count"]),
        "neun_mask_objects": int(qc_rows[3]["object_count"]) if neun_available else np.nan,
        "accepted_cfos_roi_objects": int(accepted_cfos_rois),
        "duplicate_cfos_roi_to_nucleus_mappings": int(accepted_cfos_rois - len(cfos_nuclei)),
        **counts,
        "cfos_over_dapi": safe_ratio(len(cfos_nuclei), len(dapi_labels), "c-FOS/DAPI"),
        "npy_over_dapi": safe_ratio(len(npy_nuclei), len(dapi_labels), "NPY/DAPI"),
        "cfos_npy_over_npy": safe_ratio(len(double_nuclei), len(npy_nuclei), "c-FOS and NPY/NPY"),
        "double_over_dapi": safe_ratio(len(double_nuclei), len(dapi_labels), "c-FOS and NPY/DAPI"),
        "cfos_neun_over_neun": (
            safe_ratio(len(cfos_neun_nuclei), len(neun_nuclei), "c-FOS and NeuN/NeuN")
            if neun_nuclei is not None and cfos_neun_nuclei is not None
            else np.nan
        ),
        "cfos_min_roi_area_px": CFOS_MIN_ROI_AREA_PX,
        "cfos_dapi_min_overlap_fraction": CFOS_DAPI_OVERLAP_FRACTION,
        "npy_min_overlap_px": NPY_MIN_OVERLAP_PX,
        "npy_dapi_min_overlap_fraction": NPY_DAPI_OVERLAP_FRACTION,
        "neun_min_overlap_px": NEUN_MIN_OVERLAP_PX,
        "neun_dapi_min_overlap_fraction": NEUN_DAPI_OVERLAP_FRACTION,
        "acquisition_batch": sample.acquisition_batch,
        "animal_id_status": sample.animal_id_status,
        "neun_available": bool(neun_available),
        "neun_included": bool(neun_available),
        "neun_analysis_set": "NeuN_complete_case" if neun_available else "excluded_failed_channel",
        "neun_exclusion_reason": "" if neun_available else "Water3 C2 failed; display-only MIP is not quantifiable",
        "neun_display_only_tiff": relative_to_paper(neun_tif_path) if not neun_available else "",
        "neun_display_only_sha256": neun_display_sha256,
        "batch_caveat": (
            "Water3 is a 2026 Leica SP8 40X RGB16 21-plane z-maximum projection; "
            "legacy animals are 2025 RGB8 exports from other acquisitions. "
            "Object-count ratios reduce field-size effects but do not remove acquisition/model confounding."
        ),
    }

    expected = ARCHIVED_FIRST_EIGHT_COUNTS.get(sample.animal_id)
    if expected is None:
        row["archived_reproduction_status"] = "not_applicable_new_water3"
        for key in counts:
            row[f"archived_expected_{key}"] = np.nan
            row[f"archived_delta_{key}"] = np.nan
    else:
        matches = []
        for key, value in counts.items():
            expected_value = int(expected[key])
            row[f"archived_expected_{key}"] = expected_value
            row[f"archived_delta_{key}"] = int(value - expected_value)
            matches.append(value == expected_value)
        row["archived_reproduction_status"] = "exact_match" if all(matches) else "mismatch"
    return row, qc_rows


def iqr_flags(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return np.zeros(values.size, dtype=bool)
    q1, q3 = np.percentile(values, [25.0, 75.0])
    spread = q3 - q1
    return (values < q1 - IQR_MULTIPLIER * spread) | (values > q3 + IQR_MULTIPLIER * spread)


def add_iqr_flags(animal_values: pd.DataFrame) -> pd.DataFrame:
    result = animal_values.copy()
    for endpoint, _label in ENDPOINTS:
        flag_column = f"{endpoint}_iqr_outlier"
        result[flag_column] = False
        for condition in CONDITIONS:
            selected = result["condition"].eq(condition) & result[endpoint].notna()
            values = result.loc[selected, endpoint].to_numpy(float)
            result.loc[selected, flag_column] = iqr_flags(values)
    return result


def one_way_anova(groups: Sequence[np.ndarray]) -> dict[str, Any]:
    arrays = [np.asarray(values, dtype=float) for values in groups]
    if any(values.size < 2 or not np.isfinite(values).all() for values in arrays):
        raise ValueError("One-way ANOVA requires at least two finite animals per condition")
    result = f_oneway(*arrays)
    total = int(sum(values.size for values in arrays))
    return {
        "statistic": float(result.statistic),
        "p_value_raw": float(result.pvalue),
        "df_between": int(len(arrays) - 1),
        "df_within": int(total - len(arrays)),
        "extreme_labelings": np.nan,
        "enumerated_labelings": np.nan,
        "permutation_exhaustive": False,
    }


def exact_mwu(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
    first = np.asarray(left, dtype=float)
    second = np.asarray(right, dtype=float)
    if first.size < 2 or second.size < 2 or not np.isfinite(np.concatenate((first, second))).all():
        raise ValueError("Exact MWU requires at least two finite animals in both conditions")
    pooled = np.concatenate((first, second))
    ranks = rankdata(pooled, method="average")
    n_first = int(first.size)
    expected_rank_sum = n_first * (pooled.size + 1.0) / 2.0
    observed_distance = abs(float(np.sum(ranks[:n_first])) - expected_rank_sum)
    total_labelings = math.comb(int(pooled.size), n_first)
    extreme = sum(
        abs(float(np.sum(ranks[list(selection)])) - expected_rank_sum)
        >= observed_distance - 1e-12
        for selection in combinations(range(int(pooled.size)), n_first)
    )
    return {
        "statistic": float(mannwhitneyu(first, second, alternative="two-sided").statistic),
        "p_value_raw": float(extreme / total_labelings),
        "extreme_labelings": int(extreme),
        "enumerated_labelings": int(total_labelings),
        "permutation_exhaustive": True,
    }


def cliffs_delta_right_minus_left(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    differences = right[:, None] - left[None, :]
    return float((np.count_nonzero(differences > 0) - np.count_nonzero(differences < 0)) / differences.size)


def flagged_animals(
    data: pd.DataFrame,
    endpoint: str,
    conditions: Iterable[str],
) -> tuple[int, str]:
    selected = data["condition"].isin(tuple(conditions)) & data[f"{endpoint}_iqr_outlier"].astype(bool)
    identifiers = data.loc[selected, "animal_id"].astype(str).tolist()
    return len(identifiers), ";".join(identifiers)


def compute_exact_statistics(animal_values: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for endpoint, endpoint_label in ENDPOINTS:
        groups = {
            condition: animal_values.loc[
                animal_values["condition"].eq(condition) & animal_values[endpoint].notna(), endpoint
            ].to_numpy(float)
            for condition in CONDITIONS
        }
        expected_sizes = (2, 3, 3) if endpoint == "cfos_neun_over_neun" else (3, 3, 3)
        observed_sizes = tuple(int(groups[condition].size) for condition in CONDITIONS)
        if observed_sizes != expected_sizes:
            raise RuntimeError(
                f"{endpoint}: expected finite group sizes {expected_sizes}, observed {observed_sizes}"
            )
        if any(not np.isfinite(values).all() for values in groups.values()):
            raise RuntimeError(f"{endpoint}: non-finite animal endpoint")
        analysis_set = "NeuN_complete_case" if endpoint == "cfos_neun_over_neun" else "all_nine_animals"

        for condition in CONDITIONS:
            values = groups[condition]
            flag_count, flag_ids = flagged_animals(animal_values, endpoint, (condition,))
            rows.append({
                "endpoint": endpoint,
                "endpoint_label": endpoint_label,
                "test": "descriptive",
                "comparison": condition,
                "n": int(values.size),
                "n_left": np.nan,
                "n_right": np.nan,
                "mean": float(np.mean(values)),
                "sd": float(np.std(values, ddof=1)),
                "statistic": np.nan,
                "p_value_raw": np.nan,
                "df_between": np.nan,
                "df_within": np.nan,
                "effect_size": np.nan,
                "effect_size_definition": "",
                "p_value_method": "",
                "extreme_labelings": np.nan,
                "enumerated_labelings": np.nan,
                "permutation_exhaustive": np.nan,
                "experimental_unit": "animal",
                "iqr_flag_count": flag_count,
                "iqr_flagged_animals": flag_ids,
                "iqr_rule": "within-condition Tukey 1.5-IQR",
                "outliers_excluded": False,
                "analysis_set": analysis_set,
            })

        global_result = one_way_anova([groups[condition] for condition in CONDITIONS])
        flag_count, flag_ids = flagged_animals(animal_values, endpoint, CONDITIONS)
        rows.append({
            "endpoint": endpoint,
            "endpoint_label": endpoint_label,
            "test": "one_way_ANOVA",
            "comparison": "Water|Sucrose|Allulose",
            "n": int(sum(observed_sizes)),
            "n_left": np.nan,
            "n_right": np.nan,
            "mean": np.nan,
            "sd": np.nan,
            **global_result,
            "effect_size": np.nan,
            "effect_size_definition": "",
            "p_value_method": "parametric one-way ANOVA using scipy.stats.f_oneway; raw unadjusted p value",
            "experimental_unit": "animal",
            "iqr_flag_count": flag_count,
            "iqr_flagged_animals": flag_ids,
            "iqr_rule": "within-condition Tukey 1.5-IQR",
            "outliers_excluded": False,
            "analysis_set": analysis_set,
        })

        for left_name, right_name in PAIR_ORDER:
            result = exact_mwu(groups[left_name], groups[right_name])
            expected_pair_labelings = math.comb(
                int(groups[left_name].size + groups[right_name].size), int(groups[left_name].size)
            )
            if int(result["enumerated_labelings"]) != expected_pair_labelings:
                raise RuntimeError(
                    f"{endpoint}/{left_name}-{right_name}: pair receipt is not {expected_pair_labelings}"
                )
            flag_count, flag_ids = flagged_animals(animal_values, endpoint, (left_name, right_name))
            rows.append({
                "endpoint": endpoint,
                "endpoint_label": endpoint_label,
                "test": "Mann_Whitney_U_two_sided",
                "comparison": f"{left_name}|{right_name}",
                "n": int(groups[left_name].size + groups[right_name].size),
                "n_left": int(groups[left_name].size),
                "n_right": int(groups[right_name].size),
                "mean": np.nan,
                "sd": np.nan,
                "df_between": np.nan,
                "df_within": np.nan,
                **result,
                "effect_size": cliffs_delta_right_minus_left(groups[left_name], groups[right_name]),
                "effect_size_definition": f"Cliff's delta ({right_name} minus {left_name})",
                "p_value_method": "exhaustive two-sided centered rank-sum/U independent-animal label permutation",
                "experimental_unit": "animal",
                "iqr_flag_count": flag_count,
                "iqr_flagged_animals": flag_ids,
                "iqr_rule": "within-condition Tukey 1.5-IQR",
                "outliers_excluded": False,
                "analysis_set": analysis_set,
            })

    result = pd.DataFrame(rows)
    pairwise = result["test"].eq("Mann_Whitney_U_two_sided")
    if not result.loc[pairwise, "permutation_exhaustive"].astype(bool).all():
        raise RuntimeError("Every MWU row must have an exhaustive receipt")
    if result.loc[result["test"].eq("one_way_ANOVA"), "permutation_exhaustive"].astype(bool).any():
        raise RuntimeError("Parametric ANOVA rows must not be labeled permutation-exhaustive")
    return result


def build_manifest(raw_root: Path, qc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for record in qc.itertuples(index=False):
        for role, path_value, sha, size, used in (
            ("cellpose_segmentation", record.segmentation_path, record.segmentation_sha256, record.segmentation_size_bytes, True),
            ("channel_tiff_companion", record.tiff_path, record.tiff_sha256, record.tiff_size_bytes, False),
        ):
            rows.append({
                "animal_id": record.animal_id,
                "sample_id": record.sample_id,
                "condition": record.condition,
                "channel": record.channel,
                "role": role,
                "used_by_numeric_analysis": used,
                "path": path_value,
                "size_bytes": int(size),
                "sha256": sha,
            })
    display_neun = raw_root / "water" / "WATER_NPY3_S01_NeuN_display_only.tif"
    display_receipt = raw_root / "water" / "WATER_NPY3_S01_NeuN_display_only_provenance.json"
    for path, role in (
        (display_neun, "display_only_failed_channel"),
        (display_receipt, "display_only_failure_receipt"),
    ):
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Missing Water3 NeuN display provenance: {path}")
        rows.append({
            "animal_id": "WATER_NPY3",
            "sample_id": "water-npy3-s01",
            "condition": "Water",
            "channel": "NeuN",
            "role": role,
            "used_by_numeric_analysis": False,
            "path": relative_to_paper(path),
            "size_bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
        })
    metadata_paths = (
        raw_root / "apotome_2025_07_28_npy_samplesheet.csv",
        raw_root / "water" / "apotome_2025_07_28_npy_water_samplesheet.csv",
        raw_root / "water" / "WATER_NPY3_S01_extraction_provenance.json",
        HERE.with_name("00_extract_water3_neun_display.py"),
        HERE,
    )
    for path in metadata_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Missing provenance companion: {path}")
        rows.append({
            "animal_id": "",
            "sample_id": "",
            "condition": "",
            "channel": "",
            "role": "analysis_code" if path == HERE else "metadata_companion",
            "used_by_numeric_analysis": path == HERE,
            "path": relative_to_paper(path),
            "size_bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
        })
    return pd.DataFrame(rows).sort_values(["role", "condition", "animal_id", "channel", "path"], kind="mergesort")


def methods_text(values: pd.DataFrame, stats: pd.DataFrame) -> str:
    pair_receipts = sorted(
        stats.loc[stats["test"].eq("Mann_Whitney_U_two_sided"), "enumerated_labelings"]
        .dropna().astype(int).unique().tolist()
    )
    archived_matches = int(values["archived_reproduction_status"].eq("exact_match").sum())
    return f"""Figure 3 whole-field c-FOS/NPY/NeuN analysis

Analysis population
-------------------
Nine explicitly mapped biological animals were used: Water n=3, Sucrose n=3,
and Allulose n=3. The c-FOS/NPY endpoints use all nine. The NeuN endpoint is an
explicit complete-case analysis using Water n=2, Sucrose n=3, and Allulose n=3;
Water3 remains NaN because its C2 NeuN channel failed. Pixels, Cellpose objects,
and nuclei are not statistical replicates.

Phenotype calls
---------------
DAPI, c-FOS, and NPY/GFP label images were loaded from Cellpose dictionary NPY
files. A c-FOS ROI was retained when its area was at least {CFOS_MIN_ROI_AREA_PX}
pixels and at least {CFOS_DAPI_OVERLAP_FRACTION:.2f} of its area overlapped its
best-overlapping DAPI nucleus. Multiple retained c-FOS ROIs mapping to one DAPI
nucleus yield one c-FOS-positive nucleus. A DAPI nucleus was NPY-positive when
at least {NPY_MIN_OVERLAP_PX} pixel and at least {NPY_DAPI_OVERLAP_FRACTION:.2f}
of its area overlapped the union of NPY segmentation masks. Double-positive
nuclei are the intersection of the c-FOS-positive and NPY-positive DAPI sets.
The same marker rule calls NeuN-positive DAPI nuclei using at least
{NEUN_MIN_OVERLAP_PX} NeuN-mask pixel and at least
{NEUN_DAPI_OVERLAP_FRACTION:.2f} of DAPI-nucleus area, without dilation.

Endpoints
---------
cfos_over_dapi = c-FOS-positive nuclei / DAPI nuclei.
npy_over_dapi = NPY-positive nuclei / DAPI nuclei.
cfos_npy_over_npy = double-positive nuclei / NPY-positive nuclei.
double_over_dapi = double-positive nuclei / DAPI nuclei (secondary).
cfos_neun_over_neun = c-FOS/NeuN double-positive nuclei / NeuN-positive nuclei
(NeuN-complete subset only).

Statistics
----------
Condition order is Water, Sucrose, Allulose. Bars may display animal mean plus
or minus SD, but every finite animal value must also be shown. Global tests are
ordinary parametric one-way ANOVA using scipy.stats.f_oneway with raw,
unadjusted p values. Pairwise W-S, W-A, and S-A values are raw exact two-sided
centered rank-sum/Mann-Whitney U p values; receipts enumerate
{','.join(map(str, pair_receipts))} possible allocations, depending on n=2/3 or
n=3/3. No plus-one correction or multiple-testing adjustment is applied.
Cliff's delta is right minus left. Tukey 1.5-IQR flags are descriptive only;
all finite values are retained. No cage analysis was attempted.

Water3 NeuN failure handling
----------------------------
Water3 C2 was losslessly exported from ordinary merged series 12 as a native
uint16 maximum projection for QC/display only. It has no accepted segmentation
and is explicitly ineligible for numeric analysis. Its NeuN count and ratio are
NaN, never zero. The complete-case endpoint uses only the eight animals with
source-bound, manually curated legacy NeuN masks.

Regression/provenance checks
----------------------------
All {archived_matches} legacy animals exactly reproduced embedded audited v10.1
whole-field DAPI, c-FOS, NPY, NeuN, and intersection totals without opening or
depending on retired/Trash outputs. Every quantitative mask and companion TIFF
was present and shape-compatible. SHA-256 hashes are in source_manifest.csv.

Interpretive limitation
-----------------------
WATER_NPY3 is a provisional animal ID and must be replaced by its true ID before
final inference. Water3 was acquired in 2026 on a Leica SP8 at 40X as a 21-plane
RGB16 z-maximum projection and its DAPI mask used model dapi_Keyence4. The legacy
animals are 2025 RGB8 exports from other acquisitions and generally used model
dapi_channel. Ratios reduce field-size differences but do not eliminate this
acquisition/segmentation batch confounding. The analysis is therefore explicitly
whole-field and must not be described as ARC, ME, or VMN regional quantification.
"""


def write_csv_atomic(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def main() -> int:
    args = build_parser().parse_args()
    validate_inventory()
    raw_root = args.raw_root.expanduser().resolve()
    if not raw_root.is_dir():
        raise FileNotFoundError(f"Figure 3 raw root does not exist: {raw_root}")
    output = prepare_output_dir(args.output_dir, args.force)

    animal_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for sample in SAMPLES:
        animal_row, sample_qc = quantify_sample(raw_root, sample)
        animal_rows.append(animal_row)
        qc_rows.extend(sample_qc)

    values = pd.DataFrame(animal_rows)
    condition_type = pd.CategoricalDtype(CONDITIONS, ordered=True)
    values["condition"] = values["condition"].astype(condition_type)
    values = values.sort_values(["condition", "animal_id"], kind="mergesort").reset_index(drop=True)
    values["condition"] = values["condition"].astype(str)
    values = add_iqr_flags(values)
    for endpoint, _label in ENDPOINTS:
        values[f"{endpoint}_percent"] = values[endpoint] * 100.0

    legacy = values.loc[values["animal_id"].isin(ARCHIVED_FIRST_EIGHT_COUNTS)]
    if len(legacy) != 8 or not legacy["archived_reproduction_status"].eq("exact_match").all():
        mismatched = legacy.loc[~legacy["archived_reproduction_status"].eq("exact_match"), "animal_id"].tolist()
        raise RuntimeError(f"Legacy whole-field regression check failed: {mismatched}")

    qc = pd.DataFrame(qc_rows).sort_values(["condition", "animal_id", "channel"], kind="mergesort")
    stats = compute_exact_statistics(values)
    receipt_columns = [
        "endpoint", "endpoint_label", "test", "comparison", "statistic",
        "p_value_raw", "extreme_labelings", "enumerated_labelings",
        "permutation_exhaustive", "p_value_method", "experimental_unit",
        "effect_size", "effect_size_definition", "outliers_excluded",
        "df_between", "df_within", "n", "n_left", "n_right", "analysis_set",
    ]
    receipts = stats.loc[
        stats["test"].isin(("one_way_ANOVA", "Mann_Whitney_U_two_sided")),
        receipt_columns,
    ].copy()
    manifest = build_manifest(raw_root, qc)

    write_csv_atomic(values, output / "per_animal_cfos_npy.csv")
    write_csv_atomic(stats, output / "anova_mwu_statistics.csv")
    write_csv_atomic(receipts, output / "statistical_test_receipts.csv")
    write_csv_atomic(qc, output / "segmentation_qc.csv")
    write_csv_atomic(manifest, output / "source_manifest.csv")
    readme_tmp = output / "README.txt.tmp"
    readme_tmp.write_text(methods_text(values, stats), encoding="utf-8")
    readme_tmp.replace(output / "README.txt")

    print(f"[OK] Figure 3 animal values: {len(values)} (3/condition)")
    print(f"[OK] Legacy regression: {len(legacy)}/8 exact matches")
    anova_rows = int(stats["test"].eq("one_way_ANOVA").sum())
    mwu_receipts = stats.loc[
        stats["test"].eq("Mann_Whitney_U_two_sided"),
        ["endpoint", "comparison", "enumerated_labelings"],
    ]
    denominators = sorted({int(value) for value in mwu_receipts["enumerated_labelings"]})
    print(f"[OK] Statistics: {len(stats)} rows; ANOVA={anova_rows}, exact-MWU denominators={denominators}")
    print(f"[OK] Outputs: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
