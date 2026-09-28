#!/usr/bin/env python3
"""Export the failed Water3 C2 channel as a provenance-bound display-only MIP.

The ordinary ``TileScan 1 Merged`` scene (Bio-Formats series 12), channel C2,
contains pixels but did not yield a biologically valid NeuN stain.  This script
preserves its native uint16 maximum projection for transparent QC/display.  It
never creates a segmentation and marks the output as ineligible for numerical
analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence
from uuid import uuid4

import numpy as np
import tifffile
from scipy.ndimage import binary_dilation


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

WATER_DIR = PAPER_ROOT / "Fig3/raw/legacy_experiment_2025_07_28/water"
SOURCE_LIF_NAME = "cFos-NeuN-Npy-Ctrl-Agua-40X.lif"
SOURCE_LIFEXT_NAME = "cFos-NeuN-Npy-Ctrl-Agua-40X.lifext"
OUTPUT_NAME = "WATER_NPY3_S01_NeuN_display_only.tif"
RECEIPT_NAME = "WATER_NPY3_S01_NeuN_display_only_provenance.json"
DAPI_SEG_NAME = "WATER_NPY3_S01_dapi_seg.npy"

SOURCE_LIF_SHA256 = "db4cbdb72c8e8af6a5ab264d1526083d4bb70c4ece54e9c9ee59a8c90d7671ab"
SOURCE_LIFEXT_SHA256 = "1359dd14c59ad5d558312db0a6ff71bbf1ecaf2660e126fabba039903eb22fc7"
DAPI_SEG_SHA256 = "58be1372ada83b0f799c27a32037e95e8ee28b4eeaaba92825b545fc91abad87"
EXPECTED_STACK_SHAPE = (21, 1419, 1898)
SERIES_INDEX = 12
CHANNEL_INDEX = 2
UM_PER_PX_X = 0.5687380073800737
UM_PER_PX_Y = 0.5687377997179125
Z_STEP_UM = 0.346255


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, role: str) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing {role}: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(f"{role} SHA256 mismatch: expected {expected}, observed {observed}")


def strict_cellpose_mask(path: Path) -> np.ndarray:
    payload = np.load(path, allow_pickle=True)
    if not (isinstance(payload, np.ndarray) and payload.dtype == object and payload.shape == ()):
        raise RuntimeError(f"DAPI source is not a scalar Cellpose dictionary: {path}")
    record = payload.item()
    if not isinstance(record, dict) or "masks" not in record:
        raise RuntimeError(f"DAPI source lacks Cellpose masks: {path}")
    labels = np.asarray(record["masks"])
    if labels.shape != EXPECTED_STACK_SHAPE[1:] or not np.issubdtype(labels.dtype, np.integer):
        raise RuntimeError(f"Unexpected DAPI mask shape/dtype: {labels.shape}/{labels.dtype}")
    return labels.astype(np.int32, copy=False)


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def export_display_only(
    *,
    water_dir: Path = WATER_DIR,
    bfconvert: Path,
) -> tuple[Path, Path]:
    water_dir = water_dir.expanduser().resolve()
    if water_dir.is_symlink() or WATER_DIR.resolve() != water_dir:
        raise RuntimeError(f"Output must be the canonical Figure 3 Water raw folder: {water_dir}")
    source_lif = water_dir / SOURCE_LIF_NAME
    source_lifext = water_dir / SOURCE_LIFEXT_NAME
    dapi_seg = water_dir / DAPI_SEG_NAME
    require_hash(source_lif, SOURCE_LIF_SHA256, "Water3 LIF")
    require_hash(source_lifext, SOURCE_LIFEXT_SHA256, "Water3 LIFEXT")
    require_hash(dapi_seg, DAPI_SEG_SHA256, "Water3 DAPI segmentation")
    if not bfconvert.is_file() or not os.access(bfconvert, os.X_OK):
        raise FileNotFoundError(f"Bio-Formats bfconvert is unavailable: {bfconvert}")

    dapi_labels = strict_cellpose_mask(dapi_seg)
    with tempfile.TemporaryDirectory(prefix="figure3-water3-neun-") as temporary_dir:
        temporary_root = Path(temporary_dir)
        stack_path = temporary_root / "water3_c2_stack.ome.tif"
        subprocess.run(
            [
                str(bfconvert), "-series", str(SERIES_INDEX),
                "-channel", str(CHANNEL_INDEX), "-compression", "zlib",
                str(source_lif), str(stack_path),
            ],
            check=True,
        )
        stack = np.asarray(tifffile.imread(stack_path))
        if stack.shape != EXPECTED_STACK_SHAPE or stack.dtype != np.uint16:
            raise RuntimeError(
                f"Unexpected C2 extraction: {stack.shape}/{stack.dtype}; "
                f"expected {EXPECTED_STACK_SHAPE}/uint16"
            )
        mip = np.max(stack, axis=0).astype(np.uint16, copy=False)
        output = water_dir / OUTPUT_NAME
        staged_output = water_dir / f".{OUTPUT_NAME}.{uuid4().hex}.tmp"
        try:
            tifffile.imwrite(
                staged_output,
                mip,
                photometric="minisblack",
                compression="zlib",
                metadata=None,
            )
            check = np.asarray(tifffile.imread(staged_output))
            if check.dtype != np.uint16 or not np.array_equal(check, mip):
                raise RuntimeError("Display-only NeuN TIFF did not round-trip pixel-exactly")
            os.replace(staged_output, output)
        finally:
            if staged_output.exists():
                staged_output.unlink()

    dapi_positive = dapi_labels > 0
    ring = binary_dilation(dapi_positive, iterations=3) & ~dapi_positive
    inside_median = float(np.median(mip[dapi_positive]))
    outside_median = float(np.median(mip[~dapi_positive]))
    ring_median = float(np.median(mip[ring]))
    receipt = water_dir / RECEIPT_NAME
    payload: dict[str, Any] = {
        "schema_version": "figure3-water3-neun-display-only-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": "WATER_NPY3_S01",
        "animal_id": "WATER_NPY3",
        "animal_id_status": "provisional; replace with true animal ID before inference",
        "condition": "Water",
        "source": {
            "lif_path": str(source_lif),
            "lif_sha256": SOURCE_LIF_SHA256,
            "lifext_path": str(source_lifext),
            "lifext_sha256": SOURCE_LIFEXT_SHA256,
            "bioformats_series_index": SERIES_INDEX,
            "series_name": "TileScan 1 Merged",
            "channel_index": CHANNEL_INDEX,
            "nominal_marker": "NeuN",
            "stack_shape_z_y_x": list(EXPECTED_STACK_SHAPE),
            "dtype": "uint16",
        },
        "processing": {
            "operation": "maximum projection across all 21 source Z planes",
            "resampling": "none",
            "intensity_rescaling": False,
            "denoising": False,
            "segmentation": False,
        },
        "calibration": {
            "physical_size_x_um_per_px": UM_PER_PX_X,
            "physical_size_y_um_per_px": UM_PER_PX_Y,
            "z_step_um": Z_STEP_UM,
        },
        "output": {
            "path": str(output),
            "sha256": sha256_file(output),
            "bytes": output.stat().st_size,
            "shape_y_x": list(mip.shape),
            "dtype": str(mip.dtype),
            "minimum": int(mip.min()),
            "maximum": int(mip.max()),
            "mean": float(mip.mean()),
            "saturated_65535_pixels": int(np.count_nonzero(mip == np.iinfo(np.uint16).max)),
        },
        "failure_qc": {
            "dapi_segmentation_path": str(dapi_seg),
            "dapi_segmentation_sha256": DAPI_SEG_SHA256,
            "dapi_nuclei": int(np.unique(dapi_labels[dapi_labels > 0]).size),
            "c2_median_inside_dapi_masks": inside_median,
            "c2_median_outside_dapi_masks": outside_median,
            "c2_median_local_dapi_ring": ring_median,
            "assessment": (
                "Failed/unquantifiable NeuN channel: the C2 MIP is dominated by diffuse/background "
                "signal and DAPI nuclei generally lack the bright soma-like NeuN morphology present "
                "in the eight legacy NeuN fields."
            ),
        },
        "eligibility": {
            "display_qc_only": True,
            "numeric_analysis": False,
            "cellpose_input": False,
            "segmentation_available": False,
            "missing_value_policy": "Water3 NeuN values must remain NaN/missing, never zero",
        },
        "generator": {
            "path": str(HERE),
            "sha256": sha256_file(HERE),
            "bfconvert_path": str(bfconvert),
        },
    }
    atomic_write_json(receipt, payload)
    return output, receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--water-dir", type=Path, default=WATER_DIR)
    parser.add_argument(
        "--bfconvert",
        type=Path,
        default=Path("/home/server/anaconda3/envs/microscopy_extract/bin/bfconvert"),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output, receipt = export_display_only(water_dir=args.water_dir, bfconvert=args.bfconvert)
    print(f"display_only_neun_tiff\t{output}")
    print(f"provenance_receipt\t{receipt}")
    print("numeric_analysis\tFalse")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
