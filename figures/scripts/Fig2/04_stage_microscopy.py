#!/usr/bin/env python3
"""Stage Figure 2 panel-C native display assets with the hand-made-figure framing.

This reuses Fig2/03_reconstruct_microscopy_engine.py unchanged as the
single source of truth for CZI reading, scene-0 stitching, calibration gates and
SIFT/RANSAC field validation.  Only three display decisions are overridden, each
of which is a framing or black-point choice and none of which touches
quantification:

  * Water is oriented 270 deg clockwise instead of 90 deg, so the median
    eminence sits at the bottom of the field as in the hand-made Figure 2.  This
    is a rigid rotation of the same scene-0 mosaic - no mirroring - and the crop
    centre is remapped accordingly.
  * Sucrose and Allulose crop centres move down the mosaic so the whole median
    eminence falls inside the common 600 x 1200 um field instead of being clipped
    at the bottom edge.
  * The c-FOS display black point moves from the 0.5th percentile to the field
    median, uniformly for every condition.  The Water c-FOS channel carries a
    high background floor that the 0.5th-percentile stretch rendered as an
    all-over red wash; clipping at the 35th percentile removes it while
    leaving the Sucrose and Allulose puncta intact.
  * The accepted per-condition raw low/high limits are frozen before the tile
    registration repair, so registration cannot change displayed intensity.

The common 600 x 1200 um physical field, the DAPI stretch, the pseudocolours,
gamma, and every calibration and identity gate are unchanged.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
# The paper root is found by walking up to the directory that holds README.txt
# and scripts/setup, so this works from Fig2/ and from its scripts/Fig2 mirror
# alike. Deriving it as HERE.parent.parent silently resolved to Paper/scripts
# when run from the mirror, and every raw path under it then pointed at
# scripts/Fig2/raw_data, which does not exist.
PAPER = next((path for path in HERE.parents
              if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                  and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
FIGURE = PAPER / "Fig2"
RAW = FIGURE / "raw_data" / "Apotome"
RECONSTRUCTOR = FIGURE / "03_reconstruct_microscopy_engine.py"

# condition -> (raw czi, display reference tiff, processed czi reference)
SOURCES = {
    "Water": (
        RAW / "24_06_2025/analysis/E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_RAW_001.czi",
        RAW / "24_06_2025/analysis/E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001-0001_dapi.tif",
        RAW / "24_06_2025/analysis/E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001.czi",
    ),
    "Sucrose": (
        RAW / "Aug_2025_NPY/selected/E7FR7-1H-sacarosa-NPY488-NEUN594-CFOS694-DAPI-RAW.czi",
        RAW / "Aug_2025_NPY/selected/E7FR7-1H-sacarosa-DAPI.tif",
        None,
    ),
    "Allulose": (
        RAW / "Aug_2025_NPY/selected/E8FR6-4H-alulosa-NPY488-NEUN594-CFOS694-DAPI-02-RAW.czi",
        RAW / "Aug_2025_NPY/selected/E8FR6-4H-alulosa-DAPI.tif",
        None,
    ),
}

# condition -> (clockwise quarter turns, crop centre x, crop centre y)
FRAMING = {
    "Water": (3, 2992, 6850),
    "Sucrose": (0, 2076, 7100),
    "Allulose": (0, 2080, 5377),
}

# Display black point per marker, as a percentile of the cropped field.
LOW_PERCENTILE_BY_MARKER = {"DAPI": 0.5, "cFOS": 35.0}

# Exact accepted display limits from the pre-registration Figure 2C stage. They
# are frozen so the registration repair cannot silently change color intensity,
# contrast, gamma, saturation, or brightness.
DISPLAY_LIMITS = {
    ("Water", "DAPI"): (531.0, 15567.479999999981),
    ("Water", "cFOS"): (3855.0, 8315.0),
    ("Sucrose", "DAPI"): (831.0, 11258.448000000091),
    ("Sucrose", "cFOS"): (2147.0, 6288.586000000592),
    ("Allulose", "DAPI"): (711.0, 9652.70399999991),
    ("Allulose", "cFOS"): (2007.0, 7842.4260000002105),
}


def load_reconstructor():
    spec = importlib.util.spec_from_file_location("fig2_manmade_native", RECONSTRUCTOR)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {RECONSTRUCTOR}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def patched_linear_display(module):
    """Freeze accepted intensity limits while changing only tile registration."""

    call_order = iter(
        (condition, marker)
        for condition in ("Water", "Sucrose", "Allulose")
        for marker in ("DAPI", "cFOS")
    )

    def linear_display(channel: np.ndarray, marker: str):
        condition, expected_marker = next(call_order)
        if marker != expected_marker:
            raise RuntimeError(
                f"Unexpected display call order: {condition} {marker}"
            )
        low_percentile = LOW_PERCENTILE_BY_MARKER[marker]
        low, high = DISPLAY_LIMITS[(condition, marker)]
        if not high > low:
            raise RuntimeError(f"Degenerate display range for {marker}: {low}, {high}")
        stride = max(
            1,
            int(np.ceil(max(channel.shape) / module.DISPLAY_MAX_DIMENSION_PX)),
        )
        sampled = channel[::stride, ::stride].astype(np.float32)
        normalized = np.clip((sampled - low) / (high - low), 0.0, 1.0)
        return normalized, {
            "low_percentile": low_percentile,
            "high_percentile": module.HIGH_PERCENTILE[marker],
            "raw_low": low,
            "raw_high": high,
            "percentile_sampling_stride_px": module.PERCENTILE_SAMPLE_STRIDE,
            "display_sampling_stride_px": stride,
            "transform": "clip((raw-low)/(high-low),0,1)",
            "gamma": 1.0,
            "clahe": False,
            "spatial_background_subtraction": False,
            "display_limits_source": (
                "frozen_from_pre_registration_accepted_Figure2C; "
                "registration-only comparison"
            ),
        }

    return linear_display


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Fresh absolute native_czi output directory")
    args = parser.parse_args()
    outdir = args.outdir.expanduser()
    if not outdir.is_absolute():
        raise ValueError("Output directory must be absolute")
    if outdir.exists() or outdir.is_symlink():
        raise FileExistsError(f"Refusing to replace an existing native stage: {outdir}")

    module = load_reconstructor()
    module.FIELDS = tuple(
        dataclasses.replace(
            field,
            raw_czi=SOURCES[field.condition][0],
            display_reference_tiff=SOURCES[field.condition][1],
            processed_czi_reference=SOURCES[field.condition][2],
            rotate_clockwise_quarters=FRAMING[field.condition][0],
            crop_center_x_px=FRAMING[field.condition][1],
            crop_center_y_px=FRAMING[field.condition][2],
        )
        for field in module.FIELDS
    )
    module.linear_display = patched_linear_display(module)
    module.LOW_PERCENTILE = LOW_PERCENTILE_BY_MARKER["DAPI"]
    module.DEFAULT_OUTDIR = outdir
    # Audit only the explicitly bound native/display roles; never crawl broadly.
    module.candidate_paths = lambda: sorted(
        {
            path
            for field in module.FIELDS
            for path in (field.raw_czi, field.processed_czi_reference)
            if path is not None
        },
        key=lambda path: str(path).casefold(),
    )

    previous_argv = sys.argv
    try:
        sys.argv = [str(RECONSTRUCTOR), "--outdir", str(outdir)]
        return int(module.main())
    finally:
        sys.argv = previous_argv


if __name__ == "__main__":
    raise SystemExit(main())
