#!/usr/bin/env python3
"""Reconstruct Figure 2C directly from the native 16-bit CZI acquisitions.

This is a display-only workflow. It deliberately does not alter or feed the
animal-level c-FOS/DAPI quantification.  Scene 0 tiles are read individually,
placed at the stage coordinates stored in each CZI, and cropped to the same
600 x 1200 micrometre field for all three representative animals.  Display
contrast is a documented linear percentile stretch (gamma 1; no CLAHE and no
spatial filtering).

Run with the ``apotome_quant`` environment because it contains libCZI:

    /home/server/anaconda3/envs/apotome_quant/bin/python \
      Fig2/03_reconstruct_microscopy_engine.py (this file)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

from aicspylibczi import CziFile
import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import tifffile


HERE = Path(__file__).resolve()
PAPER = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")
ROOT = PAPER.parent
DEFAULT_OUTDIR = PAPER / "analyses" / "Fig2" / "results" / "fig2c_raw_reconstruction"

COMMON_WIDTH_UM = 600.0
COMMON_HEIGHT_UM = 1200.0
LOW_PERCENTILE = 0.5
HIGH_PERCENTILE = {"DAPI": 99.8, "cFOS": 99.95}
PERCENTILE_SAMPLE_STRIDE = 8
DISPLAY_MAX_DIMENSION_PX = 2200
SCALEBAR_UM = 200.0
DAPI_RGB = np.array([0.02, 0.58, 1.00], dtype=np.float32)
CFOS_RGB = np.array([1.00, 0.02, 0.10], dtype=np.float32)

# Seam-derived native-pixel translations. The same integer placement is used
# for DAPI and c-FOS; tiles are never interpolated or blended.
TILE_REGISTRATION = {
    "Water": {"column_dy_px": 82, "row_dx_px": -66},
    "Sucrose": {"row_dx_px": -66},
    "Allulose": {"row_dx_px": -66},
}


@dataclass(frozen=True)
class FieldSpec:
    condition: str
    animal: str
    raw_czi: Path
    dapi_channel: int
    cfos_channel: int
    cfos_assignment_evidence: str
    rotate_clockwise_quarters: int
    crop_center_x_px: int
    crop_center_y_px: int
    expected_scene0_tiles: int
    expected_oriented_shape_yx: tuple[int, int]
    expected_pixel_size_um: float
    display_reference_tiff: Path
    processed_czi_reference: Path | None


FIELDS = (
    FieldSpec(
        condition="Water",
        animal="E10_FR1-1",
        raw_czi=(
            ROOT / "24_06_2025" / "analysis"
            / "E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_RAW_001.czi"
        ),
        dapi_channel=0,
        cfos_channel=2,
        cfos_assignment_evidence=(
            "CZI filename DAPI_POMC488_CFOS647 plus metadata Channel:2 "
            "Alexa Fluor 647 (Ex 653 nm; Em 668 nm)"
        ),
        rotate_clockwise_quarters=1,
        # Centre of the established legacy-matched field mapped into the
        # clockwise-oriented RAW scene-0 mosaic by SIFT/RANSAC.
        crop_center_x_px=3248,
        crop_center_y_px=6078,
        expected_scene0_tiles=6,
        expected_oriented_shape_yx=(12493, 6240),
        expected_pixel_size_um=0.13486980359192787,
        display_reference_tiff=(
            ROOT / "24_06_2025" / "analysis"
            / "E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001-0001_dapi.tif"
        ),
        processed_czi_reference=(
            ROOT / "24_06_2025" / "analysis"
            / "E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001.czi"
        ),
    ),
    FieldSpec(
        condition="Sucrose",
        animal="E7_FR7-1H",
        raw_czi=(
            ROOT / "Aug_2025_NPY" / "selected"
            / "E7FR7-1H-sacarosa-NPY488-NEUN594-CFOS694-DAPI-RAW.czi"
        ),
        dapi_channel=0,
        cfos_channel=3,
        cfos_assignment_evidence=(
            "CZI filename NPY488_NEUN594_CFOS694_DAPI plus metadata Channel:3 "
            "Alexa Fluor 647 (Ex 653 nm; Em 668 nm)"
        ),
        rotate_clockwise_quarters=0,
        crop_center_x_px=2076,
        crop_center_y_px=5998,
        expected_scene0_tiles=4,
        expected_oriented_shape_yx=(12481, 4165),
        expected_pixel_size_um=0.16405451188059882,
        display_reference_tiff=(
            ROOT / "Aug_2025_NPY" / "selected" / "E7FR7-1H-sacarosa-DAPI.tif"
        ),
        processed_czi_reference=None,
    ),
    FieldSpec(
        condition="Allulose",
        animal="E8_FR6-4H",
        raw_czi=(
            ROOT / "Aug_2025_NPY" / "selected"
            / "E8FR6-4H-alulosa-NPY488-NEUN594-CFOS694-DAPI-02-RAW.czi"
        ),
        dapi_channel=0,
        cfos_channel=3,
        cfos_assignment_evidence=(
            "CZI filename NPY488_NEUN594_CFOS694_DAPI plus metadata Channel:3 "
            "Alexa Fluor 647 (Ex 653 nm; Em 668 nm)"
        ),
        rotate_clockwise_quarters=0,
        crop_center_x_px=2080,
        crop_center_y_px=4677,
        expected_scene0_tiles=3,
        expected_oriented_shape_yx=(9360, 4165),
        expected_pixel_size_um=0.15372185638803682,
        display_reference_tiff=(
            ROOT / "Aug_2025_NPY" / "selected" / "E8FR6-4H-alulosa-DAPI.tif"
        ),
        processed_czi_reference=None,
    ),
)


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local_name(element: object) -> str:
    return str(getattr(element, "tag", "")).split("}")[-1]


def child_text(element: object, name: str) -> str | None:
    for child in element:
        if local_name(child) == name:
            return (child.text or "").strip()
    return None


def czi_metadata(czi: CziFile) -> dict[str, object]:
    root = czi.meta
    distances: dict[str, float] = {}
    for element in root.iter():
        if local_name(element) != "Distance":
            continue
        axis = element.attrib.get("Id")
        value = child_text(element, "Value")
        if axis in {"X", "Y", "Z"} and value:
            distances[axis] = float(value)
    if not {"X", "Y"}.issubset(distances):
        raise RuntimeError("CZI lacks X/Y scaling metadata")
    # Zeiss stores Scaling/Distance values in metres even when the display unit
    # field says micrometres.
    pixel_size_x_um = distances["X"] * 1_000_000.0
    pixel_size_y_um = distances["Y"] * 1_000_000.0
    if not np.isclose(pixel_size_x_um, pixel_size_y_um, rtol=0, atol=1e-12):
        raise RuntimeError(f"Anisotropic CZI pixels: {pixel_size_x_um}, {pixel_size_y_um}")

    channels: dict[int, dict[str, object]] = {}
    for element in root.iter():
        if local_name(element) != "Channel":
            continue
        match = re.fullmatch(r"Channel:(\d+)", element.attrib.get("Id", ""))
        if not match or child_text(element, "ExcitationWavelength") is None:
            continue
        index = int(match.group(1))
        channels[index] = {
            "index": index,
            "name": element.attrib.get("Name", ""),
            "fluor": child_text(element, "Fluor"),
            "excitation_nm": float(child_text(element, "ExcitationWavelength") or "nan"),
            "emission_nm": float(child_text(element, "EmissionWavelength") or "nan"),
            "exposure_ns": int(child_text(element, "ExposureTime") or "0"),
            "pixel_type": child_text(element, "PixelType"),
            "component_bit_count": int(child_text(element, "ComponentBitCount") or "0"),
        }
    return {
        "dims": czi.dims,
        "dims_shape": czi.get_dims_shape(),
        "is_mosaic": bool(czi.is_mosaic()),
        "pixel_size_um": pixel_size_x_um,
        "channels": channels,
    }


def scene0_geometry(czi: CziFile, channel: int) -> tuple[tuple[int, int, int, int], list[dict[str, int]]]:
    scene_boxes = czi.get_all_mosaic_scene_bounding_boxes()
    if 0 not in scene_boxes:
        raise RuntimeError("Expected CZI scene 0")
    box = scene_boxes[0]
    scene = (int(box.x), int(box.y), int(box.w), int(box.h))
    tiles: list[dict[str, int]] = []
    for info, bbox in czi.get_all_mosaic_tile_bounding_boxes(C=channel, S=0).items():
        coords = dict(info.dimension_coordinates)
        if int(coords.get("S", -1)) != 0 or int(coords.get("C", -1)) != channel:
            continue
        tiles.append({
            "m": int(coords["M"]),
            "x": int(bbox.x),
            "y": int(bbox.y),
            "w": int(bbox.w),
            "h": int(bbox.h),
        })
    tiles.sort(key=lambda row: row["m"])
    return scene, tiles


def condition_from_scene_geometry(width: int, height: int, tile_count: int) -> str:
    """Bind the three audited scene-0 geometries to their treatment fields."""
    mapping = {
        (12493, 6240, 6): "Water",
        (4165, 12481, 4): "Sucrose",
        (4165, 9360, 3): "Allulose",
    }
    key = (width, height, tile_count)
    if key not in mapping:
        raise RuntimeError(f"Unrecognized Figure 2C scene geometry: {key}")
    return mapping[key]


def read_scene0_mosaic(czi: CziFile, channel: int) -> tuple[np.ndarray, dict[str, object]]:
    """Register native scene-0 tiles by audited integer translations.

    Each complete uint16 tile is copied bit-exactly. The translated Water grid
    leaves two small no-data junctions; only those coordinates are filled from
    the same RAW tile pixels at their nominal acquisition ownership. No pixel
    is averaged, interpolated, blended, or intensity transformed here.
    """
    (scene_x, scene_y, width, height), tiles = scene0_geometry(czi, channel)
    condition = condition_from_scene_geometry(width, height, len(tiles))
    mosaic = np.zeros((height, width), dtype=np.uint16)
    covered = np.zeros((height, width), dtype=np.bool_)
    tile_images: dict[int, np.ndarray] = {}
    placements: list[dict[str, int]] = []

    for tile in tiles:
        array, dims = czi.read_image(C=channel, S=0, M=tile["m"])
        image = np.squeeze(array)
        if image.ndim != 2 or image.dtype != np.uint16:
            raise RuntimeError(
                f"Expected one Gray16 plane for C={channel}, M={tile['m']}; "
                f"found {image.shape} {image.dtype} ({dims})"
            )
        if image.shape != (tile["h"], tile["w"]):
            raise RuntimeError(f"Tile shape/metadata mismatch for M={tile['m']}")
        tile_images[tile["m"]] = image

        if condition == "Water":
            row, column = divmod(tile["m"], 3)
        else:
            row, column = tile["m"], 0

        # CZI metadata has one-pixel rounding differences between neighboring
        # boxes. Normalize those into the zero-overlap acquisition grid first.
        grid_x = column * tile["w"]
        grid_y = row * tile["h"]
        grid_x1 = grid_x + tile["w"]
        grid_y1 = grid_y + tile["h"]
        if grid_x1 > width or grid_y1 > height:
            raise RuntimeError(f"Ideal tile grid falls outside {condition} scene")
        occupied = covered[grid_y:grid_y1, grid_x:grid_x1]
        if occupied.any():
            raise RuntimeError(f"Ideal grid overlaps for {condition} M={tile['m']}")
        mosaic[grid_y:grid_y1, grid_x:grid_x1] = image
        occupied[:] = True

        target_x = grid_x + row * TILE_REGISTRATION[condition]["row_dx_px"]
        target_y = grid_y + column * TILE_REGISTRATION[condition].get("column_dy_px", 0)
        placements.append({
            "m": tile["m"],
            "source_x_px": tile["x"] - scene_x,
            "source_y_px": tile["y"] - scene_y,
            "registered_x_px": target_x,
            "registered_y_px": target_y,
            "dx_px": target_x - (tile["x"] - scene_x),
            "dy_px": target_y - (tile["y"] - scene_y),
        })

    registered = np.zeros_like(mosaic)
    registered_covered = np.zeros_like(covered)
    clipped_source_pixels = 0
    for tile, placement in zip(tiles, placements, strict=True):
        image = tile_images[tile["m"]]
        target_x = placement["registered_x_px"]
        target_y = placement["registered_y_px"]
        canvas_x0 = max(0, target_x)
        canvas_y0 = max(0, target_y)
        canvas_x1 = min(width, target_x + tile["w"])
        canvas_y1 = min(height, target_y + tile["h"])
        source_x0 = canvas_x0 - target_x
        source_y0 = canvas_y0 - target_y
        source_x1 = source_x0 + canvas_x1 - canvas_x0
        source_y1 = source_y0 + canvas_y1 - canvas_y0
        destination = registered[canvas_y0:canvas_y1, canvas_x0:canvas_x1]
        occupied = registered_covered[canvas_y0:canvas_y1, canvas_x0:canvas_x1]
        if occupied.any():
            raise RuntimeError(f"Registered tiles overlap for {condition} M={tile['m']}")
        destination[:] = image[source_y0:source_y1, source_x0:source_x1]
        occupied[:] = True
        clipped_source_pixels += int(
            image.size - (source_y1 - source_y0) * (source_x1 - source_x0)
        )

    # Fill only translated-grid no-data coordinates from exact nominal RAW
    # ownership. This affects 10,824 pixels in the accepted Water crop and zero
    # pixels in the Sucrose and Allulose crops.
    fallback = (~registered_covered) & covered
    registered[fallback] = mosaic[fallback]
    registered_covered[fallback] = True

    spec = next(field for field in FIELDS if field.condition == condition)
    oriented_coverage = orient(registered_covered, spec.rotate_clockwise_quarters)
    oriented_fallback = orient(fallback, spec.rotate_clockwise_quarters)
    x0, y0, x1, y1 = centered_crop_bounds(spec, spec.expected_pixel_size_um)
    crop_uncovered = int((~oriented_coverage[y0:y1, x0:x1]).sum())
    crop_fallback = int(oriented_fallback[y0:y1, x0:x1].sum())
    if crop_uncovered:
        raise RuntimeError(
            f"Registered {condition} crop has {crop_uncovered} uncovered pixels"
        )

    return registered, {
        "scene_bbox_absolute_xywh": [scene_x, scene_y, width, height],
        "tile_count": len(tiles),
        "tiles": tiles,
        "registered_placements": placements,
        "tile_overlap_pixels_averaged": 0,
        "uncovered_scene_pixels": int((~registered_covered).sum()),
        "crop_uncovered_pixels": crop_uncovered,
        "crop_nominal_raw_fallback_pixels": crop_fallback,
        "clipped_source_pixels_outside_scene_canvas": clipped_source_pixels,
        "registration": TILE_REGISTRATION[condition],
        "registration_order": (
            "direct per-tile translations; nominal RAW ownership only at no-data junctions"
        ),
        "stitching": (
            "scene-0 native uint16 tiles placed by seam-derived integer translations; "
            "zero-overlap junction holes take exact nominal-grid RAW pixels; no "
            "interpolation, no blending, and the same transform for DAPI/c-FOS"
        ),
    }


def orient(image: np.ndarray, clockwise_quarters: int) -> np.ndarray:
    quarters = clockwise_quarters % 4
    return np.rot90(image, k=(-quarters) % 4) if quarters else image


def centered_crop_bounds(spec: FieldSpec, pixel_size_um: float) -> tuple[int, int, int, int]:
    width_px = int(round(COMMON_WIDTH_UM / pixel_size_um))
    height_px = int(round(COMMON_HEIGHT_UM / pixel_size_um))
    x0 = spec.crop_center_x_px - width_px // 2
    y0 = spec.crop_center_y_px - height_px // 2
    return x0, y0, x0 + width_px, y0 + height_px


def linear_display(channel: np.ndarray, marker: str) -> tuple[np.ndarray, dict[str, float]]:
    sample = channel[::PERCENTILE_SAMPLE_STRIDE, ::PERCENTILE_SAMPLE_STRIDE]
    finite = sample[np.isfinite(sample)]
    if not finite.size:
        raise RuntimeError(f"No display pixels for {marker}")
    low = float(np.percentile(finite, LOW_PERCENTILE))
    high = float(np.percentile(finite, HIGH_PERCENTILE[marker]))
    if not high > low:
        raise RuntimeError(f"Degenerate display range for {marker}: {low}, {high}")
    stride = max(1, int(math.ceil(max(channel.shape) / DISPLAY_MAX_DIMENSION_PX)))
    sampled = channel[::stride, ::stride].astype(np.float32)
    normalized = np.clip((sampled - low) / (high - low), 0.0, 1.0)
    return normalized, {
        "low_percentile": LOW_PERCENTILE,
        "high_percentile": HIGH_PERCENTILE[marker],
        "raw_low": low,
        "raw_high": high,
        "percentile_sampling_stride_px": PERCENTILE_SAMPLE_STRIDE,
        "display_sampling_stride_px": stride,
        "transform": "clip((raw-low)/(high-low),0,1)",
        "gamma": 1.0,
        "clahe": False,
        "spatial_background_subtraction": False,
    }


def rgb_u8(signal: np.ndarray, color: np.ndarray) -> np.ndarray:
    return np.rint(np.clip(signal[..., None] * color, 0.0, 1.0) * 255).astype(np.uint8)


def save_raw_tiff(path: Path, array: np.ndarray, pixel_size_um: float) -> None:
    tifffile.imwrite(
        path,
        array,
        imagej=True,
        resolution=(1.0 / pixel_size_um, 1.0 / pixel_size_um),
        metadata={"unit": "micron", "axes": "YX"},
        compression="deflate",
        predictor=True,
    )


def save_display_pair(stem: Path, rgb: np.ndarray, display_pixel_size_um: float) -> list[Path]:
    png = stem.with_suffix(".png")
    tiff = stem.with_suffix(".tif")
    Image.fromarray(rgb, mode="RGB").save(png, compress_level=9, optimize=False)
    tifffile.imwrite(
        tiff,
        rgb,
        photometric="rgb",
        resolution=(1.0 / display_pixel_size_um, 1.0 / display_pixel_size_um),
        metadata={"unit": "micron", "axes": "YXS"},
        compression="deflate",
        predictor=True,
    )
    return [png, tiff]


def normalize_u8(image: np.ndarray) -> np.ndarray:
    finite = image[np.isfinite(image)]
    low, high = np.percentile(finite, [1.0, 99.8])
    return np.rint(np.clip((image.astype(np.float32) - low) / (high - low), 0, 1) * 255).astype(np.uint8)


def reference_signal(path: Path, stride: int = 10) -> np.ndarray:
    source = tifffile.memmap(path)
    view = np.asarray(source[::stride, ::stride])
    if view.ndim == 3:
        view = view.max(axis=2)
    return normalize_u8(view)


def feature_match_validation(raw_oriented: np.ndarray, reference: Path) -> dict[str, object]:
    raw_preview = normalize_u8(raw_oriented[::10, ::10])
    reference_preview = reference_signal(reference, 10)
    sift = cv2.SIFT_create(nfeatures=12_000, contrastThreshold=0.01)
    key_a, descriptors_a = sift.detectAndCompute(raw_preview, None)
    key_b, descriptors_b = sift.detectAndCompute(reference_preview, None)
    if descriptors_a is None or descriptors_b is None:
        raise RuntimeError(f"SIFT found no descriptors for {reference}")
    matches = cv2.BFMatcher().knnMatch(descriptors_a, descriptors_b, k=2)
    good = [first for first, second in matches if first.distance < 0.70 * second.distance]
    if len(good) < 8:
        raise RuntimeError(f"Insufficient native-to-reference matches: {len(good)}")
    source = np.float32([key_a[match.queryIdx].pt for match in good])
    target = np.float32([key_b[match.trainIdx].pt for match in good])
    homography, inlier_mask = cv2.findHomography(source, target, cv2.RANSAC, 3.0)
    if homography is None or inlier_mask is None:
        raise RuntimeError("Native-to-reference RANSAC failed")
    inliers = int(inlier_mask.sum())
    if inliers < 75:
        raise RuntimeError(f"Only {inliers} native-to-reference inliers")
    return {
        "reference_path": str(reference),
        "reference_sha256": sha256(reference),
        "raw_preview_shape_yx": list(raw_preview.shape),
        "reference_preview_shape_yx": list(reference_preview.shape),
        "sift_keypoints_raw": len(key_a),
        "sift_keypoints_reference": len(key_b),
        "ratio_test_matches": len(good),
        "ransac_inliers": inliers,
        "ransac_inlier_ratio": inliers / len(good),
        "homography_raw_preview_to_reference_preview": homography.tolist(),
        "status": "verified_same_field",
        "role": "validation only; reference pixels are not used in reconstruction",
    }


def candidate_paths() -> list[Path]:
    pattern = re.compile(r"(E10.?FR1.?1|E7.?FR7.?1H|E8.?FR6.?4H)", re.IGNORECASE)
    roots = [
        ROOT / "24_06_2025",
        ROOT / "2025.06.24",
        ROOT / "Aug_2025_NPY",
        ROOT / "2025.07.28_NPY",
        ROOT / "2025.07.29",
    ]
    paths: list[Path] = []
    for search_root in roots:
        if not search_root.is_dir():
            continue
        paths.extend(
            path for path in search_root.rglob("*.czi")
            if pattern.search(path.name)
        )
    return sorted(set(paths), key=lambda path: str(path).casefold())


def candidate_status(path: Path) -> tuple[str, str]:
    if any(path == field.raw_czi for field in FIELDS):
        return "selected_native_source", "exact representative field in selected analysis and legacy Panel C"
    if any(field.processed_czi_reference is not None and path == field.processed_czi_reference for field in FIELDS):
        return "validation_only_processed_czi", "same Water field; processed CZI used only to validate lineage/orientation"
    lower = path.name.casefold()
    if "stiched" in lower or "stitched" in lower:
        return "excluded_other_processed_field", "processed/stitch derivative is not the selected native representative source"
    if "raw" in lower:
        return "excluded_other_raw_field_or_copy", "different slide/field suffix or redundant download copy"
    return "excluded_alternate_acquisition", "different assay, magnification, slide, or field identifier"


def audit_candidates(target_hashes: dict[Path, str], hash_all: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in candidate_paths():
        status, reason = candidate_status(path)
        digest = target_hashes.get(path)
        if digest is None and hash_all:
            digest = sha256(path)
        try:
            czi = CziFile(path)
            shape = czi.get_dims_shape()
            dims = czi.dims
            mosaic = bool(czi.is_mosaic())
            probe = "opened"
            error = ""
        except Exception as exc:  # retained explicitly in the audit
            shape = []
            dims = ""
            mosaic = False
            probe = "open_failed"
            error = f"{type(exc).__name__}: {exc}"
        rows.append({
            "path": str(path),
            "relative_path": str(path.relative_to(ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": digest or "not_computed_use_--hash-all-candidates",
            "dims": dims,
            "dims_shape_json": json.dumps(shape, sort_keys=True),
            "is_mosaic": mosaic,
            "probe_status": probe,
            "probe_error": error,
            "selection_status": status,
            "exclusion_or_selection_reason": reason,
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_vector_scalebar(ax: plt.Axes, width_um: float, height_um: float) -> None:
    x0 = width_um * 0.055
    y0 = height_um * 0.955
    ax.plot([x0, x0 + SCALEBAR_UM], [y0, y0], color="black", lw=5.0, solid_capstyle="butt", zorder=9)
    ax.plot([x0, x0 + SCALEBAR_UM], [y0, y0], color="white", lw=2.6, solid_capstyle="butt", zorder=10)
    ax.text(
        x0 + SCALEBAR_UM / 2, y0 - height_um * 0.025, f"{int(SCALEBAR_UM)} µm",
        color="white", fontsize=7.5, fontweight="bold", ha="center", va="bottom", zorder=10,
        path_effects=[],
    )


def make_qc_montage(
    outdir: Path,
    displays: dict[str, dict[str, np.ndarray]],
    geometries: dict[str, dict[str, float]],
) -> list[Path]:
    fig, axes = plt.subplots(3, 3, figsize=(8.6, 13.0), constrained_layout=True)
    for row, field in enumerate(FIELDS):
        geometry = geometries[field.condition]
        width_um = geometry["field_width_um"]
        height_um = geometry["field_height_um"]
        for column, marker in enumerate(("merge", "DAPI", "cFOS")):
            ax = axes[row, column]
            ax.imshow(
                displays[field.condition][marker],
                extent=[0, width_um, height_um, 0],
                interpolation="nearest",
                aspect="equal",
            )
            ax.set_xlim(0, width_um)
            ax.set_ylim(height_um, 0)
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                ax.set_title({"merge": "DAPI + c-FOS", "DAPI": "DAPI", "cFOS": "c-FOS"}[marker], fontweight="bold")
            if column == 0:
                ax.set_ylabel(f"{field.condition}\n{field.animal}", fontweight="bold")
            add_vector_scalebar(ax, width_um, height_um)
    fig.suptitle(
        "Figure 2C · native RAW CZI reconstruction\n"
        "common 600 × 1200 µm field · linear percentile display · gamma 1 · no CLAHE",
        fontsize=12, fontweight="bold",
    )
    png = outdir / "QC_Figure2C_raw_reconstruction.png"
    pdf = outdir / "QC_Figure2C_raw_reconstruction.pdf"
    fig.savefig(png, dpi=220, facecolor="white", metadata={"Software": "Fig2/03_reconstruct_microscopy_engine.py"})
    fig.savefig(
        pdf, dpi=220, facecolor="white",
        metadata={"Title": "Figure 2C native RAW CZI reconstruction", "CreationDate": None, "ModDate": None},
    )
    plt.close(fig)
    return [png, pdf]


def output_record(path: Path, role: str, condition: str = "") -> dict[str, object]:
    record: dict[str, object] = {
        "condition": condition,
        "role": role,
        "path": str(path),
        "relative_path": str(path.relative_to(DEFAULT_OUTDIR.parent)) if DEFAULT_OUTDIR.parent in path.parents else str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }
    if path.suffix.casefold() in {".png", ".tif", ".tiff"}:
        if path.suffix.casefold() == ".png":
            with Image.open(path) as image:
                record["shape"] = f"{image.height}x{image.width}x{len(image.getbands())}"
                record["dtype"] = "uint8"
        else:
            with tifffile.TiffFile(path) as tif:
                page = tif.pages[0]
                record["shape"] = "x".join(map(str, page.shape))
                record["dtype"] = str(page.dtype)
    return record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--hash-all-candidates", action="store_true",
        help="SHA-256 every relevant native candidate/copy, not only selected and validation sources.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    for field in FIELDS:
        for path in (field.raw_czi, field.display_reference_tiff):
            if not path.is_file():
                raise FileNotFoundError(path)
        if field.processed_czi_reference is not None and not field.processed_czi_reference.is_file():
            raise FileNotFoundError(field.processed_czi_reference)

    field_manifest: list[dict[str, object]] = []
    source_manifest: list[dict[str, object]] = []
    output_paths: list[tuple[Path, str, str]] = []
    validations: dict[str, dict[str, object]] = {}
    displays: dict[str, dict[str, np.ndarray]] = {}
    geometries: dict[str, dict[str, float]] = {}
    target_hashes: dict[Path, str] = {}

    for field in FIELDS:
        raw_hash = sha256(field.raw_czi)
        target_hashes[field.raw_czi] = raw_hash
        czi = CziFile(field.raw_czi)
        metadata = czi_metadata(czi)
        if not metadata["is_mosaic"]:
            raise RuntimeError(f"Selected source is not a mosaic CZI: {field.raw_czi}")
        pixel_size_um = float(metadata["pixel_size_um"])
        if not np.isclose(pixel_size_um, field.expected_pixel_size_um, rtol=0, atol=1e-12):
            raise RuntimeError(
                f"{field.condition} pixel calibration changed: {pixel_size_um} "
                f"!= {field.expected_pixel_size_um}"
            )
        channels = metadata["channels"]
        if field.dapi_channel not in channels or field.cfos_channel not in channels:
            raise RuntimeError(f"Missing selected channel in {field.raw_czi}")
        if channels[field.dapi_channel]["name"] != "DAPI":
            raise RuntimeError(f"Channel {field.dapi_channel} is not DAPI")
        if "647" not in str(channels[field.cfos_channel]["name"]):
            raise RuntimeError(f"Selected c-FOS channel is not the 647 channel")

        bounds = centered_crop_bounds(field, pixel_size_um)
        channel_crops: dict[str, np.ndarray] = {}
        stitch_records: dict[str, dict[str, object]] = {}
        for marker, channel_index in (("DAPI", field.dapi_channel), ("cFOS", field.cfos_channel)):
            mosaic, stitch = read_scene0_mosaic(czi, channel_index)
            oriented = orient(mosaic, field.rotate_clockwise_quarters)
            if oriented.shape != field.expected_oriented_shape_yx:
                raise RuntimeError(
                    f"{field.condition} oriented shape changed: {oriented.shape} "
                    f"!= {field.expected_oriented_shape_yx}"
                )
            if stitch["tile_count"] != field.expected_scene0_tiles:
                raise RuntimeError(
                    f"{field.condition} tile count changed: {stitch['tile_count']} "
                    f"!= {field.expected_scene0_tiles}"
                )
            x0, y0, x1, y1 = bounds
            if x0 < 0 or y0 < 0 or x1 > oriented.shape[1] or y1 > oriented.shape[0]:
                raise RuntimeError(f"Common field outside {field.condition} mosaic: {bounds}")
            channel_crops[marker] = np.ascontiguousarray(oriented[y0:y1, x0:x1])
            stitch_records[marker] = stitch
            if marker == "DAPI":
                validations[field.condition] = feature_match_validation(
                    oriented, field.display_reference_tiff
                )
            del oriented, mosaic

        # Channel-specific CZI metadata can report slightly different nominal
        # origins. Scientific identity is the registered target coordinate used
        # to copy each tile, which must be exactly the same for both channels.
        placement_keys = ("m", "registered_x_px", "registered_y_px")
        dapi_placements = [
            tuple(record[key] for key in placement_keys)
            for record in stitch_records["DAPI"]["registered_placements"]
        ]
        cfos_placements = [
            tuple(record[key] for key in placement_keys)
            for record in stitch_records["cFOS"]["registered_placements"]
        ]
        if dapi_placements != cfos_placements:
            raise RuntimeError(
                f"DAPI/c-FOS registration transforms differ for {field.condition}"
            )

        if channel_crops["DAPI"].shape != channel_crops["cFOS"].shape:
            raise RuntimeError(f"Paired raw channel shape mismatch for {field.condition}")
        height_px, width_px = channel_crops["DAPI"].shape
        width_um = width_px * pixel_size_um
        height_um = height_px * pixel_size_um
        if abs(width_um - COMMON_WIDTH_UM) > pixel_size_um / 2 + 1e-9:
            raise RuntimeError(f"Common width calibration failed for {field.condition}")
        if abs(height_um - COMMON_HEIGHT_UM) > pixel_size_um / 2 + 1e-9:
            raise RuntimeError(f"Common height calibration failed for {field.condition}")

        raw_dapi = outdir / f"{field.condition}_DAPI_native16.tif"
        raw_cfos = outdir / f"{field.condition}_cFOS_native16.tif"
        save_raw_tiff(raw_dapi, channel_crops["DAPI"], pixel_size_um)
        save_raw_tiff(raw_cfos, channel_crops["cFOS"], pixel_size_um)
        output_paths.extend([
            (raw_dapi, "native_scene0_common_crop_DAPI", field.condition),
            (raw_cfos, "native_scene0_common_crop_cFOS", field.condition),
        ])

        dapi_signal, dapi_params = linear_display(channel_crops["DAPI"], "DAPI")
        cfos_signal, cfos_params = linear_display(channel_crops["cFOS"], "cFOS")
        if dapi_signal.shape != cfos_signal.shape:
            raise RuntimeError(f"Display channel shape mismatch for {field.condition}")
        display_stride = int(dapi_params["display_sampling_stride_px"])
        if display_stride != int(cfos_params["display_sampling_stride_px"]):
            raise RuntimeError(f"Display strides differ for {field.condition}")
        dapi_rgb = rgb_u8(dapi_signal, DAPI_RGB)
        cfos_rgb = rgb_u8(cfos_signal, CFOS_RGB)
        merge_rgb = np.clip(
            dapi_rgb.astype(np.uint16) + cfos_rgb.astype(np.uint16), 0, 255
        ).astype(np.uint8)
        displays[field.condition] = {"DAPI": dapi_rgb, "cFOS": cfos_rgb, "merge": merge_rgb}

        display_pixel_size_um = pixel_size_um * display_stride
        for marker, image in (("DAPI", dapi_rgb), ("cFOS", cfos_rgb), ("merge", merge_rgb)):
            for path in save_display_pair(
                outdir / f"{field.condition}_{marker}_display", image, display_pixel_size_um
            ):
                output_paths.append((path, f"lossless_display_{marker}", field.condition))

        geometry = {
            "field_width_px": width_px,
            "field_height_px": height_px,
            "field_width_um": width_um,
            "field_height_um": height_um,
            "native_pixel_size_um": pixel_size_um,
            "display_sampling_stride_px": display_stride,
            "display_pixel_size_um": display_pixel_size_um,
            "scalebar_um": SCALEBAR_UM,
            "scalebar_native_px": SCALEBAR_UM / pixel_size_um,
            "scalebar_display_px": SCALEBAR_UM / display_pixel_size_um,
        }
        geometries[field.condition] = geometry
        field_manifest.append({
            "condition": field.condition,
            "animal": field.animal,
            "source_path": str(field.raw_czi),
            "source_sha256": raw_hash,
            "source_size_bytes": field.raw_czi.stat().st_size,
            "source_dims": metadata["dims"],
            "source_dims_shape_json": json.dumps(metadata["dims_shape"], sort_keys=True),
            "selected_scene": 0,
            "excluded_scene": 1,
            "scene_selection_basis": "scene 0 is the established multi-tile legacy field; scene 1 is a single tile and is not overlaid",
            "z_projection": "identity; selected scene has one optical plane (no Z dimension)",
            "dapi_channel_index": field.dapi_channel,
            "dapi_channel_metadata_json": json.dumps(channels[field.dapi_channel], sort_keys=True),
            "cfos_channel_index": field.cfos_channel,
            "cfos_channel_metadata_json": json.dumps(channels[field.cfos_channel], sort_keys=True),
            "cfos_assignment_evidence": field.cfos_assignment_evidence,
            "scene0_tile_count": field.expected_scene0_tiles,
            "rotation_clockwise_degrees": 90 * field.rotate_clockwise_quarters,
            "oriented_shape_yx": "x".join(map(str, field.expected_oriented_shape_yx)),
            "crop_x0_px": bounds[0],
            "crop_y0_px": bounds[1],
            "crop_x1_exclusive_px": bounds[2],
            "crop_y1_exclusive_px": bounds[3],
            "dapi_stitch_record_json": json.dumps(stitch_records["DAPI"], sort_keys=True),
            "cfos_stitch_record_json": json.dumps(stitch_records["cFOS"], sort_keys=True),
            "same_integer_registration_dapi_cfos": True,
            "registration_resampling": False,
            "registration_blending": False,
            "dapi_crop_nominal_raw_fallback_pixels": stitch_records["DAPI"]["crop_nominal_raw_fallback_pixels"],
            "cfos_crop_nominal_raw_fallback_pixels": stitch_records["cFOS"]["crop_nominal_raw_fallback_pixels"],
            **geometry,
            "dapi_display_parameters_json": json.dumps(dapi_params, sort_keys=True),
            "cfos_display_parameters_json": json.dumps(cfos_params, sort_keys=True),
            "display_processing": "per-image linear percentile stretch; additive pseudocolor merge; gamma=1; no CLAHE; no spatial filter",
            "quantification_use": "none_display_only",
            "reference_match_status": validations[field.condition]["status"],
            "reference_match_inliers": validations[field.condition]["ransac_inliers"],
        })
        source_manifest.append({
            "condition": field.condition,
            "role": "selected_native_RAW_CZI",
            "path": str(field.raw_czi),
            "sha256": raw_hash,
        })
        source_manifest.append({
            "condition": field.condition,
            "role": "display_reference_validation_only",
            "path": str(field.display_reference_tiff),
            "sha256": validations[field.condition]["reference_sha256"],
        })
        if field.processed_czi_reference is not None:
            processed_hash = sha256(field.processed_czi_reference)
            target_hashes[field.processed_czi_reference] = processed_hash
            source_manifest.append({
                "condition": field.condition,
                "role": "processed_CZI_lineage_validation_only",
                "path": str(field.processed_czi_reference),
                "sha256": processed_hash,
            })

    qc_paths = make_qc_montage(outdir, displays, geometries)
    output_paths.extend((path, "QC_montage_with_vector_200um_bars", "") for path in qc_paths)

    candidate_rows = audit_candidates(target_hashes, args.hash_all_candidates)
    selected_rows = [row for row in candidate_rows if row["selection_status"] == "selected_native_source"]
    if len(selected_rows) != 3:
        raise RuntimeError(f"Candidate audit found {len(selected_rows)} selected native sources, expected 3")
    write_csv(outdir / "candidate_container_audit.csv", candidate_rows)
    write_csv(outdir / "field_reconstruction_manifest.csv", field_manifest)
    write_csv(outdir / "source_hash_manifest.csv", source_manifest)

    parameters = {
        "script": str(HERE),
        "script_sha256": sha256(HERE),
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "common_field_um": [COMMON_WIDTH_UM, COMMON_HEIGHT_UM],
        "linear_percentiles": {
            "low": LOW_PERCENTILE,
            "DAPI_high": HIGH_PERCENTILE["DAPI"],
            "cFOS_high": HIGH_PERCENTILE["cFOS"],
        },
        "gamma": 1.0,
        "clahe": False,
        "spatial_filtering": False,
        "tile_stitching": (
            "seam-derived integer per-tile translations; identical DAPI/c-FOS "
            "transform; exact nominal RAW ownership only at no-data junctions; "
            "no interpolation, blending, or overlap averaging"
        ),
        "display_colors_rgb": {"DAPI": DAPI_RGB.tolist(), "cFOS": CFOS_RGB.tolist()},
        "scale_bar_um": SCALEBAR_UM,
        "statistics_dependency": "none; display reconstruction is isolated from quantitative analysis",
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "opencv": cv2.__version__,
        "tifffile": tifffile.__version__,
    }
    (outdir / "parameters.json").write_text(
        json.dumps(parameters, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    validation = {
        "status": "passed",
        "conditions": validations,
        "common_physical_field_validated": True,
        "native_calibration_validated": True,
        "channel_identity_validated": True,
        "scene0_only_validated": True,
        "seam_registration_validated": True,
        "same_integer_transform_for_dapi_cfos": True,
        "registration_resampling": False,
        "registration_blending": False,
        "vector_scalebar_um": SCALEBAR_UM,
        "scale_bar_calculation": "length_px = 200 um / calibrated pixel_size_um",
        "display_only_not_used_for_statistics": True,
    }
    (outdir / "validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    output_records = [output_record(path, role, condition) for path, role, condition in output_paths]
    for path, role in (
        (outdir / "candidate_container_audit.csv", "candidate_audit"),
        (outdir / "field_reconstruction_manifest.csv", "field_manifest"),
        (outdir / "source_hash_manifest.csv", "source_hash_manifest"),
        (outdir / "parameters.json", "parameters"),
        (outdir / "validation.json", "validation"),
    ):
        output_records.append(output_record(path, role))
    write_csv(outdir / "output_hash_manifest.csv", output_records)

    candidate_summary: dict[str, int] = {}
    for row in candidate_rows:
        key = str(row["selection_status"])
        candidate_summary[key] = candidate_summary.get(key, 0) + 1
    report = [
        "# Figure 2C native-CZI reconstruction",
        "",
        "Panel C is reconstructed from the original 16-bit RAW CZI scene-0 tiles; no PPTX, legacy PDF, QC PNG, or 8-bit display TIFF contributes pixels.",
        "Tiles use audited seam-derived integer translations applied identically to DAPI and c-FOS, with no interpolation or blending; only translated-grid no-data junctions use the exact RAW pixel at its nominal acquisition ownership.",
        "The earlier TIFFs are used only for SIFT/RANSAC field/orientation validation.",
        "",
        "## Display policy",
        "",
        f"- Common field: {COMMON_WIDTH_UM:.0f} × {COMMON_HEIGHT_UM:.0f} µm for every condition.",
        f"- Linear stretch: P{LOW_PERCENTILE:g} to P{HIGH_PERCENTILE['DAPI']:g} (DAPI) and P{HIGH_PERCENTILE['cFOS']:g} (c-FOS).",
        "- Gamma 1; no CLAHE; no denoising, sharpening, local equalization, or spatial background subtraction.",
        "- Display images are not used for cell counts or statistical analysis.",
        f"- The QC PDF contains calibrated vector {SCALEBAR_UM:.0f}-µm bars.",
        "",
        "## Native fields",
        "",
    ]
    for field in FIELDS:
        record = next(row for row in field_manifest if row["condition"] == field.condition)
        report.extend([
            f"- **{field.condition} ({field.animal})**: scene 0, {record['scene0_tile_count']} tiles, "
            f"DAPI C{record['dapi_channel_index']}, c-FOS C{record['cfos_channel_index']}, "
            f"{record['native_pixel_size_um']:.9f} µm/px, {record['rotation_clockwise_degrees']}° clockwise.",
        ])
    report.extend([
        "",
        "## Candidate audit",
        "",
        *(f"- {key}: {value}" for key, value in sorted(candidate_summary.items())),
        "",
        "See `candidate_container_audit.csv` for every reviewed native candidate and its exclusion reason; "
        "`source_hash_manifest.csv`, `field_reconstruction_manifest.csv`, `parameters.json`, and "
        "`output_hash_manifest.csv` make the reconstruction auditable.",
        "",
    ])
    (outdir / "README.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Reconstructed Figure 2C from native RAW CZI data: {outdir}")
    print(f"Candidate CZI containers audited: {len(candidate_rows)}")
    for field in FIELDS:
        print(
            f"  {field.condition}: {geometries[field.condition]['field_width_um']:.3f} × "
            f"{geometries[field.condition]['field_height_um']:.3f} µm; "
            f"SIFT inliers={validations[field.condition]['ransac_inliers']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
