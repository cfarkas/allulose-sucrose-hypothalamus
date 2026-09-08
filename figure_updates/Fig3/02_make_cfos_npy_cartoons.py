#!/usr/bin/env python3
"""Render matched Water/Allulose Figure 3 marker-positive-nuclei cartoons.

Panels B and C use the same three-view grammar and physical geometry. Native
DAPI/c-FOS/NPY Cellpose labels are classified before deterministic registered
display transforms. Both panels are nearest-neighbor registered into one exact physical common frame
using fixed homographies derived from their Panel A views. Black lumen boundaries
are display annotations only: both retain hash-bound raw-DAPI base components,
Water adds a deterministic superior completion from registered DAPI-label
clearance, and the final lumen is supplied by an accepted, source-bound HIL mask.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from uuid import uuid4

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Polygon, Rectangle
from scipy.ndimage import (binary_fill_holes, distance_transform_edt, gaussian_filter,
                           gaussian_filter1d, label)
from skimage.measure import regionprops
from skimage.morphology import remove_small_objects

SAMPLE_ID = "water-npy2-0001"
SOURCE_STEM = "WATER_NPY2-0001"
ANIMAL_ID = "WATER_NPY2"
CONDITION = "Water"
NATIVE_SHAPE_YX = (3749, 3768)
SOURCE_BOUNDS_XYXY = (188, 1223, 1270, 2903)
DISPLAY_SHAPE_YX = (1082, 1680)
DISPLAY_ROT90_K = 1
# Homogeneous output-display (x,y,1) -> native-source (x,y,1).
DISPLAY_TO_NATIVE_XY = ((0, -1, 1269), (1, 0, 1223), (0, 0, 1))
PX_PER_UM = 1.768592
UM_PER_PX = 1.0 / PX_PER_UM

CHANNEL_COLORS = {"DAPI": "#4f86ff", "NPY": "#00e85e", "c-FOS": "#ff37d4"}
TISSUE_FILL = "#e8edf2"
TISSUE_EDGE = "#354052"
NEGATIVE_NUCLEUS_COLOR = "#9eb9ef"
DUAL_COLOR = "#f2a900"
CFOS_MIN_ROI_AREA_PX = 5
# Gaussian sigma, in boundary points, used only to smooth the ventricular
# display outline. The traced boundary is one point per pixel step, so this is
# large enough to remove the per-pixel threshold staircase and small enough to
# keep the third ventricle's own shape.
LUMEN_DISPLAY_SMOOTHING_POINTS = 12.0
VENTRICLE_HIL_SCHEMA = "fig3_cartoon_ventricle_hil_v1"
VENTRICLE_HIL_DEFAULT_RELATIVE = Path("Fig3/hil_review/ventricle_cartoon_20260827_v1")
CFOS_ROI_TO_NUCLEUS_OVERLAP = 0.40
NPY_MIN_MARKER_PIXELS = 1
NPY_MIN_NUCLEUS_FRACTION = 0.05

# Display-only ventricular-lumen annotation. All coordinates are in the exact
# 1082x1680 Water2 registered display view, after the integer crop/rotation.
LUMEN_SEARCH_BOUNDS_XYXY = (600, 400, 1080, 850)
LUMEN_WINDOW_RADIUS_PX = 32
LUMEN_WINDOW_SIZE_PX = 2 * LUMEN_WINDOW_RADIUS_PX + 1
LUMEN_MAX_DAPI_PIXELS_PER_WINDOW = 84
LUMEN_MIN_COMPONENT_AREA_PX = 4_000
EXPECTED_LUMEN_AREA_PX = 10_626
EXPECTED_LUMEN_BBOX_YXYX = (482, 766, 706, 922)
EXPECTED_LUMEN_CENTROID_YX = (609.0617353660832, 844.0355731225296)
EXPECTED_LUMEN_BINARY_SHA256 = "e3e37c56ec5c0c33614a02209ae0facd973ff777d7f6a574c487963f51ca05f7"

EXPECTED_SOURCE_SHA256 = {
    "dapi_tiff": "d5d40182e3d4bf396ef1277815f851d111aa90bfba3d6b2bf7083080c9657044",
    "npy_tiff": "aa1fd314b875032d648f0eb068842353537c7fe00b2ba7e0bf33a99e665ef6e7",
    "cfos_tiff": "d52c02775b283a12379b579c061d13a478543bb0c21f65342d5efefd2347f03c",
    "dapi_seg": "fb77993b5d928e6f9fc19043e076554729db5721203f1780feba9675d2f12c0b",
    "npy_seg": "5af0388dea2894f3d20b7794ae424bb5f66a9dd682d794744478fabcd560819d",
    "cfos_seg": "b6e82a6c85966dba752dca04bb8e21d9a05027b44c1b3fbf465cbfc5dcfea85d",
}
EXPECTED_INSTANCE_COUNTS = {"DAPI": 11608, "NPY": 177, "c-FOS": 822}
EXPECTED_ANALYSIS_COUNTS = {
    "accepted_cfos_rois": 774,
    "duplicate_cfos_to_nucleus_mappings": 1,
    "cfos_positive_nuclei": 773,
    "npy_positive_nuclei": 208,
    "dual_positive_nuclei": 57,
}
PANEL_B_STEM = "Figure3_Panel_B_DAPI_NPY_cFOS_processed_masks"
PANEL_C_STEM = "Figure3_Panel_C_marker_positive_nuclei"
PROVENANCE_STEM = "Figure3_cFOS_NPY_cartoon_provenance"


@dataclass(frozen=True)
class SourceFile:
    role: str
    path: Path
    sha256: str
    bytes: int


@dataclass(frozen=True)
class CartoonData:
    root: Path
    dapi_labels: np.ndarray
    npy_labels: np.ndarray
    cfos_labels: np.ndarray
    sources: Mapping[str, SourceFile]
    internal_filenames: Mapping[str, str]


@dataclass(frozen=True)
class AnalysisStatus:
    accepted_cfos_roi_ids: tuple[int, ...]
    cfos_roi_to_nucleus: Mapping[int, int]
    cfos_positive_nucleus_ids: frozenset[int]
    npy_positive_nucleus_ids: frozenset[int]
    dual_positive_nucleus_ids: frozenset[int]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def binary_sha256(mask: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(mask, dtype=np.uint8).tobytes()).hexdigest()


def _assert_source(path: Path, role: str) -> SourceFile:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {role}: {path}")
    observed = sha256_file(path)
    expected = EXPECTED_SOURCE_SHA256[role]
    if observed != expected:
        raise RuntimeError(f"{role} SHA256 mismatch: expected {expected}, got {observed}")
    return SourceFile(role, path.resolve(), observed, path.stat().st_size)


def _validate_tiff(path: Path, role: str) -> None:
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) != 1:
            raise RuntimeError(f"{role} must be a single-page TIFF")
        page = tif.pages[0]
        if tuple(page.shape) != (*NATIVE_SHAPE_YX, 3) or page.dtype != np.dtype("uint8"):
            raise RuntimeError(
                f"{role} is {page.shape}/{page.dtype}; expected RGB uint8 at {NATIVE_SHAPE_YX}"
            )
        for tag_name in ("XResolution", "YResolution"):
            tag = page.tags.get(tag_name)
            if tag is None:
                raise RuntimeError(f"{role} has no {tag_name} tag")
            observed = float(tag.value[0]) / float(tag.value[1])
            if abs(observed - PX_PER_UM) > 1e-12:
                raise RuntimeError(f"{role} {tag_name} changed to {observed}")


def _load_cellpose(path: Path, role: str, expected_instances: int) -> tuple[np.ndarray, str]:
    _assert_source(path, role)
    loaded = np.load(path, allow_pickle=True)
    if not (isinstance(loaded, np.ndarray) and loaded.dtype == object and loaded.shape == ()):
        raise RuntimeError(f"{role} is not a scalar Cellpose object dictionary")
    payload = loaded.item()
    if not isinstance(payload, dict) or "masks" not in payload:
        raise RuntimeError(f"{role} has no Cellpose masks entry")
    labels = np.asarray(payload["masks"])
    if tuple(labels.shape) != NATIVE_SHAPE_YX:
        raise RuntimeError(f"{role} shape is {labels.shape}; expected {NATIVE_SHAPE_YX}")
    if not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0):
        raise RuntimeError(f"{role} labels must be nonnegative integers")
    positive = np.unique(labels)
    positive = positive[positive > 0]
    if not np.array_equal(positive, np.arange(1, positive.size + 1, dtype=positive.dtype)):
        raise RuntimeError(f"{role} labels are not contiguous")
    if positive.size != expected_instances:
        raise RuntimeError(f"{role} has {positive.size} instances; expected {expected_instances}")
    internal = str(payload.get("filename", "")).strip()
    if not internal:
        raise RuntimeError(f"{role} has no internal source filename")
    return labels, internal


def load_water2_inputs(root: Path | str | None = None) -> CartoonData:
    project_root = Path(root).expanduser().resolve() if root else next(path for path in Path(__file__).resolve().parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
    source_dir = project_root / "Fig3/raw/legacy_experiment_2025_07_28/water"
    paths = {
        "dapi_tiff": source_dir / f"{SOURCE_STEM}_dapi.tif",
        "npy_tiff": source_dir / f"{SOURCE_STEM}_GFP.tif",
        "cfos_tiff": source_dir / f"{SOURCE_STEM}_fos.tif",
        "dapi_seg": source_dir / f"{SOURCE_STEM}_dapi_seg.npy",
        "npy_seg": source_dir / f"{SOURCE_STEM}_GFP_seg.npy",
        "cfos_seg": source_dir / f"{SOURCE_STEM}_fos_seg.npy",
    }
    sources = {role: _assert_source(path, role) for role, path in paths.items()}
    for role in ("dapi_tiff", "npy_tiff", "cfos_tiff"):
        _validate_tiff(paths[role], role)
    dapi, dapi_internal = _load_cellpose(
        paths["dapi_seg"], "dapi_seg", EXPECTED_INSTANCE_COUNTS["DAPI"]
    )
    npy, npy_internal = _load_cellpose(
        paths["npy_seg"], "npy_seg", EXPECTED_INSTANCE_COUNTS["NPY"]
    )
    cfos, cfos_internal = _load_cellpose(
        paths["cfos_seg"], "cfos_seg", EXPECTED_INSTANCE_COUNTS["c-FOS"]
    )
    if not (dapi.shape == npy.shape == cfos.shape == NATIVE_SHAPE_YX):
        raise RuntimeError("Water2 masks do not share exact native coordinates")
    return CartoonData(
        project_root,
        dapi,
        npy,
        cfos,
        sources,
        {"DAPI": dapi_internal, "NPY": npy_internal, "c-FOS": cfos_internal},
    )


def classify_analysis_status(
    dapi_labels: np.ndarray,
    cfos_labels: np.ndarray,
    npy_labels: np.ndarray,
) -> AnalysisStatus:
    """Classify the complete native fields before any display operation."""
    dapi = np.asarray(dapi_labels)
    cfos = np.asarray(cfos_labels)
    npy = np.asarray(npy_labels)
    if not (dapi.shape == cfos.shape == npy.shape == NATIVE_SHAPE_YX):
        raise RuntimeError("Classification inputs are not complete native fields")
    accepted: list[int] = []
    mapping: dict[int, int] = {}
    for prop in regionprops(cfos):
        if int(prop.area) < CFOS_MIN_ROI_AREA_PX:
            continue
        values = dapi[prop.coords[:, 0], prop.coords[:, 1]]
        counts = np.bincount(values, minlength=EXPECTED_INSTANCE_COUNTS["DAPI"] + 1)
        counts[0] = 0
        best = int(np.argmax(counts))
        fraction = float(counts[best]) / float(prop.area) if best else 0.0
        if best and fraction >= CFOS_ROI_TO_NUCLEUS_OVERLAP:
            accepted.append(int(prop.label))
            mapping[int(prop.label)] = best
    npy_mask = npy > 0
    npy_positive: set[int] = set()
    for prop in regionprops(dapi):
        overlap = int(npy_mask[prop.coords[:, 0], prop.coords[:, 1]].sum())
        fraction = overlap / float(max(1, prop.area))
        if overlap >= NPY_MIN_MARKER_PIXELS and fraction >= NPY_MIN_NUCLEUS_FRACTION:
            npy_positive.add(int(prop.label))
    cfos_positive = set(mapping.values())
    dual = cfos_positive & npy_positive
    observed = {
        "accepted_cfos_rois": len(accepted),
        "duplicate_cfos_to_nucleus_mappings": len(accepted) - len(cfos_positive),
        "cfos_positive_nuclei": len(cfos_positive),
        "npy_positive_nuclei": len(npy_positive),
        "dual_positive_nuclei": len(dual),
    }
    if observed != EXPECTED_ANALYSIS_COUNTS:
        raise RuntimeError(
            f"Water2 whole-field counts changed: expected {EXPECTED_ANALYSIS_COUNTS}, got {observed}"
        )
    return AnalysisStatus(
        tuple(sorted(accepted)),
        dict(sorted(mapping.items())),
        frozenset(cfos_positive),
        frozenset(npy_positive),
        frozenset(dual),
    )


def display_transform(array: np.ndarray) -> np.ndarray:
    """Apply native[y=1223:2903,x=188:1270], then np.rot90(k=1)."""
    values = np.asarray(array)
    if tuple(values.shape[:2]) != NATIVE_SHAPE_YX:
        raise RuntimeError(f"Display input starts with {values.shape[:2]}; expected {NATIVE_SHAPE_YX}")
    x0, y0, x1, y1 = SOURCE_BOUNDS_XYXY
    displayed = np.rot90(values[y0:y1, x0:x1, ...], k=DISPLAY_ROT90_K, axes=(0, 1))
    if tuple(displayed.shape[:2]) != DISPLAY_SHAPE_YX:
        raise RuntimeError(f"Display output is {displayed.shape[:2]}; expected {DISPLAY_SHAPE_YX}")
    return displayed


def _tissue_envelope(display_dapi: np.ndarray) -> np.ndarray:
    positive = display_dapi > 0
    density = gaussian_filter(positive.astype(np.float32), sigma=27.05, mode="constant")
    tissue = remove_small_objects(density > 0.02, min_size=3636)
    if not np.any(tissue) or np.all(tissue):
        raise RuntimeError("Could not derive a bounded DAPI-supported display envelope")
    return np.asarray(tissue, dtype=bool)


def _integral_window_sum(binary: np.ndarray, radius: int) -> np.ndarray:
    """Return exact centered square-window sums without floating-point filtering."""
    values = np.asarray(binary, dtype=np.int64)
    if values.ndim != 2:
        raise RuntimeError("Lumen occupancy input must be one registered 2D mask")
    radius = int(radius)
    if radius < 1:
        raise RuntimeError("Lumen occupancy radius must be positive")
    width = 2 * radius + 1
    padded = np.pad(values, ((radius, radius), (radius, radius)), mode="constant")
    integral = np.pad(padded, ((1, 0), (1, 0)), mode="constant")
    integral = integral.cumsum(axis=0, dtype=np.int64).cumsum(axis=1, dtype=np.int64)
    sums = (
        integral[width:, width:]
        - integral[:-width, width:]
        - integral[width:, :-width]
        + integral[:-width, :-width]
    )
    if sums.shape != values.shape:
        raise RuntimeError(f"Lumen occupancy output shape changed: {sums.shape} != {values.shape}")
    return sums


def ventricular_lumen_mask(display_dapi: np.ndarray) -> np.ndarray:
    """Derive the audited Water2 central lumen as a display-only annotation.

    The rule uses only the hash-bound registered DAPI Cellpose mask: a 65x65
    integer occupancy window, a fixed central search box, and 8-connectivity.
    It is not an anatomical analysis ROI and is never used in classification or
    a denominator. Exact geometry checks make source/crop drift fail closed.
    """
    dapi = np.asarray(display_dapi)
    if dapi.shape != DISPLAY_SHAPE_YX:
        raise RuntimeError(f"Lumen DAPI shape is {dapi.shape}; expected {DISPLAY_SHAPE_YX}")
    occupancy = _integral_window_sum(dapi > 0, LUMEN_WINDOW_RADIUS_PX)
    x0, y0, x1, y1 = LUMEN_SEARCH_BOUNDS_XYXY
    search = np.zeros(DISPLAY_SHAPE_YX, dtype=bool)
    search[y0:y1, x0:x1] = True
    candidate = search & (occupancy <= LUMEN_MAX_DAPI_PIXELS_PER_WINDOW)
    components, component_count = label(candidate, structure=np.ones((3, 3), dtype=np.uint8))
    accepted: list[np.ndarray] = []
    for component_id in range(1, int(component_count) + 1):
        component = components == component_id
        ys, xs = np.where(component)
        if ys.size < LUMEN_MIN_COMPONENT_AREA_PX:
            continue
        touches_search_edge = bool(
            np.any(ys == y0) or np.any(ys == y1 - 1)
            or np.any(xs == x0) or np.any(xs == x1 - 1)
        )
        if not touches_search_edge:
            accepted.append(component)
    if len(accepted) != 1:
        raise RuntimeError(
            "Expected one enclosed DAPI-sparse Water2 lumen component; "
            f"found {len(accepted)} from {component_count} candidates"
        )
    lumen = np.asarray(accepted[0], dtype=bool)
    ys, xs = np.where(lumen)
    bbox = (int(ys.min()), int(xs.min()), int(ys.max()) + 1, int(xs.max()) + 1)
    centroid = (float(ys.mean()), float(xs.mean()))
    observed_hash = binary_sha256(lumen)
    if int(lumen.sum()) != EXPECTED_LUMEN_AREA_PX or bbox != EXPECTED_LUMEN_BBOX_YXYX:
        raise RuntimeError(
            "Water2 lumen geometry changed: "
            f"area={int(lumen.sum())}, bbox={bbox}; expected "
            f"area={EXPECTED_LUMEN_AREA_PX}, bbox={EXPECTED_LUMEN_BBOX_YXYX}"
        )
    if not np.allclose(centroid, EXPECTED_LUMEN_CENTROID_YX, atol=1e-9, rtol=0.0):
        raise RuntimeError(
            f"Water2 lumen centroid changed: {centroid} != {EXPECTED_LUMEN_CENTROID_YX}"
        )
    if observed_hash != EXPECTED_LUMEN_BINARY_SHA256:
        raise RuntimeError(
            f"Water2 lumen binary hash changed: {observed_hash} != {EXPECTED_LUMEN_BINARY_SHA256}"
        )
    return lumen


def _display_layers(
    data: CartoonData,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dapi = display_transform(data.dapi_labels)
    npy = display_transform(data.npy_labels)
    cfos = display_transform(data.cfos_labels)
    lumen = ventricular_lumen_mask(dapi)
    return dapi, npy, cfos, _tissue_envelope(dapi), lumen


def _visible_ids(labels: np.ndarray) -> frozenset[int]:
    return frozenset(int(value) for value in np.unique(labels) if int(value) > 0)


def _selected_mask(labels: np.ndarray, ids: Iterable[int]) -> np.ndarray:
    ids_array = np.asarray(sorted(set(int(value) for value in ids)), dtype=np.int32)
    if ids_array.size and (ids_array.min() < 1 or ids_array.max() > EXPECTED_INSTANCE_COUNTS["DAPI"]):
        raise RuntimeError("Selected nucleus ID is outside the native DAPI range")
    lookup = np.zeros(EXPECTED_INSTANCE_COUNTS["DAPI"] + 1, dtype=bool)
    lookup[ids_array] = True
    return lookup[labels]


def _visible_status_count(display_dapi: np.ndarray, ids: Iterable[int]) -> int:
    return len(_visible_ids(display_dapi) & frozenset(int(value) for value in ids))


def _draw_background(ax: plt.Axes, tissue: np.ndarray, lumen: np.ndarray) -> None:
    if tissue.shape != lumen.shape:
        raise RuntimeError("Cartoon tissue and lumen masks do not share registered coordinates")
    rgb = np.ones((*tissue.shape, 3), dtype=np.float32)
    # Preserve the original tissue mask for background fill. Only the gray
    # display outline is hole-filled, eliminating internal ring paths.
    rgb[tissue & ~lumen] = np.asarray(to_rgb(TISSUE_FILL), dtype=np.float32)
    ax.imshow(rgb, interpolation="nearest")
    tissue_outline = np.asarray(binary_fill_holes(tissue), dtype=bool)
    if _contour_path_counts(tissue_outline) != {"external": 1, "tree": 1}:
        raise RuntimeError("Tissue display outline must contain one external path and no internal path")
    ax.contour(tissue_outline, [0.5], colors=["white"], linewidths=3.5, alpha=0.95)
    ax.contour(tissue_outline, [0.5], colors=[TISSUE_EDGE], linewidths=1.7)


def _smoothed_boundary_polygon(
    mask: np.ndarray, smoothing_points: float = LUMEN_DISPLAY_SMOOTHING_POINTS,
) -> np.ndarray:
    """Trace a single-path display mask and return a smoothed closed polygon.

    The lumen mask comes from a per-pixel DAPI intensity threshold, so its raw
    0.5 contour follows every pixel step and reads as a tortuous line rather
    than an ependymal border, most visibly along the superior third ventricle.
    This traces the one external path and smooths its coordinates with a
    periodic Gaussian, which follows the same boundary without the staircase.

    Smoothing the traced path rather than blurring the mask is deliberate:
    blurring erodes thin structures and pinches the completed superior corridor
    of the Water lumen into more than one component, whereas a closed path
    smoothed periodically is still exactly one closed path. Coordinates are
    clipped back to the frame so a boundary that reaches an edge stays flush
    with it.

    This is display only. The hash-bound source and corrected lumen geometries,
    every analysis ROI and every denominator are computed before this boundary
    smoothing and are untouched by it.
    """
    binary = np.ascontiguousarray(np.asarray(mask, dtype=np.uint8))
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if len(contours) != 1:
        raise RuntimeError(
            f"Ventricular display boundary must trace exactly one path, got {len(contours)}"
        )
    path = np.asarray(contours[0], dtype=np.float64).reshape(-1, 2)
    if path.shape[0] < 3 * int(round(smoothing_points)):
        raise RuntimeError("Ventricular display boundary is too short to smooth")
    smoothed = np.column_stack([
        gaussian_filter1d(path[:, axis], sigma=float(smoothing_points), mode="wrap")
        for axis in (0, 1)
    ])
    height, width = binary.shape
    smoothed[:, 0] = np.clip(smoothed[:, 0], 0.0, float(width - 1))
    smoothed[:, 1] = np.clip(smoothed[:, 1], 0.0, float(height - 1))
    return smoothed


def _draw_lumen_outline(ax: plt.Axes, lumen: np.ndarray) -> None:
    """Make the lumen clean white, then draw exactly one black outer path."""
    lumen_outline = np.asarray(binary_fill_holes(lumen), dtype=bool)
    if _contour_path_counts(lumen_outline) != {"external": 1, "tree": 1}:
        raise RuntimeError("Ventricular display boundary must contain exactly one outer path")
    polygon = _smoothed_boundary_polygon(lumen_outline)
    # This final opaque layer covers DAPI objects and any earlier tissue-outline
    # pixels inside the display-only lumen. It deliberately does not alter the
    # hash-bound lumen geometry used to generate the single outer contour.
    ax.add_patch(Polygon(
        polygon, closed=True, facecolor="white", edgecolor="black",
        linewidth=1.7, joinstyle="round", zorder=28,
    ))


def _add_scalebar(ax: plt.Axes, shape: tuple[int, int]) -> None:
    height, width = shape
    length_px = 100.0 / UM_PER_PX
    x0 = width - max(30.0, 0.026 * width) - length_px
    y0 = height - max(25.0, 0.055 * height)
    bar_h = max(11.0, 0.017 * height)
    pad = max(5.0, 0.025 * length_px)
    ax.add_patch(Rectangle(
        (x0 - pad, y0 - bar_h - 65), length_px + 2 * pad, bar_h + 72,
        facecolor="black", edgecolor="black", linewidth=1.2, alpha=0.92, zorder=30,
    ))
    ax.add_patch(Rectangle(
        (x0, y0 - bar_h), length_px, bar_h,
        facecolor="white", edgecolor="white", linewidth=1.0, zorder=31,
    ))
    label = ax.text(
        x0 + length_px / 2, y0 - bar_h - 4, "100 µm",
        color="white", ha="center", va="bottom", fontsize=6.2,
        fontweight="bold", zorder=32,
    )
    label.set_path_effects([pe.withStroke(linewidth=3.2, foreground="black")])


def _finish_axis(ax: plt.Axes, shape: tuple[int, int], lumen: np.ndarray) -> None:
    _draw_lumen_outline(ax, lumen)
    _add_scalebar(ax, shape)
    ax.set_xlim(-0.5, shape[1] - 0.5)
    ax.set_ylim(shape[0] - 0.5, -0.5)
    ax.set_axis_off()


def _base_figure(
    letter: str,
    title: str,
    *,
    bottom: float,
) -> tuple[plt.Figure, Sequence[plt.Axes]]:
    # The panel assets contain only the cartoons and their essential labels.
    # Method prose belongs in the figure legend/provenance, not below the art.
    fig, axes = plt.subplots(1, 3, figsize=(7.44, 2.55), facecolor="white")
    fig.subplots_adjust(left=0.025, right=0.988, bottom=bottom, top=0.735, wspace=0.045)
    fig.text(0.012, 0.970, letter, fontsize=22, fontweight="bold", ha="left", va="top")
    fig.text(0.500, 0.955, title, fontsize=12.5, fontweight="bold", ha="center", va="top")
    return fig, tuple(axes)


def draw_processed_mask_cartoon(data: CartoonData, panel_letter: str = "B") -> plt.Figure:
    dapi, npy, cfos, tissue, lumen = _display_layers(data)
    fig, axes = _base_figure(
        panel_letter, "Processed masks (Water2 registered view)", bottom=0.035,
    )
    for ax, marker, labels in zip(axes, ("DAPI", "NPY", "c-FOS"), (dapi, npy, cfos)):
        _draw_background(ax, tissue, lumen)
        rgba = np.zeros((*labels.shape, 4), dtype=np.float32)
        rgba[..., :3] = to_rgb(CHANNEL_COLORS[marker])
        rgba[..., 3] = (labels > 0).astype(np.float32) * (0.78 if marker == "DAPI" else 0.90)
        ax.imshow(rgba, interpolation="nearest")
        noun = "nuclei" if marker == "DAPI" else "objects"
        ax.set_title(
            f"{marker} segmentation\ndisplay-field n={len(_visible_ids(labels)):,} {noun}",
            color=CHANNEL_COLORS[marker], fontsize=10.5, fontweight="bold", pad=3,
        )
        _finish_axis(ax, labels.shape, lumen)
    return fig


def _draw_status_axis(
    ax: plt.Axes,
    dapi: np.ndarray,
    tissue: np.ndarray,
    lumen: np.ndarray,
    layers: Sequence[tuple[Iterable[int], str]],
    title: str,
) -> None:
    _draw_background(ax, tissue, lumen)
    base = np.zeros((*dapi.shape, 4), dtype=np.float32)
    base[..., :3] = to_rgb(NEGATIVE_NUCLEUS_COLOR)
    base[..., 3] = (dapi > 0).astype(np.float32) * 0.25
    ax.imshow(base, interpolation="nearest")
    for ids, color in layers:
        rgba = np.zeros((*dapi.shape, 4), dtype=np.float32)
        rgba[..., :3] = to_rgb(color)
        rgba[..., 3] = _selected_mask(dapi, ids).astype(np.float32) * 0.94
        ax.imshow(rgba, interpolation="nearest")
    ax.set_title(title, fontsize=10.5, fontweight="bold", color="#172033", pad=3)
    _finish_axis(ax, dapi.shape, lumen)


def draw_analysis_status_cartoon(
    data: CartoonData,
    status: AnalysisStatus,
    panel_letter: str = "C",
) -> plt.Figure:
    dapi, _, _, tissue, lumen = _display_layers(data)
    fig, axes = _base_figure(
        panel_letter, "Marker-positive DAPI nuclei (Water2 registered view)", bottom=0.120,
    )
    npy_n = _visible_status_count(dapi, status.npy_positive_nucleus_ids)
    cfos_n = _visible_status_count(dapi, status.cfos_positive_nucleus_ids)
    dual_n = _visible_status_count(dapi, status.dual_positive_nucleus_ids)
    _draw_status_axis(
        axes[0], dapi, tissue, lumen,
        ((status.npy_positive_nucleus_ids, CHANNEL_COLORS["NPY"]),),
        f"NPY-positive nuclei\ndisplay-field n={npy_n:,}",
    )
    _draw_status_axis(
        axes[1], dapi, tissue, lumen,
        ((status.cfos_positive_nucleus_ids, CHANNEL_COLORS["c-FOS"]),),
        f"c-FOS-positive nuclei\ndisplay-field n={cfos_n:,}",
    )
    npy_only = status.npy_positive_nucleus_ids - status.dual_positive_nucleus_ids
    cfos_only = status.cfos_positive_nucleus_ids - status.dual_positive_nucleus_ids
    _draw_status_axis(
        axes[2], dapi, tissue, lumen,
        ((npy_only, CHANNEL_COLORS["NPY"]),
         (cfos_only, CHANNEL_COLORS["c-FOS"]),
         (status.dual_positive_nucleus_ids, DUAL_COLOR)),
        f"Combined display-field status\nNPY {npy_n:,} · cFOS {cfos_n:,} · dual {dual_n:,}",
    )
    fig.legend(
        handles=(
            Patch(facecolor=NEGATIVE_NUCLEUS_COLOR, alpha=0.45, label="Other DAPI nuclei"),
            Patch(facecolor=CHANNEL_COLORS["NPY"], label="NPY-positive nucleus"),
            Patch(facecolor=CHANNEL_COLORS["c-FOS"], label="c-FOS-positive nucleus"),
            Patch(facecolor=DUAL_COLOR, label="Dual-positive nucleus"),
        ),
        loc="lower center", bbox_to_anchor=(0.50, 0.002), ncol=4,
        frameon=False, fontsize=7.0, handlelength=1.4, columnspacing=1.4,
    )
    return fig


def _atomic_save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.{uuid4().hex}.tmp{path.suffix}")
    try:
        fig.savefig(
            temporary, format=path.suffix.lstrip("."), dpi=dpi,
            facecolor="white", edgecolor="none",
            metadata={
                "Title": path.stem,
                "Subject": "Figure 3 matched Water/Allulose c-FOS/NPY marker-positive-nuclei cartoon",
                "Creator": Path(__file__).name,
            },
        )
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _atomic_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _visible_counts(data: CartoonData, status: AnalysisStatus) -> tuple[dict[str, int], dict[str, int]]:
    dapi = display_transform(data.dapi_labels)
    mask_counts = {
        "DAPI": len(_visible_ids(dapi)),
        "NPY": len(_visible_ids(display_transform(data.npy_labels))),
        "c-FOS": len(_visible_ids(display_transform(data.cfos_labels))),
    }
    status_counts = {
        "npy_positive_nuclei": _visible_status_count(dapi, status.npy_positive_nucleus_ids),
        "cfos_positive_nuclei": _visible_status_count(dapi, status.cfos_positive_nucleus_ids),
        "dual_positive_nuclei": _visible_status_count(dapi, status.dual_positive_nucleus_ids),
    }
    return mask_counts, status_counts


def _provenance_payload(
    data: CartoonData,
    status: AnalysisStatus,
    outputs: Sequence[Path],
) -> dict[str, Any]:
    visible_masks, visible_status = _visible_counts(data, status)
    labels_by_marker = {
        "DAPI": data.dapi_labels,
        "NPY": data.npy_labels,
        "c-FOS": data.cfos_labels,
    }
    whole_status = {
        "accepted_cfos_rois": len(status.accepted_cfos_roi_ids),
        "duplicate_cfos_to_nucleus_mappings": len(status.accepted_cfos_roi_ids) - len(status.cfos_positive_nucleus_ids),
        "cfos_positive_nuclei": len(status.cfos_positive_nucleus_ids),
        "npy_positive_nuclei": len(status.npy_positive_nucleus_ids),
        "dual_positive_nuclei": len(status.dual_positive_nucleus_ids),
    }
    display_dapi = display_transform(data.dapi_labels)
    lumen = ventricular_lumen_mask(display_dapi)
    lumen_ys, lumen_xs = np.where(lumen)
    return {
        "schema_version": 3,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {
            "path": _relative(Path(__file__), data.root),
            "sha256": sha256_file(Path(__file__)),
        },
        "sample": {
            "sample_id": SAMPLE_ID,
            "animal_id": ANIMAL_ID,
            "animal_id_status": "confirmed",
            "condition": CONDITION,
            "native_shape_yx": list(NATIVE_SHAPE_YX),
            "physical_size_xy_um_per_px": [UM_PER_PX, UM_PER_PX],
            "scalebar_um": 100.0,
        },
        "sources": {
            role: {"path": _relative(source.path, data.root), "sha256": source.sha256, "bytes": source.bytes}
            for role, source in data.sources.items()
        },
        "processed_masks": {
            marker: {
                "whole_field_instances": int(labels.max()),
                "whole_field_positive_pixels": int(np.count_nonzero(labels)),
                "whole_field_binary_sha256": binary_sha256(labels > 0),
                "display_field_instances": visible_masks[marker],
                "internal_source_filename": data.internal_filenames[marker],
                "semantics": "DAPI nuclear reference" if marker == "DAPI" else "Cellpose marker-object mask; not a cell call",
            }
            for marker, labels in labels_by_marker.items()
        },
        "analysis": {
            "classification_coordinates": "complete native label fields before display transform",
            "whole_field_counts": whole_status,
            "whole_field_counts_validated": whole_status == EXPECTED_ANALYSIS_COUNTS,
            "display_field_counts": visible_status,
            "cfos_rule": {
                "minimum_roi_area_px": CFOS_MIN_ROI_AREA_PX,
                "minimum_best_roi_fraction_overlapping_one_dapi_nucleus": CFOS_ROI_TO_NUCLEUS_OVERLAP,
            },
            "npy_rule": {
                "minimum_marker_pixels": NPY_MIN_MARKER_PIXELS,
                "minimum_marker_fraction_of_nucleus": NPY_MIN_NUCLEUS_FRACTION,
            },
        },
        "ventricular_lumen_outline": {
            "role": "display-only DAPI-supported annotation; never an analysis ROI or denominator",
            "source": "hash-bound registered Water2 DAPI Cellpose mask",
            "search_bounds_xyxy_half_open": list(LUMEN_SEARCH_BOUNDS_XYXY),
            "window_shape_px": [LUMEN_WINDOW_SIZE_PX, LUMEN_WINDOW_SIZE_PX],
            "maximum_dapi_positive_pixels_per_window": LUMEN_MAX_DAPI_PIXELS_PER_WINDOW,
            "maximum_dapi_positive_fraction": (
                LUMEN_MAX_DAPI_PIXELS_PER_WINDOW / float(LUMEN_WINDOW_SIZE_PX ** 2)
            ),
            "component_connectivity": 8,
            "minimum_component_area_px": LUMEN_MIN_COMPONENT_AREA_PX,
            "selected_area_px": int(lumen.sum()),
            "selected_bbox_yxyx_half_open": [
                int(lumen_ys.min()), int(lumen_xs.min()),
                int(lumen_ys.max()) + 1, int(lumen_xs.max()) + 1,
            ],
            "selected_centroid_yx": [float(lumen_ys.mean()), float(lumen_xs.mean())],
            "selected_binary_sha256": binary_sha256(lumen),
            "render": "opaque white lumen interior followed by one black outer contour",
            "display_boundary_smoothing": {
                "applied": True,
                "method": "periodic 1-D Gaussian over the single traced external path",
                "sigma_boundary_points": LUMEN_DISPLAY_SMOOTHING_POINTS,
                "why": (
                    "the lumen mask is a per-pixel DAPI intensity threshold, so its raw "
                    "0.5 contour is a staircase that reads as a tortuous line rather than "
                    "an ependymal border"
                ),
                "analysis_role": (
                    "display annotation only; the unsmoothed mask defines every analysis "
                    "ROI, denominator and hash-bound geometry"
                ),
            },
        },
        "display_transform": {
            "native_source_bounds_xyxy_half_open": list(SOURCE_BOUNDS_XYXY),
            "operations_in_order": [
                "native[y=1223:2903, x=188:1270]",
                "numpy.rot90(k=1, axes=(0, 1))",
            ],
            "rotation": "90 degrees counter-clockwise",
            "output_shape_yx": list(DISPLAY_SHAPE_YX),
            "display_to_native_xy_homogeneous_matrix": [list(row) for row in DISPLAY_TO_NATIVE_XY],
            "interpolation": "none; integer array slicing and index permutation only",
            "same_transform_for_all_masks": True,
            "output_field_um_wh": [DISPLAY_SHAPE_YX[1] * UM_PER_PX, DISPLAY_SHAPE_YX[0] * UM_PER_PX],
        },
        "display_semantics": {
            "panel_B": "display-field support of DAPI/NPY/c-FOS Cellpose masks",
            "panel_C": "display-field DAPI instances carrying whole-field c-FOS/NPY classifications",
            "tissue_envelope": "original registered DAPI-supported mask controls fill; hole-filled mask supplies one external gray display path only",
            "ventricular_lumen_outline": (
                "opaque clean-white interior and one black display path; hash-bound raw-DAPI base components; Water alone adds "
                "registered DAPI-label-clearance superior completion; not quantified"
            ),
            "colors": CHANNEL_COLORS | {"dual": DUAL_COLOR, "tissue": TISSUE_FILL, "lumen_outline": "black"},
        },
        "outputs": [
            {"path": _relative(path, data.root), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in outputs
        ],
    }


def _write_csv(path: Path, data: CartoonData, status: AnalysisStatus, outputs: Sequence[Path]) -> None:
    visible_masks, _ = _visible_counts(data, status)
    mask_by_role = {
        "dapi_seg": ("DAPI", data.dapi_labels),
        "npy_seg": ("NPY", data.npy_labels),
        "cfos_seg": ("c-FOS", data.cfos_labels),
    }
    rows: list[dict[str, Any]] = []
    for role, source in data.sources.items():
        item = mask_by_role.get(role)
        marker = item[0] if item else None
        labels = item[1] if item else None
        rows.append({
            "role": role,
            "path": _relative(source.path, data.root),
            "sha256": source.sha256,
            "bytes": source.bytes,
            "native_shape_yx": "3749x3768",
            "whole_field_instances": int(labels.max()) if labels is not None else "",
            "display_field_instances": visible_masks[marker] if marker else "",
            "notes": "native label source; whole-field classification" if labels is not None else "native TIFF companion",
        })
    display_dapi = display_transform(data.dapi_labels)
    lumen = ventricular_lumen_mask(display_dapi)
    rows.append({
        "role": "ventricular_lumen_display_annotation",
        "path": "",
        "sha256": binary_sha256(lumen),
        "bytes": "",
        "native_shape_yx": "1082x1680 registered display",
        "whole_field_instances": 1,
        "display_field_instances": 1,
        "notes": (
            "display only; 65x65 exact DAPI occupancy <=84; "
            "search xyxy=600,400,1080,850; area=10626; bbox yxyx=482,766,706,922"
        ),
    })
    for output in outputs:
        rows.append({
            "role": "rendered_panel",
            "path": _relative(output, data.root),
            "sha256": sha256_file(output),
            "bytes": output.stat().st_size,
            "native_shape_yx": "",
            "whole_field_instances": "",
            "display_field_instances": "",
            "notes": "1x3 output; 1082x1680-pixel registered field; integer crop plus CCW rotation",
        })
    fields = (
        "role", "path", "sha256", "bytes", "native_shape_yx",
        "whole_field_instances", "display_field_instances", "notes",
    )
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def render_all(
    root: Path | str | None = None,
    output_dir: Path | str | None = None,
    dpi: int = 300,
) -> Mapping[str, Path]:
    data = load_water2_inputs(root)
    status = classify_analysis_status(data.dapi_labels, data.cfos_labels, data.npy_labels)
    outdir = Path(output_dir).expanduser().resolve() if output_dir else data.root / "analyses/Fig3/results/final_run/panels"
    outdir.mkdir(parents=True, exist_ok=True)
    output = {
        "panel_b_png": outdir / f"{PANEL_B_STEM}.png",
        "panel_b_pdf": outdir / f"{PANEL_B_STEM}.pdf",
        "panel_c_png": outdir / f"{PANEL_C_STEM}.png",
        "panel_c_pdf": outdir / f"{PANEL_C_STEM}.pdf",
        "provenance_json": outdir / f"{PROVENANCE_STEM}.json",
        "provenance_csv": outdir / f"{PROVENANCE_STEM}.csv",
    }
    panel_b = draw_processed_mask_cartoon(data)
    try:
        _atomic_save(panel_b, output["panel_b_png"], dpi)
        _atomic_save(panel_b, output["panel_b_pdf"], dpi)
    finally:
        plt.close(panel_b)
    panel_c = draw_analysis_status_cartoon(data, status)
    try:
        _atomic_save(panel_c, output["panel_c_png"], dpi)
        _atomic_save(panel_c, output["panel_c_pdf"], dpi)
    finally:
        plt.close(panel_c)
    panel_paths = tuple(output[key] for key in ("panel_b_png", "panel_b_pdf", "panel_c_png", "panel_c_pdf"))
    _atomic_text(output["provenance_json"], json.dumps(_provenance_payload(data, status, panel_paths), indent=2, sort_keys=True) + "\n")
    _write_csv(output["provenance_csv"], data, status, panel_paths)
    return output



# ---------------------------------------------------------------------------
# Matched-condition Panel B/C implementation (schema v6).
# The legacy Water-only helpers above remain import-compatible, but render_all
# below intentionally replaces their old processed-mask/status output pairing.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MatchedSampleSpec:
    sample_id: str
    animal_id: str
    condition: str
    source_dir: str
    filenames: Mapping[str, str]
    native_shape_yx: tuple[int, int]
    output_shape_yx: tuple[int, int]
    display_to_native_xy: tuple[tuple[float, float, float], ...]
    transform_kind: str
    um_per_px: float
    tiff_px_per_um: float
    expected_sha256: Mapping[str, str]
    expected_instances: Mapping[str, int]
    expected_analysis_counts: Mapping[str, int]
    lumen: Mapping[str, Any]


WATER_CARTOON_SPEC = MatchedSampleSpec(
    sample_id="water-npy2-0001", animal_id="WATER_NPY2", condition="Water",
    source_dir="water",
    filenames={
        "dapi_tiff": "WATER_NPY2-0001_dapi.tif",
        "npy_tiff": "WATER_NPY2-0001_GFP.tif",
        "cfos_tiff": "WATER_NPY2-0001_fos.tif",
        "dapi_seg": "WATER_NPY2-0001_dapi_seg.npy",
        "npy_seg": "WATER_NPY2-0001_GFP_seg.npy",
        "cfos_seg": "WATER_NPY2-0001_fos_seg.npy",
    },
    native_shape_yx=(3749, 3768), output_shape_yx=(866, 1374),
    display_to_native_xy=(
        (0.0, -1.1747212143924286, 1236.6542858319217),
        (1.1747212143924286, 0.0, 1255.9665257124016),
        (0.0, 0.0, 1.0),
    ),
    transform_kind="fixed Panel-A-derived homography into the shared physical common frame",
    um_per_px=0.664212669961432, tiff_px_per_um=1.768592,
    expected_sha256=dict(EXPECTED_SOURCE_SHA256),
    expected_instances=dict(EXPECTED_INSTANCE_COUNTS),
    expected_analysis_counts=dict(EXPECTED_ANALYSIS_COUNTS),
    lumen={
        "expected_display_percentiles": (0.15234375, 191.72410644531215),
        "expected_component_id": 12,
        "expected_area_px": 29609, "expected_bbox_yxyx": (219, 537, 606, 802),
        "expected_centroid_yx": (474.75271032456345, 686.9537302847107),
        "expected_binary_sha256": "4bed23a808318227fb8b9f2e9a09c54e2bc37287383a2c3005f94bf40aba5b82",
        "superior_completion": {
            "corridor_x_half_open": (604, 769), "apex_xy": (708, 219),
            "expected_centerline_top_xy": (698, 0),
            "expected_centerline_x_range": (661, 708),
            "expected_centerline_total_variation_px": 106,
            "expected_centerline_sha256_i32_le": "44372b1db1143fa5c1bf90278b663ba22e9ab2152698a7e8f0b5f8b683001c67",
            "expected_extension_area_px": 5618,
            "expected_extension_bbox_yxyx": (0, 645, 222, 712),
            "expected_extension_binary_sha256": "f733aa88ee10a1960a38024b0e1537ac4a53c787f7e5f951287f4916b7068c2a",
            "expected_final_area_px": 35541,
            "expected_final_bbox_yxyx": (0, 537, 606, 802),
            "expected_final_centroid_yx": (415.9878168875383, 685.8594862271742),
            "expected_final_binary_sha256": "1e462a1232813f3ed4b36f6893d12b114473f1ff3cf7c089bfccf363e2e70063",
            "expected_final_top_row_x_inclusive": (692, 704),
        },
    },
)

ALLULOSE_CARTOON_SPEC = MatchedSampleSpec(
    sample_id="e8-fr6-4-alulosa", animal_id="E8_FR6-4", condition="Allulose",
    source_dir=".",
    filenames={
        "dapi_tiff": "E8-FR6-4_ALULOSA_dapi.tif",
        "npy_tiff": "E8-FR6-4_ALULOSA_gfp.tif",
        "cfos_tiff": "E8-FR6-4_ALULOSA_fos.tif",
        "dapi_seg": "E8-FR6-4_ALULOSA_dapi_seg.npy",
        "npy_seg": "E8-FR6-4_ALULOSA_gfp_seg.npy",
        "cfos_seg": "E8-FR6-4_ALULOSA_fos_seg.npy",
    },
    native_shape_yx=(946, 1382), output_shape_yx=(866, 1374),
    display_to_native_xy=(
        (0.991402, -0.125351, 59.858489),
        (-0.125762, -0.991309, 949.166712),
        (-0.000001, 0.0, 1.0),
    ),
    transform_kind="fixed Panel-A homography; cv2 INTER_NEAREST for labels",
    um_per_px=0.664212669961432, tiff_px_per_um=1.505541,
    expected_sha256={
        "dapi_tiff": "59318a5c319d2d72c5aa3bc672742852149ea9038cdd143b0df2b5c5e9fe1f6b",
        "npy_tiff": "4c21c8f06844f488757ce9c8e9b82b33fc7e430c6b9fe8d3e45f0a6bac897da1",
        "cfos_tiff": "a545e444fa819068b6ff0b1ab7d7488df543cc15f3f172126e03620a048de3c7",
        "dapi_seg": "f3ba80789acb7a2d1b79241236fef7fbee11d66a3f3db7c8f959e8c077965314",
        "npy_seg": "81978bcca94d0e74e6cdd5d53291283aa54e6ab9288d4531072ac4cf889c0d14",
        "cfos_seg": "b0df57d75a7307980a81cb34d280044e3cf6ea7b54fe82b12899838577a51dca",
    },
    expected_instances={"DAPI": 3372, "NPY": 265, "c-FOS": 328},
    expected_analysis_counts={
        "accepted_cfos_rois": 306,
        "duplicate_cfos_to_nucleus_mappings": 1,
        "cfos_positive_nuclei": 305,
        "npy_positive_nuclei": 494,
        "dual_positive_nuclei": 131,
    },
    lumen={
        "expected_display_percentiles": (1.5751953125, 255.0),
        "expected_component_id": 5,
        "expected_area_px": 62387, "expected_bbox_yxyx": (0, 621, 653, 861),
        "expected_centroid_yx": (377.91698591052625, 705.0579127061727),
        "expected_binary_sha256": "c4521d0f26dc6c1eb160343ef6f91a1516d16c53f1e65f6b8336e8977490b973",
    },
)

PANEL_B_STEM = "Figure3_Panel_B_Water_marker_positive_nuclei"
PANEL_C_STEM = "Figure3_Panel_C_Allulose_marker_positive_nuclei"


@dataclass(frozen=True)
class MatchedCartoonData:
    root: Path
    spec: MatchedSampleSpec
    dapi_labels: np.ndarray
    npy_labels: np.ndarray
    cfos_labels: np.ndarray
    sources: Mapping[str, SourceFile]
    internal_filenames: Mapping[str, str]


def _matched_source(path: Path, role: str, spec: MatchedSampleSpec) -> SourceFile:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {spec.condition} {role}: {path}")
    observed = sha256_file(path)
    expected = str(spec.expected_sha256[role])
    if observed != expected:
        raise RuntimeError(
            f"{spec.condition} {role} SHA256 mismatch: expected {expected}, got {observed}"
        )
    return SourceFile(role, path.resolve(), observed, path.stat().st_size)


def _matched_validate_tiff(path: Path, role: str, spec: MatchedSampleSpec) -> None:
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) != 1:
            raise RuntimeError(f"{spec.condition} {role} must be a single-page TIFF")
        page = tif.pages[0]
        if tuple(page.shape) != (*spec.native_shape_yx, 3) or page.dtype != np.dtype("uint8"):
            raise RuntimeError(
                f"{spec.condition} {role} is {page.shape}/{page.dtype}; expected "
                f"RGB uint8 at {spec.native_shape_yx}"
            )
        expected_resolution = float(spec.tiff_px_per_um)
        for tag_name in ("XResolution", "YResolution"):
            tag = page.tags.get(tag_name)
            if tag is None:
                raise RuntimeError(f"{spec.condition} {role} has no {tag_name} tag")
            observed = float(tag.value[0]) / float(tag.value[1])
            if abs(observed - expected_resolution) > 1e-6:
                raise RuntimeError(
                    f"{spec.condition} {role} {tag_name} changed: {observed} != {expected_resolution}"
                )


def _matched_load_labels(
    path: Path, role: str, marker: str, spec: MatchedSampleSpec,
) -> tuple[np.ndarray, str]:
    _matched_source(path, role, spec)
    loaded = np.load(path, allow_pickle=True)
    if not (isinstance(loaded, np.ndarray) and loaded.dtype == object and loaded.shape == ()):
        raise RuntimeError(f"{spec.condition} {role} is not a scalar Cellpose dictionary")
    payload = loaded.item()
    if not isinstance(payload, dict) or "masks" not in payload:
        raise RuntimeError(f"{spec.condition} {role} has no masks entry")
    labels = np.asarray(payload["masks"])
    if tuple(labels.shape) != spec.native_shape_yx:
        raise RuntimeError(
            f"{spec.condition} {role} shape changed: {labels.shape} != {spec.native_shape_yx}"
        )
    if not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0):
        raise RuntimeError(f"{spec.condition} {role} labels must be nonnegative integers")
    positive = np.unique(labels[labels > 0])
    if not np.array_equal(positive, np.arange(1, positive.size + 1, dtype=positive.dtype)):
        raise RuntimeError(f"{spec.condition} {role} labels are not contiguous")
    expected_instances = int(spec.expected_instances[marker])
    if positive.size != expected_instances:
        raise RuntimeError(
            f"{spec.condition} {role} has {positive.size} instances; expected {expected_instances}"
        )
    internal = str(payload.get("filename", "")).strip()
    if not internal:
        raise RuntimeError(f"{spec.condition} {role} has no internal source filename")
    return labels, internal


def _matched_load_inputs(
    spec: MatchedSampleSpec, root: Path | str | None = None,
) -> MatchedCartoonData:
    project_root = Path(root).expanduser().resolve() if root else next(path for path in Path(__file__).resolve().parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
    base = project_root / "Fig3/raw/legacy_experiment_2025_07_28" / spec.source_dir
    paths = {role: base / filename for role, filename in spec.filenames.items()}
    sources = {role: _matched_source(path, role, spec) for role, path in paths.items()}
    for role in ("dapi_tiff", "npy_tiff", "cfos_tiff"):
        _matched_validate_tiff(paths[role], role, spec)
    dapi, dapi_internal = _matched_load_labels(paths["dapi_seg"], "dapi_seg", "DAPI", spec)
    npy, npy_internal = _matched_load_labels(paths["npy_seg"], "npy_seg", "NPY", spec)
    cfos, cfos_internal = _matched_load_labels(paths["cfos_seg"], "cfos_seg", "c-FOS", spec)
    if not (dapi.shape == npy.shape == cfos.shape == spec.native_shape_yx):
        raise RuntimeError(f"{spec.condition} masks do not share exact native coordinates")
    return MatchedCartoonData(
        project_root, spec, dapi, npy, cfos, sources,
        {"DAPI": dapi_internal, "NPY": npy_internal, "c-FOS": cfos_internal},
    )


def load_allulose_inputs(root: Path | str | None = None) -> MatchedCartoonData:
    """Load the hash-bound Allulose E8_FR6-4 representative and masks."""
    return _matched_load_inputs(ALLULOSE_CARTOON_SPEC, root)


def _matched_classify(data: MatchedCartoonData) -> AnalysisStatus:
    dapi, cfos, npy = data.dapi_labels, data.cfos_labels, data.npy_labels
    accepted: list[int] = []
    mapping: dict[int, int] = {}
    maximum_dapi = int(dapi.max())
    for prop in regionprops(cfos):
        if int(prop.area) < CFOS_MIN_ROI_AREA_PX:
            continue
        values = dapi[prop.coords[:, 0], prop.coords[:, 1]]
        counts = np.bincount(values.astype(np.int64), minlength=maximum_dapi + 1)
        counts[0] = 0
        best = int(np.argmax(counts)) if counts.sum() else 0
        fraction = float(counts[best]) / float(prop.area) if best else 0.0
        if best and fraction >= CFOS_ROI_TO_NUCLEUS_OVERLAP:
            accepted.append(int(prop.label))
            mapping[int(prop.label)] = best
    npy_mask = npy > 0
    npy_positive: set[int] = set()
    for prop in regionprops(dapi):
        overlap = int(npy_mask[prop.coords[:, 0], prop.coords[:, 1]].sum())
        fraction = overlap / float(max(1, prop.area))
        if overlap >= NPY_MIN_MARKER_PIXELS and fraction >= NPY_MIN_NUCLEUS_FRACTION:
            npy_positive.add(int(prop.label))
    cfos_positive = set(mapping.values())
    dual = cfos_positive & npy_positive
    observed = {
        "accepted_cfos_rois": len(accepted),
        "duplicate_cfos_to_nucleus_mappings": len(accepted) - len(cfos_positive),
        "cfos_positive_nuclei": len(cfos_positive),
        "npy_positive_nuclei": len(npy_positive),
        "dual_positive_nuclei": len(dual),
    }
    if observed != dict(data.spec.expected_analysis_counts):
        raise RuntimeError(
            f"{data.spec.condition} whole-field counts changed: expected "
            f"{dict(data.spec.expected_analysis_counts)}, got {observed}"
        )
    return AnalysisStatus(
        tuple(sorted(accepted)), dict(sorted(mapping.items())),
        frozenset(cfos_positive), frozenset(npy_positive), frozenset(dual),
    )


def _matched_display_transform(array: np.ndarray, spec: MatchedSampleSpec) -> np.ndarray:
    values = np.asarray(array)
    if tuple(values.shape[:2]) != spec.native_shape_yx or values.ndim != 2:
        raise RuntimeError(
            f"{spec.condition} display input is {values.shape}; expected one field at {spec.native_shape_yx}"
        )
    height, width = spec.output_shape_yx
    matrix = np.asarray(spec.display_to_native_xy, dtype=np.float64)
    displayed = cv2.warpPerspective(
        values.astype(np.float32), matrix, (width, height),
        flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP,
        borderMode=cv2.BORDER_CONSTANT, borderValue=0,
    )
    displayed = np.rint(displayed).astype(np.int32)
    if tuple(displayed.shape) != spec.output_shape_yx or int(displayed.max()) > int(values.max()):
        raise RuntimeError(f"{spec.condition} nearest-neighbor label registration failed")
    return displayed


def _contour_path_counts(mask: np.ndarray) -> dict[str, int]:
    binary = np.ascontiguousarray(np.asarray(mask, dtype=np.uint8))
    external, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    tree, _ = cv2.findContours(binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    return {"external": len(external), "tree": len(tree)}


def _water_superior_lumen_completion(
    base_lumen: np.ndarray, displayed_dapi_labels: np.ndarray, rule: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Complete the Water lumen to y=0 from registered DAPI-label clearance."""
    base = np.asarray(base_lumen, dtype=bool)
    dapi = np.asarray(displayed_dapi_labels)
    if base.shape != dapi.shape:
        raise RuntimeError("Water lumen and registered DAPI labels do not share coordinates")
    settings = rule["superior_completion"]
    x0, x1 = (int(value) for value in settings["corridor_x_half_open"])
    apex_x, apex_y = (int(value) for value in settings["apex_xy"])
    ys, xs = np.where(base)
    top = int(ys.min())
    top_xs = xs[ys == top]
    if top != apex_y or top_xs.size != 1 or int(top_xs[0]) != apex_x:
        raise RuntimeError("Water base lumen no longer has the unique preregistered superior apex")

    clearance = distance_transform_edt(dapi == 0)
    corridor_width = x1 - x0
    row_count = apex_y + 1
    score = np.full((row_count, corridor_width), -np.inf, dtype=np.float64)
    predecessor = np.full((row_count, corridor_width), -1, dtype=np.int32)
    score[0, :] = clearance[0, x0:x1]
    for y in range(1, row_count):
        for local_x in range(corridor_width):
            previous_start = max(0, local_x - 1)
            previous_stop = min(corridor_width, local_x + 2)
            previous_values = score[y - 1, previous_start:previous_stop]
            # np.argmax returns the first maximum, implementing the declared
            # leftmost tie-break because candidates are ordered left-to-right.
            previous_x = previous_start + int(np.argmax(previous_values))
            predecessor[y, local_x] = previous_x
            score[y, local_x] = score[y - 1, previous_x] + clearance[y, x0 + local_x]

    centerline_x = np.empty(row_count, dtype=np.int32)
    local_x = apex_x - x0
    for y in range(apex_y, -1, -1):
        centerline_x[y] = x0 + local_x
        if y:
            local_x = int(predecessor[y, local_x])
    centerline_sha = hashlib.sha256(
        np.ascontiguousarray(centerline_x.astype("<i4")).tobytes()
    ).hexdigest()
    expected_top_x, expected_top_y = settings["expected_centerline_top_xy"]
    if (int(centerline_x[0]), 0) != (int(expected_top_x), int(expected_top_y)):
        raise RuntimeError("Water superior-completion centerline top endpoint changed")
    if (int(centerline_x.min()), int(centerline_x.max())) != tuple(settings["expected_centerline_x_range"]):
        raise RuntimeError("Water superior-completion centerline x range changed")
    total_variation = int(np.abs(np.diff(centerline_x)).sum())
    if total_variation != int(settings["expected_centerline_total_variation_px"]):
        raise RuntimeError("Water superior-completion centerline variation changed")
    if centerline_sha != str(settings["expected_centerline_sha256_i32_le"]):
        raise RuntimeError("Water superior-completion centerline receipt changed")

    extension_u8 = np.zeros(base.shape, dtype=np.uint8)
    radii: list[int] = []
    for y, x in enumerate(centerline_x.tolist()):
        radius = max(1, int(np.floor(clearance[y, x])) - 1)
        radii.append(radius)
        cv2.circle(extension_u8, (int(x), int(y)), radius, 1, thickness=-1)
    extension = extension_u8.astype(bool)
    ext_y, ext_x = np.where(extension)
    extension_bbox = [int(ext_y.min()), int(ext_x.min()), int(ext_y.max()) + 1, int(ext_x.max()) + 1]
    extension_sha = binary_sha256(extension)
    if int(extension.sum()) != int(settings["expected_extension_area_px"]):
        raise RuntimeError("Water superior-completion extension area changed")
    if extension_bbox != list(settings["expected_extension_bbox_yxyx"]):
        raise RuntimeError("Water superior-completion extension bounds changed")
    if extension_sha != str(settings["expected_extension_binary_sha256"]):
        raise RuntimeError("Water superior-completion extension binary changed")

    union = base | extension
    final = np.asarray(binary_fill_holes(union), dtype=bool)
    final_y, final_x = np.where(final)
    final_bbox = [int(final_y.min()), int(final_x.min()), int(final_y.max()) + 1, int(final_x.max()) + 1]
    final_centroid = [float(final_y.mean()), float(final_x.mean())]
    final_sha = binary_sha256(final)
    top_row_x = np.where(final[0])[0]
    top_span = [int(top_row_x.min()), int(top_row_x.max())]
    path_counts = _contour_path_counts(final)
    if int(final.sum()) != int(settings["expected_final_area_px"]):
        raise RuntimeError("Water completed lumen area changed")
    if final_bbox != list(settings["expected_final_bbox_yxyx"]):
        raise RuntimeError("Water completed lumen bounds changed")
    if not np.allclose(final_centroid, settings["expected_final_centroid_yx"], atol=1e-9, rtol=0.0):
        raise RuntimeError("Water completed lumen centroid changed")
    if final_sha != str(settings["expected_final_binary_sha256"]):
        raise RuntimeError("Water completed lumen binary changed")
    if top_span != list(settings["expected_final_top_row_x_inclusive"]):
        raise RuntimeError("Water completed lumen top span changed")
    if path_counts != {"external": 1, "tree": 1}:
        raise RuntimeError(f"Water completed lumen must have one outer path: {path_counts}")

    receipt = {
        "applied": True,
        "analysis_role": "display annotation only; never an analysis ROI or denominator",
        "distance_map_source": "Euclidean distance transform computed from the registered DAPI Cellpose-label complement",
        "centerline_rule": (
            "top-to-fixed-apex dynamic programming maximizing cumulative EDT; predecessor x-1/x/x+1; "
            "leftmost argmax tie-break; path is not required to remain exclusively in the label complement"
        ),
        "corridor_x_half_open": [x0, x1],
        "fixed_apex_xy": [apex_x, apex_y],
        "centerline_top_xy": [int(centerline_x[0]), 0],
        "centerline_x_range": [int(centerline_x.min()), int(centerline_x.max())],
        "centerline_total_variation_px": total_variation,
        "centerline_zero_clearance_pixel_count": int(np.count_nonzero(clearance[np.arange(row_count), centerline_x] == 0)),
        "centerline_x_sha256_i32_le": centerline_sha,
        "tube_radius_rule": "max(1, floor(EDT) - 1)",
        "tube_radius_min_median_max_px": [int(min(radii)), float(np.median(radii)), int(max(radii))],
        "extension_area_px": int(extension.sum()),
        "extension_bbox_yxyx_half_open": extension_bbox,
        "extension_binary_sha256": extension_sha,
        "base_extension_overlap_px": int(np.count_nonzero(base & extension)),
        "union_pre_fill_area_px": int(union.sum()),
        "holes_filled_px": int(final.sum() - union.sum()),
        "final_area_px": int(final.sum()),
        "final_bbox_yxyx_half_open": final_bbox,
        "final_centroid_yx": final_centroid,
        "final_top_row_x_inclusive": top_span,
        "final_binary_sha256": final_sha,
        "final_connected_components_8": int(label(final, structure=np.ones((3, 3), dtype=np.uint8))[1]),
        "final_contour_path_counts": path_counts,
    }
    return final, receipt


def _matched_lumen_source_result(
    data: MatchedCartoonData,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return the hash-bound raw-DAPI ventricular proposal and receipt."""
    spec = data.spec
    raw = np.asarray(tifffile.imread(data.sources["dapi_tiff"].path))
    if raw.shape != (*spec.native_shape_yx, 3) or raw.dtype != np.uint8:
        raise RuntimeError(f"{spec.condition} raw DAPI TIFF geometry changed: {raw.shape}/{raw.dtype}")
    scalar = np.max(raw[..., :3], axis=-1).astype(np.float32)
    height, width = spec.output_shape_yx
    matrix = np.asarray(spec.display_to_native_xy, dtype=np.float64)
    displayed = cv2.warpPerspective(
        scalar, matrix, (width, height),
        flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP,
        borderMode=cv2.BORDER_CONSTANT, borderValue=0,
    )
    positive = displayed[np.isfinite(displayed) & (displayed > 0)]
    if positive.size < 100:
        raise RuntimeError(f"{spec.condition} raw DAPI common-frame population is too small")
    low, high = (float(value) for value in np.percentile(positive, [0.5, 99.5]))
    expected_low, expected_high = (float(value) for value in spec.lumen["expected_display_percentiles"])
    if not np.allclose((low, high), (expected_low, expected_high), atol=1e-10, rtol=0.0):
        raise RuntimeError(
            f"{spec.condition} raw-DAPI display percentiles changed: {(low, high)} != "
            f"{(expected_low, expected_high)}"
        )
    normalized = np.clip((displayed - low) / (high - low), 0.0, 1.0)
    sigma_px = 1.0 / float(spec.um_per_px)
    smooth = gaussian_filter(normalized, sigma=sigma_px, mode="nearest")
    x0, x1 = int(np.floor(0.30 * width)), int(np.floor(0.70 * width))
    y0, y1 = 0, int(np.floor(0.88 * height))
    roi = np.zeros((height, width), dtype=bool)
    roi[y0:y1, x0:x1] = True
    candidate = roi & (smooth <= 0.08)
    components, component_count = label(candidate, structure=np.ones((3, 3), dtype=np.uint8))
    sx0, sx1 = int(np.floor(0.44 * width)), int(np.floor(0.56 * width))
    sy0, sy1 = int(np.floor(0.45 * height)), int(np.floor(0.70 * height))
    accepted: list[tuple[int, np.ndarray]] = []
    for component_id in range(1, int(component_count) + 1):
        component = components == component_id
        ys, xs = np.where(component)
        if ys.size < 2000:
            continue
        intersects_seed = bool(np.any(component[sy0:sy1, sx0:sx1]))
        touches_forbidden_edge = bool(
            np.any(xs == x0) or np.any(xs == x1 - 1) or np.any(ys == y1 - 1)
        )
        if intersects_seed and not touches_forbidden_edge:
            accepted.append((component_id, component))
    if not accepted:
        raise RuntimeError(f"No accepted {spec.condition} raw-DAPI lumen component")
    accepted.sort(key=lambda item: int(item[1].sum()), reverse=True)
    if len(accepted) > 1 and int(accepted[0][1].sum()) == int(accepted[1][1].sum()):
        raise RuntimeError(f"{spec.condition} lumen largest-component selection is tied")
    component_id, base = accepted[0]
    base = np.asarray(base, dtype=bool)
    ys, xs = np.where(base)
    base_bbox = [int(ys.min()), int(xs.min()), int(ys.max()) + 1, int(xs.max()) + 1]
    base_centroid = [float(ys.mean()), float(xs.mean())]
    base_sha = binary_sha256(base)
    if component_id != int(spec.lumen["expected_component_id"]):
        raise RuntimeError(f"{spec.condition} lumen component ID changed")
    if int(base.sum()) != int(spec.lumen["expected_area_px"]) or base_bbox != list(spec.lumen["expected_bbox_yxyx"]):
        raise RuntimeError(f"{spec.condition} base lumen geometry changed")
    if not np.allclose(base_centroid, spec.lumen["expected_centroid_yx"], atol=1e-9, rtol=0.0):
        raise RuntimeError(f"{spec.condition} base lumen centroid changed")
    if base_sha != str(spec.lumen["expected_binary_sha256"]):
        raise RuntimeError(f"{spec.condition} base lumen binary changed")

    base_receipt = {
        "source": "registered raw DAPI TIFF, normalized/smoothed low-intensity component",
        "component_id": int(component_id), "area_px": int(base.sum()),
        "bbox_yxyx_half_open": base_bbox, "centroid_yx": base_centroid,
        "binary_sha256": base_sha, "contour_path_counts_before_hole_fill": _contour_path_counts(base),
    }
    if spec.condition == "Water":
        displayed_dapi = _matched_display_transform(data.dapi_labels, spec)
        final, completion = _water_superior_lumen_completion(base, displayed_dapi, spec.lumen)
    else:
        final = np.asarray(binary_fill_holes(base), dtype=bool)
        completion = {
            "applied": False,
            "reason": "Allulose base raw-DAPI component already reaches the superior frame edge",
            "analysis_role": "display annotation only; never an analysis ROI or denominator",
        }
    source_final = np.asarray(final, dtype=bool)
    source_y, source_x = np.where(source_final)
    source_path_counts = _contour_path_counts(source_final)
    if source_path_counts != {"external": 1, "tree": 1}:
        raise RuntimeError(
            f"{spec.condition} source lumen must have exactly one contour path: {source_path_counts}"
        )
    receipt = {
        "raw_dapi_base_component": base_receipt,
        "superior_completion": completion,
        "hole_fill_before_contouring": True,
        "source_final_after_completion_and_hole_fill": {
            "area_px": int(source_final.sum()),
            "bbox_yxyx_half_open": [
                int(source_y.min()), int(source_x.min()),
                int(source_y.max()) + 1, int(source_x.max()) + 1,
            ],
            "centroid_yx": [float(source_y.mean()), float(source_x.mean())],
            "binary_sha256": binary_sha256(source_final),
            "connected_components_8": int(label(
                source_final, structure=np.ones((3, 3), dtype=np.uint8)
            )[1]),
            "contour_path_counts": source_path_counts,
        },
    }
    return source_final, receipt


def _resolve_ventricle_hil_dir(
    root: Path, ventricle_hil_dir: Path | str | None,
) -> Path:
    hil_dir = (
        Path(ventricle_hil_dir).expanduser().resolve()
        if ventricle_hil_dir is not None
        else (root / VENTRICLE_HIL_DEFAULT_RELATIVE).resolve()
    )
    if not hil_dir.is_dir() or hil_dir.is_symlink():
        raise RuntimeError(
            f"Accepted Figure 3 ventricle HIL directory is missing/non-normal: {hil_dir}. "
            "Run Fig3/02a_review_cfos_npy_ventricles.py and accept Water and Allulose first."
        )
    return hil_dir


def _load_hil_lumen_mask(
    data: MatchedCartoonData,
    source_lumen: np.ndarray,
    ventricle_hil_dir: Path | str | None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load and verify one accepted HIL lumen polygon and its raster mask."""
    hil_dir = _resolve_ventricle_hil_dir(data.root, ventricle_hil_dir)
    slug = data.spec.condition.lower()
    receipt_path = hil_dir / f"{slug}_ventricle_receipt.json"
    mask_path = hil_dir / f"{slug}_ventricle_mask.png"
    for path in (receipt_path, mask_path):
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"Missing/non-normal accepted {data.spec.condition} HIL file: {path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") != VENTRICLE_HIL_SCHEMA or receipt.get("status") != "accepted":
        raise RuntimeError(f"{data.spec.condition} HIL receipt is not accepted schema v1")
    expected_identity = {
        "condition": data.spec.condition,
        "sample_id": data.spec.sample_id,
        "animal_id": data.spec.animal_id,
        "coordinate_frame_yx": list(data.spec.output_shape_yx),
    }
    observed_identity = {key: receipt.get(key) for key in expected_identity}
    if observed_identity != expected_identity:
        raise RuntimeError(
            f"{data.spec.condition} HIL identity/geometry mismatch: {observed_identity}"
        )
    reviewer = str(receipt.get("reviewer", "")).strip()
    session = str(receipt.get("session", "")).strip()
    accepted_at = str(receipt.get("accepted_at_utc", "")).strip()
    try:
        accepted_time = datetime.fromisoformat(accepted_at.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RuntimeError(f"{data.spec.condition} HIL accepted_at_utc is invalid") from exc
    if not reviewer or not session or accepted_time.tzinfo is None:
        raise RuntimeError(f"{data.spec.condition} HIL receipt lacks review identity/timezone")

    source_binding = receipt.get("source_binding", {})
    expected_binding = {
        "dapi_tiff_sha256": data.sources["dapi_tiff"].sha256,
        "registered_dapi_binary_sha256": binary_sha256(
            _matched_display_transform(data.dapi_labels, data.spec) > 0
        ),
        "raw_dapi_lumen_source_binary_sha256": binary_sha256(source_lumen),
    }
    observed_binding = {key: source_binding.get(key) for key in expected_binding}
    if observed_binding != expected_binding:
        raise RuntimeError(
            f"{data.spec.condition} HIL source binding mismatch: {observed_binding}"
        )

    mask_u8 = cv2.imread(str(mask_path), cv2.IMREAD_UNCHANGED)
    if mask_u8 is None or mask_u8.ndim != 2 or tuple(mask_u8.shape) != data.spec.output_shape_yx:
        raise RuntimeError(f"{data.spec.condition} HIL mask geometry is invalid")
    unique = set(int(value) for value in np.unique(mask_u8))
    if not unique.issubset({0, 255}) or 255 not in unique:
        raise RuntimeError(f"{data.spec.condition} HIL mask must be nonempty binary 0/255 PNG")
    mask = mask_u8 == 255
    if not np.array_equal(mask, binary_fill_holes(mask)):
        raise RuntimeError(f"{data.spec.condition} HIL lumen contains an internal hole")
    component_count = int(label(mask, structure=np.ones((3, 3), dtype=np.uint8))[1])
    path_counts = _contour_path_counts(mask)
    if component_count != 1 or path_counts != {"external": 1, "tree": 1}:
        raise RuntimeError(
            f"{data.spec.condition} HIL lumen must be exactly one filled polygon: "
            f"components={component_count}, contours={path_counts}"
        )

    vertices = np.asarray(receipt.get("polygon_vertices_xy", []), dtype=np.float64)
    if vertices.ndim != 2 or vertices.shape[1:] != (2,) or vertices.shape[0] < 3:
        raise RuntimeError(f"{data.spec.condition} HIL polygon vertices are malformed")
    if not np.all(np.isfinite(vertices)):
        raise RuntimeError(f"{data.spec.condition} HIL polygon vertices are non-finite")
    raster_vertices = np.rint(vertices).astype(np.int32)
    height, width = data.spec.output_shape_yx
    if (
        np.any(raster_vertices[:, 0] < 0) or np.any(raster_vertices[:, 0] >= width)
        or np.any(raster_vertices[:, 1] < 0) or np.any(raster_vertices[:, 1] >= height)
    ):
        raise RuntimeError(f"{data.spec.condition} HIL polygon leaves the common frame")
    rerasterized = np.zeros(data.spec.output_shape_yx, dtype=np.uint8)
    cv2.fillPoly(rerasterized, [raster_vertices.reshape(-1, 1, 2)], 1)
    if not np.array_equal(rerasterized > 0, mask):
        raise RuntimeError(f"{data.spec.condition} HIL mask does not match its recorded polygon")

    ys, xs = np.where(mask)
    geometry = {
        "area_px": int(mask.sum()),
        "bbox_yxyx_half_open": [
            int(ys.min()), int(xs.min()), int(ys.max()) + 1, int(xs.max()) + 1,
        ],
        "centroid_yx": [float(ys.mean()), float(xs.mean())],
        "binary_sha256": binary_sha256(mask),
        "connected_components_8": component_count,
        "contour_path_counts": path_counts,
    }
    recorded_geometry = receipt.get("mask_geometry", {})
    for key, value in geometry.items():
        observed = recorded_geometry.get(key)
        if key == "centroid_yx":
            if not np.allclose(observed, value, atol=1e-9, rtol=0.0):
                raise RuntimeError(f"{data.spec.condition} HIL centroid receipt mismatch")
        elif observed != value:
            raise RuntimeError(f"{data.spec.condition} HIL geometry receipt mismatch: {key}")
    if receipt.get("mask_filename") != mask_path.name:
        raise RuntimeError(f"{data.spec.condition} HIL mask filename receipt mismatch")
    mask_file_sha = sha256_file(mask_path)
    if receipt.get("mask_file_sha256") != mask_file_sha:
        raise RuntimeError(f"{data.spec.condition} HIL mask file hash mismatch")
    vertices_hash = hashlib.sha256(
        np.ascontiguousarray(raster_vertices, dtype="<i4").tobytes()
    ).hexdigest()
    if receipt.get("raster_vertices_sha256_i32_le") != vertices_hash:
        raise RuntimeError(f"{data.spec.condition} HIL vertex hash mismatch")
    return mask, {
        "method": "human-in-the-loop closed polygon in the exact registered common frame",
        "analysis_role": "display annotation only; never an analysis ROI or denominator",
        "receipt_path": _relative(receipt_path, data.root),
        "receipt_sha256": sha256_file(receipt_path),
        "mask_path": _relative(mask_path, data.root),
        "mask_file_sha256": mask_file_sha,
        "reviewer": reviewer,
        "session": session,
        "accepted_at_utc": accepted_at,
        "notes": str(receipt.get("notes", "")),
        "polygon_vertex_count": int(vertices.shape[0]),
        "raster_vertices_sha256_i32_le": vertices_hash,
        "source_binding": expected_binding,
        "mask_geometry": geometry,
    }


def _matched_lumen_result(
    data: MatchedCartoonData,
    ventricle_hil_dir: Path | str | None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return the accepted, fail-closed HIL ventricular mask and receipt."""
    source_lumen, source_receipt = _matched_lumen_source_result(data)
    final, hil_receipt = _load_hil_lumen_mask(data, source_lumen, ventricle_hil_dir)
    ys, xs = np.where(final)
    paths = _contour_path_counts(final)
    return final, {
        **source_receipt,
        "human_in_loop_annotation": hil_receipt,
        "final_area_px": int(final.sum()),
        "final_bbox_yxyx_half_open": [
            int(ys.min()), int(xs.min()), int(ys.max()) + 1, int(xs.max()) + 1,
        ],
        "final_centroid_yx": [float(ys.mean()), float(xs.mean())],
        "final_binary_sha256": binary_sha256(final),
        "final_connected_components_8": int(label(
            final, structure=np.ones((3, 3), dtype=np.uint8)
        )[1]),
        "final_contour_path_counts": paths,
    }


def _matched_lumen_mask(
    data: MatchedCartoonData, ventricle_hil_dir: Path | str | None,
) -> np.ndarray:
    return _matched_lumen_result(data, ventricle_hil_dir)[0]

def _matched_layers(
    data: MatchedCartoonData, ventricle_hil_dir: Path | str | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dapi = _matched_display_transform(data.dapi_labels, data.spec)
    npy = _matched_display_transform(data.npy_labels, data.spec)
    cfos = _matched_display_transform(data.cfos_labels, data.spec)
    tissue = _tissue_envelope(dapi)
    lumen = _matched_lumen_mask(data, ventricle_hil_dir)
    return dapi, npy, cfos, tissue, lumen


def _matched_selected_mask(labels: np.ndarray, ids: Iterable[int]) -> np.ndarray:
    ids_array = np.asarray(sorted(set(int(value) for value in ids)), dtype=np.int64)
    maximum = int(labels.max())
    if ids_array.size and (ids_array.min() < 1 or ids_array.max() > maximum):
        raise RuntimeError("Selected nucleus ID is outside the native DAPI range")
    lookup = np.zeros(maximum + 1, dtype=bool)
    lookup[ids_array] = True
    return lookup[labels]


def _matched_add_scalebar(ax: plt.Axes, shape: tuple[int, int], um_per_px: float) -> None:
    height, width = shape
    length_px = 100.0 / float(um_per_px)
    x0 = width - max(30.0, 0.026 * width) - length_px
    y0 = height - max(25.0, 0.055 * height)
    bar_h = max(9.0, 0.015 * height)
    pad = max(5.0, 0.025 * length_px)
    ax.add_patch(Rectangle(
        (x0 - pad, y0 - bar_h - 55), length_px + 2 * pad, bar_h + 62,
        facecolor="black", edgecolor="black", linewidth=1.1, alpha=0.92, zorder=30,
    ))
    ax.add_patch(Rectangle(
        (x0, y0 - bar_h), length_px, bar_h,
        facecolor="white", edgecolor="white", linewidth=1.0, zorder=31,
    ))
    label_artist = ax.text(
        x0 + length_px / 2, y0 - bar_h - 4, "100 µm",
        color="white", ha="center", va="bottom", fontsize=6.2,
        fontweight="bold", zorder=32,
    )
    label_artist.set_path_effects([pe.withStroke(linewidth=3.2, foreground="black")])


def _matched_finish_axis(
    ax: plt.Axes, shape: tuple[int, int], lumen: np.ndarray, um_per_px: float,
) -> None:
    _draw_lumen_outline(ax, lumen)
    _matched_add_scalebar(ax, shape, um_per_px)
    ax.set_xlim(-0.5, shape[1] - 0.5)
    ax.set_ylim(shape[0] - 0.5, -0.5)
    ax.set_axis_off()


def _matched_status_axis(
    ax: plt.Axes, dapi: np.ndarray, tissue: np.ndarray, lumen: np.ndarray,
    layers: Sequence[tuple[Iterable[int], str]], title: str, um_per_px: float,
) -> None:
    _draw_background(ax, tissue, lumen)
    base = np.zeros((*dapi.shape, 4), dtype=np.float32)
    base[..., :3] = to_rgb(NEGATIVE_NUCLEUS_COLOR)
    base[..., 3] = (dapi > 0).astype(np.float32) * 0.25
    ax.imshow(base, interpolation="nearest")
    for ids, color in layers:
        rgba = np.zeros((*dapi.shape, 4), dtype=np.float32)
        rgba[..., :3] = to_rgb(color)
        rgba[..., 3] = _matched_selected_mask(dapi, ids).astype(np.float32) * 0.94
        ax.imshow(rgba, interpolation="nearest")
    ax.set_title(title, fontsize=10.5, fontweight="bold", color="#172033", pad=3)
    _matched_finish_axis(ax, dapi.shape, lumen, um_per_px)


def _matched_status_cartoon(
    data: MatchedCartoonData, status: AnalysisStatus, panel_letter: str,
    ventricle_hil_dir: Path | str | None,
) -> plt.Figure:
    dapi, _, _, tissue, lumen = _matched_layers(data, ventricle_hil_dir)
    fig, axes = _base_figure(
        panel_letter,
        f"{data.spec.condition} — marker-positive DAPI nuclei",
        bottom=0.120,
    )
    visible = _visible_ids(dapi)
    npy_n = len(visible & status.npy_positive_nucleus_ids)
    cfos_n = len(visible & status.cfos_positive_nucleus_ids)
    dual_n = len(visible & status.dual_positive_nucleus_ids)
    _matched_status_axis(
        axes[0], dapi, tissue, lumen,
        ((status.npy_positive_nucleus_ids, CHANNEL_COLORS["NPY"]),),
        "NPY-positive nuclei", data.spec.um_per_px,
    )
    _matched_status_axis(
        axes[1], dapi, tissue, lumen,
        ((status.cfos_positive_nucleus_ids, CHANNEL_COLORS["c-FOS"]),),
        "c-FOS-positive nuclei", data.spec.um_per_px,
    )
    npy_only = status.npy_positive_nucleus_ids - status.dual_positive_nucleus_ids
    cfos_only = status.cfos_positive_nucleus_ids - status.dual_positive_nucleus_ids
    _matched_status_axis(
        axes[2], dapi, tissue, lumen,
        ((npy_only, CHANNEL_COLORS["NPY"]),
         (cfos_only, CHANNEL_COLORS["c-FOS"]),
         (status.dual_positive_nucleus_ids, DUAL_COLOR)),
        "Combined status", data.spec.um_per_px,
    )
    fig.legend(
        handles=(
            Patch(facecolor=NEGATIVE_NUCLEUS_COLOR, alpha=0.45, label="Other DAPI nuclei"),
            Patch(facecolor=CHANNEL_COLORS["NPY"], label="NPY-positive nucleus"),
            Patch(facecolor=CHANNEL_COLORS["c-FOS"], label="c-FOS-positive nucleus"),
            Patch(facecolor=DUAL_COLOR, label="Dual-positive nucleus"),
        ),
        loc="lower center", bbox_to_anchor=(0.50, 0.002), ncol=4,
        frameon=False, fontsize=7.0, handlelength=1.4, columnspacing=1.4,
    )
    return fig


def _matched_visible_counts(
    data: MatchedCartoonData, status: AnalysisStatus,
) -> tuple[dict[str, int], dict[str, int]]:
    displayed = {
        "DAPI": _matched_display_transform(data.dapi_labels, data.spec),
        "NPY": _matched_display_transform(data.npy_labels, data.spec),
        "c-FOS": _matched_display_transform(data.cfos_labels, data.spec),
    }
    visible = _visible_ids(displayed["DAPI"])
    return (
        {marker: len(_visible_ids(labels)) for marker, labels in displayed.items()},
        {
            "npy_positive_nuclei": len(visible & status.npy_positive_nucleus_ids),
            "cfos_positive_nuclei": len(visible & status.cfos_positive_nucleus_ids),
            "dual_positive_nuclei": len(visible & status.dual_positive_nucleus_ids),
        },
    )


def _matched_tissue_outline_receipt(data: MatchedCartoonData) -> dict[str, Any]:
    displayed_dapi = _matched_display_transform(data.dapi_labels, data.spec)
    tissue = _tissue_envelope(displayed_dapi)
    outline = np.asarray(binary_fill_holes(tissue), dtype=bool)
    paths = _contour_path_counts(outline)
    expected = {
        "Water": {
            "original_area_px": 1130258, "outline_area_px": 1132291, "holes_filled_px": 2033,
            "outline_binary_sha256": "37492cf3f736f56f810488e0b2a497e89b2c8877976af21ae9e94e18ee5e0d9b",
        },
        "Allulose": {
            "original_area_px": 1096130, "outline_area_px": 1100602, "holes_filled_px": 4472,
            "outline_binary_sha256": "d8053c8bffd6fb5e7b44d9623e9cc786c9c14beed4f8bece1073910a3b0b660a",
        },
    }[data.spec.condition]
    observed = {
        "original_area_px": int(tissue.sum()),
        "outline_area_px": int(outline.sum()),
        "holes_filled_px": int(outline.sum() - tissue.sum()),
        "outline_binary_sha256": binary_sha256(outline),
    }
    if observed != expected or paths != {"external": 1, "tree": 1}:
        raise RuntimeError(f"{data.spec.condition} external tissue-outline receipt changed: {observed}/{paths}")
    return {
        "fill_mask_policy": "original registered DAPI-supported tissue mask controls background fill",
        "outline_mask_policy": "scipy.ndimage.binary_fill_holes(original tissue mask) used for display contour only",
        **observed,
        "contour_path_counts": paths,
        "internal_gray_contour_path_count": 0,
        "analysis_role": "display only; never an anatomical analysis ROI or denominator",
    }


def _matched_lumen_receipt(
    data: MatchedCartoonData, ventricle_hil_dir: Path | str | None,
) -> dict[str, Any]:
    lumen, algorithm = _matched_lumen_result(data, ventricle_hil_dir)
    ys, xs = np.where(lumen)
    rule = data.spec.lumen
    return {
        "role": "display-only ventricular annotation; never an analysis ROI or denominator",
        "base_source": f"hash-bound registered {data.spec.condition} raw DAPI TIFF",
        "common_frame_shape_yx": list(data.spec.output_shape_yx),
        "common_frame_um_per_px": data.spec.um_per_px,
        "raw_display_percentiles_0_5_99_5": list(rule["expected_display_percentiles"]),
        "gaussian_sigma_px": 1.0 / data.spec.um_per_px,
        "base_threshold_normalized_intensity_lte": 0.08,
        "search_bounds_xyxy_half_open": [
            int(np.floor(0.30 * data.spec.output_shape_yx[1])), 0,
            int(np.floor(0.70 * data.spec.output_shape_yx[1])),
            int(np.floor(0.88 * data.spec.output_shape_yx[0])),
        ],
        "base_component_selection": "unique largest seed-intersecting component; top touch allowed; ROI lateral/bottom touch rejected",
        "raw_dapi_base_component": algorithm["raw_dapi_base_component"],
        "superior_completion": algorithm["superior_completion"],
        "source_final_after_completion_and_hole_fill": algorithm[
            "source_final_after_completion_and_hole_fill"
        ],
        "human_in_loop_annotation": algorithm["human_in_loop_annotation"],
        "hole_fill_before_contouring": True,
        "selected_area_px": int(lumen.sum()),
        "selected_bbox_yxyx_half_open": [int(ys.min()), int(xs.min()), int(ys.max()) + 1, int(xs.max()) + 1],
        "selected_centroid_yx": [float(ys.mean()), float(xs.mean())],
        "selected_binary_sha256": binary_sha256(lumen),
        "connected_components_8": algorithm["final_connected_components_8"],
        "contour_path_counts": algorithm["final_contour_path_counts"],
        "render": (
            "opaque white corrected midline interior followed by exactly one black outer boundary; "
            "no internal subdivision"
        ),
        "interior_cleanup": {
            "opaque_white_fill": True,
            "covers_dapi_objects_and_prior_outline_pixels": True,
            "white_backing_stroke": False,
        },
    }

def _matched_sample_receipt(
    data: MatchedCartoonData, status: AnalysisStatus,
    ventricle_hil_dir: Path | str | None,
) -> dict[str, Any]:
    visible_masks, visible_status = _matched_visible_counts(data, status)
    labels = {"DAPI": data.dapi_labels, "NPY": data.npy_labels, "c-FOS": data.cfos_labels}
    whole_status = {
        "accepted_cfos_rois": len(status.accepted_cfos_roi_ids),
        "duplicate_cfos_to_nucleus_mappings": len(status.accepted_cfos_roi_ids) - len(status.cfos_positive_nucleus_ids),
        "cfos_positive_nuclei": len(status.cfos_positive_nucleus_ids),
        "npy_positive_nuclei": len(status.npy_positive_nucleus_ids),
        "dual_positive_nuclei": len(status.dual_positive_nucleus_ids),
    }
    return {
        "sample_id": data.spec.sample_id,
        "animal_id": data.spec.animal_id,
        "condition": data.spec.condition,
        "native_shape_yx": list(data.spec.native_shape_yx),
        "physical_size_xy_um_per_px": [data.spec.um_per_px, data.spec.um_per_px],
        "scalebar_um": 100.0,
        "sources": {
            role: {
                "path": _relative(source.path, data.root), "sha256": source.sha256,
                "bytes": source.bytes,
            }
            for role, source in data.sources.items()
        },
        "processed_masks": {
            marker: {
                "whole_field_instances": int(mask.max()),
                "whole_field_positive_pixels": int(np.count_nonzero(mask)),
                "whole_field_binary_sha256": binary_sha256(mask > 0),
                "display_field_instances": visible_masks[marker],
                "internal_source_filename": data.internal_filenames[marker],
                "semantics": "DAPI nuclear reference" if marker == "DAPI" else "Cellpose marker-object mask; not a cell call",
            }
            for marker, mask in labels.items()
        },
        "analysis": {
            "classification_coordinates": "complete native label fields before display transform",
            "whole_field_counts": whole_status,
            "whole_field_counts_validated": whole_status == dict(data.spec.expected_analysis_counts),
            "display_field_counts": visible_status,
            "cfos_rule": {
                "minimum_roi_area_px": CFOS_MIN_ROI_AREA_PX,
                "minimum_best_roi_fraction_overlapping_one_dapi_nucleus": CFOS_ROI_TO_NUCLEUS_OVERLAP,
            },
            "npy_rule": {
                "minimum_marker_pixels": NPY_MIN_MARKER_PIXELS,
                "minimum_marker_fraction_of_nucleus": NPY_MIN_NUCLEUS_FRACTION,
            },
        },
        "display_transform": {
            "operation": data.spec.transform_kind,
            "output_shape_yx": list(data.spec.output_shape_yx),
            "display_to_native_xy_homogeneous_matrix": [list(row) for row in data.spec.display_to_native_xy],
            "label_interpolation": "nearest neighbor only",
            "same_transform_for_all_masks": True,
            "output_field_um_wh": [
                data.spec.output_shape_yx[1] * data.spec.um_per_px,
                data.spec.output_shape_yx[0] * data.spec.um_per_px,
            ],
        },
        "external_tissue_outline": _matched_tissue_outline_receipt(data),
        "ventricular_lumen_outline": _matched_lumen_receipt(data, ventricle_hil_dir),
    }


def _matched_provenance(
    datasets: Sequence[tuple[str, MatchedCartoonData, AnalysisStatus]], outputs: Sequence[Path],
    ventricle_hil_dir: Path | str | None,
) -> dict[str, Any]:
    root = datasets[0][1].root
    return {
        "schema_version": 9,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {"path": _relative(Path(__file__), root), "sha256": sha256_file(Path(__file__))},
        "panel_map": {
            "B": "Water marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
            "C": "Allulose marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
        },
        "matched_panel_grammar": {
            "views_in_order": ["NPY-positive DAPI nuclei", "c-FOS-positive DAPI nuclei", "combined/dual status"],
            "same_physical_canvas": True,
            "same_orientation_policy": "fixed Panel-A-derived transforms into one 866x1374 physical common frame",
            "same_colors": True,
            "same_title_and_scale_bar_grammar": True,
            "panel_title_x_figure_fraction": 0.5,
            "panel_title_horizontal_alignment": "center",
            "in_panel_method_footnotes": False,
            "in_panel_status_counts": False,
            "view_title_text": ["NPY-positive nuclei", "c-FOS-positive nuclei", "Combined status"],
            "count_location": "Figure 3 caption, panel legends, and provenance only",
            "method_location": "Figure 3 caption and panel legends",
        },
        "samples": {
            data.spec.condition: _matched_sample_receipt(data, status, ventricle_hil_dir)
            for _, data, status in datasets
        },
        "display_semantics": {
            "tissue_envelope": (
                "original registered DAPI-supported mask controls fill; binary_fill_holes supplies exactly one "
                "external gray display path and zero internal gray paths; not an anatomical ROI or denominator"
            ),
            "ventricular_lumen_outline": (
                "opaque clean-white interior and one black display path; hash-bound raw-DAPI source components; Water alone adds "
                "registered DAPI-label-clearance superior completion; both final lumen polygons are accepted in the dedicated "
                "HIL reviewer and bound to their source images, raster masks, reviewers and timestamps; not quantified"
            ),
            "colors": CHANNEL_COLORS | {"dual": DUAL_COLOR, "tissue": TISSUE_FILL, "lumen_outline": "black"},
        },
        "outputs": [
            {"path": _relative(path, root), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in outputs
        ],
    }


def _matched_write_csv(
    path: Path, datasets: Sequence[tuple[str, MatchedCartoonData, AnalysisStatus]], outputs: Sequence[Path],
    ventricle_hil_dir: Path | str | None,
) -> None:
    rows: list[dict[str, Any]] = []
    for panel, data, status in datasets:
        visible_masks, visible_status = _matched_visible_counts(data, status)
        labels_by_role = {
            "dapi_seg": ("DAPI", data.dapi_labels),
            "npy_seg": ("NPY", data.npy_labels),
            "cfos_seg": ("c-FOS", data.cfos_labels),
        }
        for role, source in data.sources.items():
            marker_labels = labels_by_role.get(role)
            marker = marker_labels[0] if marker_labels else None
            labels = marker_labels[1] if marker_labels else None
            rows.append({
                "panel": panel, "condition": data.spec.condition, "sample_id": data.spec.sample_id,
                "role": role, "path": _relative(source.path, data.root),
                "sha256": source.sha256, "bytes": source.bytes,
                "native_shape_yx": "x".join(map(str, data.spec.native_shape_yx)),
                "whole_field_instances": int(labels.max()) if labels is not None else "",
                "display_field_instances": visible_masks[marker] if marker else "",
                "notes": "native label source; whole-field classification" if labels is not None else "native TIFF companion",
            })
        tissue_outline = _matched_tissue_outline_receipt(data)
        rows.append({
            "panel": panel, "condition": data.spec.condition, "sample_id": data.spec.sample_id,
            "role": "external_tissue_display_outline", "path": "",
            "sha256": tissue_outline["outline_binary_sha256"], "bytes": "",
            "native_shape_yx": "x".join(map(str, data.spec.output_shape_yx)),
            "whole_field_instances": 1, "display_field_instances": 1,
            "notes": json.dumps(tissue_outline, sort_keys=True),
        })
        lumen = _matched_lumen_receipt(data, ventricle_hil_dir)
        rows.append({
            "panel": panel, "condition": data.spec.condition, "sample_id": data.spec.sample_id,
            "role": "ventricular_lumen_display_annotation", "path": "",
            "sha256": lumen["selected_binary_sha256"], "bytes": "",
            "native_shape_yx": "x".join(map(str, data.spec.output_shape_yx)),
            "whole_field_instances": 1, "display_field_instances": 1,
            "notes": json.dumps(lumen, sort_keys=True),
        })
        rows.append({
            "panel": panel, "condition": data.spec.condition, "sample_id": data.spec.sample_id,
            "role": "display_field_analysis_counts", "path": "", "sha256": "", "bytes": "",
            "native_shape_yx": "x".join(map(str, data.spec.output_shape_yx)),
            "whole_field_instances": "", "display_field_instances": visible_status["dual_positive_nuclei"],
            "notes": json.dumps(visible_status, sort_keys=True),
        })
    for output in outputs:
        rows.append({
            "panel": "B/C", "condition": "matched Water/Allulose", "sample_id": "",
            "role": "rendered_panel", "path": _relative(output, datasets[0][1].root),
            "sha256": sha256_file(output), "bytes": output.stat().st_size,
            "native_shape_yx": "", "whole_field_instances": "", "display_field_instances": "",
            "notes": "matched 1x3 status output; method prose intentionally excluded from panel art",
        })
    fields = (
        "panel", "condition", "sample_id", "role", "path", "sha256", "bytes",
        "native_shape_yx", "whole_field_instances", "display_field_instances", "notes",
    )
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def render_all(
    root: Path | str | None = None,
    output_dir: Path | str | None = None,
    dpi: int = 600,
    ventricle_hil_dir: Path | str | None = None,
) -> Mapping[str, Path]:
    water = _matched_load_inputs(WATER_CARTOON_SPEC, root)
    allulose = _matched_load_inputs(ALLULOSE_CARTOON_SPEC, root)
    water_status = _matched_classify(water)
    allulose_status = _matched_classify(allulose)
    outdir = Path(output_dir).expanduser().resolve() if output_dir else water.root / "analyses/Fig3/results/final_run/panels"
    outdir.mkdir(parents=True, exist_ok=True)
    output = {
        "panel_b_png": outdir / f"{PANEL_B_STEM}.png",
        "panel_b_pdf": outdir / f"{PANEL_B_STEM}.pdf",
        "panel_c_png": outdir / f"{PANEL_C_STEM}.png",
        "panel_c_pdf": outdir / f"{PANEL_C_STEM}.pdf",
        "provenance_json": outdir / f"{PROVENANCE_STEM}.json",
        "provenance_csv": outdir / f"{PROVENANCE_STEM}.csv",
    }
    datasets = (("B", water, water_status), ("C", allulose, allulose_status))
    for panel, data, status in datasets:
        fig = _matched_status_cartoon(data, status, panel, ventricle_hil_dir)
        try:
            for suffix in ("png", "pdf"):
                _atomic_save(fig, output[f"panel_{panel.lower()}_{suffix}"], dpi)
        finally:
            plt.close(fig)
    panel_paths = tuple(output[key] for key in ("panel_b_png", "panel_b_pdf", "panel_c_png", "panel_c_pdf"))
    _atomic_text(
        output["provenance_json"],
        json.dumps(
            _matched_provenance(datasets, panel_paths, ventricle_hil_dir),
            indent=2,
            sort_keys=True,
        ) + "\n",
    )
    _matched_write_csv(
        output["provenance_csv"], datasets, panel_paths, ventricle_hil_dir
    )
    return output

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument(
        "--ventricle-hil-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing accepted Water/Allulose ventricle HIL masks "
            "and receipts; defaults to Fig3/hil_review/ventricle_cartoon_20260827_v1."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.dpi < 300:
        raise ValueError("Publication cartoon DPI must be at least 300")
    for role, path in render_all(
        args.root, args.output_dir, args.dpi, args.ventricle_hil_dir
    ).items():
        print(f"{role}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
