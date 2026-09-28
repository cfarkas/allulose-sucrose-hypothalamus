#!/usr/bin/env python3
"""Render Supplementary Figure S5 from the reviewed July 2026 ACTH/CLIP analysis.

Supplementary Figure S5 includes eight WT animals and male NPY-transgenic NPY-M
(Water, Sucrose and Allulose n=3 each), as confirmed by the author on 2026-09-08. All publication graphics and bilingual
isolated panels are fixed at 600 dpi. ACTH/CLIP Cellpose objects smaller than
the acquisition-specific DAPI nuclear-area fifth percentile are omitted from
cell calls and from the ROI-gated microscopy panel; source masks are never modified.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import shutil
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import tifffile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.ndimage import binary_fill_holes, gaussian_filter
from skimage.segmentation import find_boundaries
from PIL import Image, ImageDraw

Image.MAX_IMAGE_PIXELS = None


DPI = 600
PANEL_SIZE = (6.2, 3.75)
SPATIAL_PANEL_SIZE = (8.1, 6.25)
N_RADIAL_SHELLS = 6
CONDITIONS = ("Water", "Sucrose", "Allulose")
CONDITION_COLORS = {
    "Water": "#9bd7e8",
    "Sucrose": "#e33c3c",
    "Allulose": "#36c875",
}
MEAN_LINE_COLORS = {
    "Water": "#2f91b8",
    "Sucrose": "#c51b1f",
    "Allulose": "#159447",
}
CONDITION_ES = {"Water": "Agua", "Sucrose": "Sacarosa", "Allulose": "Alulosa"}
TINTS = {
    "DAPI": np.asarray((0.22, 0.45, 1.00), dtype=np.float32),
    "cFOS": np.asarray((1.00, 0.05, 0.72), dtype=np.float32),
    "ACTH_CLIP": np.asarray((1.00, 0.58, 0.08), dtype=np.float32),
}
REGION_COLORS = {"ARC": "#00adb5", "ME": "#e67e22", "VMN": "#be3eae"}
SOURCE_TABLE_NAMES = (
    "annotation_status.csv",
    "per_analysis_unit_region_summary.csv",
    "per_animal_region_summary.csv",
    "per_image_region_summary.csv",
    "per_nucleus_assignments.csv",
    "marker_roi_assignments.csv",
    "acth_clip_dapi_size_filter_audit.csv",
    "count_invariants.csv",
    "endpoint_values_long.csv",
    "statistics_all_values.csv",
    "analysis_readme.txt",
)
PANEL_NAMES = {
    "A": ("DAPI", "DAPI"),
    "B": ("c-FOS", "c-FOS"),
    "C": ("ACTH/CLIP", "ACTH/CLIP"),
    "D": ("Cellpose and HIL anatomy", "Cellpose y anatomía HIL"),
    "E": ("Merged channels", "Canales combinados"),
    "F": ("ACTH/CLIP ROI microscopy", "Microscopía de ROIs ACTH/CLIP"),
    "G": ("Total c-FOS activation", "Activación total de c-FOS"),
    "H": ("ACTH/CLIP-cell activation", "Activación de células ACTH/CLIP"),
    "I": ("Spatial c-FOS occurrence", "Ocurrencia espacial de c-FOS"),
    "J": ("Spatial double-positive occurrence", "Ocurrencia espacial doble positiva"),
}
EXPECTED_REVIEW_MIP_PIXEL_SHA256 = {
    "2026-07-09_FR4-1_Allulose": "9b57da6d433cac9d20ae884a41a65f8d562dd418d0f6f5f6a7993dcbd004aec9",
    "2026-07-09_FR5-1_Allulose": "a62b70496992dd96280b343b5921b043493bae2e722a0f6a4473db097d97c328",
    "2026-07-09_FR5-2_Allulose": "947123b99140e98bb83845e269c390cc8f36fea5e029f4d76330f57e9c604b8d",
    "2026-07-09_FR5-3_Sucrose": "2883e699f5a96e1fc501d803035536e123648bd3411c8dcd97f58d5c34a0e536",
    "2026-07-09_FR5-4_Sucrose": "6e2e64fafb224c7ae504dd9b629032c25a009a4ebe380ac37f76f658fcb41cfc",
    "2026-07-10_AGUA-M_Water": "f1bd18a6835c05c6ba7787fb4b614252864eab0228711b5d50526347b4691df7",
    "2026-07-10_FR4-3_Sucrose": "39a76fe442d192f37773c07770d538d3e7eb903f5f9d87aff147bed7773f6758",
    "2026-07-10_FR5-5_Water": "89a70d2ead3b5974c9804202259d020290e604244fe350a5b67959247a01f1d5",
    "2026-07-10_NPY-M_Water": "c5629347d9333cd5e62f77b8458e59d239fbf0a1558e724f6ee8375c36ffc66f",
}


def expected_review_inventory(count: int) -> set[str]:
    full = set(EXPECTED_REVIEW_MIP_PIXEL_SHA256)
    if count == 9:
        return full
    if count == 8:
        return full - {"2026-07-10_NPY-M_Water"}
    raise ValueError(f"Unsupported S5 cohort size: {count}")



def sha256_file(path: Path, block_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def resolve_work_path(workdir: Path, value: Any) -> Path:
    """Resolve a work-tree path after staging or repository relocation.

    Fresh builds record paths relative to workdir. Historical accepted
    analyses recorded absolute paths inside their temporary build directory;
    those are safely re-rooted only when their suffix begins at a known
    Supplementary Figure S5 work-tree directory.
    """
    raw = Path(str(value)).expanduser()
    candidates: list[Path] = []
    if raw.is_absolute():
        candidates.append(raw)
    else:
        candidates.extend((workdir / raw, raw))
    for marker in ("prepared", "annotations", "region_masks", "results"):
        if marker in raw.parts:
            candidates.append(
                workdir.joinpath(*raw.parts[raw.parts.index(marker) :])
            )
            break
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Supplementary Figure S5 work-tree file is missing: {value!r}; "
        f"resolved relative to {workdir}"
    )


def validate_spatial_tissue_hil_source(
    directory: Path,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, Any]]:
    """Validate frozen I/J HIL without the excluded generated analyses tree."""
    lexical = directory.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_dir():
        raise FileNotFoundError(
            f"Finalized spatial tissue HIL directory is missing or linked: {lexical}"
        )
    root = lexical.resolve()
    receipt_path = root / "SPATIAL_TISSUE_HIL_RECEIPT.json"
    manifest_path = root / "SPATIAL_TISSUE_HIL_MANIFEST.csv"
    if any(path.is_symlink() or not path.is_file() for path in (receipt_path, manifest_path)):
        raise FileNotFoundError(
            f"Finalized spatial tissue HIL receipt/manifest missing or linked: {root}"
        )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    schema = str(receipt.get("schema", ""))
    common_valid = (
        receipt.get("status") == "human_accepted"
        and receipt.get("panels") == ["I", "J"]
        and receipt.get("channels_used") == ["DAPI"]
        and receipt.get("accepted_acquisitions") in {8, 9}
        and receipt.get("pending_acquisitions") == 0
        and receipt.get("raw_images_modified") is False
        and receipt.get("raw_segmentations_modified") is False
        and receipt.get("manifest") == manifest_path.name
    )
    legacy_valid = (
        schema in {"figs4_spatial_tissue_hil_v1", "figs5_spatial_tissue_hil_v1"}
        and receipt.get("ventricular_contour_available") is False
        and receipt.get("one_outer_polygon_per_acquisition") is True
    )
    current_valid = (
        schema == "figs5_spatial_3v_hil_v2"
        and receipt.get("ventricular_contour_available") is True
        and receipt.get("one_ventricle_polygon_per_acquisition") is True
        and receipt.get("used_as_3v_exclusion_mask") is True
        and receipt.get("used_for_tissue_geometry") is False
    )
    if not common_valid or not (legacy_valid or current_valid):
        raise ValueError("Spatial tissue HIL receipt violates the S5 I/J contract")
    if sha256_file(manifest_path) != str(receipt.get("manifest_sha256", "")):
        raise ValueError("Spatial tissue HIL manifest hash mismatch")

    manifest = pd.read_csv(manifest_path)
    if int(receipt.get("accepted_acquisitions", -1)) != len(manifest):
        raise ValueError("Spatial HIL receipt and manifest counts differ")
    mask_path_column = "tissue_mask_path" if legacy_valid else "ventricle_mask_path"
    mask_hash_column = "tissue_mask_sha256" if legacy_valid else "ventricle_mask_sha256"
    pixel_count_column = "tissue_pixels" if legacy_valid else "ventricle_pixels"
    required = {
        "acquisition_id", "sample_id", "animal_id", "condition", "decision",
        "annotation_path", "annotation_sha256", mask_path_column, mask_hash_column,
        "image_height", "image_width", pixel_count_column, "channels_used",
    }
    if required - set(manifest.columns):
        raise ValueError("Spatial tissue HIL manifest columns differ")
    acquisitions = manifest["acquisition_id"].astype(str)
    if (
        len(manifest) not in {8, 9}
        or acquisitions.duplicated().any()
        or set(acquisitions) != expected_review_inventory(len(manifest))
        or manifest["condition"].astype(str).value_counts().to_dict()
        != {"Allulose": 3, "Sucrose": 3, "Water": len(manifest) - 6}
        or (manifest["decision"].astype(str) != "include").any()
        or (manifest["channels_used"].astype(str) != "DAPI only").any()
    ):
        raise ValueError("Spatial tissue HIL source inventory differs")

    def source_member(value: Any, label: str) -> Path:
        relative = Path(str(value))
        if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
            raise ValueError(f"Unsafe spatial HIL {label} path: {value!r}")
        target = root.joinpath(*relative.parts)
        if (
            target.is_symlink()
            or not target.is_file()
            or not target.resolve().is_relative_to(root)
        ):
            raise FileNotFoundError(f"Spatial HIL {label} is missing or unsafe: {target}")
        return target

    masks: dict[str, np.ndarray] = {}
    contours: dict[str, np.ndarray] = {}
    for _, row in manifest.iterrows():
        acquisition = str(row["acquisition_id"])
        annotation = source_member(row["annotation_path"], "annotation")
        mask_path = source_member(row[mask_path_column], "mask")
        if sha256_file(annotation) != str(row["annotation_sha256"]):
            raise ValueError(f"{acquisition}: spatial HIL annotation hash mismatch")
        if sha256_file(mask_path) != str(row[mask_hash_column]):
            raise ValueError(f"{acquisition}: spatial HIL mask hash mismatch")
        record = json.loads(annotation.read_text(encoding="utf-8"))
        semantics_valid = (
            legacy_valid
            and record.get("region") == "OUTER_TISSUE"
            and record.get("ventricular_contour_available") is False
        ) or (
            current_valid
            and record.get("region") == "3V"
            and record.get("ventricular_contour_available") is True
            and record.get("used_for_tissue_geometry") is False
            and record.get("used_as_3v_exclusion_mask") is True
        )
        if (
            record.get("acquisition_id") != acquisition
            or record.get("sample_id") != str(row["sample_id"])
            or record.get("animal_id") != str(row["animal_id"])
            or record.get("condition") != str(row["condition"])
            or record.get("decision") != "include"
            or record.get("channels_used") != ["DAPI"]
            or not semantics_valid
        ):
            raise ValueError(f"{acquisition}: invalid spatial HIL annotation contract")
        points = np.asarray(record.get("points", []), dtype=float)
        width = int(record.get("image_width", 0))
        height = int(record.get("image_height", 0))
        if (
            points.ndim != 2
            or points.shape[0] < 3
            or points.shape[1] != 2
            or not np.isfinite(points).all()
            or width != int(row["image_width"])
            or height != int(row["image_height"])
            or width <= 0
            or height <= 0
            or float(points[:, 0].min()) < -1.0
            or float(points[:, 1].min()) < -1.0
            or float(points[:, 0].max()) > width + 1.0
            or float(points[:, 1].max()) > height + 1.0
        ):
            raise ValueError(f"{acquisition}: invalid spatial HIL point contour")
        mask = np.squeeze(np.asarray(tifffile.imread(str(mask_path)), dtype=np.uint8))
        if (
            mask.ndim != 2
            or mask.shape != (height, width)
            or not np.any(mask)
            or int(np.count_nonzero(mask)) != int(row[pixel_count_column])
        ):
            raise ValueError(f"{acquisition}: spatial HIL mask shape/content differs")
        masks[acquisition] = mask > 0
        contours[acquisition] = points
    return masks, contours, receipt


def load_spatial_tissue_hil(
    directory: Path, sheet: pd.DataFrame, workdir: Path
) -> tuple[
    dict[str, np.ndarray], dict[str, np.ndarray], dict[str, Any]
]:
    """Hash-validate legacy or current DAPI HIL and return 3V masks/contours.

    The finalized 2026-08-28 v1 receipt is immutable and historically labels
    these polygons as OUTER_TISSUE. A visual/raw-DAPI audit established that
    the polygons follow the 3V cavity. Their raster masks are therefore used
    only as negative/exclusion masks for I/J spatial maps and analyses, never
    as positive tissue masks or as crop/normalization geometry. Point contours
    remain the visible 3V outlines.
    """
    receipt_path = directory / "SPATIAL_TISSUE_HIL_RECEIPT.json"
    manifest_path = directory / "SPATIAL_TISSUE_HIL_MANIFEST.csv"
    if not receipt_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(
            f"Finalized spatial tissue HIL receipt/manifest missing: {directory}"
        )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    schema = str(receipt.get("schema", ""))
    common_valid = (
        receipt.get("status") == "human_accepted"
        and receipt.get("panels") == ["I", "J"]
        and receipt.get("channels_used") == ["DAPI"]
    )
    legacy_valid = (
        schema in {"figs4_spatial_tissue_hil_v1", "figs5_spatial_tissue_hil_v1"}
        and receipt.get("ventricular_contour_available") is False
        and receipt.get("one_outer_polygon_per_acquisition") is True
    )
    current_valid = (
        schema == "figs5_spatial_3v_hil_v2"
        and receipt.get("ventricular_contour_available") is True
        and receipt.get("one_ventricle_polygon_per_acquisition") is True
        and receipt.get("used_as_3v_exclusion_mask") is True
        and receipt.get("used_for_tissue_geometry") is False
    )
    if not common_valid or not (legacy_valid or current_valid):
        raise ValueError(
            "Spatial tissue HIL receipt violates the S5 I/J contract"
        )
    if sha256_file(manifest_path) != str(
        receipt.get("manifest_sha256", "")
    ):
        raise ValueError("Spatial tissue HIL manifest hash mismatch")
    manifest = pd.read_csv(manifest_path)
    if int(receipt.get("accepted_acquisitions", -1)) != len(manifest):
        raise ValueError("Spatial HIL receipt and manifest counts differ")
    expected = set(sheet["acquisition_id"].astype(str))
    if expected != expected_review_inventory(len(sheet)):
        raise ValueError("Canonical reviewed DAPI MIP pixel-hash inventory differs from S5 sheet")
    observed = set(manifest["acquisition_id"].astype(str))
    if observed != expected or len(manifest) != len(sheet):
        raise ValueError(
            f"Spatial tissue HIL acquisition mismatch: "
            f"expected={sorted(expected)} observed={sorted(observed)}"
        )
    sheet_by_id = sheet.copy()
    sheet_by_id.index = sheet_by_id["acquisition_id"].astype(str)
    masks: dict[str, np.ndarray] = {}
    contours: dict[str, np.ndarray] = {}
    mask_path_column = (
        "tissue_mask_path" if legacy_valid else "ventricle_mask_path"
    )
    mask_hash_column = (
        "tissue_mask_sha256" if legacy_valid else "ventricle_mask_sha256"
    )
    for _, record in manifest.iterrows():
        acquisition = str(record["acquisition_id"])
        annotation = directory / str(record["annotation_path"])
        mask_path = directory / str(record[mask_path_column])
        if sha256_file(annotation) != str(record["annotation_sha256"]):
            raise ValueError(
                f"{acquisition}: spatial HIL annotation hash mismatch"
            )
        if sha256_file(mask_path) != str(record[mask_hash_column]):
            raise ValueError(
                f"{acquisition}: spatial HIL contour-mask hash mismatch"
            )
        annotation_record = json.loads(
            annotation.read_text(encoding="utf-8")
        )
        legacy_annotation = (
            legacy_valid
            and annotation_record.get("region") == "OUTER_TISSUE"
            and annotation_record.get(
                "ventricular_contour_available"
            ) is False
        )
        current_annotation = (
            current_valid
            and annotation_record.get("region") == "3V"
            and annotation_record.get(
                "ventricular_contour_available"
            ) is True
            and annotation_record.get("used_for_tissue_geometry") is False
            and annotation_record.get("used_as_3v_exclusion_mask") is True
        )
        if (
            annotation_record.get("channels_used") != ["DAPI"]
            or not (legacy_annotation or current_annotation)
        ):
            raise ValueError(
                f"{acquisition}: invalid spatial HIL annotation semantics"
            )
        points = np.asarray(annotation_record.get("points", []), dtype=float)
        width = int(annotation_record.get("image_width", 0))
        height = int(annotation_record.get("image_height", 0))
        if (
            points.ndim != 2
            or points.shape[0] < 3
            or points.shape[1] != 2
            or not np.isfinite(points).all()
            or width <= 0
            or height <= 0
            or float(points[:, 0].min()) < -1.0
            or float(points[:, 1].min()) < -1.0
            or float(points[:, 0].max()) > width + 1.0
            or float(points[:, 1].max()) > height + 1.0
        ):
            raise ValueError(
                f"{acquisition}: invalid spatial HIL point contour"
            )
        mask = np.squeeze(
            np.asarray(tifffile.imread(str(mask_path)), dtype=np.uint8)
        )
        dapi_path = resolve_work_path(
            workdir, sheet_by_id.loc[acquisition, "dapi_mip"]
        )
        dapi_seg_path = resolve_work_path(
            workdir, sheet_by_id.loc[acquisition, "dapi_seg"]
        )
        dapi = load_image(dapi_path)
        if sha256_file(dapi_path) != str(record["source_dapi_mip_sha256"]):
            actual_pixel_hash = hashlib.sha256(
                np.ascontiguousarray(dapi).tobytes()
            ).hexdigest()
            expected_pixel_hash = EXPECTED_REVIEW_MIP_PIXEL_SHA256.get(acquisition)
            if actual_pixel_hash != expected_pixel_hash:
                raise ValueError(
                    f"{acquisition}: reviewed DAPI MIP file and pixel hashes mismatch"
                )
        if sha256_file(dapi_seg_path) != str(
            record["source_dapi_segmentation_sha256"]
        ):
            raise ValueError(
                f"{acquisition}: reviewed DAPI segmentation hash mismatch"
            )
        if mask.shape != dapi.shape or not np.any(mask):
            raise ValueError(
                f"{acquisition}: spatial HIL mask shape/content mismatch"
            )
        if (height, width) != dapi.shape:
            raise ValueError(
                f"{acquisition}: annotation/DAPI dimensions differ"
            )
        masks[acquisition] = mask > 0
        contours[acquisition] = points
    return masks, contours, receipt


def bool_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t"}


def load_image(path: Path) -> np.ndarray:
    array = np.squeeze(np.asarray(tifffile.imread(str(path))))
    while array.ndim > 2:
        array = np.max(array, axis=0)
    if array.ndim != 2:
        raise ValueError(f"Expected 2-D image after projection: {path} {array.shape}")
    return array


def load_segmentation(path: Path) -> np.ndarray:
    raw = np.load(str(path), allow_pickle=True)
    if isinstance(raw, np.ndarray) and raw.dtype == object and raw.shape == ():
        raw = raw.item()
    if isinstance(raw, dict):
        raw = raw["masks"]
    array = np.squeeze(np.asarray(raw))
    while array.ndim > 2:
        array = array[0]
    if array.ndim != 2:
        raise ValueError(f"Expected 2-D segmentation: {path} {array.shape}")
    return array.astype(np.int32, copy=False)


def normalize_channel(array: np.ndarray, low: float, high: float, gamma: float) -> np.ndarray:
    values = np.asarray(array, dtype=np.float32)
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.zeros(values.shape, dtype=np.float32)
    lo, hi = np.percentile(finite, (low, high))
    if not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    scaled = np.clip((values - lo) / (hi - lo), 0.0, 1.0)
    return np.power(scaled, gamma, dtype=np.float32)


def tinted(channel: np.ndarray, marker: str) -> np.ndarray:
    return np.clip(channel[..., None] * TINTS[marker], 0.0, 1.0)


def robust_crop(mask: np.ndarray, pixel_size_um: float, padding_um: float = 130.0) -> tuple[int, int, int, int]:
    yy, xx = np.nonzero(mask)
    if not yy.size:
        raise ValueError("Representative HIL region mask is empty")
    padding = max(20, int(round(padding_um / pixel_size_um)))
    y0 = max(0, int(yy.min()) - padding)
    y1 = min(mask.shape[0], int(yy.max()) + padding + 1)
    x0 = max(0, int(xx.min()) - padding)
    x1 = min(mask.shape[1], int(xx.max()) + padding + 1)
    target_aspect = 1.42
    height = y1 - y0
    width = x1 - x0
    if width / max(height, 1) < target_aspect:
        wanted = int(round(height * target_aspect))
        extra = wanted - width
        x0 = max(0, x0 - extra // 2)
        x1 = min(mask.shape[1], x1 + extra - extra // 2)
    else:
        wanted = int(round(width / target_aspect))
        extra = wanted - height
        y0 = max(0, y0 - extra // 2)
        y1 = min(mask.shape[0], y1 + extra - extra // 2)
    return y0, y1, x0, x1


def condition_label(condition: str, language: str) -> str:
    return CONDITION_ES.get(condition, condition) if language == "es" else condition


def subtitle(context: dict[str, Any], language: str) -> str:
    condition = condition_label(context["condition"], language)
    if language == "es":
        return f"{context['sample_id']} · {condition} · ayuno 16 h · primera exposición"
    return f"{context['sample_id']} · {condition} · 16 h fast · first exposure"


def add_figure_letter(fig: plt.Figure, letter: str, color: str = "black") -> None:
    """Place every panel letter at one fixed figure-relative coordinate."""
    fig.text(
        0.018, 0.925, letter, ha="left", va="top",
        fontsize=20, fontweight="bold", color=color, zorder=1000,
    )


def add_panel_letter(ax: plt.Axes, letter: str, color: str = "white") -> None:
    add_figure_letter(ax.figure, letter, color)


def add_scale_bar(
    ax: plt.Axes,
    shape: tuple[int, int],
    pixel_size_um: float,
    length_um: float = 100.0,
    *,
    y_fraction: float | None = None,
) -> None:
    pixels = int(round(length_um / pixel_size_um))
    if pixels <= 0 or pixels >= shape[1] * 0.6:
        return
    x1 = shape[1] - max(18, int(shape[1] * 0.04))
    x0 = x1 - pixels
    y = (
        shape[0] - max(18, int(shape[0] * 0.06))
        if y_fraction is None
        else shape[0] * y_fraction
    )
    ax.plot([x0, x1], [y, y], color="white", linewidth=4.2, solid_capstyle="butt", zorder=15)
    ax.text((x0 + x1) / 2, y - max(8, shape[0] * 0.018), f"{length_um:g} µm",
            color="white", ha="center", va="bottom", fontsize=9, fontweight="bold", zorder=15)


def save_figure(fig: plt.Figure, png: Path, pdf: Path) -> None:
    fig.savefig(png, dpi=DPI, facecolor=fig.get_facecolor())
    fig.savefig(pdf, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)


def choose_representative(sheet: pd.DataFrame, summary: pd.DataFrame) -> pd.Series:
    allulose = summary[
        (summary["condition"] == "Allulose")
        & summary["region_present"].map(bool_value)
    ].copy()
    grouped = allulose.groupby("acquisition_id", as_index=False).agg(
        n_double=("n_double", "sum"),
        n_total_acth_clip=("n_total_acth_clip", "sum"),
    )
    grouped["activation"] = grouped["n_double"] / grouped["n_total_acth_clip"].replace(0, np.nan)
    finite = grouped["activation"].dropna()
    if finite.empty:
        candidates = sorted(sheet.loc[sheet["condition"] == "Allulose", "acquisition_id"].astype(str))
        chosen = candidates[0]
    else:
        median = float(finite.median())
        grouped["distance"] = (grouped["activation"] - median).abs()
        chosen = str(grouped.sort_values(["distance", "acquisition_id"]).iloc[0]["acquisition_id"])
    row = sheet[sheet["acquisition_id"].astype(str) == chosen]
    if len(row) != 1:
        raise ValueError(f"Representative acquisition not unique: {chosen}")
    return row.iloc[0]


def load_context(workdir: Path) -> dict[str, Any]:
    sheet_path = workdir / "analysis_samplesheet.csv"
    result_dir = workdir / "results"
    required = [
        sheet_path,
        result_dir / "annotation_status.csv",
        result_dir / "per_nucleus_assignments.csv",
        result_dir / "per_analysis_unit_region_summary.csv",
        result_dir / "marker_roi_assignments.csv",
        result_dir / "statistics_all_values.csv",
        result_dir / "acth_clip_dapi_size_filter_audit.csv",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Run reviewed Supplementary Figure S5 quantification first; missing: {missing}")

    sheet = pd.read_csv(sheet_path)
    if len(sheet) not in {8, 9} or sheet["animal_id"].nunique() != len(sheet):
        raise ValueError("S5 requires eight original animals or those animals plus author-confirmed NPY-M")
    if set(sheet["acquisition_id"]) != expected_review_inventory(len(sheet)):
        raise ValueError("S5 acquisition inventory differs from the reviewed cohort")
    npy = sheet["acquisition_id"].eq("2026-07-10_NPY-M_Water")
    if not sheet.loc[~npy, "genotype"].eq("WT").all():
        raise ValueError("The eight original S5 animals must retain their WT assignment")
    if npy.any():
        row = sheet.loc[npy].iloc[0]
        if row["genotype"] != "NPY-Tg" or row["condition"] != "Water":
            raise ValueError("NPY-M genotype/condition differs from the author's confirmation")
    counts = sheet["condition"].value_counts().to_dict()
    if counts != {"Allulose": 3, "Sucrose": 3, "Water": len(sheet) - 6}:
        raise ValueError(f"Unexpected Supplementary Figure S5 cohort: {counts}")

    status = pd.read_csv(result_dir / "annotation_status.csv")
    if status["status"].astype(str).eq("pending").any():
        raise ValueError("Supplementary Figure S5 HIL review still has pending acquisitions")
    nuclei = pd.read_csv(result_dir / "per_nucleus_assignments.csv")
    summary = pd.read_csv(result_dir / "per_analysis_unit_region_summary.csv")
    rois = pd.read_csv(result_dir / "marker_roi_assignments.csv")
    statistics = pd.read_csv(result_dir / "statistics_all_values.csv")
    size_audit = pd.read_csv(result_dir / "acth_clip_dapi_size_filter_audit.csv")
    representative = choose_representative(sheet, summary)
    acquisition = str(representative["acquisition_id"])

    dapi_image = load_image(resolve_work_path(workdir, representative["dapi_mip"]))
    cfos_image = load_image(resolve_work_path(workdir, representative["cfos_mip"]))
    acth_image = load_image(resolve_work_path(workdir, representative["acth_clip_mip"]))
    dapi_seg = load_segmentation(resolve_work_path(workdir, representative["dapi_seg"]))
    cfos_seg = load_segmentation(resolve_work_path(workdir, representative["cfos_seg"]))
    acth_seg = load_segmentation(resolve_work_path(workdir, representative["acth_clip_seg"]))
    shapes = {array.shape for array in (dapi_image, cfos_image, acth_image, dapi_seg, cfos_seg, acth_seg)}
    if len(shapes) != 1:
        raise ValueError(f"Representative image/mask shape mismatch: {shapes}")

    region_path = workdir / "region_masks" / f"{acquisition}_regions.tif"
    if not region_path.is_file():
        raise FileNotFoundError(region_path)
    region_mask = np.squeeze(np.asarray(tifffile.imread(str(region_path)), dtype=np.uint8))
    if region_mask.shape != dapi_image.shape:
        raise ValueError(f"Region mask shape mismatch: {region_mask.shape} vs {dapi_image.shape}")
    pixel_size = float(representative["pixel_size_x_um"])
    bounds = robust_crop(region_mask > 0, pixel_size)

    dapi_norm = normalize_channel(dapi_image, 1.0, 99.8, 0.72)
    cfos_norm = normalize_channel(cfos_image, 1.0, 99.75, 0.78)
    acth_norm = normalize_channel(acth_image, 1.0, 99.75, 0.82)
    channels = {
        "DAPI": tinted(dapi_norm, "DAPI"),
        "cFOS": tinted(cfos_norm, "cFOS"),
        "ACTH_CLIP": tinted(acth_norm, "ACTH_CLIP"),
    }
    merged = np.clip(channels["DAPI"] + channels["cFOS"] + channels["ACTH_CLIP"], 0.0, 1.0)

    rep_rois = rois[
        (rois["acquisition_id"].astype(str) == acquisition)
        & (rois["marker"].astype(str) == "ACTH_CLIP")
    ].copy()
    accepted_labels = rep_rois.loc[
        rep_rois["accepted"].map(bool_value), "roi_label"
    ].astype(int).to_numpy()
    retained_mask = np.isin(acth_seg, accepted_labels)
    gated_acth = acth_norm * retained_mask.astype(np.float32)
    gated_merge = np.clip(
        channels["DAPI"] * 0.38 + channels["cFOS"] * 0.58
        + tinted(gated_acth, "ACTH_CLIP"), 0.0, 1.0,
    )

    return {
        "workdir": workdir,
        "result_dir": result_dir,
        "sheet": sheet,
        "status": status,
        "nuclei": nuclei,
        "summary": summary,
        "rois": rois,
        "statistics": statistics,
        "size_audit": size_audit,
        "representative": representative,
        "acquisition_id": acquisition,
        "sample_id": str(representative["sample_id"]),
        "condition": str(representative["condition"]),
        "pixel_size_um": pixel_size,
        "bounds": bounds,
        "dapi_image": dapi_image,
        "cfos_image": cfos_image,
        "acth_image": acth_image,
        "dapi_seg": dapi_seg,
        "cfos_seg": cfos_seg,
        "acth_seg": acth_seg,
        "region_mask": region_mask,
        "channels": channels,
        "merged": merged,
        "retained_acth_mask": retained_mask,
        "gated_merge": gated_merge,
        "accepted_acth_labels": accepted_labels,
    }


def crop(array: np.ndarray, bounds: tuple[int, int, int, int]) -> np.ndarray:
    y0, y1, x0, x1 = bounds
    return np.asarray(array[y0:y1, x0:x1])


def image_panel(context: dict[str, Any], letter: str, marker: str, language: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    fig.patch.set_facecolor("black")
    image = crop(context["channels"][marker], context["bounds"])
    ax.imshow(image, interpolation="nearest")
    add_panel_letter(ax, letter)
    title = PANEL_NAMES[letter][1 if language == "es" else 0]
    ax.set_title(f"{title}\n{subtitle(context, language)}", color="white", fontsize=12, fontweight="bold", pad=7)
    add_scale_bar(ax, image.shape[:2], context["pixel_size_um"])
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.87, bottom=0.01)
    return fig


def draw_region_contours(ax: plt.Axes, mask: np.ndarray, dashed: bool = True) -> None:
    for code, name in ((1, "ARC"), (2, "ME"), (3, "VMN")):
        binary = mask == code
        if binary.any():
            ax.contour(
                binary.astype(float), levels=[0.5], colors=[REGION_COLORS[name]],
                linewidths=1.35, linestyles="--" if dashed else "-",
            )


def cartoon_panel(context: dict[str, Any], language: str) -> plt.Figure:
    y0, y1, x0, x1 = context["bounds"]
    dapi = context["dapi_seg"][y0:y1, x0:x1]
    regions = context["region_mask"][y0:y1, x0:x1]
    local = context["nuclei"]
    local = local[local["acquisition_id"].astype(str) == context["acquisition_id"]]
    rgb = np.zeros((*dapi.shape, 3), dtype=np.float32)
    rgb[dapi > 0] = (0.12, 0.22, 0.42)
    phenotype_colors = {
        "cFOS_only": (0.95, 0.08, 0.68),
        "ACTH_CLIP_only": (1.00, 0.58, 0.08),
        "cFOS_ACTH_CLIP": (1.00, 1.00, 1.00),
    }
    for phenotype, color in phenotype_colors.items():
        labels = local.loc[local["phenotype"].astype(str) == phenotype, "nucleus_label"].astype(int).to_numpy()
        if labels.size:
            rgb[np.isin(dapi, labels)] = color
    boundary = find_boundaries(dapi, mode="inner")
    rgb[boundary & (dapi > 0)] = np.maximum(rgb[boundary & (dapi > 0)], 0.58)

    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    fig.patch.set_facecolor("black")
    ax.imshow(rgb, interpolation="nearest")
    draw_region_contours(ax, regions)
    add_scale_bar(
        ax, rgb.shape[:2], context["pixel_size_um"], y_fraction=0.80,
    )
    add_panel_letter(ax, "D")
    ax.set_title(PANEL_NAMES["D"][1 if language == "es" else 0], color="white", fontsize=12, fontweight="bold", pad=7)
    labels = [
        ("DAPI", (0.22, 0.45, 1.0)),
        ("c-FOS", phenotype_colors["cFOS_only"]),
        ("ACTH/CLIP", phenotype_colors["ACTH_CLIP_only"]),
        ("c-FOS+ ACTH/CLIP+", phenotype_colors["cFOS_ACTH_CLIP"]),
    ]
    handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=color,
                      markeredgecolor="white", markersize=7, label=label) for label, color in labels]
    handles.extend(Line2D([0], [0], color=REGION_COLORS[name], linestyle="--",
                          linewidth=2.3, label=name) for name in ("ARC", "ME", "VMN"))
    ax.legend(handles=handles, loc="lower left", ncol=3, frameon=False, labelcolor="white", fontsize=6.6)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.91, bottom=0.01)
    return fig


def choose_inset_nucleus(context: dict[str, Any]) -> tuple[float, float, str]:
    cells = context["nuclei"]
    cells = cells[cells["acquisition_id"].astype(str) == context["acquisition_id"]].copy()
    preferred = cells[cells["phenotype"].astype(str) == "cFOS_ACTH_CLIP"].copy()
    label = "c-FOS+ ACTH/CLIP+"
    if preferred.empty:
        preferred = cells[cells["is_acth_clip"].map(bool_value)].copy()
        label = "ACTH/CLIP+"
    if preferred.empty:
        y0, y1, x0, x1 = context["bounds"]
        return (y0 + y1) / 2, (x0 + x1) / 2, "review field"
    y0, y1, x0, x1 = context["bounds"]
    cy, cx = (y0 + y1) / 2, (x0 + x1) / 2
    preferred["distance"] = np.hypot(preferred["centroid_y"] - cy, preferred["centroid_x"] - cx)
    row = preferred.sort_values(["distance", "nucleus_label"]).iloc[0]
    return float(row["centroid_y"]), float(row["centroid_x"]), label


def merged_panel(
    context: dict[str, Any], language: str, roi_gated: bool
) -> plt.Figure:
    letter = "F" if roi_gated else "E"
    full = context["gated_merge"] if roi_gated else context["merged"]
    main = crop(full, context["bounds"])
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    fig.patch.set_facecolor("black")
    ax.imshow(main, interpolation="nearest")
    add_panel_letter(ax, letter)
    title = PANEL_NAMES[letter][1 if language == "es" else 0]
    ax.set_title(f"{title}\n{subtitle(context, language)}", color="white", fontsize=11.2, fontweight="bold", pad=7)
    add_scale_bar(ax, main.shape[:2], context["pixel_size_um"])

    cy, cx, cell_label = choose_inset_nucleus(context)
    half = max(24, int(round(48.0 / context["pixel_size_um"])))
    ylo = max(0, int(round(cy)) - half)
    yhi = min(full.shape[0], int(round(cy)) + half)
    xlo = max(0, int(round(cx)) - half)
    xhi = min(full.shape[1], int(round(cx)) + half)
    detail_image = full[ylo:yhi, xlo:xhi]
    y0, y1, x0, x1 = context["bounds"]
    rect = Rectangle(
        (xlo - x0, ylo - y0), xhi - xlo, yhi - ylo,
        fill=False, edgecolor="white", linewidth=3.2,
        linestyle=(0, (4, 2)), joinstyle="miter",
    )
    ax.add_patch(rect)
    inset = ax.inset_axes([0.035, 0.055, 0.34, 0.34])
    inset.imshow(detail_image, interpolation="nearest")
    inset.set_title(cell_label, color="white", fontsize=7.5, fontweight="bold", pad=2)
    inset.set_facecolor("black")
    add_scale_bar(
        inset, detail_image.shape[:2], context["pixel_size_um"], 20.0,
    )
    for spine in inset.spines.values():
        spine.set_edgecolor("white")
        spine.set_linewidth(3.0)
    inset.set_xticks([])
    inset.set_yticks([])
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.87, bottom=0.01)
    return fig


def format_probability(value: Any) -> str:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if not np.isfinite(numeric):
        return "NA"
    if numeric < 0.001:
        return "<0.001"
    return f"{numeric:.3f}"


def endpoint_statistics_text(
    context: dict[str, Any], endpoint: str, region: str, language: str
) -> str:
    stats = context["statistics"].copy()
    stats = stats[
        (stats["endpoint"].astype(str) == endpoint)
        & (stats["region"].astype(str) == region)
    ]

    def selected(test: str, group1: str, group2: str | None = None) -> pd.Series | None:
        rows = stats[
            (stats["test"].astype(str) == test)
            & (stats["group1"].astype(str) == group1)
        ]
        if group2 is not None:
            rows = rows[rows["group2"].astype(str) == group2]
        return None if rows.empty else rows.iloc[0]

    anova = selected("ANOVA", "ALL")
    kw = selected("Kruskal-Wallis", "ALL")
    pair = selected("Mann-Whitney U", "Sucrose", "Allulose")
    lines = [
        f"ANOVA p={format_probability(np.nan if anova is None else anova.get('p_value'))}",
        f"KW p={format_probability(np.nan if kw is None else kw.get('p_value'))}",
    ]
    pair_valid = (
        pair is not None
        and not bool_value(pair.get("insufficient_n"))
        and np.isfinite(pd.to_numeric(pd.Series([pair.get("p_value")]), errors="coerce").iloc[0])
    )
    if pair_valid:
        q_text = format_probability(pair.get("q_value_bh"))
        delta = pd.to_numeric(pd.Series([pair.get("cliffs_delta")]), errors="coerce").iloc[0]
        delta_text = "NA" if not np.isfinite(delta) else f"{delta:.2f}"
        lines.append(f"MWU S–A: qBH={q_text}")
        lines.append(f"Cliff δ={delta_text}")
    else:
        lines.append("MWU: n insuficiente" if language == "es" else "MWU: insufficient n")
    return "\n".join(lines)


def endpoint_panel(context: dict[str, Any], letter: str, endpoint: str, language: str) -> plt.Figure:
    data = context["summary"].copy()
    data = data[data["region_present"].map(bool_value)]
    # Panel H is an ARC-specific ACTH/CLIP-cell endpoint. Keep G as the
    # two-region overview, but centre one ARC graph in H.
    regions = ("ARC",) if letter == "H" else ("ME", "ARC")
    fig, raw_axes = plt.subplots(1, len(regions), figsize=PANEL_SIZE, sharey=False)
    axes = [raw_axes] if len(regions) == 1 else list(raw_axes)
    fig.patch.set_facecolor("white")
    rng = np.random.default_rng(105 if letter == "G" else 106)
    for ax, region in zip(axes, regions):
        subset = data[data["region"].astype(str) == region]
        local_values = pd.to_numeric(subset[endpoint], errors="coerce").dropna()
        local_peak = float(local_values.max()) if not local_values.empty else 0.05
        local_peak = max(local_peak, 0.025)
        ylimit = local_peak * 2.12
        for index, condition in enumerate(CONDITIONS):
            values = pd.to_numeric(
                subset.loc[subset["condition"].astype(str) == condition, endpoint],
                errors="coerce",
            ).dropna().to_numpy(float)
            if values.size:
                mean = float(np.mean(values))
                sd = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
                ax.bar(
                    index, mean, yerr=sd, width=0.58,
                    color=CONDITION_COLORS[condition], edgecolor="#111111",
                    linewidth=1.7, capsize=4,
                    error_kw={"elinewidth": 1.6, "capthick": 1.6}, zorder=2,
                )
                jitter = rng.uniform(-0.065, 0.065, len(values))
                ax.scatter(
                    index + jitter, values, s=36,
                    facecolor=CONDITION_COLORS[condition], edgecolor="#111111",
                    linewidth=1.0, zorder=4,
                )
                bar_top = max(float(np.max(values)), mean + sd)
                ax.text(
                    index, bar_top + ylimit * 0.025, f"n={values.size}",
                    ha="center", va="bottom", fontsize=9.1, fontweight="bold",
                )
            else:
                ax.text(index, ylimit * 0.025, "n=0", ha="center", va="bottom",
                        fontsize=9.1, fontweight="bold", color="#555555")
        ax.text(
            0.045, 0.955, endpoint_statistics_text(context, endpoint, region, language),
            transform=ax.transAxes, ha="left", va="top", fontsize=9.35,
            fontweight="bold", linespacing=1.22,
            bbox={"boxstyle": "round,pad=0.42", "facecolor": "white",
                  "edgecolor": "#111111", "linewidth": 1.45,
                  "alpha": 0.96}, zorder=8,
        )
        ax.set_title(region, fontsize=11.8, fontweight="bold", pad=6)
        ax.set_xticks(range(len(CONDITIONS)))
        ax.set_xticklabels(
            [condition_label(c, language) for c in CONDITIONS],
            rotation=18, ha="right", fontsize=9.0, fontweight="bold",
        )
        ax.set_ylim(0, ylimit)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#111111")
            spine.set_linewidth(1.65)
        ax.tick_params(width=1.35, length=3.5, labelsize=8.7, color="#111111")
        ax.grid(axis="y", color="#d7dde0", linewidth=0.7, alpha=0.78)
        ax.set_axisbelow(True)
    ylabel_en = {
        "total_cfos_over_dapi": "c-FOS+ nuclei / DAPI nuclei",
        "double_over_total_acth_clip": "c-FOS+ ACTH/CLIP+ / ACTH/CLIP+",
    }
    ylabel_es = {
        "total_cfos_over_dapi": "núcleos c-FOS+ / núcleos DAPI",
        "double_over_total_acth_clip": "c-FOS+ ACTH/CLIP+ / ACTH/CLIP+",
    }
    axes[0].set_ylabel(
        (ylabel_es if language == "es" else ylabel_en)[endpoint],
        fontsize=9.6, fontweight="bold",
    )
    title = PANEL_NAMES[letter][1 if language == "es" else 0]
    fig.suptitle(title, fontsize=12.8, fontweight="bold", y=0.99)
    add_figure_letter(fig, letter, "black")
    note = (
        "Media ± DE; una adquisición por animal."
        if language == "es"
        else "Mean ± SD; one acquisition per animal."
    )
    fig.text(0.985, 0.018, note, ha="right", va="bottom",
             fontsize=7.25, color="#303844")
    if len(regions) == 1:
        fig.subplots_adjust(left=0.25, right=0.75, top=0.79, bottom=0.245)
    else:
        fig.subplots_adjust(left=0.13, right=0.985, top=0.79, bottom=0.245, wspace=0.27)
    return fig


def covariance_normalized_radial_shells(
    points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build six covariance-normalized shells containing equal DAPI counts."""
    centered = points - np.mean(points, axis=0, keepdims=True)
    covariance = np.cov(centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if not np.isfinite(eigenvalues).all() or float(np.max(eigenvalues)) <= 0:
        raise ValueError("Degenerate DAPI centroid covariance")
    eigenvalues = np.maximum(eigenvalues, float(np.max(eigenvalues)) * 1e-8)
    whitened = centered @ eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues))
    radius = np.sqrt(np.sum(whitened * whitened, axis=1))
    finite_edges = np.quantile(radius, np.linspace(0.0, 1.0, N_RADIAL_SHELLS + 1))
    if np.any(np.diff(finite_edges) <= 0):
        raise ValueError("Repeated DAPI radial quantiles prevent six-shell construction")
    assignment_edges = finite_edges.copy()
    assignment_edges[0] = -np.inf
    assignment_edges[-1] = np.inf
    shell_index = np.searchsorted(assignment_edges[1:-1], radius, side="right")
    if int(shell_index.min()) != 0 or int(shell_index.max()) != N_RADIAL_SHELLS - 1:
        raise RuntimeError("Radial shell indexing did not cover all six shells")
    return shell_index.astype(np.int16), finite_edges, eigenvalues


def spatial_pseudo_f(matrix: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    grand = np.mean(matrix, axis=0)
    total_ss = float(np.sum((matrix - grand) ** 2))
    between_ss = 0.0
    for group_index in range(len(CONDITIONS)):
        group = matrix[labels == group_index]
        if group.shape[0] == 0:
            raise ValueError("A spatial permutation contains an empty condition")
        between_ss += float(group.shape[0] * np.sum((np.mean(group, axis=0) - grand) ** 2))
    within_ss = total_ss - between_ss
    pseudo_f = (
        math.inf if within_ss <= 0
        else (between_ss / (len(CONDITIONS) - 1))
        / (within_ss / (matrix.shape[0] - len(CONDITIONS)))
    )
    r_squared = between_ss / total_ss if total_ss > 0 else math.nan
    return float(pseudo_f), float(r_squared)


def enumerate_spatial_allocations(group_sizes: tuple[int, int, int]):
    total = int(sum(group_sizes))
    indices = tuple(range(total))
    for first in itertools.combinations(indices, group_sizes[0]):
        remaining = tuple(index for index in indices if index not in first)
        for second in itertools.combinations(remaining, group_sizes[1]):
            labels = np.full(total, 2, dtype=np.int8)
            labels[list(first)] = 0
            labels[list(second)] = 1
            yield labels


def benjamini_hochberg(raw_p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(raw_p_values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = ranked * values.size / np.arange(1, values.size + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def nuclei_inside_binary_mask(
    frame: pd.DataFrame, mask: np.ndarray, acquisition: str
) -> np.ndarray:
    """Return strict centroid membership in a native-resolution binary mask."""
    if frame.empty:
        return np.zeros(0, dtype=bool)
    binary = np.asarray(mask, dtype=bool)
    if binary.ndim != 2 or not binary.any():
        raise ValueError(f"{acquisition}: invalid 3V exclusion mask")
    x_float = pd.to_numeric(
        frame["centroid_x"], errors="coerce"
    ).to_numpy(dtype=float)
    y_float = pd.to_numeric(
        frame["centroid_y"], errors="coerce"
    ).to_numpy(dtype=float)
    finite = np.isfinite(x_float) & np.isfinite(y_float)
    x_index = np.zeros(len(frame), dtype=np.int64)
    y_index = np.zeros(len(frame), dtype=np.int64)
    x_index[finite] = np.rint(x_float[finite]).astype(np.int64)
    y_index[finite] = np.rint(y_float[finite]).astype(np.int64)
    valid = (
        finite
        & (x_index >= 0)
        & (x_index < binary.shape[1])
        & (y_index >= 0)
        & (y_index < binary.shape[0])
    )
    if not valid.all():
        raise ValueError(
            f"{acquisition}: one or more nucleus centroids fall outside "
            "the reviewed 3V mask dimensions"
        )
    return binary[y_index, x_index]


def rasterize_normalized_polygon(
    points: np.ndarray, shape: tuple[int, int]
) -> np.ndarray:
    """Rasterize a normalized closed polygon for zero-density display masking."""
    height, width = map(int, shape)
    polygon = np.asarray(points, dtype=float)
    if (
        height <= 0
        or width <= 0
        or polygon.ndim != 2
        or polygon.shape[0] < 3
        or polygon.shape[1] != 2
        or not np.isfinite(polygon).all()
    ):
        raise ValueError("Invalid normalized 3V polygon")
    polygon = np.clip(polygon, 0.0, 1.0)
    canvas = Image.new("1", (width, height), 0)
    ImageDraw.Draw(canvas).polygon(
        [
            (
                float(x) * (width - 1),
                float(y) * (height - 1),
            )
            for x, y in polygon
        ],
        fill=1,
    )
    return np.asarray(canvas, dtype=bool)


def spatial_radial_analysis(
    context: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute 3V-excluded animal profiles and exhaustive label tests."""
    nuclei = context["nuclei"].copy()
    nuclei = nuclei[nuclei["region"].astype(str).isin({"ARC", "ME"})]
    ventricle_masks = context.get("spatial_ventricle_masks", {})
    expected = set(context["sheet"]["acquisition_id"].astype(str))
    if set(ventricle_masks) != expected:
        raise ValueError(
            "Per-acquisition accepted 3V masks are required for I/J analysis"
        )
    feature_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    endpoint_masks = {
        "cfos": lambda frame: frame["is_cfos"].map(bool_value).to_numpy(bool),
        "double": lambda frame: (
            frame["is_cfos"].map(bool_value).to_numpy(bool)
            & frame["is_acth_clip"].map(bool_value).to_numpy(bool)
        ),
    }
    endpoint_labels = {
        "cfos": "c-FOS-positive nuclei / DAPI nuclei",
        "double": "c-FOS and ACTH/CLIP double-positive nuclei / DAPI nuclei",
    }
    for _, sample in context["sheet"].sort_values("acquisition_id").iterrows():
        acquisition = str(sample["acquisition_id"])
        before = nuclei[
            nuclei["acquisition_id"].astype(str) == acquisition
        ].copy()
        before = before.sort_values("nucleus_label", kind="mergesort")
        inside_3v = nuclei_inside_binary_mask(
            before, ventricle_masks[acquisition], acquisition
        )
        before_cfos = endpoint_masks["cfos"](before)
        before_double = endpoint_masks["double"](before)
        section = before.loc[~inside_3v].copy()
        if len(section) < N_RADIAL_SHELLS:
            raise RuntimeError(f"Fewer nuclei than radial shells: {acquisition}")
        points = section[["centroid_y", "centroid_x"]].to_numpy(dtype=np.float64)
        shell_index, edges, eigenvalues = covariance_normalized_radial_shells(points)
        shell_counts = [int(np.count_nonzero(shell_index == index))
                        for index in range(N_RADIAL_SHELLS)]
        if max(shell_counts) - min(shell_counts) > 1:
            raise RuntimeError(f"Unequal-DAPI shells for {acquisition}: {shell_counts}")
        audit_rows.append({
            "acquisition_id": acquisition,
            "animal_id": str(sample["animal_id"]),
            "condition": str(sample["condition"]),
            "arc_me_dapi_nuclei_before_3v_exclusion": int(len(before)),
            "dapi_nuclei_excluded_inside_3v": int(inside_3v.sum()),
            "cfos_nuclei_excluded_inside_3v": int(
                np.count_nonzero(inside_3v & before_cfos)
            ),
            "double_nuclei_excluded_inside_3v": int(
                np.count_nonzero(inside_3v & before_double)
            ),
            "arc_me_dapi_nuclei_after_3v_exclusion": int(len(section)),
            "arc_me_dapi_nuclei": int(len(section)),
            "ventricle_exclusion_applied": True,
            "ventricle_exclusion_source": (
                "accepted per-acquisition DAPI-only HIL 3V mask"
            ),
            "shell_min_dapi_nuclei": min(shell_counts),
            "shell_max_dapi_nuclei": max(shell_counts),
            "whitening_eigenvalue_small": float(np.min(eigenvalues)),
            "whitening_eigenvalue_large": float(np.max(eigenvalues)),
            "innermost_whitened_edge": float(edges[0]),
            "outermost_whitened_edge": float(edges[-1]),
        })
        for endpoint, mask_builder in endpoint_masks.items():
            positive_mask = mask_builder(section)
            for shell_zero in range(N_RADIAL_SHELLS):
                selected = shell_index == shell_zero
                denominator = int(np.count_nonzero(selected))
                positive = int(np.count_nonzero(positive_mask & selected))
                occurrence = positive / denominator
                feature_rows.append({
                    "endpoint": endpoint,
                    "endpoint_label": endpoint_labels[endpoint],
                    "acquisition_id": acquisition,
                    "animal_id": str(sample["animal_id"]),
                    "condition": str(sample["condition"]),
                    "radial_shell": shell_zero + 1,
                    "shell_dapi_nuclei": denominator,
                    "shell_positive_nuclei": positive,
                    "occurrence": occurrence,
                    "occurrence_percent": 100.0 * occurrence,
                    "arcsin_sqrt_occurrence": float(np.arcsin(np.sqrt(occurrence))),
                    "experimental_unit": "biological animal (one acquisition per animal)",
                    "ventricle_exclusion_applied": True,
                    "ventricle_exclusion_source": (
                        "accepted per-acquisition DAPI-only HIL 3V mask"
                    ),
                })
    features = pd.DataFrame(feature_rows).sort_values(
        ["endpoint", "condition", "animal_id", "radial_shell"]
    ).reset_index(drop=True)
    condition_to_index = {condition: index for index, condition in enumerate(CONDITIONS)}
    statistic_rows: list[dict[str, Any]] = []
    permutation_rows: list[dict[str, Any]] = []
    for endpoint in ("cfos", "double"):
        selected = features[features["endpoint"].eq(endpoint)]
        animal_ids = sorted(selected["animal_id"].astype(str).unique())
        animal_condition = (
            selected[["animal_id", "condition"]].drop_duplicates()
            .set_index("animal_id")["condition"].astype(str).to_dict()
        )
        observed = np.asarray(
            [condition_to_index[animal_condition[animal]] for animal in animal_ids],
            dtype=np.int8,
        )
        group_sizes = tuple(int(np.count_nonzero(observed == index)) for index in range(3))
        if min(group_sizes) < 2:
            raise RuntimeError(f"Spatial endpoint {endpoint} lacks two animals per condition")
        pivot = selected.pivot(
            index="animal_id", columns="radial_shell", values="arcsin_sqrt_occurrence"
        ).reindex(animal_ids)
        if list(pivot.columns) != list(range(1, N_RADIAL_SHELLS + 1)) or pivot.isna().any().any():
            raise RuntimeError(f"Incomplete six-shell feature matrix for {endpoint}")
        matrix = pivot.to_numpy(dtype=np.float64)
        observed_f, observed_r2 = spatial_pseudo_f(matrix, observed)
        expected = math.comb(len(animal_ids), group_sizes[0]) * math.comb(
            len(animal_ids) - group_sizes[0], group_sizes[1]
        )
        extreme = 0
        enumerated = 0
        observed_seen = 0
        for allocation_index, labels in enumerate(
            enumerate_spatial_allocations(group_sizes), start=1
        ):
            pseudo_f, r_squared = spatial_pseudo_f(matrix, labels)
            is_observed = bool(np.array_equal(labels, observed))
            is_extreme = bool(pseudo_f >= observed_f - 1e-12)
            observed_seen += int(is_observed)
            extreme += int(is_extreme)
            enumerated += 1
            permutation_rows.append({
                "endpoint": endpoint,
                "allocation_index": allocation_index,
                "water_animals": ";".join(a for a, label in zip(animal_ids, labels) if label == 0),
                "sucrose_animals": ";".join(a for a, label in zip(animal_ids, labels) if label == 1),
                "allulose_animals": ";".join(a for a, label in zip(animal_ids, labels) if label == 2),
                "pseudo_f": pseudo_f,
                "r_squared": r_squared,
                "is_observed_allocation": is_observed,
                "is_as_or_more_extreme": is_extreme,
            })
        if enumerated != expected or observed_seen != 1:
            raise RuntimeError(
                f"Spatial enumeration invariant failed for {endpoint}: "
                f"{enumerated}/{expected}, observed={observed_seen}"
            )
        statistic_rows.append({
            "endpoint": endpoint,
            "endpoint_label": selected["endpoint_label"].iloc[0],
            "test": "exact animal-label PERMANOVA",
            "comparison": "Water|Sucrose|Allulose",
            "n_animals": len(animal_ids),
            "n_water": group_sizes[0],
            "n_sucrose": group_sizes[1],
            "n_allulose": group_sizes[2],
            "feature_count": N_RADIAL_SHELLS,
            "feature_transform": "arcsin(sqrt(shell positive nuclei / shell DAPI nuclei))",
            "distance": "Euclidean on six-shell animal vectors",
            "pseudo_f": observed_f,
            "r_squared": observed_r2,
            "p_value_exact": float(extreme / enumerated),
            "extreme_labelings": extreme,
            "enumerated_labelings": enumerated,
            "permutation_exhaustive": True,
            "experimental_unit": "biological animal",
            "cells_as_independent_replicates": False,
            "analysis_tier": "exploratory",
            "ventricle_exclusion": (
                "accepted per-acquisition DAPI-only HIL 3V mask applied "
                "before shell construction and endpoint calculation"
            ),
        })
    statistics = pd.DataFrame(statistic_rows)
    statistics["bh_q_value_two_endpoint_family"] = benjamini_hochberg(
        statistics["p_value_exact"].to_numpy(dtype=float)
    )
    statistics["multiple_testing_family"] = (
        "two declared S5 ARC/ME spatial profiles: c-FOS/DAPI and "
        "c-FOS+ACTH/CLIP/DAPI"
    )
    audit = pd.DataFrame(audit_rows).sort_values("acquisition_id").reset_index(drop=True)
    permutations = pd.DataFrame(permutation_rows)
    return features, statistics, audit, permutations


def spatial_grids(
    context: dict[str, Any], phenotype: str, bins_y: int = 34, bins_x: int = 50
) -> tuple[
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    list[dict[str, Any]],
]:
    sheet = context["sheet"].copy()
    sheet.index = sheet["acquisition_id"].astype(str)
    nuclei = context["nuclei"].copy()
    if phenotype == "cfos":
        nuclei = nuclei[nuclei["is_cfos"].map(bool_value)]
    elif phenotype == "double":
        nuclei = nuclei[
            nuclei["is_cfos"].map(bool_value) & nuclei["is_acth_clip"].map(bool_value)
        ]
    else:
        raise ValueError(phenotype)
    nuclei = nuclei[nuclei["region"].astype(str).isin({"ARC", "ME"})]
    ventricle_masks = context.get("spatial_ventricle_masks", {})
    expected = set(sheet.index.astype(str))
    if set(ventricle_masks) != expected:
        raise ValueError(
            "Per-acquisition accepted 3V masks are required for I/J maps"
        )
    by_condition: dict[str, list[np.ndarray]] = {
        condition: [] for condition in CONDITIONS
    }
    by_condition_tissue: dict[str, list[np.ndarray]] = {
        condition: [] for condition in CONDITIONS
    }
    long_rows: list[dict[str, Any]] = []
    for acquisition, row in sheet.iterrows():
        region_path = (
            context["workdir"]
            / "region_masks"
            / f"{acquisition}_regions.tif"
        )
        tissue_mask = np.squeeze(
            np.asarray(tifffile.imread(str(region_path)), dtype=np.uint8)
        ) > 0
        ventricle_mask = np.asarray(
            ventricle_masks[str(acquisition)], dtype=bool
        )
        if ventricle_mask.shape != tissue_mask.shape:
            raise ValueError(
                f"{acquisition}: 3V/anatomy mask shape mismatch"
            )
        yy, xx = np.nonzero(tissue_mask)
        if not yy.size:
            raise ValueError(f"{acquisition}: empty accepted anatomy HIL")
        x_span = max(1.0, float(xx.max() - xx.min()))
        y_span = max(1.0, float(yy.max() - yy.min()))
        subset = nuclei[
            nuclei["acquisition_id"].astype(str) == str(acquisition)
        ].copy()
        inside_3v = nuclei_inside_binary_mask(
            subset, ventricle_mask, str(acquisition)
        )
        subset = subset.loc[~inside_3v].copy()
        x_norm = (
            subset["centroid_x"].to_numpy(float) - float(xx.min())
        ) / x_span
        y_norm = (
            subset["centroid_y"].to_numpy(float) - float(yy.min())
        ) / y_span
        valid = (x_norm >= 0) & (x_norm <= 1) & (y_norm >= 0) & (y_norm <= 1)
        hist, _, _ = np.histogram2d(
            y_norm[valid], x_norm[valid], bins=(bins_y, bins_x),
            range=((0, 1), (0, 1)),
        )
        tissue_hist, _, _ = np.histogram2d(
            (yy.astype(float) - float(yy.min())) / y_span,
            (xx.astype(float) - float(xx.min())) / x_span,
            bins=(bins_y, bins_x), range=((0, 1), (0, 1)),
        )
        allowed_mask = tissue_mask & ~ventricle_mask
        allowed_y, allowed_x = np.nonzero(allowed_mask)
        allowed_hist, _, _ = np.histogram2d(
            (allowed_y.astype(float) - float(yy.min())) / y_span,
            (allowed_x.astype(float) - float(xx.min())) / x_span,
            bins=(bins_y, bins_x), range=((0, 1), (0, 1)),
        )
        ventricle_y, ventricle_x = np.nonzero(ventricle_mask)
        ventricle_hist, _, _ = np.histogram2d(
            (ventricle_y.astype(float) - float(yy.min())) / y_span,
            (ventricle_x.astype(float) - float(xx.min())) / x_span,
            bins=(bins_y, bins_x), range=((0, 1), (0, 1)),
        )
        tissue = tissue_hist > 0
        ventricle = ventricle_hist > 0
        analysis_domain = (allowed_hist > 0) & ~ventricle
        occurrence = ((hist > 0) & analysis_domain).astype(np.float32)
        smoothed = gaussian_filter(occurrence, sigma=1.25, mode="nearest")
        smoothed[~analysis_domain] = 0.0
        condition = str(row["condition"])
        by_condition[condition].append(smoothed)
        by_condition_tissue[condition].append(tissue.astype(np.float32))
        boundary_source = (
            "outer envelope of accepted DAPI-only ARC/ME/VMN anatomy HIL"
        )
        for iy in range(bins_y):
            for ix in range(bins_x):
                long_rows.append({
                    "phenotype": phenotype,
                    "acquisition_id": acquisition,
                    "animal_id": row["animal_id"],
                    "condition": condition,
                    "y_bin": iy,
                    "x_bin": ix,
                    "binary_occurrence": float(occurrence[iy, ix]),
                    "smoothed_occurrence": float(smoothed[iy, ix]),
                    "tissue_occupancy": float(tissue[iy, ix]),
                    "analysis_domain_occupancy": float(
                        analysis_domain[iy, ix]
                    ),
                    "ventricle_occupancy": float(ventricle[iy, ix]),
                    "ventricle_exclusion_applied": True,
                    "ventricle_exclusion_source": (
                        "accepted per-acquisition DAPI-only HIL 3V mask"
                    ),
                    "tissue_boundary_source": boundary_source,
                })
    mean_grids = {
        condition: np.mean(by_condition[condition], axis=0)
        if by_condition[condition]
        else np.zeros((bins_y, bins_x), dtype=np.float32)
        for condition in CONDITIONS
    }
    mean_tissue = {
        condition: np.mean(by_condition_tissue[condition], axis=0)
        if by_condition_tissue[condition]
        else np.zeros((bins_y, bins_x), dtype=np.float32)
        for condition in CONDITIONS
    }
    display_ventricle = {
        condition: rasterize_normalized_polygon(
            context["condition_ventricle_contours"][condition]["points"],
            (bins_y, bins_x),
        )
        for condition in CONDITIONS
    }
    for condition in CONDITIONS:
        mean_grids[condition][display_ventricle[condition]] = 0.0
    return mean_grids, mean_tissue, display_ventricle, long_rows


def outer_tissue_envelope(tissue_occupancy: np.ndarray) -> np.ndarray:
    """Fill holes without changing the HIL-reviewed outer tissue contour."""
    support = np.asarray(tissue_occupancy >= 0.5, dtype=bool)
    if not support.any():
        support = np.asarray(tissue_occupancy > 0, dtype=bool)
    return np.asarray(binary_fill_holes(support), dtype=bool) if support.any() else support


def select_condition_ventricle_contours(
    context: dict[str, Any], raster_size: int = 128
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Normalize reviewed 3V points and choose one clean medoid per condition.

    Tissue coordinates always come from the accepted DAPI-only ARC/ME/VMN
    anatomy mask. The 3V points are normalized into that same frame for the
    visible outline and condition-level zero-density display mask. Native
    per-acquisition 3V raster masks exclude nuclei from spatial inference.
    A single observed contour is selected per condition instead of averaging
    outlines, which avoids synthetic branches when fields differ in orientation.
    """
    reviewed = context.get("spatial_ventricle_contours", {})
    expected = set(context["sheet"]["acquisition_id"].astype(str))
    if set(reviewed) != expected:
        raise ValueError(
            "A finalized 3V HIL contour is required for every S5 acquisition"
        )
    sheet = context["sheet"].copy()
    sheet.index = sheet["acquisition_id"].astype(str)
    candidates: dict[str, list[dict[str, Any]]] = {
        condition: [] for condition in CONDITIONS
    }
    for acquisition in sorted(expected):
        region_path = (
            context["workdir"]
            / "region_masks"
            / f"{acquisition}_regions.tif"
        )
        region_mask = np.squeeze(
            np.asarray(tifffile.imread(str(region_path)), dtype=np.uint8)
        ) > 0
        yy, xx = np.nonzero(region_mask)
        if not yy.size:
            raise ValueError(
                f"{acquisition}: accepted anatomy HIL has no tissue pixels"
            )
        x_span = max(1.0, float(xx.max() - xx.min()))
        y_span = max(1.0, float(yy.max() - yy.min()))
        points_raw = np.asarray(reviewed[acquisition], dtype=float)
        normalized_raw = np.column_stack(
            (
                (points_raw[:, 0] - float(xx.min())) / x_span,
                (points_raw[:, 1] - float(yy.min())) / y_span,
            )
        )
        normalized = np.clip(normalized_raw, 0.0, 1.0)
        retained = np.ones(len(normalized), dtype=bool)
        if len(normalized) > 1:
            retained[1:] = (
                np.linalg.norm(np.diff(normalized, axis=0), axis=1) > 1e-9
            )
        normalized = normalized[retained]
        normalized_raw = normalized_raw[retained]
        points_raw = points_raw[retained]
        if len(normalized) < 3:
            raise ValueError(
                f"{acquisition}: normalized 3V contour has fewer than 3 points"
            )
        canvas = Image.new("1", (raster_size, raster_size), 0)
        ImageDraw.Draw(canvas).polygon(
            [
                (
                    float(x) * (raster_size - 1),
                    float(y) * (raster_size - 1),
                )
                for x, y in normalized
            ],
            fill=1,
        )
        condition = str(sheet.loc[acquisition, "condition"])
        candidates[condition].append(
            {
                "acquisition_id": acquisition,
                "condition": condition,
                "points": normalized,
                "points_raw": points_raw,
                "points_normalized_unclipped": normalized_raw,
                "raster": np.asarray(canvas, dtype=bool),
            }
        )
    selected: dict[str, dict[str, Any]] = {}
    source_rows: list[dict[str, Any]] = []
    for condition in CONDITIONS:
        group = candidates[condition]
        if not group:
            raise ValueError(f"No reviewed 3V contour for {condition}")
        for candidate in group:
            distances: list[float] = []
            for other in group:
                union = np.count_nonzero(candidate["raster"] | other["raster"])
                intersection = np.count_nonzero(
                    candidate["raster"] & other["raster"]
                )
                distances.append(
                    0.0 if union == 0 else 1.0 - intersection / union
                )
            candidate["medoid_distance_sum"] = float(sum(distances))
        winner = min(
            group,
            key=lambda item: (
                item["medoid_distance_sum"],
                item["acquisition_id"],
            ),
        )
        selected[condition] = {
            "acquisition_id": winner["acquisition_id"],
            "points": winner["points"],
            "selection_method": (
                "minimum summed pairwise Jaccard distance among observed "
                "condition contours; acquisition ID breaks exact ties"
            ),
        }
        for candidate in group:
            for point_index, (raw, unclipped, clipped) in enumerate(
                zip(
                    candidate["points_raw"],
                    candidate["points_normalized_unclipped"],
                    candidate["points"],
                ),
                start=1,
            ):
                source_rows.append(
                    {
                        "condition": condition,
                        "acquisition_id": candidate["acquisition_id"],
                        "selected_as_condition_medoid": (
                            candidate["acquisition_id"]
                            == winner["acquisition_id"]
                        ),
                        "medoid_distance_sum": candidate[
                            "medoid_distance_sum"
                        ],
                        "point_index": point_index,
                        "x_raw_px": float(raw[0]),
                        "y_raw_px": float(raw[1]),
                        "x_normalized_unclipped": float(unclipped[0]),
                        "y_normalized_unclipped": float(unclipped[1]),
                        "x_normalized_display": float(clipped[0]),
                        "y_normalized_display": float(clipped[1]),
                        "normalization_geometry": (
                            "accepted DAPI-only ARC/ME/VMN anatomy HIL bbox"
                        ),
                        "geometry_role": (
                            "3V analysis exclusion plus display overlay; "
                            "never outer-tissue crop or normalization"
                        ),
                    }
                )
    return selected, source_rows


def spatial_endpoint_labels(phenotype: str, language: str) -> tuple[str, str]:
    if phenotype == "cfos":
        if language == "es":
            return "Ocurrencia radial de c-FOS", "núcleos c-FOS⁺ / núcleos DAPI (%)"
        return "Radial c-FOS occurrence", "c-FOS⁺ nuclei / DAPI nuclei (%)"
    if language == "es":
        return (
            "Ocurrencia radial doble positiva",
            "c-FOS⁺ ACTH/CLIP⁺ / núcleos DAPI (%)",
        )
    return (
        "Radial double-positive occurrence",
        "c-FOS⁺ ACTH/CLIP⁺ / DAPI nuclei (%)",
    )


def heatmap_panel(
    context: dict[str, Any], letter: str, phenotype: str, language: str
) -> tuple[plt.Figure, list[dict[str, Any]]]:
    grids, tissue_grids, ventricle_grids, rows = spatial_grids(context, phenotype)
    maximum = max(float(np.max(grid)) for grid in grids.values())
    maximum = maximum if maximum > 0 else 1.0
    fig = plt.figure(figsize=SPATIAL_PANEL_SIZE, facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=(1.25, 1.0), hspace=0.56)
    top = outer[0].subgridspec(1, 4, width_ratios=(1.0, 1.0, 1.0, 0.055), wspace=0.12)
    map_axes = [fig.add_subplot(top[0, index]) for index in range(3)]
    map_color_axis = fig.add_subplot(top[0, 3])
    cmap = plt.get_cmap(
        "magma" if phenotype == "cfos" else "viridis"
    ).copy()
    cmap.set_bad(cmap(0.0))
    shown = None
    for ax, condition in zip(map_axes, CONDITIONS):
        display_grid = np.ma.array(
            grids[condition], mask=ventricle_grids[condition]
        )
        shown = ax.imshow(
            display_grid, cmap=cmap, vmin=0, vmax=maximum,
            origin="upper", aspect="auto", interpolation="nearest",
        )
        envelope = outer_tissue_envelope(tissue_grids[condition])
        if envelope.any() and not envelope.all():
            ax.contour(
                envelope.astype(float), levels=[0.5], colors=["#c8c8c8"],
                linewidths=1.15, linestyles=[(0, (2, 2))], zorder=7,
            )
        ventricle = context["condition_ventricle_contours"][condition]
        ventricle_points = np.asarray(ventricle["points"], dtype=float)
        ventricle_points = np.vstack(
            (ventricle_points, ventricle_points[0])
        )
        ax.plot(
            ventricle_points[:, 0] * (grids[condition].shape[1] - 1),
            ventricle_points[:, 1] * (grids[condition].shape[0] - 1),
            color="white", linewidth=2.45, linestyle=(0, (5, 2.4)),
            solid_capstyle="round", dash_capstyle="round", zorder=9,
        )
        n_units = int((context["sheet"]["condition"].astype(str) == condition).sum())
        ax.set_title(
            f"{condition_label(condition, language)}\n={n_units}",
            fontsize=9.4, fontweight="bold", pad=5,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
            spine.set_color("#222222")
    map_label = (
        "mean acquisition occurrence"
        if language == "en"
        else "ocurrencia media por adquisición"
    )
    map_colorbar = fig.colorbar(shown, cax=map_color_axis)
    map_colorbar.set_label(map_label, fontsize=7.7, labelpad=5)
    map_colorbar.ax.tick_params(labelsize=7.2, width=0.8)
    ventricle_label = (
        "3V excluded"
        if language == "en"
        else "3V excluido"
    )
    tissue_label = (
        "tissue edge" if language == "en" else "borde tisular"
    )
    map_axes[0].legend(
        handles=[
            Line2D(
                [0], [0], color="white", linewidth=2.45,
                linestyle=(0, (5, 2.4)), label=ventricle_label,
            ),
            Line2D(
                [0], [0], color="#c8c8c8", linewidth=1.15,
                linestyle=(0, (2, 2)), label=tissue_label,
            ),
        ],
        loc="lower left", frameon=True, facecolor="black", edgecolor="white",
        framealpha=0.72, labelcolor="white", fontsize=5.8, handlelength=1.8,
        ncol=2, columnspacing=0.75, handletextpad=0.45, borderpad=0.35,
    )

    bottom = outer[1].subgridspec(
        1, 3, width_ratios=(2.28, 1.42, 0.055), wspace=0.34
    )
    profile_axis = fig.add_subplot(bottom[0, 0])
    means_axis = fig.add_subplot(bottom[0, 1])
    means_color_axis = fig.add_subplot(bottom[0, 2])
    endpoint = "cfos" if phenotype == "cfos" else "double"
    features = context["spatial_radial_features"]
    selected = features[features["endpoint"].astype(str) == endpoint]
    statistic = context["spatial_radial_statistics"]
    statistic = statistic[statistic["endpoint"].astype(str) == endpoint].iloc[0]
    x_values = np.arange(1, N_RADIAL_SHELLS + 1, dtype=float)
    for condition in CONDITIONS:
        condition_rows = selected[selected["condition"].astype(str) == condition]
        profiles: list[np.ndarray] = []
        for animal_id in condition_rows["animal_id"].drop_duplicates().astype(str):
            profile = (
                condition_rows[condition_rows["animal_id"].astype(str) == animal_id]
                .sort_values("radial_shell")["occurrence_percent"].to_numpy(float)
            )
            if profile.size != N_RADIAL_SHELLS:
                raise RuntimeError(f"Incomplete S5 spatial profile: {animal_id}/{endpoint}")
            profiles.append(profile)
            profile_axis.plot(
                x_values, profile, color=MEAN_LINE_COLORS[condition], alpha=0.23,
                linewidth=0.95, marker="o", markersize=2.8,
                markeredgewidth=0.0, zorder=2,
            )
        mean_profile = np.mean(np.stack(profiles), axis=0)
        profile_axis.plot(
            x_values, mean_profile, color=MEAN_LINE_COLORS[condition],
            linewidth=2.45, marker="o", markersize=5.0,
            markerfacecolor=CONDITION_COLORS[condition],
            markeredgecolor="#172126", markeredgewidth=0.65, zorder=4,
        )
    ymax = float(selected["occurrence_percent"].max())
    profile_axis.set_ylim(0.0, max(1.0, ymax * 1.24))
    profile_axis.set_xlim(0.72, N_RADIAL_SHELLS + 0.28)
    profile_axis.set_xticks(x_values)
    inner = "interior" if language == "es" else "inner"
    outer_label = "exterior" if language == "es" else "outer"
    profile_axis.set_xticklabels([f"1\n{inner}", "2", "3", "4", "5", f"6\n{outer_label}"], fontsize=7.4)
    profile_title, profile_ylabel = spatial_endpoint_labels(phenotype, language)
    profile_axis.set_title(profile_title, fontsize=9.5, fontweight="bold", loc="left", pad=7)
    profile_axis.set_xlabel(
        "Capa radial de igual densidad DAPI" if language == "es"
        else "Equal-DAPI-density radial shell",
        fontsize=8.0, labelpad=2,
    )
    profile_axis.set_ylabel(profile_ylabel, fontsize=8.0, labelpad=3)
    p_value = float(statistic["p_value_exact"])
    q_value = float(statistic["bh_q_value_two_endpoint_family"])
    stats_label = (
        f"PERMANOVA exacta p={format_probability(p_value)}\nBH q={format_probability(q_value)}"
        if language == "es"
        else f"Exact PERMANOVA p={format_probability(p_value)}\nBH q={format_probability(q_value)}"
    )
    profile_axis.text(
        0.97, 0.95, stats_label, transform=profile_axis.transAxes,
        ha="right", va="top", fontsize=8.55, fontweight="bold",
        bbox={"boxstyle": "round,pad=0.36", "facecolor": "white",
              "edgecolor": "#26353b", "linewidth": 1.05, "alpha": 0.94},
        zorder=8,
    )
    profile_axis.grid(axis="y", color="#d9dee1", linewidth=0.65, alpha=0.9)
    profile_axis.set_axisbelow(True)
    profile_axis.spines[["top", "right"]].set_visible(False)
    profile_axis.tick_params(width=0.85, length=3.0, labelsize=7.4, color="#30393d")
    legend_handles = [
        Line2D([0], [0], color=MEAN_LINE_COLORS[condition], marker="o",
               markerfacecolor=CONDITION_COLORS[condition], markeredgecolor="#172126",
               linewidth=2.0, markersize=4.5, label=condition_label(condition, language))
        for condition in CONDITIONS
    ]
    profile_axis.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.31),
        ncol=3, frameon=False, fontsize=7.2, columnspacing=1.2, handlelength=1.5,
    )

    condition_means = np.vstack([
        selected[selected["condition"].astype(str) == condition]
        .groupby("radial_shell", sort=True)["occurrence_percent"].mean()
        .reindex(range(1, N_RADIAL_SHELLS + 1)).to_numpy(float)
        for condition in CONDITIONS
    ])
    heat_image = means_axis.imshow(
        condition_means, cmap="inferno", vmin=0.0,
        vmax=max(float(np.max(condition_means)), 1e-12),
        aspect="auto", interpolation="nearest",
    )
    means_axis.set_title(
        "Medias por condición" if language == "es" else "Condition means",
        fontsize=9.2, fontweight="bold", pad=7,
    )
    means_axis.set_xticks(np.arange(N_RADIAL_SHELLS))
    means_axis.set_xticklabels(range(1, N_RADIAL_SHELLS + 1), fontsize=7.2)
    means_axis.set_xlabel(
        "capa radial" if language == "es" else "radial shell", fontsize=7.6, labelpad=2
    )
    means_axis.set_yticks(np.arange(len(CONDITIONS)))
    short_labels = ("A", "S", "Al") if language == "es" else ("W", "S", "A")
    means_axis.set_yticklabels(short_labels, fontsize=7.7, fontweight="bold")
    for tick, condition in zip(means_axis.get_yticklabels(), CONDITIONS):
        tick.set_color(MEAN_LINE_COLORS[condition])
    for row_index in range(len(CONDITIONS)):
        for column_index in range(N_RADIAL_SHELLS):
            value = float(condition_means[row_index, column_index])
            red, green, blue, _ = heat_image.cmap(heat_image.norm(value))
            luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
            means_axis.text(
                column_index, row_index, f"{value:.1f}", ha="center", va="center",
                fontsize=6.55, fontweight="bold",
                color="#172126" if luminance >= 0.56 else "white",
            )
    means_axis.set_xticks(np.arange(-0.5, N_RADIAL_SHELLS, 1), minor=True)
    means_axis.set_yticks(np.arange(-0.5, len(CONDITIONS), 1), minor=True)
    means_axis.grid(which="minor", color="white", linewidth=0.55, alpha=0.72)
    means_axis.tick_params(which="minor", bottom=False, left=False)
    means_axis.tick_params(which="major", length=0)
    for spine in means_axis.spines.values():
        spine.set_color("#45565d")
        spine.set_linewidth(0.85)
    means_colorbar = fig.colorbar(heat_image, cax=means_color_axis)
    means_colorbar.set_label("%", fontsize=8.0)
    means_colorbar.ax.tick_params(labelsize=7.0, width=0.8)

    title = PANEL_NAMES[letter][1 if language == "es" else 0]
    fig.suptitle(title, fontsize=13.2, fontweight="bold", y=0.99)
    add_figure_letter(fig, letter, "black")
    coordinate_label = (
        "normalized anatomy-HIL tissue coordinates · 3V excluded"
        if language == "en"
        else "coordenadas tisulares HIL anatómicas normalizadas · 3V excluido"
    )
    fig.text(0.445, 0.505, coordinate_label, ha="center", fontsize=7.5)
    fig.subplots_adjust(left=0.075, right=0.965, top=0.87, bottom=0.125)
    return fig, rows


def render_panels(context: dict[str, Any], output_dir: Path, language: str) -> tuple[dict[str, Path], list[dict[str, Any]]]:
    suffix = "_spanish" if language == "es" else ""
    panel_dir = output_dir / "panels"
    panel_dir.mkdir(parents=True, exist_ok=True)
    builders: dict[str, Callable[[], plt.Figure]] = {
        "A": lambda: image_panel(context, "A", "DAPI", language),
        "B": lambda: image_panel(context, "B", "cFOS", language),
        "C": lambda: image_panel(context, "C", "ACTH_CLIP", language),
        "D": lambda: cartoon_panel(context, language),
        "E": lambda: merged_panel(context, language, False),
        "F": lambda: merged_panel(context, language, True),
        "G": lambda: endpoint_panel(context, "G", "total_cfos_over_dapi", language),
        "H": lambda: endpoint_panel(context, "H", "double_over_total_acth_clip", language),
    }
    outputs: dict[str, Path] = {}
    spatial_rows: list[dict[str, Any]] = []
    for letter in "ABCDEFGH":
        stem = f"Figure_S5_panel_{letter}{suffix}"
        png = panel_dir / f"{stem}.png"
        pdf = panel_dir / f"{stem}.pdf"
        save_figure(builders[letter](), png, pdf)
        outputs[letter] = png
    for letter, phenotype in (("I", "cfos"), ("J", "double")):
        fig, rows = heatmap_panel(context, letter, phenotype, language)
        stem = f"Figure_S5_panel_{letter}{suffix}"
        png = panel_dir / f"{stem}.png"
        pdf = panel_dir / f"{stem}.pdf"
        save_figure(fig, png, pdf)
        outputs[letter] = png
        if language == "en":
            spatial_rows.extend(rows)
    return outputs, spatial_rows


def compose_master(panel_paths: dict[str, Path], output_dir: Path, language: str) -> tuple[Path, Path]:
    """Assemble 2×3 microscopy, vertically stacked G/H, then wide I/J."""
    if language != "en":
        raise ValueError(
            "S5 multipanel masters are English-only; Spanish is isolated panels only"
        )
    suffix = "_es" if language == "es" else ""
    fig = plt.figure(figsize=(16.5, 15.35))
    fig.patch.set_facecolor("white")
    outer = fig.add_gridspec(2, 1, height_ratios=(1.37, 1.0), hspace=0.018)
    top = outer[0].subgridspec(
        6, 3, width_ratios=(1.0, 1.0, 1.45), hspace=0.025, wspace=0.025
    )
    bottom = outer[1].subgridspec(1, 2, wspace=0.025)
    positions = {
        "A": top[0:2, 0], "B": top[0:2, 1],
        "C": top[2:4, 0], "D": top[2:4, 1],
        "E": top[4:6, 0], "F": top[4:6, 1],
        "G": top[0:3, 2], "H": top[3:6, 2],
        "I": bottom[0, 0], "J": bottom[0, 1],
    }
    for letter in ("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"):
        ax = fig.add_subplot(positions[letter])
        image = np.asarray(Image.open(panel_paths[letter]).convert("RGB"))
        ax.imshow(image)
        ax.axis("off")
    if language == "es":
        title = "Ayuno de 16 h y primera exposición · c-FOS y ACTH/CLIP"
    else:
        title = "16 h fasting and first-ever exposure · c-FOS and ACTH/CLIP"
    fig.suptitle(title, fontsize=19, fontweight="bold", y=0.996)
    fig.subplots_adjust(
        left=0.007, right=0.993, top=0.969, bottom=0.006,
    )
    png = output_dir / f"Figure_S5{suffix}.png"
    pdf = output_dir / f"Figure_S5{suffix}.pdf"
    save_figure(fig, png, pdf)
    return png, pdf


def write_legends(context: dict[str, Any], output_dir: Path) -> None:
    size_audit = context["size_audit"]
    rejected = int(pd.to_numeric(size_audit["n_rejected_likely_noncell_by_size"], errors="coerce").sum())
    raw = int(pd.to_numeric(size_audit["n_acth_clip_rois_raw"], errors="coerce").sum())
    spatial = context["spatial_radial_statistics"].set_index("endpoint")
    cfos_p = format_probability(spatial.loc["cfos", "p_value_exact"])
    cfos_q = format_probability(spatial.loc["cfos", "bh_q_value_two_endpoint_family"])
    double_p = format_probability(spatial.loc["double", "p_value_exact"])
    double_q = format_probability(spatial.loc["double", "bh_q_value_two_endpoint_family"])
    allocations = int(spatial.loc["cfos", "enumerated_labelings"])
    exclusion = context["spatial_3v_exclusion_audit"]
    dapi_excluded = int(
        exclusion["dapi_nuclei_excluded_inside_3v"].sum()
    )
    cfos_excluded = int(
        exclusion["cfos_nuclei_excluded_inside_3v"].sum()
    )
    double_excluded = int(
        exclusion["double_nuclei_excluded_inside_3v"].sum()
    )
    nine = len(context["sheet"]) == 9
    cohort_en = ("Nine animals were analyzed (Water, Sucrose and Allulose n=3 each): eight WT animals and one male NPY-transgenic animal (NPY-M) in the Water group. ACTH/CLIP immunoreactivity was used as a POMC-related signal; the NPY-GFP channel of NPY-M was not used for ACTH/CLIP or c-FOS classification." if nine else "Eight WT animals were analyzed: Water n=2, Sucrose n=3 and Allulose n=3.")
    cohort_es = ("Se analizaron nueve animales (Agua, Sacarosa y Alulosa n=3 por grupo): ocho WT y un macho transgénico NPY (NPY-M) en Agua. La inmunorreactividad ACTH/CLIP se utilizó como señal relacionada con POMC; el canal NPY-GFP de NPY-M no se utilizó para clasificar ACTH/CLIP ni c-FOS." if nine else "Se analizaron ocho animales WT: Agua n=2, Sacarosa n=3 y Alulosa n=3.")
    water_en = ("All three pairwise contrasts, including Water, are provided in the source table. Regional q values use the Benjamini-Hochberg family of all finite tests in statistics_all_values.csv." if nine else "Water contrasts are descriptive because Water has n=2.")
    water_es = ("Los tres contrastes por pares, incluidos los de Agua, se entregan en la tabla fuente. Los valores q regionales corresponden al ajuste Benjamini-Hochberg de todas las pruebas con p finita de statistics_all_values.csv." if nine else "Los contrastes con Agua son descriptivos porque Agua tiene n=2.")
    english = f"""Supplementary Figure S5. POMC-related neuronal responses to oral solutions after a 16-hour fast.

{cohort_en} The experiment tested whether first exposure to water, sucrose or allulose preferentially recruited satiety-associated cells after food deprivation. ACTH/CLIP was detected with antibody F-3 (Santa Cruz Biotechnology, sc-373878; 1:100). A–C, DAPI, c-FOS and ACTH/CLIP from one representative allulose field ({context['sample_id']}). D, nuclear classifications and DAPI-guided ARC, median eminence (ME) and ventromedial nucleus boundaries. E, merged channels and a magnified cell. F, ACTH/CLIP signal within accepted cellular regions. G, c-FOS-positive/DAPI-positive fractions in ME and ARC. H, c-FOS-positive fractions among ACTH/CLIP-positive cells in ARC. Points represent animals and bars show mean ± SD. Regional ANOVA and Kruskal–Wallis p values and the displayed sucrose–allulose Mann–Whitney q and Cliff’s delta are reported in the panels. All pairwise comparisons and the Benjamini–Hochberg family of regional tests are provided in the source tables.

I–J, spatial occurrence of c-FOS-positive and c-FOS/ACTH/CLIP double-positive nuclei, respectively. Light dashed lines mark the external tissue envelope; white dashed lines show a representative human-delineated third ventricle for each condition. Each animal’s ventricular mask excludes intraventricular centroids before analysis ({dapi_excluded} DAPI, {cfos_excluded} c-FOS-positive and {double_excluded} double-positive objects excluded). Six covariance-normalized shells contain approximately equal DAPI counts. Thin curves represent animals; thick curves and heat maps show condition means. Arcsine-square-root-transformed profiles were compared by exact animal-label PERMANOVA ({allocations} allocations); BH adjustment covers the two spatial endpoints. c-FOS: p={cfos_p}, q={cfos_q}; double-positive: p={double_p}, q={double_q}. ACTH/CLIP regions below the acquisition-specific fifth percentile of DAPI nuclear area were excluded as likely fragments ({rejected}/{raw}). The animal is the independent unit. Scale bars are shown within microscopy panels.
"""
    spanish = f"""Figura suplementaria S5. Respuesta de células relacionadas con POMC a soluciones orales después de 16 horas de ayuno.

{cohort_es} El experimento evaluó si la primera ingesta de agua, sacarosa o alulosa reclutaba preferentemente células vinculadas a la saciedad después de la privación de alimento. Se utilizó ACTH/CLIP F-3 (Santa Cruz Biotechnology, sc-373878; 1:100). A–C, DAPI, c-FOS y ACTH/CLIP de un campo representativo de alulosa ({context['sample_id']}). D, clasificación nuclear y límites de ARC, eminencia media (ME) y núcleo ventromedial, delimitados sobre DAPI. E, canales combinados y ampliación celular. F, señal ACTH/CLIP en regiones celulares aceptadas. G, fracción c-FOS/DAPI en ME y ARC. H, fracción c-FOS positiva entre células ACTH/CLIP positivas del ARC. Los puntos representan animales; las barras, media ± DE. Los paneles muestran p de ANOVA y Kruskal–Wallis, q de Mann–Whitney para sacarosa–alulosa y delta de Cliff. Las tablas fuente contienen todas las comparaciones por pares y la familia regional ajustada por Benjamini–Hochberg.

I–J, ocurrencia espacial de núcleos c-FOS positivos y dobles c-FOS/ACTH/CLIP. Las líneas discontinuas claras delimitan el tejido; las blancas muestran un tercer ventrículo representativo delineado manualmente para cada condición. La máscara ventricular de cada animal excluyó centroides intraventriculares antes del análisis ({dapi_excluded} objetos DAPI, {cfos_excluded} c-FOS positivos y {double_excluded} dobles). Se definieron seis anillos normalizados por covarianza con cantidades DAPI aproximadamente iguales. Las líneas finas representan animales y las gruesas y mapas de calor, medias grupales. Los perfiles transformados mediante arcoseno-raíz cuadrada se compararon con PERMANOVA exacta a nivel animal ({allocations} asignaciones), con ajuste BH para dos variables. c-FOS: p={cfos_p}, q={cfos_q}; doble positividad: p={double_p}, q={double_q}. Se excluyeron regiones ACTH/CLIP menores que el percentil 5 del área nuclear DAPI de su adquisición como posibles fragmentos ({rejected}/{raw}). La unidad independiente es el animal. Las barras de escala se indican en las micrografías.
"""
    (output_dir / "Figure_S5_LEGEND.txt").write_text(english, encoding="utf-8")
    (output_dir / "Figure_S5_LEGEND_spanish.txt").write_text(spanish, encoding="utf-8")


def validate_outputs(output_dir: Path, paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty Supplementary Figure S5 output: {path}")
        record = {
            "path": str(path),
            "relative_path": str(path.relative_to(output_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        if path.suffix.lower() == ".png":
            with Image.open(path) as image:
                dpi = image.info.get("dpi", (0.0, 0.0))
                record["width_px"] = image.width
                record["height_px"] = image.height
                record["dpi_x"] = float(dpi[0])
                record["dpi_y"] = float(dpi[1])
                if min(float(dpi[0]), float(dpi[1])) < 599.0:
                    raise RuntimeError(f"PNG is below 600 dpi: {path} {dpi}")
        rows.append(record)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Render an English S5 master plus bilingual isolated A-J panels at 600 dpi.")
    parser.add_argument("--workdir", type=Path, default=Path("analyses/FigS4/july_cfos_acth_clip_work_20260828_v1"))
    parser.add_argument("--output-dir", type=Path, default=Path("analyses/FigS5/figure/acth_clip_candidate_20260828_v1"))
    parser.add_argument(
        "--spatial-hil-dir", type=Path, default=None,
        help=("Finalized dedicated DAPI review directory. Its reviewed "
              "3V masks exclude nuclei and density from I/J spatial analyses "
              "while remaining excluded from outer-tissue crop and "
              "normalization geometry; "
              "I/J use the accepted DAPI anatomy-HIL outer envelope."),
    )
    parser.add_argument(
        "--validate-spatial-hil-source-only",
        action="store_true",
        help=(
            "Hash-validate the finalized public DAPI-only I/J HIL inputs "
            "without requiring a generated analysis work tree, then exit."
        ),
    )
    parser.add_argument(
        "--validate-spatial-hil-only",
        action="store_true",
        help=(
            "Hash-validate the finalized DAPI-only I/J HIL against the "
            "analysis sheet, print a short receipt, and exit without rendering."
        ),
    )
    parser.add_argument("--dpi", type=int, default=DPI)
    args = parser.parse_args()
    if args.dpi != DPI:
        raise SystemExit(f"Supplementary Figure S5 publication graphics are fixed at {DPI} dpi")
    workdir = args.workdir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if args.validate_spatial_hil_source_only and args.validate_spatial_hil_only:
        parser.error("choose only one spatial HIL validation mode")
    if args.validate_spatial_hil_source_only:
        if args.spatial_hil_dir is None:
            raise SystemExit(
                "--validate-spatial-hil-source-only requires --spatial-hil-dir"
            )
        spatial_hil_dir = args.spatial_hil_dir.expanduser().resolve()
        spatial_masks, spatial_contours, spatial_receipt = (
            validate_spatial_tissue_hil_source(spatial_hil_dir)
        )
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "validation": "source_only_public_inputs",
                    "review_status": spatial_receipt["status"],
                    "review_schema": spatial_receipt["schema"],
                    "channels_used": spatial_receipt["channels_used"],
                    "panels": spatial_receipt["panels"],
                    "accepted_acquisitions": len(spatial_masks),
                    "point_contours": len(spatial_contours),
                    "review_dir": str(spatial_hil_dir),
                },
                indent=2,
            )
        )
        return 0
    if args.validate_spatial_hil_only:
        if args.spatial_hil_dir is None:
            raise SystemExit(
                "--validate-spatial-hil-only requires --spatial-hil-dir"
            )
        validation_context = load_context(workdir)
        spatial_hil_dir = args.spatial_hil_dir.expanduser().resolve()
        spatial_masks, spatial_contours, spatial_receipt = (
            load_spatial_tissue_hil(
                spatial_hil_dir, validation_context["sheet"], workdir
            )
        )
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "review_status": spatial_receipt["status"],
                    "review_schema": spatial_receipt["schema"],
                    "channels_used": spatial_receipt["channels_used"],
                    "panels": spatial_receipt["panels"],
                    "receipt_recorded_ventricular_contour_available": spatial_receipt[
                        "ventricular_contour_available"
                    ],
                    "audit_reclassification": (
                        "reviewed v1 contours follow the 3V cavity"
                    ),
                    "geometry_use": (
                        "3V exclusion mask plus display outline; excluded from "
                        "outer-tissue crop and normalization geometry"
                    ),
                    "accepted_acquisitions": len(spatial_masks),
                    "point_contours": len(spatial_contours),
                    "review_dir": str(spatial_hil_dir),
                },
                indent=2,
            )
        )
        return 0
    if args.spatial_hil_dir is None:
        raise SystemExit(
            "Production S5 rendering requires --spatial-hil-dir so I/J "
            "include hash-validated 3V overlays"
        )
    if output_dir.exists():
        raise FileExistsError(f"Choose an absent Supplementary Figure S5 output directory: {output_dir}")
    output_dir.mkdir(parents=True)

    script_path = Path(__file__).resolve()
    paper_root = next(
        parent for parent in script_path.parents
        if (
            (parent / "README.txt").is_file()
            and (parent / "scripts" / "setup").is_dir()
            and (parent / "FigS5").is_dir()
        )
    )
    result_dir = workdir / "results"
    source_dir = output_dir / "source_data"
    legend_dir = output_dir / "legends"
    provenance_dir = output_dir / "provenance"
    for directory in (source_dir, legend_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=False)

    context = load_context(workdir)
    spatial_hil_dir = args.spatial_hil_dir.expanduser().resolve()
    spatial_masks, spatial_contours, spatial_receipt = (
        load_spatial_tissue_hil(spatial_hil_dir, context["sheet"], workdir)
    )
    context["spatial_ventricle_masks"] = spatial_masks
    context["spatial_ventricle_contours"] = spatial_contours
    context["spatial_tissue_hil_dir"] = spatial_hil_dir
    context["spatial_tissue_hil_receipt"] = spatial_receipt
    context["spatial_tissue_boundary_source"] = (
        "outer envelope of accepted DAPI-only ARC/ME/VMN anatomy HIL"
    )
    context["spatial_review_geometry_use"] = (
        "native reviewed 3V masks exclude nuclei/density from I/J; condition "
        "medoids are visible outlines; outer-tissue crop/normalization unchanged"
    )
    (
        context["condition_ventricle_contours"],
        context["spatial_ventricle_contour_rows"],
    ) = select_condition_ventricle_contours(context)
    radial_features, radial_statistics, shell_audit, permutation_receipts = spatial_radial_analysis(context)
    context["spatial_radial_features"] = radial_features
    context["spatial_radial_statistics"] = radial_statistics

    exclusion_columns = [
        "acquisition_id",
        "animal_id",
        "condition",
        "arc_me_dapi_nuclei_before_3v_exclusion",
        "dapi_nuclei_excluded_inside_3v",
        "cfos_nuclei_excluded_inside_3v",
        "double_nuclei_excluded_inside_3v",
        "arc_me_dapi_nuclei_after_3v_exclusion",
        "ventricle_exclusion_source",
    ]
    exclusion_audit = shell_audit.loc[:, exclusion_columns].copy()
    context["spatial_3v_exclusion_audit"] = exclusion_audit
    english_panels, spatial_rows = render_panels(context, output_dir, "en")
    render_panels(context, output_dir, "es")
    master_png, master_pdf = compose_master(english_panels, output_dir, "en")
    pd.DataFrame(spatial_rows).to_csv(source_dir / "spatial_occurrence_source.csv", index=False)
    pd.DataFrame(context["spatial_ventricle_contour_rows"]).to_csv(
        source_dir / "spatial_ventricle_contour_source.csv", index=False
    )
    radial_features.to_csv(source_dir / "spatial_radial_shell_values.csv", index=False)
    radial_statistics.to_csv(source_dir / "spatial_exact_permanova.csv", index=False)
    shell_audit.to_csv(source_dir / "spatial_radial_shell_audit.csv", index=False)
    exclusion_audit.to_csv(
        source_dir / "spatial_3v_exclusion_audit.csv", index=False
    )
    permutation_receipts.to_csv(source_dir / "spatial_exact_permutation_receipts.csv", index=False)
    write_legends(context, legend_dir)

    for name in SOURCE_TABLE_NAMES:
        source = result_dir / name
        if not source.is_file():
            raise FileNotFoundError(f"Required Supplementary Figure S5 source table is missing: {source}")
        shutil.copy2(source, source_dir / name)

    all_outputs = sorted(output_dir.rglob("*.png")) + sorted(output_dir.rglob("*.pdf"))
    manifest_rows = validate_outputs(output_dir, all_outputs)
    output_manifest = provenance_dir / "output_manifest.csv"
    pd.DataFrame(manifest_rows).to_csv(output_manifest, index=False)

    correction_path = None
    if context["spatial_tissue_hil_dir"] is not None:
        correction_path = (
            provenance_dir / "SPATIAL_BOUNDARY_SEMANTICS_CORRECTION.json"
        )
        correction_record = {
            "schema": "figs5_spatial_boundary_semantics_correction_v3",
            "figure": "S5",
            "finding": (
                "The dedicated DAPI-review polygons follow the 3V cavity, "
                "not the outer tissue boundary."
            ),
            "cause_of_rejected_render": (
                "The polygons were incorrectly consumed as filled tissue "
                "masks, changing the normalization box and drawing a cavity "
                "contour through the tissue."
            ),
            "corrective_action": (
                "Preserve and hash-validate the immutable review; reclassify "
                "each raster polygon as a negative 3V mask; exclude in-mask "
                "nuclear centroids before I/J shell construction and endpoint "
                "calculation; force displayed density to zero inside the "
                "condition-medoid 3V; use only accepted DAPI ARC/ME/VMN anatomy "
                "for the outer-tissue crop and normalization box."
            ),
            "raw_images_modified": False,
            "raw_segmentations_modified": False,
            "review_annotations_modified": False,
            "legacy_raster_masks_used_for_tissue_geometry": False,
            "reviewed_raster_masks_used_as_3v_exclusion": True,
            "ventricular_contour_drawn": True,
            "ventricular_contour_role": (
                "condition-medoid visible outline and zero-density display "
                "mask; native per-acquisition masks exclude nuclei from I/J"
            ),
            "ventricular_contour_selection": (
                "minimum summed pairwise Jaccard distance among observed "
                "condition contours; acquisition ID breaks exact ties"
            ),
            "selected_acquisitions_by_condition": {
                condition: context["condition_ventricle_contours"][
                    condition
                ]["acquisition_id"]
                for condition in CONDITIONS
            },
            "review_dir": str(context["spatial_tissue_hil_dir"]),
            "review_schema": context["spatial_tissue_hil_receipt"]["schema"],
            "review_receipt_sha256": sha256_file(
                context["spatial_tissue_hil_dir"]
                / "SPATIAL_TISSUE_HIL_RECEIPT.json"
            ),
            "ventricle_source_table": (
                "source_data/spatial_ventricle_contour_source.csv"
            ),
            "ventricle_source_table_sha256": sha256_file(
                source_dir / "spatial_ventricle_contour_source.csv"
            ),

            "ventricle_exclusion_audit": (
                "source_data/spatial_3v_exclusion_audit.csv"
            ),
            "ventricle_exclusion_audit_sha256": sha256_file(
                source_dir / "spatial_3v_exclusion_audit.csv"
            ),
            "excluded_totals": {
                "dapi_nuclei": int(
                    exclusion_audit[
                        "dapi_nuclei_excluded_inside_3v"
                    ].sum()
                ),
                "cfos_nuclei": int(
                    exclusion_audit[
                        "cfos_nuclei_excluded_inside_3v"
                    ].sum()
                ),
                "double_nuclei": int(
                    exclusion_audit[
                        "double_nuclei_excluded_inside_3v"
                    ].sum()
                ),
            },
        }
        correction_path.write_text(
            json.dumps(correction_record, indent=2) + "\n", encoding="utf-8"
        )
    if len(context["sheet"]) == 9:
        sensitivity_dir = paper_root / "FigS5/provenance/inclusion_NPY_M_20260908/wt_only_reference"
        sensitivity_receipt = json.loads((sensitivity_dir / "REFERENCE_RECEIPT.json").read_text())
        for filename, expected_hash in sensitivity_receipt["files_sha256"].items():
            source = sensitivity_dir / filename
            if sha256_file(source) != expected_hash:
                raise ValueError(f"WT-only sensitivity reference changed: {source}")
            shutil.copy2(source, source_dir / filename)
    ancillary_paths = sorted(source_dir.iterdir()) + sorted(legend_dir.iterdir())
    if correction_path is not None:
        ancillary_paths.append(correction_path)
    ancillary = [
        {
            "relative_path": str(path.relative_to(output_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in ancillary_paths if path.is_file()
    ]
    analysis_script = paper_root / "FigS5" / "02a_july_cfos_acth_clip_human_in_loop_s5.py"
    renderer_script = Path(__file__).resolve()
    anatomy_dir = "human_final_including_NPY_M_20260908" if len(context["sheet"]) == 9 else "human_final_20260828_v1"
    hil_receipt = paper_root / "FigS5" / "hil_review" / anatomy_dir / "HIL_RECEIPT.json"
    if not hil_receipt.is_file():
        raise FileNotFoundError(f"Accepted Supplementary Figure S5 HIL receipt is missing: {hil_receipt}")
    receipt = {
        "schema_version": 10,
        "figure": "S5",
        "source_dataset": "FigS5",
        "dpi": DPI,
        "workdir": str(workdir),
        "output_dir": str(output_dir),
        "experiment": "16 h fasting; first-ever solution exposure",
        "antibody": "ACTH/CLIP",
        "cohort": {"WT": {"Water": 2, "Sucrose": 3, "Allulose": 3}, "NPY-Tg": {"Water": 1} if len(context["sheet"]) == 9 else {}, "excluded": [] if len(context["sheet"]) == 9 else ["NPY-M"], "counts": context["sheet"]["condition"].value_counts().to_dict()},
        "representative_acquisition": context["acquisition_id"],
        "representative_sample": context["sample_id"],
        "representative_condition": context["condition"],
        "representative_rule": "Allulose animal closest to the Allulose median ACTH/CLIP activation fraction",
        "acth_clip_size_rule": "reject ROI area below acquisition-specific DAPI nuclear-area fifth percentile",
        "raw_cellpose_masks_modified": False,
        "panel_layout": {
            "microscopy_2x3": [["A", "B"], ["C", "D"], ["E", "F"]],
            "right_statistics_vertical": ["G", "H"],
            "bottom_spatial": ["I", "J"],
        },
        "panel_f_title_shortened": True,
        "panel_h_regions": ["ARC"],
        "panel_f_roi_gated_natural_intensity": True,
        "master_title_omits_figure_number": True,
        "panel_letters_use_fixed_figure_coordinates": True,
        "spatial_colorbars_use_dedicated_axes": True,
        "inset_delimiter_linewidth_pt": 3.2,
        "heatmap_tissue_boundary": (
            "condition-consensus outer envelope of accepted DAPI-only "
            "ARC/ME/VMN anatomy masks; holes filled; thin light dashed line"
        ),
        "heatmap_ventricle_overlay": (
            "condition-medoid accepted HIL 3V point contour; bold white "
            "dashed line; density masked to zero inside"
        ),
        "spatial_tissue_boundary_source": context["spatial_tissue_boundary_source"],
        "spatial_dedicated_review_geometry_used": False,
        "spatial_dedicated_review_used_for_tissue_geometry": False,
        "spatial_dedicated_review_used_for_ventricle_overlay": True,
        "spatial_dedicated_review_used_for_ventricle_exclusion": True,
        "spatial_dedicated_review_geometry_disposition": context[
            "spatial_review_geometry_use"
        ],
        "spatial_ventricle_selected_acquisitions": {
            condition: context["condition_ventricle_contours"][
                condition
            ]["acquisition_id"]
            for condition in CONDITIONS
        },
        "spatial_ventricle_selection_method": (
            "minimum summed pairwise Jaccard distance among observed "
            "condition contours; acquisition ID breaks exact ties"
        ),
        "spatial_ventricle_source_table": (
            "source_data/spatial_ventricle_contour_source.csv"
        ),
        "spatial_ventricle_source_table_sha256": sha256_file(
            source_dir / "spatial_ventricle_contour_source.csv"
        ),

        "spatial_ventricle_exclusion_audit": (
            "source_data/spatial_3v_exclusion_audit.csv"
        ),
        "spatial_ventricle_exclusion_audit_sha256": sha256_file(
            source_dir / "spatial_3v_exclusion_audit.csv"
        ),
        "spatial_ventricle_excluded_totals": {
            "dapi_nuclei": int(
                exclusion_audit[
                    "dapi_nuclei_excluded_inside_3v"
                ].sum()
            ),
            "cfos_nuclei": int(
                exclusion_audit[
                    "cfos_nuclei_excluded_inside_3v"
                ].sum()
            ),
            "double_nuclei": int(
                exclusion_audit[
                    "double_nuclei_excluded_inside_3v"
                ].sum()
            ),
        },
        "spatial_tissue_hil_schema": context[
            "spatial_tissue_hil_receipt"
        ]["schema"],
        "spatial_boundary_semantics_correction": (
            None if correction_path is None
            else str(correction_path)
        ),
        "spatial_boundary_semantics_correction_sha256": (
            None if correction_path is None
            else sha256_file(correction_path)
        ),
        "spatial_tissue_hil_dir": (
            None if context["spatial_tissue_hil_dir"] is None
            else str(context["spatial_tissue_hil_dir"])
        ),
        "spatial_tissue_hil_receipt_sha256": (
            None if context["spatial_tissue_hil_dir"] is None
            else sha256_file(
                context["spatial_tissue_hil_dir"]
                / "SPATIAL_TISSUE_HIL_RECEIPT.json"
            )
        ),
        "spatial_analysis": (
            "accepted per-acquisition 3V exclusion before six "
            "covariance-normalized equal-DAPI-density shells; arcsin-sqrt "
            f"proportions; exhaustive animal-label PERMANOVA across {int(radial_statistics.iloc[0]['enumerated_labelings'])} "
            "allocations; BH across two endpoints"
        ),
        "spatial_statistics": radial_statistics.to_dict(orient="records"),
        "statistics_annotation": "G (ME/ARC) and H (ARC only): Kruskal-Wallis and one-way ANOVA; Sucrose-vs-Allulose Mann-Whitney only when n>=3/group, BH q and Cliff delta. I/J: exact PERMANOVA and two-endpoint BH q",
        "analysis_script": str(analysis_script),
        "analysis_script_sha256": sha256_file(analysis_script),
        "renderer_script": str(renderer_script),
        "renderer_script_sha256": sha256_file(renderer_script),
        "accepted_hil_receipt": str(hil_receipt),
        "accepted_hil_receipt_sha256": sha256_file(hil_receipt),
        "english_panels": 10,
        "spanish_panels": 10,
        "master_language": "English only",
        "spanish_multipanel_master": False,
        "masters": [str(master_png), str(master_pdf)],
        "outputs": manifest_rows,
        "ancillary_files": ancillary,
    }
    (provenance_dir / "BUILD_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[DONE] Supplementary Figure S5 candidate: {master_png}")
    print("[DONE] Multipanel master: English only; Spanish output is isolated panels")
    print(f"[DONE] 20 bilingual isolated panels at {DPI} dpi: {output_dir / 'panels'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
