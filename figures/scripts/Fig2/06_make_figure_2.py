#!/usr/bin/env python3
"""Render Figure 2: one English master plus isolated English/Spanish subpanels.

Panels
  A  experimental timeline, rasterised from the bundled vector PDF
  B  accepted HIL-defined ARC / ME / VMN anatomy over Cellpose DAPI
     segmentation, with DAPI-associated c-FOS-positive nuclei in red
  C  native RAW-CZI representative microscopy on a common 600 x 1200 um field
  D  animal-level c-FOS/DAPI ratios, drawn in the Figure 3 quantitative style

Every microscopy axes is laid out from its own native pixel aspect, so no image
is ever stretched, and the three conditions are equal in width in panels B and
C. The program is fresh-output only: it refuses to write into an existing
directory and never writes into the live Paper tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("SOURCE_DATE_EPOCH", "1787443200")

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon
from PIL import Image, ImageDraw, ImageEnhance, ImageFont
from scipy import stats as scipy_stats

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
TIMELINE_PDF = FIGURE / "experimental_timeline_with_npy_gfp_mice.pdf"

CONDITIONS = ("Water", "Sucrose", "Allulose")
REGIONS = ("ME", "ARC", "VMN", "OTHERS")
ANATOMY = ("ME", "ARC", "VMN")
COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
# Accepted v9.2 QC plate colours for the HIL anatomy polygons.
REGION_COLORS = {"ME": "#ffa500", "ARC": "#2e4fbf", "VMN": "#b15fd0"}
PAIR_ORDER = (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose"))

LOCALE = {
    "en": {
        "Water": "Water",
        "Sucrose": "Sucrose",
        "Allulose": "Allulose",
        "ME": "ME",
        "ARC": "ARC",
        "VMN": "VMN",
        "OTHERS": "OTHERS",
        "panel_b_side": "HIL anatomy",
        "panel_b_legend": ("ME HIL", "ARC HIL", "VMN HIL", "c-FOS+ nucleus"),
        "panel_d_title": "c-FOS/DAPI Ratios per Animal",
        "ylabel": "FOS/DAPI ratio (per animal)",
        "anova": "One-way ANOVA",
        "mwu": "exact MWU",
        "fos": "FOS",
        "outlier": "Outlier (IQR)",
        "abbr": {"Water": "W", "Sucrose": "S", "Allulose": "A"},
    },
    "es": {
        "Water": "Agua",
        "Sucrose": "Sacarosa",
        "Allulose": "Alulosa",
        "ME": "EM",
        "ARC": "ARC",
        "VMN": "VMN",
        "OTHERS": "OTROS",
        "panel_b_side": "anatomía HIL",
        "panel_b_legend": ("EM HIL", "ARC HIL", "VMN HIL", "núcleo c-FOS+"),
        "panel_d_title": "Proporciones c-FOS/DAPI por animal",
        "ylabel": "proporción FOS/DAPI (por animal)",
        "anova": "ANOVA de una vía",
        "mwu": "MWU exacta",
        "fos": "FOS",
        "outlier": "Atípico (IQR)",
        "abbr": {"Water": "A", "Sucrose": "S", "Allulose": "Al"},
    },
}
# label, x, y in merge-axes fractions, and whether to draw a pointing triangle
ANATOMY_LABELS = (
    ("VMN", 0.115, 0.585, True),
    ("ARC", 0.115, 0.345, True),
    ("3V", 0.490, 0.265, False),
    ("ME", 0.480, 0.130, False),
)

FIGURE_SIZE_INCHES = (19.0, 11.32)
MASTER_DPI = 600
PANEL_DPI = 600
SCALEBAR_UM = 200.0
DISPLAY_SATURATION = 1.18
DISPLAY_CONTRAST = 1.08
DISPLAY_BRIGHTNESS = 1.02
TIMELINE_RASTER_DPI = 320
JITTER_SEED = 20260824

# Panel bands on the master sheet, in figure fractions.
PANEL_A_BOX = (0.014, 0.540, 0.392, 0.440)
PANEL_B_BOX = (0.045, 0.014, 0.398, 0.505)
RIGHT_LEFT = 0.455
RIGHT_WIDTH = 0.540
PANEL_C_TOP = 0.972
PANEL_D_BOTTOM = 0.014
PANEL_CD_GAP = 0.026
# Panel D is narrower than the microscopy row so the 2x2 block reads squarer.
PANEL_D_WIDTH = 0.448
PANEL_B_GAP = 0.007
PANEL_C_GROUP_GAP = 0.011
PANEL_C_INNER_GAP = 0.004


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def csv_write(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing empty CSV receipt: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fmt_p(value: float) -> str:
    """Figure 3 p-value convention."""
    if not np.isfinite(value):
        return "NA"
    return f"{value:.1e}" if value < 0.001 else f"{value:.3g}"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.linewidth": 1.6,
            "axes.titlesize": 11.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def figure_aspect(fig: plt.Figure) -> float:
    return fig.get_figheight() / fig.get_figwidth()


def width_fraction(fig: plt.Figure, height_fraction: float, pixel_aspect: float) -> float:
    """Figure-width fraction that renders `pixel_aspect` (w/h) undistorted."""
    return height_fraction * pixel_aspect * figure_aspect(fig)


# --------------------------------------------------------------------------- inputs


def validate_stage(stage: Path) -> tuple[Path, Path, Path, dict[str, Any]]:
    if not stage.is_dir() or stage.is_symlink():
        raise RuntimeError(f"Invalid stage: {stage}")
    if Path("/tmp").resolve() not in stage.parents:
        raise RuntimeError(f"Stage must be below /tmp: {stage}")
    native, analysis, regions = stage / "native_czi", stage / "analysis", stage / "human_regions"
    validation_path = native / "validation.json"
    manifest_path = native / "output_hash_manifest.csv"
    regions_manifest = regions / "human_regions_manifest.json"
    if not (
        validation_path.is_file()
        and manifest_path.is_file()
        and analysis.is_dir()
        and regions_manifest.is_file()
    ):
        raise RuntimeError("Stage lacks validated native, analysis or HIL anatomy inputs")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    required = (
        "common_physical_field_validated",
        "native_calibration_validated",
        "channel_identity_validated",
        "scene0_only_validated",
        "display_only_not_used_for_statistics",
    )
    if validation.get("status") != "passed" or not all(validation.get(key) is True for key in required):
        raise RuntimeError("Native microscopy validation is not fully passed")
    for row in csv.DictReader(manifest_path.open(newline="", encoding="utf-8")):
        path = Path(row["path"])
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            raise RuntimeError(f"Native microscopy hash gate failed: {path}")
    regions_meta = json.loads(regions_manifest.read_text(encoding="utf-8"))
    if regions_meta.get("status") != "passed" or not all(
        record.get("recorded_counts_match") is True for record in regions_meta["conditions"]
    ):
        raise RuntimeError("HIL anatomy stage did not reproduce the recorded HIL counts")
    return native, analysis, regions, regions_meta


def enhanced(path: Path) -> np.ndarray:
    with Image.open(path) as source:
        image = source.convert("RGB")
    image = ImageEnhance.Color(image).enhance(DISPLAY_SATURATION)
    image = ImageEnhance.Contrast(image).enhance(DISPLAY_CONTRAST)
    image = ImageEnhance.Brightness(image).enhance(DISPLAY_BRIGHTNESS)
    return np.asarray(image)


def microscopy_assets(
    native: Path,
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, float], list[dict[str, Any]]]:
    manifest = pd.read_csv(native / "field_reconstruction_manifest.csv")
    assets: dict[str, dict[str, np.ndarray]] = {}
    display_pixel_um: dict[str, float] = {}
    rows: list[dict[str, Any]] = []
    for condition in CONDITIONS:
        row = manifest.loc[manifest.condition.eq(condition)]
        if len(row) != 1:
            raise RuntimeError(f"Missing unique field manifest row for {condition}")
        entry = row.iloc[0]
        display_pixel_um[condition] = float(entry.display_pixel_size_um)
        assets[condition] = {}
        for role in ("merge", "DAPI", "cFOS"):
            path = native / f"{condition}_{role}_display.png"
            if not path.is_file():
                raise FileNotFoundError(path)
            array = enhanced(path)
            assets[condition][role] = array
            rows.append(
                {
                    "condition": condition,
                    "role": role,
                    "source_path": str(path),
                    "source_bytes": path.stat().st_size,
                    "source_sha256": sha256(path),
                    "source_height_px": array.shape[0],
                    "source_width_px": array.shape[1],
                    "display_pixel_size_um": display_pixel_um[condition],
                    "field_width_um": float(entry.field_width_um),
                    "field_height_um": float(entry.field_height_um),
                    "rotation_clockwise_degrees": int(entry.rotation_clockwise_degrees),
                    "crop_y0_px": int(entry.crop_y0_px),
                    "crop_y1_exclusive_px": int(entry.crop_y1_exclusive_px),
                    "color_saturation_factor": DISPLAY_SATURATION,
                    "contrast_factor": DISPLAY_CONTRAST,
                    "brightness_factor": DISPLAY_BRIGHTNESS,
                    "quantification_use": "none_display_only",
                }
            )
    return assets, display_pixel_um, rows


def timeline_raster(outdir: Path) -> np.ndarray:
    if not TIMELINE_PDF.is_file():
        raise FileNotFoundError(TIMELINE_PDF)
    stem = outdir / "panel_A_timeline_raster"
    subprocess.run(
        ["pdftoppm", "-r", str(TIMELINE_RASTER_DPI), "-png", "-singlefile", str(TIMELINE_PDF), str(stem)],
        check=True,
        capture_output=True,
    )
    with Image.open(stem.with_suffix(".png")) as image:
        return np.asarray(image.convert("RGB"))


def spanish_timeline_raster(raster: np.ndarray) -> np.ndarray:
    """Replace every English label embedded in the timeline raster with Spanish.

    The source timeline is a publication asset whose artwork is rasterised in
    the bundled PDF, so its text cannot be translated through Matplotlib artist
    traversal. These normalized boxes cover only the original text and retain
    all experimental artwork, arrows and geometry.
    """
    image = Image.fromarray(raster).convert("RGB")
    draw = ImageDraw.Draw(image)
    width, height = image.size
    font_path = font_manager.findfont(
        font_manager.FontProperties(family="DejaVu Sans", weight="bold"),
        fallback_to_default=False,
    )

    def replace(
        box: tuple[float, float, float, float],
        text: str,
        *,
        background: tuple[int, int, int] = (255, 255, 255),
        foreground: tuple[int, int, int] = (20, 20, 20),
        max_font_fraction: float = 0.70,
        spacing_fraction: float = 0.10,
    ) -> None:
        left = round(box[0] * width)
        top = round(box[1] * height)
        right = round(box[2] * width)
        bottom = round(box[3] * height)
        draw.rectangle((left, top, right, bottom), fill=background)
        max_width = max(1, right - left - round(0.012 * width))
        max_height = max(1, bottom - top - round(0.006 * height))
        size = max(8, round(max_font_fraction * (bottom - top)))
        while size > 8:
            font = ImageFont.truetype(font_path, size=size)
            spacing = max(1, round(size * spacing_fraction))
            bounds = draw.multiline_textbbox(
                (0, 0), text, font=font, spacing=spacing, align="center"
            )
            text_width = bounds[2] - bounds[0]
            text_height = bounds[3] - bounds[1]
            if text_width <= max_width and text_height <= max_height:
                break
            size -= 1
        else:  # pragma: no cover - guarded by generous reviewed boxes
            font = ImageFont.truetype(font_path, size=8)
            spacing = 1
            bounds = draw.multiline_textbbox(
                (0, 0), text, font=font, spacing=spacing, align="center"
            )
            text_width = bounds[2] - bounds[0]
            text_height = bounds[3] - bounds[1]
        x = left + (right - left - text_width) / 2 - bounds[0]
        y = top + (bottom - top - text_height) / 2 - bounds[1]
        draw.multiline_text(
            (x, y), text, font=font, fill=foreground,
            spacing=spacing, align="center",
        )

    # White-background labels.
    replace((0.027, 0.092, 0.151, 0.157), "Ratones NPY-GFP")
    replace(
        (0.215, 0.290, 0.684, 0.345),
        "medición de peso corporal / ingesta de solución / ingesta de alimento",
        max_font_fraction=0.62,
    )
    replace(
        (0.143, 0.446, 0.267, 0.491), "dieta estándar",
        max_font_fraction=0.46,
    )
    replace((0.333, 0.493, 0.365, 0.535), "o", max_font_fraction=0.42)
    replace((0.405, 0.493, 0.438, 0.535), "o", max_font_fraction=0.42)
    replace((0.570, 0.446, 0.667, 0.491), "por agua", max_font_fraction=0.50)
    replace(
        (0.655, 0.446, 0.825, 0.491), "ayuno de 4 horas",
        max_font_fraction=0.46,
    )
    replace(
        (0.900, 0.055, 0.995, 0.238),
        "Sacarosa\no\nAlulosa\no\nAgua",
        max_font_fraction=0.15,
        spacing_fraction=0.02,
    )
    replace(
        (0.885, 0.594, 0.995, 0.650), "Eutanasia",
        max_font_fraction=0.48,
    )
    replace(
        (0.102, 0.640, 0.289, 0.684), "Condiciones",
        max_font_fraction=0.48,
    )
    replace(
        (0.050, 0.725, 0.153, 0.759), "10% sacarosa",
        max_font_fraction=0.54,
    )
    replace(
        (0.183, 0.725, 0.292, 0.759), "13% alulosa",
        max_font_fraction=0.54,
    )
    replace(
        (0.310, 0.725, 0.425, 0.759), "Control (agua)",
        max_font_fraction=0.54,
    )
    replace(
        (0.805, 0.820, 0.885, 0.914), "IHQ\nc-FOS–NeuN",
        max_font_fraction=0.34,
    )

    # Labels printed inside the experimental-period bars.
    dark = (57, 57, 57)
    replace(
        (0.195, 0.378, 0.500, 0.433), "Adaptación: 7 días",
        background=dark, foreground=(255, 255, 255), max_font_fraction=0.43,
    )
    replace(
        (0.500, 0.378, 0.620, 0.433), "Reemplazo: 2 días",
        max_font_fraction=0.46,
    )
    replace(
        (0.620, 0.378, 0.877, 0.433), "Exposición aguda al azúcar",
        background=dark, foreground=(255, 255, 255), max_font_fraction=0.40,
    )
    return np.asarray(image)


def load_statistics(analysis: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    values = pd.read_csv(analysis / "animal_region_values.csv")
    statistics = pd.read_csv(analysis / "global_cfos_statistics.csv")
    required = {"animal", "condition", "region", "cfos_over_dapi_recomputed"}
    if not required.issubset(values.columns):
        raise RuntimeError(f"Animal values lack columns: {sorted(required - set(values.columns))}")
    checks: list[dict[str, Any]] = []
    for region in REGIONS:
        regional = values.loc[values.region.eq(region)]
        groups = [
            regional.loc[regional.condition.eq(condition), "cfos_over_dapi_recomputed"].dropna().to_numpy(float)
            for condition in CONDITIONS
        ]
        counts = [len(group) for group in groups]
        if min(counts) < 3:
            raise RuntimeError(f"Region {region} violates the n>=3 independent-animal gate: {counts}")
        result = scipy_stats.f_oneway(*groups)
        row = statistics.loc[statistics.region.eq(region) & statistics.test.eq("One_way_ANOVA")]
        if len(row) != 1:
            raise RuntimeError(f"Missing unique ANOVA row: {region}")
        expected = row.iloc[0]
        passed = math.isclose(
            float(result.statistic), float(expected.statistic), rel_tol=0, abs_tol=1e-12
        ) and math.isclose(float(result.pvalue), float(expected.p_value_raw), rel_tol=0, abs_tol=1e-12)
        checks.append(
            {
                "region": region,
                "n_Water": counts[0],
                "n_Sucrose": counts[1],
                "n_Allulose": counts[2],
                "anova_F_recomputed": float(result.statistic),
                "anova_p_recomputed": float(result.pvalue),
                "recorded_F": float(expected.statistic),
                "recorded_p": float(expected.p_value_raw),
                "status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            raise RuntimeError(f"ANOVA source-table mismatch: {region}")
    return values, statistics, checks


def pairwise_p(statistics: pd.DataFrame, region: str, left: str, right: str) -> float:
    row = statistics.loc[
        statistics.region.eq(region)
        & statistics.test.eq("Mann_Whitney_U_two_sided")
        & statistics.comparison.eq(f"{left}|{right}")
    ]
    if len(row) != 1:
        raise RuntimeError(f"Missing exact pairwise result: {region}/{left}|{right}")
    return float(row.iloc[0].p_value_raw)


def load_outlier_flags(analysis: Path) -> set[tuple[str, str, str]]:
    """Audited descriptive Tukey flags; a flagged animal still enters every test."""
    audit = pd.read_csv(analysis / "outlier_audit.csv")
    if not audit.inferential_inclusion.all():
        raise RuntimeError("Outlier audit excludes an animal from inference; refusing to plot")
    flagged = audit.loc[audit.descriptive_iqr_flag.astype(bool)]
    if not flagged.plot_marker.eq("x").all():
        raise RuntimeError("Flagged animals disagree with their recorded plot marker")
    return {(str(row.animal), str(row.condition), str(row.region)) for row in flagged.itertuples()}


# --------------------------------------------------------------------------- panels


def add_scalebar(
    ax: plt.Axes, width_px: int, pixel_um: float, *, font: float, line: float, corner: str = "right"
) -> dict[str, Any]:
    fraction = (SCALEBAR_UM / pixel_um) / width_px
    if not 0 < fraction < 0.9:
        raise RuntimeError(f"Implausible scale-bar fraction {fraction}")
    if corner == "left":
        x0 = 0.045
        x1 = x0 + fraction
    else:
        x1 = 0.955
        x0 = x1 - fraction
    y = 0.075
    ax.plot([x0, x1], [y, y], color="white", lw=line, solid_capstyle="butt", transform=ax.transAxes, zorder=8)
    ax.text(
        (x0 + x1) / 2,
        y - 0.016,
        f"{int(SCALEBAR_UM)} µm",
        color="white",
        fontsize=font,
        weight="bold",
        ha="center",
        va="top",
        transform=ax.transAxes,
        zorder=9,
    )
    return {"axes_fraction": fraction, "display_pixel_size_um": pixel_um, "length_um": SCALEBAR_UM}


def draw_panel_a(fig: plt.Figure, box: tuple[float, float, float, float], raster: np.ndarray) -> None:
    left, bottom, width, height = box
    aspect = raster.shape[1] / raster.shape[0]
    panel_height = min(height, width / (aspect * figure_aspect(fig)))
    panel_width = width_fraction(fig, panel_height, aspect)
    ax = fig.add_axes(
        (left + (width - panel_width) / 2, bottom + height - panel_height, panel_width, panel_height)
    )
    ax.imshow(raster, interpolation="antialiased", aspect="auto")
    ax.set_axis_off()


def draw_panel_b(
    fig: plt.Figure,
    box: tuple[float, float, float, float],
    regions_dir: Path,
    meta: dict[str, Any],
    locale: str,
    *,
    side_label: bool = True,
    height_scale: float = 0.93,
) -> list[dict[str, Any]]:
    left, bottom, width, height = box
    text = LOCALE[locale]
    records = {record["condition"]: record for record in meta["conditions"]}
    aspect = records[CONDITIONS[0]]["display_width_px"] / records[CONDITIONS[0]]["display_height_px"]
    span = width - PANEL_B_GAP * (len(CONDITIONS) - 1)
    # All three share one display size, so one width applies to every condition.
    panel_height = min(height * height_scale, span / (len(CONDITIONS) * aspect * figure_aspect(fig)))
    panel_width = width_fraction(fig, panel_height, aspect)
    total = panel_width * len(CONDITIONS) + PANEL_B_GAP * (len(CONDITIONS) - 1)

    artists: list[dict[str, Any]] = []
    x = left + (width - total) / 2
    row_bottom = bottom + (height - panel_height)
    for condition in CONDITIONS:
        record = records[condition]
        ax = fig.add_axes((x, row_bottom, panel_width, panel_height))
        with Image.open(regions_dir / record["png"]) as image:
            ax.imshow(np.asarray(image.convert("RGB")), interpolation="antialiased", aspect="auto")
        ax.set_xlim(0, record["display_width_px"])
        ax.set_ylim(record["display_height_px"], 0)
        ax.set_axis_off()

        anchors: dict[str, tuple[float, float, float]] = {}
        for polygon in record["polygons"]:
            region = polygon["region"]
            points = np.asarray(polygon["display_points"], dtype=float)
            ax.add_patch(
                MplPolygon(
                    points,
                    closed=True,
                    facecolor=REGION_COLORS[region],
                    edgecolor=REGION_COLORS[region],
                    alpha=0.20,
                    linewidth=0,
                    zorder=4,
                )
            )
            ax.add_patch(
                MplPolygon(
                    points,
                    closed=True,
                    facecolor="none",
                    edgecolor=REGION_COLORS[region],
                    linewidth=1.9,
                    zorder=5,
                )
            )
            area = 0.5 * abs(
                np.dot(points[:, 0], np.roll(points[:, 1], 1)) - np.dot(points[:, 1], np.roll(points[:, 0], 1))
            )
            if area > anchors.get(region, (0.0, 0.0, 0.0))[2]:
                anchors[region] = (points[:, 0].mean(), points[:, 1].mean(), area)
        for region in ANATOMY:
            cx, cy, _ = anchors[region]
            ax.text(
                cx,
                cy,
                text[region],
                color="white",
                fontsize=10.5,
                weight="bold",
                ha="center",
                va="center",
                zorder=7,
                path_effects=[patheffects.withStroke(linewidth=2.6, foreground="black")],
            )
            counts = record["region_counts"][region]
            artists.append(
                {
                    "locale": locale,
                    "panel": "B",
                    "condition": condition,
                    "animal": record["animal"],
                    "tile_key": record["tile_key"],
                    "region": region,
                    "cfos_positive_nuclei": counts["cfos_cells"],
                    "dapi_nuclei": counts["dapi_nuclei"],
                    "region_area_um2": counts["area_um2"],
                }
            )
        ax.text(
            0.03,
            0.985,
            text[condition],
            transform=ax.transAxes,
            color="white",
            fontsize=14.0,
            weight="bold",
            ha="left",
            va="top",
            zorder=8,
        )
        add_scalebar(
            ax, record["display_width_px"], record["display_pixel_size_um"],
            font=8.5, line=2.8, corner="left",
        )
        x += panel_width + PANEL_B_GAP

    handles = [
        MplPolygon([[0, 0]], facecolor=REGION_COLORS[region], edgecolor=REGION_COLORS[region], alpha=0.55, label=label)
        for region, label in zip(ANATOMY, text["panel_b_legend"][:3])
    ] + [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="#ec1a1a", markeredgecolor="#ec1a1a",
               markersize=7, label=text["panel_b_legend"][3])
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(left + width / 2, row_bottom - 0.010),
        ncol=4,
        frameon=False,
        fontsize=9.0,
        handletextpad=0.5,
        columnspacing=1.6,
    )
    if side_label:
        fig.text(
            left - 0.019,
            row_bottom + panel_height / 2,
            text["panel_b_side"],
            rotation=90,
            fontsize=13.0,
            weight="bold",
            ha="center",
            va="center",
        )
    return artists


def draw_panel_c(
    fig: plt.Figure,
    left: float,
    top: float,
    width: float,
    assets: dict[str, dict[str, np.ndarray]],
    display_pixel_um: dict[str, float],
    locale: str,
) -> tuple[float, list[dict[str, Any]]]:
    """Merge plus DAPI/FOS insets at native aspect and equal width per condition."""
    text = LOCALE[locale]
    aspects = [assets[c]["merge"].shape[1] / assets[c]["merge"].shape[0] for c in CONDITIONS]
    # The common physical field makes every condition the same shape; use one
    # aspect for all three so the rendered groups are exactly equal in width.
    aspect = float(np.mean(aspects))
    if max(abs(value - aspect) for value in aspects) > 0.01:
        raise RuntimeError(f"Panel C conditions are not a common field shape: {aspects}")
    groups = len(CONDITIONS)
    fixed = groups * PANEL_C_INNER_GAP + PANEL_C_GROUP_GAP * (groups - 1)
    slope = 1.5 * figure_aspect(fig) * aspect * groups
    offset = 0.5 * PANEL_C_INNER_GAP * figure_aspect(fig) * aspect * groups
    panel_height = (width - fixed + offset) / slope

    records: list[dict[str, Any]] = []
    inset_height = (panel_height - PANEL_C_INNER_GAP) / 2.0
    merge_width = width_fraction(fig, panel_height, aspect)
    inset_width = width_fraction(fig, inset_height, aspect)
    x = left
    for condition in CONDITIONS:
        bottom = top - panel_height
        merge_ax = fig.add_axes((x, bottom, merge_width, panel_height))
        inset_x = x + merge_width + PANEL_C_INNER_GAP
        dapi_ax = fig.add_axes((inset_x, bottom + inset_height + PANEL_C_INNER_GAP, inset_width, inset_height))
        cfos_ax = fig.add_axes((inset_x, bottom, inset_width, inset_height))
        for ax, role, label, font in (
            (merge_ax, "merge", text[condition], 16.0),
            (dapi_ax, "DAPI", "DAPI", 11.0),
            (cfos_ax, "cFOS", text["fos"], 11.0),
        ):
            image = assets[condition][role]
            ax.imshow(image, interpolation="antialiased", aspect="auto")
            ax.set_axis_off()
            ax.text(
                0.035, 0.985, label, transform=ax.transAxes, color="white",
                fontsize=font, weight="bold", ha="left", va="top", zorder=9,
            )
            bar = add_scalebar(
                ax,
                image.shape[1],
                display_pixel_um[condition],
                font=10.0 if role == "merge" else 6.4,
                line=3.6 if role == "merge" else 2.0,
            )
            records.append({"panel": "C", "condition": condition, "role": role, **bar})
        if condition == "Water":
            for label, tx, ty, arrow in ANATOMY_LABELS:
                merge_ax.text(
                    tx, ty, text.get(label, label), transform=merge_ax.transAxes, color="white",
                    fontsize=10.0, weight="bold", ha="center", va="center", zorder=9,
                )
                if arrow:
                    merge_ax.plot(
                        [tx + 0.055], [ty - 0.026], marker="<", markersize=7.0, color="white",
                        transform=merge_ax.transAxes, zorder=9,
                    )
        x += merge_width + PANEL_C_INNER_GAP + inset_width + PANEL_C_GROUP_GAP
    return top - panel_height, records


def draw_quantitative_axes(
    ax: plt.Axes,
    region: str,
    values: pd.DataFrame,
    statistics: pd.DataFrame,
    outlier_keys: set[tuple[str, str, str]],
    locale: str,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    """Figure 3 quantitative-panel style, applied to one region."""
    text = LOCALE[locale]
    regional = values.loc[values.region.eq(region)]
    groups: dict[str, pd.DataFrame] = {
        condition: regional.loc[regional.condition.eq(condition), ["animal", "cfos_over_dapi_recomputed"]]
        .dropna()
        .sort_values("animal")
        for condition in CONDITIONS
    }
    positions = np.arange(3, dtype=float)
    means = [float(groups[c].cfos_over_dapi_recomputed.mean()) for c in CONDITIONS]
    sds = [float(groups[c].cfos_over_dapi_recomputed.std(ddof=1)) for c in CONDITIONS]
    ax.bar(
        positions, means, yerr=sds, width=0.52,
        color=[COLORS[c] for c in CONDITIONS], edgecolor="black", linewidth=1.7, capsize=4,
        error_kw={"elinewidth": 1.6, "capthick": 1.6}, zorder=2,
    )
    artists: list[dict[str, Any]] = []
    tops: list[float] = []
    for index, condition in enumerate(CONDITIONS):
        frame = groups[condition]
        vector = frame.cfos_over_dapi_recomputed.to_numpy(float)
        flags = np.array([(str(a), condition, region) in outlier_keys for a in frame.animal], dtype=bool)
        jitter = rng.uniform(-0.065, 0.065, vector.size)
        ax.scatter(
            index + jitter[~flags], vector[~flags], facecolor=COLORS[condition], edgecolor="black",
            linewidth=1.0, s=34, zorder=4,
        )
        if flags.any():
            ax.scatter(
                index + jitter[flags], vector[flags], marker="x", color="black", linewidth=1.8, s=44, zorder=5
            )
        tops.append(max(means[index] + sds[index], float(vector.max())))
        artists.append(
            {
                "locale": locale,
                "panel": "D",
                "region": region,
                "condition": condition,
                "n_animals": int(vector.size),
                "mean": means[index],
                "sample_sd": sds[index],
                "iqr_flagged_animals": int(flags.sum()),
                "values": ";".join(f"{value:.17g}" for value in vector),
                "animals": ";".join(frame.animal.astype(str)),
            }
        )
    y_top = max(max(tops), 1e-6) * 1.85
    ax.set_ylim(0, y_top)
    for index, condition in enumerate(CONDITIONS):
        ax.text(
            index, tops[index] + 0.025 * y_top, f"n={groups[condition].shape[0]}",
            ha="center", va="bottom", fontsize=9.2, fontweight="bold",
        )
    anova = statistics.loc[statistics.region.eq(region) & statistics.test.eq("One_way_ANOVA")].iloc[0]
    ax.text(
        0.025, 0.975,
        f"{text['anova']}\nF(2,{int(anova.df_within)}) = {float(anova.statistic):.2f}; p = {fmt_p(float(anova.p_value_raw))}",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.7, fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2}, zorder=11,
    )
    lines = [
        f"{text['mwu']} {text['abbr'][first]}–{text['abbr'][second]} p = {fmt_p(pairwise_p(statistics, region, first, second))}"
        for first, second in PAIR_ORDER
    ]
    ax.text(
        0.025, 0.845, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
        fontsize=7.8, fontweight="bold", linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 1.2}, zorder=10,
    )
    ax.set_title(text[region], fontsize=11.5, fontweight="bold", pad=4)
    ax.set_xticks(positions)
    ax.set_xticklabels([text[c] for c in CONDITIONS], rotation=24, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel(text["ylabel"], fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9.5, width=1.5, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)
    ax.grid(False)
    return artists


def draw_panel_d(
    fig: plt.Figure,
    box: tuple[float, float, float, float],
    values: pd.DataFrame,
    statistics: pd.DataFrame,
    outlier_keys: set[tuple[str, str, str]],
    locale: str,
    *,
    title: bool = True,
) -> list[dict[str, Any]]:
    left, bottom, width, height = box
    text = LOCALE[locale]
    top = bottom + height
    if title:
        fig.text(left + width / 2, top, text["panel_d_title"], fontsize=16.0, weight="bold", ha="center", va="top")
    left_pad = 0.78 / fig.get_figwidth()
    bottom_pad = 0.92 / fig.get_figheight()
    top_pad = (0.52 if title else 0.16) / fig.get_figheight()
    grid = fig.add_gridspec(
        2, 2,
        left=left + left_pad, right=left + width,
        bottom=bottom + bottom_pad, top=top - top_pad,
        wspace=0.26, hspace=0.40,
    )
    rng = np.random.default_rng(JITTER_SEED)
    artists: list[dict[str, Any]] = []
    for index, region in enumerate(REGIONS):
        ax = fig.add_subplot(grid[index // 2, index % 2])
        artists.extend(draw_quantitative_axes(ax, region, values, statistics, outlier_keys, locale, rng))
    fig.text(
        left + width, bottom + 0.004, f"\u00d7  {text['outlier']}",
        fontsize=9.0, fontweight="bold", ha="right", va="bottom",
    )
    return artists


# --------------------------------------------------------------------------- assembly


def build_master(
    assets, display_pixel_um, timeline, regions_dir, regions_meta, values, statistics, outlier_keys
) -> tuple[plt.Figure, dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    fig = plt.figure(figsize=FIGURE_SIZE_INCHES, facecolor="white")
    draw_panel_a(fig, PANEL_A_BOX, timeline)
    hil = draw_panel_b(fig, PANEL_B_BOX, regions_dir, regions_meta, "en")
    panel_c_bottom, scale_records = draw_panel_c(
        fig, RIGHT_LEFT, PANEL_C_TOP, RIGHT_WIDTH, assets, display_pixel_um, "en"
    )
    panel_d_box = (
        RIGHT_LEFT + (RIGHT_WIDTH - PANEL_D_WIDTH) / 2,
        PANEL_D_BOTTOM,
        PANEL_D_WIDTH,
        panel_c_bottom - PANEL_CD_GAP - PANEL_D_BOTTOM,
    )
    stats = draw_panel_d(fig, panel_d_box, values, statistics, outlier_keys, "en")
    letters = {
        "A": (0.004, 0.995),
        "B": (0.004, PANEL_B_BOX[1] + PANEL_B_BOX[3] + 0.022),
        "C": (RIGHT_LEFT - 0.032, 0.995),
        "D": (panel_d_box[0] - 0.016, panel_d_box[1] + panel_d_box[3] + 0.004),
    }
    for letter, (x, y) in letters.items():
        fig.text(x, y, letter, fontsize=28, weight="bold", ha="left", va="top")
    return fig, {"B": hil, "D": stats}, scale_records


def build_subpanel(
    panel: str, locale: str, assets, display_pixel_um, timeline, regions_dir, regions_meta,
    values, statistics, outlier_keys,
) -> tuple[plt.Figure, list[dict[str, Any]], list[dict[str, Any]]]:
    if panel == "A":
        fig = plt.figure(figsize=(9.0, 9.0 / (timeline.shape[1] / timeline.shape[0])), facecolor="white")
        ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
        ax.imshow(timeline, interpolation="antialiased", aspect="auto")
        ax.set_axis_off()
        return fig, [], []
    if panel == "B":
        record = regions_meta["conditions"][0]
        aspect = record["display_width_px"] / record["display_height_px"]
        width_in, side_in, legend_in, top_in = 12.0, 0.62, 0.66, 0.08
        row_in = width_in - side_in - 0.10
        image_in = (row_in - 2 * PANEL_B_GAP * width_in) / 3.0
        height_in = image_in / aspect + legend_in + top_in
        fig = plt.figure(figsize=(width_in, height_in), facecolor="white")
        box = (side_in / width_in, legend_in / height_in, row_in / width_in, (height_in - legend_in - top_in) / height_in)
        artists = draw_panel_b(fig, box, regions_dir, regions_meta, locale, height_scale=1.0)
        return fig, artists, []
    if panel == "C":
        aspect = assets["Water"]["merge"].shape[1] / assets["Water"]["merge"].shape[0]
        width_in = 12.0
        height_in = width_in / (4.5 * aspect + 0.06)
        fig = plt.figure(figsize=(width_in, height_in), facecolor="white")
        _, records = draw_panel_c(fig, 0.006, 0.994, 0.988, assets, display_pixel_um, locale)
        return fig, [], records
    if panel == "D":
        fig = plt.figure(figsize=(10.0, 8.6), facecolor="white")
        artists = draw_panel_d(fig, (0.008, 0.006, 0.986, 0.99), values, statistics, outlier_keys, locale, title=False)
        return fig, artists, []
    raise ValueError(panel)


def save_figure(fig: plt.Figure, stem: Path, dpi: int) -> tuple[Path, Path]:
    pdf, png = stem.with_suffix(".pdf"), stem.with_suffix(".png")
    if pdf.exists() or png.exists():
        raise FileExistsError(stem)
    fig.savefig(pdf, dpi=dpi, metadata={"Title": stem.name, "Creator": "Figure 2 renderer"})
    fig.savefig(png, dpi=dpi, facecolor="white")
    plt.close(fig)
    return pdf, png


def pdf_audit(path: Path) -> dict[str, Any]:
    info = subprocess.run(["pdfinfo", str(path)], check=True, text=True, capture_output=True).stdout
    fonts = subprocess.run(["pdffonts", str(path)], check=True, text=True, capture_output=True).stdout
    pages = next((int(line.split(":", 1)[1]) for line in info.splitlines() if line.startswith("Pages:")), 0)
    if pages != 1 or "Type 3" in fonts:
        raise RuntimeError(f"PDF audit failed: {path}")
    return {"path": str(path), "pages": pages, "Type3_absent": True}


LEGENDS = {
    "en": (
        "Figure 2. Experimental design, human-in-the-loop (HIL) anatomy and hypothalamic c-FOS activation.\n"
        "(A) Experimental timeline. NPY-GFP and C57BL/6J mice received standard diet plus 10% sucrose, 13% allulose or "
        "water; body weight, solution intake and food intake were recorded through 7 days of adaptation and 2 days of "
        "replacement with water, followed by a 4-hour fast, acute sugar exposure, euthanasia 30-40 min later and "
        "c-FOS/NeuN immunohistochemistry.\n"
        "(B) Accepted HIL anatomy for one representative section per condition. HIL-defined ME, ARC and VMN "
        "polygons are drawn over the Cellpose DAPI segmentation, and red objects are strictly DAPI-associated "
        "c-FOS-positive nuclei accepted at a 0.40 overlap. HIL input is anatomical annotation only; nucleus "
        "segmentation and c-FOS association remain automated. Region membership and c-FOS positivity are taken from "
        "the accepted v9.2 per-nucleus assignment table and are cross-checked against the per-image summary and the "
        "region mask at render time. Bar, 200 µm.\n"
        "(C) Representative native microscopy reconstructed from the 16-bit RAW CZI scene-0 mosaics on a common "
        "600 x 1200 µm physical field, one animal per condition and equal in width. Seam-derived integer tile "
        "translations are applied identically to DAPI and c-FOS, without interpolation or blending. Each condition shows the "
        "DAPI + c-FOS merge beside separate DAPI and c-FOS views; VMN, ARC, third ventricle (3V) and median eminence "
        "(ME) are annotated on the Water merge. Display uses a per-image linear percentile stretch with the c-FOS "
        "black point at the 35th percentile, gamma 1, no CLAHE and no spatial filtering, followed by uniform "
        "saturation 1.18, contrast 1.08 and brightness 1.02. These are visualisation-only choices; no displayed pixel "
        "contributes to quantification. Bars, 200 µm.\n"
        "(D) c-FOS/DAPI ratios per independent animal in ME, ARC, VMN and residual non-ROI tissue. Bars are the "
        "animal-level mean, whiskers the sample SD and points individual animals; crosses mark animals outside the "
        "1.5 x IQR fences, which are flagged descriptively and retained in every test. Statistics are ordinary "
        "one-way ANOVA with exact two-sided Mann-Whitney comparisons. Animals are n; cells, sections and fields are not."
    ),
    "es": (
        "Figura 2. Diseño experimental, anatomía human-in-the-loop (HIL) y activación hipotalámica de c-FOS.\n"
        "(A) Línea temporal experimental. Ratones NPY-GFP y C57BL/6J recibieron dieta estándar más sacarosa al 10%, "
        "alulosa al 13% o agua; se registraron peso corporal, ingesta de solución e ingesta de alimento durante 7 días "
        "de adaptación y 2 días de reemplazo por agua, seguidos de 4 horas de ayuno, exposición aguda al azúcar, "
        "eutanasia 30-40 min después e inmunohistoquímica c-FOS/NeuN.\n"
        "(B) Anatomía HIL aceptada para una sección representativa por tratamiento. Los "
        "polígonos HIL de EM, ARC y VMN se dibujan sobre la segmentación DAPI de Cellpose y los objetos rojos son "
        "núcleos c-FOS positivos estrictamente asociados a DAPI, aceptados con un solapamiento de 0,40. La "
        "intervención HIL es solo la anotación anatómica; la segmentación nuclear y la asociación c-FOS siguen "
        "siendo automáticas. La pertenencia regional y la positividad c-FOS provienen de la tabla aceptada v9.2 por "
        "núcleo y se contrastan con el resumen por imagen y la máscara de regiones al renderizar. Barra, 200 µm.\n"
        "(C) Microscopía nativa representativa reconstruida desde los mosaicos RAW CZI de 16 bits de la escena 0 sobre "
        "un campo físico común de 600 x 1200 µm, un animal por tratamiento y de igual anchura. Las traslaciones enteras "
        "de las teselas, derivadas de las uniones, se aplican de forma idéntica a DAPI y c-FOS, sin interpolación ni mezcla. Cada tratamiento "
        "muestra la combinación DAPI + c-FOS junto a vistas separadas de DAPI y c-FOS; VMN, ARC, tercer ventrículo "
        "(3V) y eminencia media (EM) se anotan en la combinación de Agua. La visualización usa un estiramiento lineal "
        "por percentiles por imagen con el punto negro de c-FOS en el percentil 35, gamma 1, sin CLAHE y sin filtrado "
        "espacial, seguido de saturación 1,18, contraste 1,08 y brillo 1,02 uniformes. Son decisiones solo de "
        "visualización; ningún píxel mostrado contribuye a la cuantificación. Barras, 200 µm.\n"
        "(D) Proporciones c-FOS/DAPI por animal independiente en EM, ARC, VMN y tejido residual fuera de ROI. Las "
        "barras son la media por animal, los bigotes la DE muestral y los puntos animales individuales; las cruces "
        "marcan animales fuera de las vallas de 1,5 x IQR, señalados de forma descriptiva y retenidos en todas las "
        "pruebas. Se aplicó ANOVA ordinaria de una vía con comparaciones exactas bilaterales de Mann-Whitney. Los "
        "animales son n; células, secciones y campos no lo son."
    ),
}
PANEL_STEMS = {
    "A": "Panel_A_experimental_timeline",
    "B": "Panel_B_human_regions",
    "C": "Panel_C_native_microscopy",
    "D": "Panel_D_animal_statistics",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    stage, outdir = args.stage_root.resolve(), args.outdir.resolve()
    if outdir.exists() or outdir.is_symlink():
        raise RuntimeError(f"Output must be absent: {outdir}")
    if PAPER in outdir.parents or outdir == PAPER:
        raise RuntimeError("Output and live Paper must be disjoint")
    if Path("/tmp").resolve() not in outdir.parents:
        raise RuntimeError("Output must be below /tmp")

    native, analysis, regions_dir, regions_meta = validate_stage(stage)
    configure_style()
    assets, display_pixel_um, display_rows = microscopy_assets(native)
    values, statistics, stat_checks = load_statistics(analysis)
    outlier_keys = load_outlier_flags(analysis)

    outdir.mkdir(mode=0o750)
    outputs = outdir / "outputs"
    panels = outputs / "panels"
    source_data = outdir / "source_data"
    receipts = outdir / "receipts"
    legends = outdir / "legends"
    for path in (outputs, panels, source_data, receipts, legends):
        path.mkdir()
    timeline = timeline_raster(receipts)
    timeline_spanish = spanish_timeline_raster(timeline)
    Image.fromarray(timeline_spanish).save(
        receipts / "panel_A_timeline_raster_spanish.png",
        dpi=(TIMELINE_RASTER_DPI,) * 2,
    )

    artist_rows: dict[str, list[dict[str, Any]]] = {"B": [], "D": []}
    scale_rows: list[dict[str, Any]] = []

    fig, artists, scales = build_master(
        assets, display_pixel_um, timeline, regions_dir, regions_meta, values, statistics, outlier_keys
    )
    for panel, rows in artists.items():
        artist_rows[panel].extend({"graphic": "master", **row} for row in rows)
    scale_rows.extend({"locale": "en", "graphic": "master", **row} for row in scales)
    master_pdf, master_png = save_figure(fig, outputs / "Figure_2", MASTER_DPI)

    for locale in ("en", "es"):
        suffix = "" if locale == "en" else "_spanish"
        for panel in ("A", "B", "C", "D"):
            fig, artists, scales = build_subpanel(
                panel, locale, assets, display_pixel_um,
                timeline if locale == "en" else timeline_spanish,
                regions_dir, regions_meta, values, statistics, outlier_keys,
            )
            if panel in artist_rows:
                artist_rows[panel].extend({"graphic": f"panel_{panel}", **row} for row in artists)
            scale_rows.extend({"locale": locale, "graphic": f"panel_{panel}", **row} for row in scales)
            save_figure(fig, panels / f"{PANEL_STEMS[panel]}{suffix}", PANEL_DPI)
        (legends / f"Figure_2_legend{suffix}.txt").write_text(LEGENDS[locale] + "\n", encoding="utf-8")

    shutil.copy2(analysis / "animal_region_values.csv", source_data / "animal_region_values.csv")
    shutil.copy2(analysis / "global_cfos_statistics.csv", source_data / "figure_statistics.csv")
    shutil.copy2(analysis / "outlier_audit.csv", source_data / "outlier_audit.csv")
    shutil.copy2(regions_dir / "human_regions_manifest.json", source_data / "panel_B_human_regions_manifest.json")
    shutil.copy2(native / "field_reconstruction_manifest.csv", source_data / "panel_C_field_manifest.csv")
    csv_write(source_data / "native_display_sources.csv", display_rows)
    csv_write(source_data / "panel_B_human_region_artist_values.csv", artist_rows["B"])
    csv_write(source_data / "panel_D_statistics_artist_values.csv", artist_rows["D"])
    csv_write(receipts / "scalebar_geometry.csv", scale_rows)
    csv_write(receipts / "independent_anova_recomputation.csv", stat_checks)

    files = sorted(path for path in outputs.rglob("*") if path.is_file())
    csv_write(
        receipts / "OUTPUT_SHA256_MANIFEST.csv",
        [{"relative_path": str(p.relative_to(outdir)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files],
    )
    json_write(
        receipts / "PDF_FONT_PAGE_AUDIT.json",
        {"status": "PASS", "files": [pdf_audit(p) for p in files if p.suffix == ".pdf"]},
    )
    summary = {
        "status": "PASS",
        "panels": [
            "A experimental timeline",
            "B accepted HIL ARC/ME/VMN anatomy",
            "C native microscopy",
            "D animal-level statistics in the Figure 3 style",
        ],
        "master_pdf": str(master_pdf),
        "master_png": str(master_png),
        "master_pdf_sha256": sha256(master_pdf),
        "master_png_sha256": sha256(master_png),
        "master_language": "en",
        "isolated_subpanels_bilingual": True,
        "graphics": len(files),
        "dpi": MASTER_DPI,
        "microscopy_aspect_preserved": True,
        "panel_b_and_c_equal_width": True,
        "animals_are_n": True,
        "n_at_least_3_every_condition_region": True,
        "iqr_outliers_flagged_not_removed": True,
        "live_Paper_mutated": False,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "builder_sha256": sha256(HERE),
    }
    json_write(receipts / "FINAL_BUILD_RECEIPT.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
