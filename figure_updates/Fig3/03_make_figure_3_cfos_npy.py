#!/usr/bin/env python3
"""Assemble publication Figure 3 from raw microscopy and animal-level data.

Panel A is reconstructed deterministically from the source TIFF channels; it
does not reuse a PowerPoint render or screenshot. Panels B and C are generated
by 02_make_cfos_npy_cartoons.py. Panels D and E show the two retained
animal-level ratios as percentages. Panels F and G show descriptive,
Cellpose-derived spatial occurrence distributions made by script 05. Omnibus
values in D/E are ordinary one-way ANOVAs from the analysis receipt; pairwise
values are raw exact two-sided Mann-Whitney U tests. The removed NPY/DAPI
endpoint remains explicitly preserved in supplemental source data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.transforms import Bbox
import cv2
import numpy as np
import pandas as pd
from PIL import Image
import tifffile
from skimage.measure import find_contours, regionprops


HERE = Path(__file__).resolve().parent
PAPER_ROOT = next(path for path in Path(__file__).resolve().parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
DEFAULT_ANALYSIS_DIR = PAPER_ROOT / "analyses" / "Fig3" / "results" / "final_run"
DEFAULT_OUTPUT_DIR = PAPER_ROOT / "Fig3"
DEFAULT_RAW_ROOT = PAPER_ROOT / "Fig3" / "raw" / "legacy_experiment_2025_07_28"
DEFAULT_SPATIAL_DIR = PAPER_ROOT / "analyses" / "Fig3" / "results" / "spatial_distribution_run"

FIGURE_NAME = "Figure_3_cFos_NPY"
# Raised from 300 on 2026-08-26 so Figure 3's PDF rasters match the other
# figures. validate_final.py enforces a 300 dpi floor, which 600 clears.
PDF_RASTER_DPI = 600
CONDITIONS = ("Water", "Sucrose", "Allulose")
COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
PAIR_ORDER = (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose"))

ENDPOINTS = (
    (
        "D", "cfos_over_dapi", "Total c-FOS activation",
        "c-FOS-positive nuclei (% DAPI)",
        ("cfos_over_dapi", "cfos_dapi_ratio", "cfos_positive_over_dapi"),
    ),
    (
        "E", "double_over_npy", "c-FOS activation of NPY nuclei",
        "c-FOS-positive NPY nuclei (% NPY)",
        ("double_over_npy", "cfos_npy_over_npy", "cfos_npy_fraction"),
    ),
)

# This endpoint was removed from the main plate at the user's request. It is
# still normalized, validated, and exported with its ANOVA/MWU receipts so the
# previously displayed result is never silently discarded.
SUPPLEMENTAL_ENDPOINTS = (
    (
        "npy_over_dapi", "NPY-positive nuclei (% DAPI)",
        ("npy_over_dapi", "npy_dapi_ratio", "npy_positive_over_dapi"),
    ),
)
ALL_ENDPOINT_SPECS = tuple(
    (endpoint, aliases) for _, endpoint, _, _, aliases in ENDPOINTS
) + tuple((endpoint, aliases) for endpoint, _, aliases in SUPPLEMENTAL_ENDPOINTS)


# Source fields and fixed raw-to-display transforms were audited against the
# author-selected reference morphology. The matrices map display coordinates
# back into the source TIFF and are inverted by cv2.warpPerspective. Water2 is
# the exact Water source field in Image 1; its integer crop and 90-degree
# counter-clockwise rotation are shared with the explanatory cartoons.
REPRESENTATIVE_FIELDS = (
    {
        "condition": "Water",
        "sample": "WATER_NPY2",
        "paths": {
            "npy": "water/WATER_NPY2-0001_GFP.tif",
            "cfos": "water/WATER_NPY2-0001_fos.tif",
        },
        "mask_paths": {
            "dapi": "water/WATER_NPY2-0001_dapi_seg.npy",
            "cfos": "water/WATER_NPY2-0001_fos_seg.npy",
            "npy": "water/WATER_NPY2-0001_GFP_seg.npy",
        },
        "expected_analysis_counts": {
            "dapi_nuclei": 11608, "cfos_positive_nuclei": 773,
            "npy_positive_nuclei": 208, "dual_positive_nuclei": 57,
            "inset_centroid_dual_positive_nuclei": 24,
            "inset_assigned_npy_roi_contours": 23,
            "inset_dual_nuclei_with_contours": 23,
        },
        "output_shape_yx": (1082, 1680),
        "display_to_raw": (
            (0.0, -1.0, 1269.0),
            (1.0, 0.0, 1223.0),
            (0.0, 0.0, 1.0),
        ),
        "um_per_px": 0.565421533061328,
        "reference_match": (
            "exact Image1 Water source; integer raw crop x=188:1270, "
            "y=1223:2903 followed by 90-degree counter-clockwise rotation"
        ),
        "mirror_zoom_horizontally_for_homologous_side": True,
    },
    {
        "condition": "Sucrose",
        "sample": "E7_FR7-5",
        "paths": {
            "npy": "E7-FR7-5_gfp.tif",
            "cfos": "E7-FR7-5_fos.tif",
        },
        "mask_paths": {
            "dapi": "E7-FR7-5_dapi_seg.npy",
            "cfos": "E7-FR7-5_fos_seg.npy",
            "npy": "E7-FR7-5_gfp_seg.npy",
        },
        "expected_analysis_counts": {
            "dapi_nuclei": 3241, "cfos_positive_nuclei": 81,
            "npy_positive_nuclei": 302, "dual_positive_nuclei": 19,
            "inset_centroid_dual_positive_nuclei": 7,
            "inset_assigned_npy_roi_contours": 6,
            "inset_dual_nuclei_with_contours": 6,
        },
        "output_shape_yx": (848, 1370),
        "display_to_raw": (
            (0.333764, -0.942664, 666.501835),
            (0.942559, 0.333778, -27.068433),
            (0.0, 0.0, 1.0),
        ),
        "um_per_px": 0.664212669961432,
        "reference_match": "exact Image1 source field; audited fixed rotation/crop",
    },
    {
        "condition": "Allulose",
        "sample": "E8_FR6-4",
        "paths": {
            "npy": "E8-FR6-4_ALULOSA_gfp.tif",
            "cfos": "E8-FR6-4_ALULOSA_fos.tif",
        },
        "mask_paths": {
            "dapi": "E8-FR6-4_ALULOSA_dapi_seg.npy",
            "cfos": "E8-FR6-4_ALULOSA_fos_seg.npy",
            "npy": "E8-FR6-4_ALULOSA_gfp_seg.npy",
        },
        "expected_analysis_counts": {
            "dapi_nuclei": 3372, "cfos_positive_nuclei": 305,
            "npy_positive_nuclei": 494, "dual_positive_nuclei": 131,
            "inset_centroid_dual_positive_nuclei": 38,
            "inset_assigned_npy_roi_contours": 36,
            "inset_dual_nuclei_with_contours": 36,
        },
        "output_shape_yx": (866, 1374),
        "display_to_raw": (
            (0.991402, -0.125351, 59.858489),
            (-0.125762, -0.991309, 949.166712),
            (-0.000001, 0.0, 1.0),
        ),
        "um_per_px": 0.664212669961432,
        "reference_match": "best registered Image1 raw source; audited fixed reflection/rotation/crop",
    },
)

MARKER_COLORS = {
    "npy": np.array([0.0, 1.0, 0.0], dtype=float),
    "cfos": np.array([1.0, 0.0, 0.0], dtype=float),
}

# Native-mask classification contract shared with scripts 01 and 02. These
# constants determine Panel A ROI contours; display intensity never determines
# a contour or enters the c-FOS/NPY cell call.
CFOS_MIN_ROI_AREA_PX = 5
CFOS_ROI_TO_NUCLEUS_OVERLAP = 0.40
NPY_MIN_MARKER_PIXELS = 1
NPY_MIN_NUCLEUS_FRACTION = 0.05
NPY_ROI_ASSIGNMENT_STRICT_MAJORITY = 0.50
# A single clean black display stroke. The 1.56-pt geometry and native
# ROI membership are unchanged; no backing/path effect is drawn.
NPY_ROI_CONTOUR_COLOR = "#000000"
NPY_ROI_CONTOUR_WIDTH_PT = 1.56
PANEL_A_TILE_BOX_ASPECT = 0.63
# Use the previously idle side and bottom margins to enlarge the two spatial
# heatmaps without changing the 16-inch plate or distorting source pixels.
SPATIAL_ROW_WIDTH_RATIOS = (0.01, 7.72, 0.09, 7.72, 0.01)
SPATIAL_PANEL_LABEL_X = 0.008
SPATIAL_PANEL_LABEL_Y = 0.995
SPATIAL_PANEL_MIN_BOTTOM_FRACTION = 0.010
SPATIAL_SOURCE_PANEL_FIGSIZE = (6.40, 3.12)
SPATIAL_SOURCE_ASPECT_ATOL = 1.0e-3


def norm_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def canonical_condition(value: object) -> str:
    token = norm_token(value)
    aliases = {
        "water": "Water", "agua": "Water", "h2o": "Water",
        "sucrose": "Sucrose", "sacarosa": "Sucrose", "sucrosa": "Sucrose",
        "allulose": "Allulose", "alulosa": "Allulose", "psicose": "Allulose",
    }
    return aliases.get(token, str(value).strip())


def find_column(frame: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lookup = {norm_token(column): str(column) for column in frame.columns}
    for candidate in candidates:
        if norm_token(candidate) in lookup:
            return lookup[norm_token(candidate)]
    return None


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_endpoint(value: object) -> str:
    token = norm_token(value)
    aliases: dict[str, str] = {}
    for endpoint, names in ALL_ENDPOINT_SPECS:
        for name in (endpoint, *names):
            aliases[norm_token(name)] = endpoint
    aliases.update({
        "cfosdapi": "cfos_over_dapi",
        "totalcfosactivationfraction": "cfos_over_dapi",
        "npydapi": "npy_over_dapi",
        "doublepositivenpyfraction": "double_over_npy",
        "cfospositivenpynucleiovernpy": "double_over_npy",
    })
    return aliases.get(token, str(value).strip())


def _ratio(frame: pd.DataFrame, numerator_names: Sequence[str], denominator_names: Sequence[str]) -> pd.Series | None:
    numerator = find_column(frame, numerator_names)
    denominator = find_column(frame, denominator_names)
    if numerator is None or denominator is None:
        return None
    num = pd.to_numeric(frame[numerator], errors="coerce")
    den = pd.to_numeric(frame[denominator], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        result = num / den
    return result.where(np.isfinite(result) & den.gt(0))


def normalize_animal_values(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept the canonical wide table and a simple endpoint/value long table."""
    condition_col = find_column(frame, ("condition", "cond", "treatment", "group"))
    animal_col = find_column(frame, ("animal", "animal_id", "animalid", "mouse", "subject"))
    if condition_col is None or animal_col is None:
        raise ValueError("per_animal_cfos_npy.csv requires animal and condition columns")
    sample_col = find_column(frame, ("sample", "sample_id", "tile_key", "field"))
    endpoint_col = find_column(frame, ("endpoint", "metric", "measure"))
    value_col = find_column(frame, ("value", "ratio", "endpoint_value", "raw_value"))

    base = pd.DataFrame({
        "animal": frame[animal_col].astype(str).str.strip(),
        "condition": frame[condition_col].map(canonical_condition),
        "sample": frame[sample_col].astype(str).str.strip() if sample_col else "",
    })
    if endpoint_col and value_col:
        long = base.copy()
        long["endpoint"] = frame[endpoint_col].map(canonical_endpoint)
        long["value"] = pd.to_numeric(frame[value_col], errors="coerce")
        target = {endpoint for endpoint, _ in ALL_ENDPOINT_SPECS}
        long = long[long["endpoint"].isin(target)].copy()
        duplicates = long.duplicated(["animal", "condition", "endpoint"], keep=False)
        if duplicates.any():
            keys = long.loc[duplicates, ["animal", "condition", "endpoint"]].drop_duplicates()
            raise ValueError(f"Duplicate animal/condition/endpoint rows:\n{keys.to_string(index=False)}")
        values = long.pivot(index=["animal", "condition"], columns="endpoint", values="value").reset_index()
        samples = long.groupby(["animal", "condition"], sort=False)["sample"].agg(
            lambda x: ";".join(sorted({item for item in x if item}))
        ).reset_index()
        values = values.merge(samples, on=["animal", "condition"], how="left", validate="one_to_one")
    else:
        values = base.copy()
        for endpoint, aliases in ALL_ENDPOINT_SPECS:
            source = find_column(frame, (endpoint, *aliases))
            if source is not None:
                values[endpoint] = pd.to_numeric(frame[source], errors="coerce")
        derivations = {
            "cfos_over_dapi": _ratio(
                frame, ("cfos_cells", "cfos_nuclei", "cfos_positive_nuclei"),
                ("dapi_nuclei", "dapi_cells", "dapi_positive_nuclei"),
            ),
            "npy_over_dapi": _ratio(
                frame, ("npy_cells", "npy_nuclei", "npy_positive_nuclei"),
                ("dapi_nuclei", "dapi_cells", "dapi_positive_nuclei"),
            ),
            "double_over_npy": _ratio(
                frame, ("cfos_npy_cells", "double_positive_nuclei", "cfos_npy_positive_nuclei"),
                ("npy_cells", "npy_nuclei", "npy_positive_nuclei"),
            ),
        }
        for endpoint, derived in derivations.items():
            if endpoint not in values and derived is not None:
                values[endpoint] = derived

    target_columns = [endpoint for endpoint, _ in ALL_ENDPOINT_SPECS]
    missing = [column for column in target_columns if column not in values]
    if missing:
        raise ValueError(f"Missing required endpoint columns: {missing}")
    values = values[values["condition"].isin(CONDITIONS)].copy()
    if set(values["condition"]) != set(CONDITIONS):
        raise ValueError(f"Expected Water/Sucrose/Allulose; found {sorted(set(values['condition']))}")
    if values["animal"].eq("").any():
        raise ValueError("Blank animal identifier in per-animal table")
    if values.duplicated(["animal", "condition"]).any():
        raise ValueError("Per-animal table is not unique by animal and condition")
    conflicts = values.groupby("animal")["condition"].nunique()
    if (conflicts > 1).any():
        raise ValueError(f"Animals assigned to multiple conditions: {conflicts[conflicts > 1].index.tolist()}")
    for endpoint in target_columns:
        values[endpoint] = pd.to_numeric(values[endpoint], errors="coerce")
        finite = values[endpoint][np.isfinite(values[endpoint])]
        if finite.empty or (finite < -1e-12).any() or (finite > 1 + 1e-9).any():
            raise ValueError(f"{endpoint} must contain finite ratios in [0, 1]")
        counts = values.groupby("condition")[endpoint].apply(lambda x: int(np.isfinite(x).sum()))
        if any(int(counts.get(condition, 0)) < 2 for condition in CONDITIONS):
            raise ValueError(f"{endpoint} requires >=2 finite animals per condition; counts={counts.to_dict()}")
    order = {condition: index for index, condition in enumerate(CONDITIONS)}
    values["_order"] = values["condition"].map(order)
    return values.sort_values(["_order", "animal"], kind="mergesort").drop(columns="_order").reset_index(drop=True)


def _p_value(row: pd.Series) -> float:
    for name in ("p_value_raw", "p_exact_raw", "p_raw", "p_value", "p"):
        column = find_column(pd.DataFrame(columns=row.index), (name,))
        if column is not None:
            value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
            if np.isfinite(value):
                value = float(value)
                if not 0 <= value <= 1:
                    raise ValueError(f"Invalid p value {value}")
                return value
    raise ValueError("Exact-statistics row has no finite raw p value")


def _row_conditions(row: pd.Series, columns: Mapping[str, str | None]) -> tuple[str, ...]:
    first = canonical_condition(row[columns["group1"]]) if columns["group1"] else ""
    second = canonical_condition(row[columns["group2"]]) if columns["group2"] else ""
    if first in CONDITIONS and second in CONDITIONS:
        return first, second
    comparison = str(row[columns["comparison"]]) if columns["comparison"] else ""
    found: list[str] = []
    for token in re.split(r"\s*(?:\||/|,|;|\bvs\.?\b|[-–—])\s*", comparison, flags=re.I):
        condition = canonical_condition(token)
        if condition in CONDITIONS and condition not in found:
            found.append(condition)
    return tuple(found)


def exact_statistics_for_endpoint(
    stats: pd.DataFrame, endpoint: str,
) -> tuple[float, float, dict[tuple[str, str], float], pd.DataFrame]:
    endpoint_col = find_column(stats, ("endpoint", "metric", "measure"))
    test_col = find_column(stats, ("test", "test_name", "analysis"))
    statistic_col = find_column(stats, ("statistic", "f_statistic", "f"))
    if endpoint_col is None or test_col is None or statistic_col is None:
        raise ValueError("exact_statistics.csv requires endpoint, test, and statistic columns")
    endpoint_rows = stats[stats[endpoint_col].map(canonical_endpoint).eq(endpoint)].copy()
    if endpoint_rows.empty:
        raise ValueError(f"No statistics found for endpoint {endpoint}")
    columns = {
        "group1": find_column(stats, ("group1", "first_group", "condition1")),
        "group2": find_column(stats, ("group2", "second_group", "condition2")),
        "comparison": find_column(stats, ("comparison", "contrast", "groups")),
    }
    omnibus_rows: list[pd.Series] = []
    stale_kw_rows: list[pd.Series] = []
    pair_rows: dict[tuple[str, str], list[pd.Series]] = {pair: [] for pair in PAIR_ORDER}
    for _, row in endpoint_rows.iterrows():
        test = norm_token(row[test_col])
        groups = _row_conditions(row, columns)
        if "anova" in test:
            omnibus_rows.append(row)
        elif "kruskal" in test:
            stale_kw_rows.append(row)
        elif "mann" in test or "whitney" in test or "ranksum" in test:
            group_set = set(groups)
            for pair in PAIR_ORDER:
                if group_set == set(pair):
                    pair_rows[pair].append(row)
                    break
    if len(omnibus_rows) != 1:
        raise ValueError(f"Expected one ordinary one-way ANOVA row for {endpoint}; found {len(omnibus_rows)}")
    if stale_kw_rows:
        raise ValueError(f"Stale Kruskal-Wallis receipt remains for {endpoint}; refusing mixed omnibus claims")
    omnibus = omnibus_rows[0]
    f_statistic = float(pd.to_numeric(pd.Series([omnibus[statistic_col]]), errors="coerce").iloc[0])
    if not np.isfinite(f_statistic) or f_statistic < 0:
        raise ValueError(f"Invalid ANOVA F statistic for {endpoint}: {f_statistic}")
    selected = [omnibus]
    pairs: dict[tuple[str, str], float] = {}
    for pair in PAIR_ORDER:
        rows = pair_rows[pair]
        if len(rows) != 1:
            raise ValueError(f"Expected one exact two-sided MWU row for {endpoint}/{pair}; found {len(rows)}")
        selected.append(rows[0])
        pairs[pair] = _p_value(rows[0])
    return _p_value(omnibus), f_statistic, pairs, pd.DataFrame(selected)


def tukey_outliers(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.zeros(values.shape, dtype=bool)
    if values.size < 4:
        return result
    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    if not np.isfinite(iqr) or iqr <= 0:
        return result
    return (values < q1 - 1.5 * iqr) | (values > q3 + 1.5 * iqr)


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.1e}" if value < 0.001 else f"{value:.3g}"


def draw_quantitative_panel(
    ax: plt.Axes,
    values: pd.DataFrame,
    endpoint: str,
    ylabel: str,
    anova_p: float,
    anova_f: float,
    pairs: Mapping[tuple[str, str], float],
    seed: int,
) -> None:
    groups = {
        condition: 100.0 * pd.to_numeric(
            values.loc[values["condition"].eq(condition), endpoint], errors="coerce"
        ).dropna().to_numpy(float)
        for condition in CONDITIONS
    }
    positions = np.arange(3, dtype=float)
    means = [float(np.mean(groups[condition])) for condition in CONDITIONS]
    sds = [float(np.std(groups[condition], ddof=1)) if groups[condition].size >= 2 else 0.0 for condition in CONDITIONS]
    ax.bar(
        positions, means, yerr=sds, width=0.52,
        color=[COLORS[condition] for condition in CONDITIONS],
        edgecolor="black", linewidth=1.7, capsize=4,
        error_kw={"elinewidth": 1.6, "capthick": 1.6}, zorder=2,
    )
    rng = np.random.default_rng(seed)
    tops: list[float] = []
    for index, condition in enumerate(CONDITIONS):
        group = groups[condition]
        outliers = tukey_outliers(group)
        jitter = rng.uniform(-0.065, 0.065, group.size)
        ax.scatter(
            index + jitter[~outliers], group[~outliers],
            facecolor=COLORS[condition], edgecolor="black", linewidth=1.0,
            s=34, zorder=4,
        )
        if outliers.any():
            ax.scatter(
                index + jitter[outliers], group[outliers], marker="x",
                color="black", linewidth=1.8, s=44, zorder=5,
            )
        tops.append(max(means[index] + sds[index], float(np.max(group))))
    y_top = max(max(tops), 1e-4) * 1.85
    ax.set_ylim(0, y_top)
    for index, condition in enumerate(CONDITIONS):
        ax.text(
            index, tops[index] + 0.025 * y_top, f"n={groups[condition].size}",
            ha="center", va="bottom", fontsize=9.2, fontweight="bold",
        )
    df_within = sum(group.size for group in groups.values()) - len(CONDITIONS)
    ax.text(
        0.025, 0.975,
        f"One-way ANOVA\nF(2,{df_within}) = {anova_f:.2f}; p = {fmt_p(anova_p)}",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.7,
        fontweight="bold", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2}, zorder=11,
    )
    lines = [
        f"exact MWU {first[0]}–{second[0]} p = {fmt_p(pairs[(first, second)])}"
        for first, second in PAIR_ORDER
    ]
    ax.text(
        0.025, 0.865, "\n".join(lines), transform=ax.transAxes,
        ha="left", va="top", fontsize=7.8, fontweight="bold", linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 1.2}, zorder=10,
    )
    ax.set_title("Whole field", fontsize=11.5, fontweight="bold", pad=4)
    ax.set_xticks(positions)
    ax.set_xticklabels(CONDITIONS, rotation=24, ha="right", fontsize=10, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9.5, width=1.5, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)
    ax.grid(False)


def _marker_plane(path: Path) -> np.ndarray:
    array = np.asarray(tifffile.imread(path))
    if array.ndim == 3 and array.shape[-1] in (3, 4):
        array = np.max(array[..., :3], axis=-1)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2-D marker image or RGB display TIFF: {path} -> {array.shape}")
    if not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"Expected integer microscopy pixels: {path} -> {array.dtype}")
    return array.astype(np.float32, copy=False)


def _warp_to_display(source: np.ndarray, spec: Mapping[str, object]) -> np.ndarray:
    height, width = (int(value) for value in spec["output_shape_yx"])
    matrix = spec["display_to_raw"]
    if matrix is None:
        if source.shape != (height, width):
            raise ValueError(
                f"Native representative field shape changed for {spec['sample']}: "
                f"{source.shape} != {(height, width)}"
            )
        return source.copy()
    display_to_raw = np.asarray(matrix, dtype=np.float64)
    if display_to_raw.shape != (3, 3) or not np.isfinite(display_to_raw).all():
        raise ValueError(f"Invalid fixed display-to-raw transform for {spec['sample']}")
    return cv2.warpPerspective(
        source,
        display_to_raw,
        (width, height),
        flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0,
    )


def _display_stretch(source: np.ndarray) -> tuple[np.ndarray, float, float]:
    finite = source[np.isfinite(source)]
    positive = finite[finite > 0]
    population = positive if positive.size >= 100 else finite
    if population.size < 100:
        raise ValueError("Microscopy plane has too few finite pixels")
    low, high = np.percentile(population, [0.5, 99.8])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        raise ValueError(f"Invalid microscopy display range: {low}, {high}")
    scaled = np.clip((source - low) / (high - low), 0.0, 1.0)
    # A fixed display gamma is applied after the per-field/channel percentile
    # stretch. The result must not be used for cross-condition intensity claims.
    scaled = np.power(scaled, 0.85, dtype=np.float32)
    return scaled.astype(np.float32, copy=False), float(low), float(high)


def _colorize(signal: np.ndarray, color: np.ndarray) -> np.ndarray:
    return np.clip(signal[..., None] * color.reshape(1, 1, 3), 0.0, 1.0)


def _choose_zoom_box(npy_signal: np.ndarray, cfos_signal: np.ndarray) -> tuple[int, int, int, int]:
    if npy_signal.shape != cfos_signal.shape:
        raise ValueError("NPY and c-FOS display fields are not registered")
    height, width = npy_signal.shape
    box_width = max(80, int(round(width * 0.22)))
    box_height = max(80, int(round(height * 0.28)))
    score = npy_signal * cfos_signal
    sigma = max(2.0, min(height, width) * 0.012)
    score = cv2.GaussianBlur(score.astype(np.float32), (0, 0), sigmaX=sigma, sigmaY=sigma)
    margin_x = box_width // 2 + 2
    margin_y = box_height // 2 + 2
    masked = np.full(score.shape, -np.inf, dtype=np.float32)
    if height > 2 * margin_y and width > 2 * margin_x:
        masked[margin_y : height - margin_y, margin_x : width - margin_x] = score[
            margin_y : height - margin_y, margin_x : width - margin_x
        ]
    if not np.isfinite(masked).any() or float(np.nanmax(masked)) <= 0:
        center_y, center_x = height // 2, width // 2
    else:
        center_y, center_x = np.unravel_index(int(np.nanargmax(masked)), masked.shape)
    x0 = int(np.clip(center_x - box_width // 2, 0, width - box_width))
    y0 = int(np.clip(center_y - box_height // 2, 0, height - box_height))
    return x0, y0, box_width, box_height


def _load_cellpose_labels(path: Path, role: str) -> tuple[np.ndarray, dict[str, object]]:
    source = require_file(path, role)
    payload = np.load(source, allow_pickle=True)
    if not (isinstance(payload, np.ndarray) and payload.dtype == object and payload.shape == ()):
        raise ValueError(f"{role} is not a scalar Cellpose dictionary: {source}")
    record = payload.item()
    if not isinstance(record, dict) or "masks" not in record:
        raise ValueError(f"{role} has no Cellpose masks entry: {source}")
    labels = np.asarray(record["masks"])
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0):
        raise ValueError(f"{role} must contain one nonnegative integer label field: {source}")
    positive = np.unique(labels[labels > 0])
    if positive.size == 0 or not np.array_equal(
        positive, np.arange(1, int(positive[-1]) + 1, dtype=positive.dtype),
    ):
        raise ValueError(f"{role} labels are empty or noncontiguous: {source}")
    return labels, {
        "role": role,
        "path": str(source),
        "sha256": sha256_file(source),
        "bytes": source.stat().st_size,
        "shape_yx": [int(value) for value in labels.shape],
        "dtype": str(labels.dtype),
        "instance_count": int(positive.size),
        "internal_source_filename": str(record.get("filename", "")),
    }


def _analysis_positive_sets(
    dapi_labels: np.ndarray,
    cfos_labels: np.ndarray,
    npy_labels: np.ndarray,
) -> tuple[frozenset[int], frozenset[int], frozenset[int], dict[str, int]]:
    dapi = np.asarray(dapi_labels)
    cfos = np.asarray(cfos_labels)
    npy = np.asarray(npy_labels)
    if dapi.shape != cfos.shape or dapi.shape != npy.shape:
        raise ValueError("DAPI, c-FOS, and NPY native label fields do not share coordinates")
    maximum_nucleus = int(dapi.max())
    accepted_mappings: list[int] = []
    for prop in regionprops(cfos):
        if int(prop.area) < CFOS_MIN_ROI_AREA_PX:
            continue
        values = dapi[prop.coords[:, 0], prop.coords[:, 1]]
        counts = np.bincount(values.astype(np.int64), minlength=maximum_nucleus + 1)
        counts[0] = 0
        best = int(np.argmax(counts)) if counts.sum() else 0
        fraction = float(counts[best]) / float(prop.area) if best else 0.0
        if best and fraction >= CFOS_ROI_TO_NUCLEUS_OVERLAP:
            accepted_mappings.append(best)
    cfos_positive = frozenset(accepted_mappings)
    marker = npy > 0
    npy_positive: set[int] = set()
    for prop in regionprops(dapi):
        overlap = int(marker[prop.coords[:, 0], prop.coords[:, 1]].sum())
        fraction = overlap / float(max(1, prop.area))
        if overlap >= NPY_MIN_MARKER_PIXELS and fraction >= NPY_MIN_NUCLEUS_FRACTION:
            npy_positive.add(int(prop.label))
    npy_positive_frozen = frozenset(npy_positive)
    dual_positive = cfos_positive & npy_positive_frozen
    counts = {
        "dapi_nuclei": maximum_nucleus,
        "accepted_cfos_roi_mappings": len(accepted_mappings),
        "duplicate_cfos_roi_mappings": len(accepted_mappings) - len(cfos_positive),
        "cfos_positive_nuclei": len(cfos_positive),
        "npy_positive_nuclei": len(npy_positive_frozen),
        "dual_positive_nuclei": len(dual_positive),
    }
    return cfos_positive, npy_positive_frozen, dual_positive, counts


def _warp_labels_to_display(labels: np.ndarray, spec: Mapping[str, object]) -> np.ndarray:
    height, width = (int(value) for value in spec["output_shape_yx"])
    matrix = spec["display_to_raw"]
    if matrix is None:
        if labels.shape != (height, width):
            raise ValueError(
                f"Native representative label shape changed for {spec['sample']}: "
                f"{labels.shape} != {(height, width)}"
            )
        return labels.astype(np.int32, copy=True)
    display_to_raw = np.asarray(matrix, dtype=np.float64)
    if display_to_raw.shape != (3, 3) or not np.isfinite(display_to_raw).all():
        raise ValueError(f"Invalid label transform for {spec['sample']}")
    # Label IDs are below 2^24, so float32 stores them exactly. INTER_NEAREST
    # preserves instance IDs under the same fixed display geometry as the TIFF.
    displayed = cv2.warpPerspective(
        labels.astype(np.float32), display_to_raw, (width, height),
        flags=cv2.INTER_NEAREST | cv2.WARP_INVERSE_MAP,
        borderMode=cv2.BORDER_CONSTANT, borderValue=0,
    )
    displayed = np.rint(displayed).astype(np.int32)
    if displayed.shape != (height, width) or int(displayed.max()) > int(labels.max()):
        raise ValueError(f"Nearest-neighbor label transform failed for {spec['sample']}")
    return displayed


def _binary_mask_sha256(mask: np.ndarray) -> str:
    values = np.ascontiguousarray(np.asarray(mask, dtype=np.uint8))
    return hashlib.sha256(values.tobytes()).hexdigest()


def _native_npy_roi_mappings(
    dapi_labels: np.ndarray,
    npy_labels: np.ndarray,
) -> list[dict[str, object]]:
    """Map each complete native NPY ROI to DAPI by strict overlap majority."""
    dapi = np.asarray(dapi_labels)
    npy = np.asarray(npy_labels)
    if dapi.shape != npy.shape:
        raise ValueError("Native DAPI and NPY labels do not share coordinates")
    maximum_nucleus = int(dapi.max())
    dapi_area = {int(prop.label): int(prop.area) for prop in regionprops(dapi)}
    mappings: list[dict[str, object]] = []
    for prop in regionprops(npy):
        roi_id = int(prop.label)
        coordinates = np.ascontiguousarray(prop.coords, dtype=np.int32)
        values = dapi[coordinates[:, 0], coordinates[:, 1]]
        counts = np.bincount(values.astype(np.int64), minlength=maximum_nucleus + 1)
        counts[0] = 0
        overlap_ids = np.flatnonzero(counts)
        overlaps = [
            {"nucleus_id": int(nucleus_id), "overlap_pixels": int(counts[nucleus_id])}
            for nucleus_id in sorted(
                (int(value) for value in overlap_ids),
                key=lambda nucleus_id: (-int(counts[nucleus_id]), nucleus_id),
            )
        ]
        total_overlap = int(counts.sum())
        best_pixels = int(counts.max()) if total_overlap else 0
        best_candidates = [
            int(nucleus_id) for nucleus_id in overlap_ids
            if int(counts[nucleus_id]) == best_pixels
        ]
        best_nucleus = best_candidates[0] if len(best_candidates) == 1 else None
        dominance = best_pixels / float(total_overlap) if total_overlap else 0.0
        if total_overlap == 0:
            status = "excluded_no_dapi_overlap"
            ambiguity = "NPY ROI has no native DAPI-label overlap"
            assigned_nucleus = None
        elif len(best_candidates) != 1:
            status = "excluded_tied_best_overlap"
            ambiguity = "two or more DAPI nuclei tie for maximum overlap"
            assigned_nucleus = None
        elif dominance <= NPY_ROI_ASSIGNMENT_STRICT_MAJORITY:
            status = "excluded_no_strict_majority"
            ambiguity = "unique best DAPI nucleus does not contain >50% of DAPI-overlapping ROI pixels"
            assigned_nucleus = None
        else:
            status = "assigned_strict_majority"
            ambiguity = "none under the preregistered strict-majority rule"
            assigned_nucleus = int(best_nucleus)
        centroid_y, centroid_x = (float(value) for value in prop.centroid)
        mappings.append({
            "npy_roi_id": roi_id,
            "status": status,
            "ambiguity": ambiguity,
            "assigned_nucleus_id": assigned_nucleus,
            "best_nucleus_id_if_unique": int(best_nucleus) if best_nucleus is not None else None,
            "best_nucleus_candidates_if_tied": best_candidates,
            "best_overlap_pixels": best_pixels,
            "total_dapi_overlap_pixels": total_overlap,
            "best_fraction_of_dapi_overlap": dominance,
            "distinct_overlapping_dapi_nuclei": len(overlaps),
            "dapi_overlaps": overlaps,
            "native_roi_area_px": int(prop.area),
            "native_roi_centroid_xy": [centroid_x, centroid_y],
            "native_roi_coordinate_sha256": hashlib.sha256(coordinates.tobytes()).hexdigest(),
            "best_overlap_fraction_of_npy_roi_area": (
                best_pixels / float(prop.area) if int(prop.area) else 0.0
            ),
            "best_overlap_fraction_of_dapi_nucleus_area": (
                best_pixels / float(dapi_area[int(best_nucleus)])
                if best_nucleus is not None else None
            ),
        })
    return mappings


def _double_positive_npy_roi_layers(
    raw_root: Path,
    spec: Mapping[str, object],
    zoom_box: tuple[int, int, int, int],
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    sample = str(spec["sample"])
    labels_by_marker: dict[str, np.ndarray] = {}
    source_receipts: dict[str, dict[str, object]] = {}
    source_paths: list[Path] = []
    for marker in ("dapi", "cfos", "npy"):
        source = require_file(
            raw_root / str(spec["mask_paths"][marker]),
            f"{sample} {marker} Cellpose mask",
        )
        labels_by_marker[marker], source_receipts[marker] = _load_cellpose_labels(
            source, f"{sample} {marker} Cellpose mask",
        )
        source_paths.append(source)
    dapi = labels_by_marker["dapi"]
    _, _, dual_positive, observed = _analysis_positive_sets(
        dapi, labels_by_marker["cfos"], labels_by_marker["npy"],
    )
    expected = dict(spec["expected_analysis_counts"])
    for key in (
        "dapi_nuclei", "cfos_positive_nuclei",
        "npy_positive_nuclei", "dual_positive_nuclei",
    ):
        if int(observed[key]) != int(expected[key]):
            raise RuntimeError(f"{sample} {key} changed: {observed[key]} != {expected[key]}")

    displayed_dapi = _warp_labels_to_display(dapi, spec)
    display_dapi_props = {int(prop.label): prop for prop in regionprops(displayed_dapi)}
    x0, y0, box_width, box_height = (int(value) for value in zoom_box)
    x1, y1 = x0 + box_width, y0 + box_height
    inset_dual_ids = sorted(
        int(nucleus_id) for nucleus_id in dual_positive
        if (
            int(nucleus_id) in display_dapi_props
            and x0 <= float(display_dapi_props[int(nucleus_id)].centroid[1]) < x1
            and y0 <= float(display_dapi_props[int(nucleus_id)].centroid[0]) < y1
        )
    )
    expected_inset = int(expected["inset_centroid_dual_positive_nuclei"])
    if len(inset_dual_ids) != expected_inset:
        raise RuntimeError(
            f"{sample} inset dual-positive centroid count changed: "
            f"{len(inset_dual_ids)} != {expected_inset}"
        )
    inset_dual_set = set(inset_dual_ids)

    mappings = _native_npy_roi_mappings(dapi, labels_by_marker["npy"])
    mapping_status_counts: dict[str, int] = {}
    for record in mappings:
        status = str(record["status"])
        mapping_status_counts[status] = mapping_status_counts.get(status, 0) + 1
    strict_inset_mappings = [
        record for record in mappings
        if (
            record["status"] == "assigned_strict_majority"
            and int(record["assigned_nucleus_id"]) in inset_dual_set
        )
    ]
    strict_inset_mappings.sort(key=lambda record: int(record["npy_roi_id"]))
    candidates_by_nucleus: dict[int, list[dict[str, object]]] = {}
    for record in strict_inset_mappings:
        nucleus_id = int(record["assigned_nucleus_id"])
        candidates_by_nucleus.setdefault(nucleus_id, []).append(record)
    selected_mappings: list[dict[str, object]] = []
    duplicate_display_exclusions: list[dict[str, object]] = []
    for nucleus_id, candidates in sorted(candidates_by_nucleus.items()):
        retained = min(
            candidates,
            key=lambda record: (
                -int(record["best_overlap_pixels"]),
                int(record["npy_roi_id"]),
            ),
        )
        selected_mappings.append(retained)
        for candidate in candidates:
            if candidate is retained:
                continue
            duplicate_display_exclusions.append({
                "assigned_dual_positive_nucleus_id": nucleus_id,
                "excluded_npy_roi_id": int(candidate["npy_roi_id"]),
                "excluded_native_overlap_pixels": int(candidate["best_overlap_pixels"]),
                "retained_npy_roi_id": int(retained["npy_roi_id"]),
                "retained_native_overlap_pixels": int(retained["best_overlap_pixels"]),
                "reason": (
                    "display-only one-ROI-per-nucleus rule retains greatest native overlap "
                    "pixel count; lowest NPY ROI ID breaks a tie"
                ),
            })
    selected_mappings.sort(key=lambda record: int(record["npy_roi_id"]))

    displayed_npy = _warp_labels_to_display(labels_by_marker["npy"], spec)
    display_npy_props = {int(prop.label): prop for prop in regionprops(displayed_npy)}
    inset_npy = displayed_npy[y0:y1, x0:x1]
    layers: list[dict[str, object]] = []
    outlined_receipts: list[dict[str, object]] = []
    for mapping in selected_mappings:
        roi_id = int(mapping["npy_roi_id"])
        full_prop = display_npy_props.get(roi_id)
        if full_prop is None:
            raise RuntimeError(f"{sample} assigned NPY ROI {roi_id} vanished after nearest transform")
        inset_mask = np.asarray(inset_npy == roi_id, dtype=bool)
        if not np.any(inset_mask):
            raise RuntimeError(f"{sample} assigned NPY ROI {roi_id} does not intersect its inset")
        full_mask = np.asarray(displayed_npy == roi_id, dtype=bool)
        full_paths = find_contours(
            np.pad(full_mask.astype(np.uint8), 1, mode="constant"),
            0.5, fully_connected="high", positive_orientation="low",
        )
        if not full_paths:
            raise RuntimeError(f"{sample} assigned NPY ROI {roi_id} has no full-frame boundary")
        visible_paths: list[np.ndarray] = []
        for path in full_paths:
            translated = np.asarray(path, dtype=float).copy()
            translated[:, 0] = translated[:, 0] - 1.0 - y0
            translated[:, 1] = translated[:, 1] - 1.0 - x0
            min_y, min_x = np.min(translated, axis=0)
            max_y, max_x = np.max(translated, axis=0)
            if not (max_x < -0.5 or min_x > box_width - 0.5 or max_y < -0.5 or min_y > box_height - 0.5):
                visible_paths.append(translated)
        if not visible_paths:
            raise RuntimeError(f"{sample} assigned NPY ROI {roi_id} has no contour intersecting its inset")
        touches_edge = bool(
            np.any(inset_mask[0, :]) or np.any(inset_mask[-1, :])
            or np.any(inset_mask[:, 0]) or np.any(inset_mask[:, -1])
        )
        display_coordinates = np.ascontiguousarray(full_prop.coords, dtype=np.int32)
        receipt = {
            "npy_roi_id": roi_id,
            "assigned_dual_positive_nucleus_id": int(mapping["assigned_nucleus_id"]),
            "native_roi_area_px": int(mapping["native_roi_area_px"]),
            "native_roi_centroid_xy": mapping["native_roi_centroid_xy"],
            "native_roi_coordinate_sha256": mapping["native_roi_coordinate_sha256"],
            "assignment_best_overlap_pixels": int(mapping["best_overlap_pixels"]),
            "assignment_total_dapi_overlap_pixels": int(mapping["total_dapi_overlap_pixels"]),
            "assignment_best_fraction_of_dapi_overlap": float(mapping["best_fraction_of_dapi_overlap"]),
            "assignment_distinct_overlapping_dapi_nuclei": int(mapping["distinct_overlapping_dapi_nuclei"]),
            "display_full_roi_area_px": int(full_prop.area),
            "display_full_roi_coordinate_sha256": hashlib.sha256(display_coordinates.tobytes()).hexdigest(),
            "display_full_frame_contour_path_count": len(full_paths),
            "display_visible_contour_path_count": len(visible_paths),
            "display_inset_roi_pixels": int(inset_mask.sum()),
            "display_inset_fraction_of_full_roi": int(inset_mask.sum()) / float(full_prop.area),
            "display_inset_binary_sha256": _binary_mask_sha256(inset_mask),
            "touches_inset_edge": touches_edge,
            "crop_edge_closure_added": False,
        }
        layers.append({
            "npy_roi_id": roi_id,
            "full_frame_contour_paths_in_inset_coordinates": visible_paths,
        })
        outlined_receipts.append(receipt)

    outlined_nucleus_ids = sorted({
        int(record["assigned_dual_positive_nucleus_id"]) for record in outlined_receipts
    })
    expected_rois = int(expected["inset_assigned_npy_roi_contours"])
    expected_nuclei = int(expected["inset_dual_nuclei_with_contours"])
    if len(layers) != expected_rois or len(outlined_nucleus_ids) != expected_nuclei:
        raise RuntimeError(
            f"{sample} outlined NPY ROI/nucleus counts changed: "
            f"{len(layers)}/{len(outlined_nucleus_ids)} != {expected_rois}/{expected_nuclei}"
        )

    multiple_mapping_candidates = [
        {
            "dual_positive_nucleus_id": nucleus_id,
            "candidate_npy_roi_ids": sorted(
                int(record["npy_roi_id"]) for record in candidates
            ),
            "candidate_roi_count": len(candidates),
        }
        for nucleus_id, candidates in sorted(candidates_by_nucleus.items())
        if len(candidates) > 1
    ]
    uncovered_ids = sorted(inset_dual_set.difference(outlined_nucleus_ids))
    uncovered: list[dict[str, object]] = []
    for nucleus_id in uncovered_ids:
        candidates = []
        for mapping in mappings:
            overlap_by_nucleus = {
                int(item["nucleus_id"]): int(item["overlap_pixels"])
                for item in mapping["dapi_overlaps"]
            }
            if nucleus_id not in overlap_by_nucleus:
                continue
            candidates.append({
                "npy_roi_id": int(mapping["npy_roi_id"]),
                "overlap_pixels_with_uncovered_nucleus": overlap_by_nucleus[nucleus_id],
                "mapping_status": mapping["status"],
                "assigned_nucleus_id": mapping["assigned_nucleus_id"],
                "best_nucleus_id_if_unique": mapping["best_nucleus_id_if_unique"],
                "best_fraction_of_dapi_overlap": mapping["best_fraction_of_dapi_overlap"],
            })
        uncovered.append({
            "dual_positive_nucleus_id": nucleus_id,
            "reason": "no NPY ROI passes strict-majority assignment to this nucleus",
            "overlapping_npy_roi_candidates": candidates,
        })
    relevant_ambiguities = [
        record for record in mappings
        if (
            record["status"] != "assigned_strict_majority"
            and any(int(item["nucleus_id"]) in inset_dual_set for item in record["dapi_overlaps"])
        )
    ]

    mapping_bytes = json.dumps(mappings, sort_keys=True, separators=(",", ":")).encode("utf-8")
    outlined_bytes = json.dumps(
        outlined_receipts, sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    receipt = {
        "schema_version": "figure3-panel-a-npy-roi-contours-v1",
        "sample": sample,
        "classification_coordinates": "complete native Cellpose label fields before display transform",
        "analysis_counts": observed,
        "npy_roi_assignment_rule": (
            "use all native NPY ROI pixels overlapping native DAPI labels; assign only when one DAPI "
            "nucleus contains a strict >50% majority of those DAPI-overlapping ROI pixels"
        ),
        "ambiguity_policy": (
            "fail closed: exclude no-DAPI-overlap ROIs, tied best overlaps, and unique-best ROIs "
            "without a strict majority; the full strict mapping remains recorded before display selection"
        ),
        "display_one_roi_per_assigned_nucleus_rule": (
            "among strict-majority inset candidates assigned to the same nucleus, retain the ROI with "
            "the largest native overlap-pixel count; lowest NPY ROI ID breaks a tie"
        ),
        "intensity_used_for_roi_assignment_or_display_selection": False,
        "classification_rules": {
            "cfos_minimum_roi_area_px": CFOS_MIN_ROI_AREA_PX,
            "cfos_minimum_best_roi_fraction_overlapping_one_dapi_nucleus": CFOS_ROI_TO_NUCLEUS_OVERLAP,
            "npy_positive_nucleus_minimum_marker_pixels": NPY_MIN_MARKER_PIXELS,
            "npy_positive_nucleus_minimum_marker_fraction": NPY_MIN_NUCLEUS_FRACTION,
            "npy_roi_to_dapi_strict_majority_threshold": NPY_ROI_ASSIGNMENT_STRICT_MAJORITY,
            "strict_inequality": True,
        },
        "whole_field_dual_positive_nucleus_count": len(dual_positive),
        "whole_field_dual_positive_nucleus_ids": sorted(int(value) for value in dual_positive),
        "zoom_box_xywh_half_open": [x0, y0, box_width, box_height],
        "inset_dual_positive_nucleus_rule": (
            "nearest-neighbor transformed DAPI-label centroid inside the half-open zoom box"
        ),
        "inset_dual_positive_nucleus_count": len(inset_dual_ids),
        "inset_dual_positive_nucleus_ids": inset_dual_ids,
        "native_npy_roi_mapping_status_counts": mapping_status_counts,
        "native_npy_roi_mapping_count": len(mappings),
        "native_npy_roi_mapping_receipt_sha256": hashlib.sha256(mapping_bytes).hexdigest(),
        "native_npy_roi_mappings": mappings,
        "strict_majority_inset_candidate_npy_roi_count": len(strict_inset_mappings),
        "strict_majority_inset_candidate_npy_roi_ids": [
            int(record["npy_roi_id"]) for record in strict_inset_mappings
        ],
        "display_duplicate_assignment_exclusions": duplicate_display_exclusions,
        "outlined_npy_roi_count": len(outlined_receipts),
        "outlined_npy_roi_ids": [int(record["npy_roi_id"]) for record in outlined_receipts],
        "outlined_dual_positive_nucleus_count": len(outlined_nucleus_ids),
        "outlined_dual_positive_nucleus_ids": outlined_nucleus_ids,
        "outlined_roi_receipt_sha256": hashlib.sha256(outlined_bytes).hexdigest(),
        "outlined_rois": outlined_receipts,
        "multiple_strict_majority_npy_roi_candidates_for_one_dual_positive_nucleus": multiple_mapping_candidates,
        "inset_dual_positive_nuclei_without_assigned_roi": uncovered,
        "relevant_fail_closed_ambiguity_exclusions": relevant_ambiguities,
        "contour_render": {
            "geometry": (
                "every connected boundary is derived on the complete nearest-neighbor registered "
                "display ROI-ID field before any inset operation"
            ),
            "crop_policy": (
                "full-frame paths are translated into inset coordinates; the rendered stroke is "
                "clipped to the exact microscopy-image extent before axes letterboxing and then by "
                "the axes, without adding an artificial closure at an inset edge"
            ),
            "full_native_roi_identity_preserved": True,
            "stroke_color_name": "black",
            "stroke_color_hex": NPY_ROI_CONTOUR_COLOR,
            "linewidth_pt": NPY_ROI_CONTOUR_WIDTH_PT,
            "white_backing_applied": False,
            "stroke_count_per_path": 1,
            "image_extent_clip_xyxy_data": [
                -0.5, -0.5, box_width - 0.5, box_height - 0.5,
            ],
            "black_letterbox_stroke_allowed": False,
            "path_geometry_modified_for_clipping": False,
        },
        "mask_sources": source_receipts,
        "display_to_raw_homogeneous_matrix": spec["display_to_raw"],
        "display_output_shape_yx": list(spec["output_shape_yx"]),
        "label_interpolation": "nearest neighbor only",
    }
    return layers, receipt, source_paths


def _draw_npy_roi_contours(
    ax: plt.Axes,
    layers: Sequence[Mapping[str, object]],
    image_shape_yx: Sequence[int],
) -> tuple[int, int]:
    image_height, image_width = (int(value) for value in image_shape_yx)
    if image_height <= 0 or image_width <= 0:
        raise ValueError(f"Invalid ROI contour image shape: {tuple(image_shape_yx)}")
    # Axes clipping alone includes the black letterbox surrounding an equal-
    # sized microscopy tile. Clip each line to the actual image rectangle so
    # an edge-touching full-frame ROI stops where microscopy pixels stop.
    image_clip = Rectangle(
        (-0.5, -0.5), image_width, image_height,
        transform=ax.transData, facecolor="none", edgecolor="none",
    )
    image_xlim = (-0.5, image_width - 0.5)
    image_ylim = (image_height - 0.5, -0.5)
    if not (
        np.allclose(ax.get_xlim(), image_xlim, atol=1e-9, rtol=0.0)
        and np.allclose(ax.get_ylim(), image_ylim, atol=1e-9, rtol=0.0)
    ):
        raise RuntimeError("ROI contours must be drawn immediately after the microscopy image")
    # Each line disables x/y scaling locally so complete paths cannot change
    # the image-set limits or the surrounding letterbox geometry.
    drawn_rois = 0
    drawn_paths = 0
    for layer in layers:
        paths = layer["full_frame_contour_paths_in_inset_coordinates"]
        if not paths:
            raise RuntimeError(
                "NPY ROI {} has no drawable inset boundary".format(layer["npy_roi_id"])
            )
        for path in paths:
            path = np.asarray(path, dtype=float)
            x = path[:, 1]
            y = path[:, 0]
            line, = ax.plot(
                x, y, color=NPY_ROI_CONTOUR_COLOR,
                linewidth=NPY_ROI_CONTOUR_WIDTH_PT,
                solid_capstyle="round", solid_joinstyle="round",
                antialiased=True, zorder=15, clip_on=True,
                scalex=False, scaley=False,
            )
            line.set_clip_path(image_clip)
            drawn_paths += 1
        drawn_rois += 1
    if not (
        np.allclose(ax.get_xlim(), image_xlim, atol=1e-9, rtol=0.0)
        and np.allclose(ax.get_ylim(), image_ylim, atol=1e-9, rtol=0.0)
    ):
        raise RuntimeError("Complete ROI paths changed the microscopy image limits")
    return drawn_rois, drawn_paths


def _add_scale_bar(
    ax: plt.Axes,
    shape_yx: tuple[int, int],
    um_per_px: float,
    *,
    length_um: float,
) -> None:
    height, width = shape_yx
    bar_px = float(length_um) / float(um_per_px)
    if not (0 < bar_px < 0.45 * width):
        raise ValueError(f"Invalid scale bar: {length_um:g} µm -> {bar_px:.1f} px for width {width}")
    x1 = width * 0.955
    x0 = x1 - bar_px
    y = height * 0.925
    ax.plot([x0, x1], [y, y], color="white", linewidth=3.0, solid_capstyle="butt", zorder=12)
    ax.text(
        (x0 + x1) / 2, y - height * 0.028, f"{length_um:g} µm",
        color="white", ha="center", va="bottom", fontsize=6.8, fontweight="bold", zorder=12,
    )


def build_raw_microscopy_panel(
    raw_root: Path,
    png_path: Path,
    pdf_path: Path,
    provenance_path: Path,
    *,
    dpi: int = 300,
) -> tuple[np.ndarray, list[Path]]:
    """Render Panel A from raw NPY and c-FOS channel TIFFs."""
    raw_root = raw_root.expanduser().resolve()
    records: list[dict[str, object]] = []
    rendered: list[dict[str, object]] = []
    source_paths: list[Path] = []
    for spec in REPRESENTATIVE_FIELDS:
        signals: dict[str, np.ndarray] = {}
        channel_ranges: dict[str, dict[str, float]] = {}
        channel_paths: dict[str, Path] = {}
        for marker_name, relative_path in spec["paths"].items():
            path = require_file(raw_root / str(relative_path), f"{spec['sample']} {marker_name} raw TIFF")
            source_paths.append(path)
            channel_paths[marker_name] = path
            raw_source = np.asarray(tifffile.imread(path))
            raw_plane = _marker_plane(path)
            display_plane = _warp_to_display(raw_plane, spec)
            signals[marker_name], low, high = _display_stretch(display_plane)
            channel_ranges[marker_name] = {"percentile_0_5": low, "percentile_99_8": high}
            records.append({
                "condition": spec["condition"],
                "sample": spec['sample'],
                "marker": marker_name,
                "path": str(path),
                "sha256": sha256_file(path),
                "source_shape_yx": list(raw_plane.shape),
                "source_dtype": str(raw_source.dtype),
                "display_percentile_0_5": low,
                "display_percentile_99_8": high,
            })
        colors = {name: _colorize(signals[name], MARKER_COLORS[name]) for name in ("npy", "cfos")}
        merge = np.clip(colors["npy"] + colors["cfos"], 0.0, 1.0)
        base_zoom_box = _choose_zoom_box(signals["npy"], signals["cfos"])
        zoom_box = base_zoom_box
        if bool(spec.get("mirror_zoom_horizontally_for_homologous_side", False)):
            base_x, base_y, base_width, base_height = base_zoom_box
            display_width = int(signals["npy"].shape[1])
            zoom_box = (display_width - (base_x + base_width), base_y, base_width, base_height)
        x0, y0, box_width, box_height = zoom_box
        roi_layers, roi_receipt, roi_sources = _double_positive_npy_roi_layers(
            raw_root, spec, zoom_box,
        )
        source_paths.extend(roi_sources)
        roi_receipt["zoom_side_selection"] = {
            "base_zoom_box_xywh_half_open": list(base_zoom_box),
            "final_zoom_box_xywh_half_open": list(zoom_box),
            "water_homologous_side_rule": (
                "Water only: exact horizontal mirror x = display_width - (base_x + width); "
                "y, width, and height unchanged to match the viewer-left Sucrose/Allulose side"
                if zoom_box != base_zoom_box else "not applied; base deterministic zoom retained"
            ),
            "display_intensity_used_for_roi_assignment": False,
        }
        rendered.append({
            "spec": spec,
            "signals": signals,
            "colors": colors,
            "merge": merge,
            "base_zoom_box": base_zoom_box,
            "zoom_box": zoom_box,
            "zoom_merge": merge[y0 : y0 + box_height, x0 : x0 + box_width],
            "zoom_npy": signals["npy"][y0 : y0 + box_height, x0 : x0 + box_width],
            "zoom_cfos": signals["cfos"][y0 : y0 + box_height, x0 : x0 + box_width],
            "npy_roi_layers": roi_layers,
            "npy_roi_receipt": roi_receipt,
            "channel_ranges": channel_ranges,
            "channel_paths": channel_paths,
        })

    outlined_roi_counts = {
        str(item["spec"]["condition"]): int(item["npy_roi_receipt"]["outlined_npy_roi_count"])
        for item in rendered
    }
    outlined_nucleus_counts = {
        str(item["spec"]["condition"]): int(
            item["npy_roi_receipt"]["outlined_dual_positive_nucleus_count"]
        )
        for item in rendered
    }
    inset_dual_counts = {
        str(item["spec"]["condition"]): int(
            item["npy_roi_receipt"]["inset_dual_positive_nucleus_count"]
        )
        for item in rendered
    }
    if inset_dual_counts != {"Water": 24, "Sucrose": 7, "Allulose": 38}:
        raise RuntimeError(f"Unexpected inset dual-positive nucleus counts: {inset_dual_counts}")
    if outlined_roi_counts != {"Water": 23, "Sucrose": 6, "Allulose": 36}:
        raise RuntimeError(f"Unexpected outlined NPY ROI counts: {outlined_roi_counts}")
    if outlined_nucleus_counts != {"Water": 23, "Sucrose": 6, "Allulose": 36}:
        raise RuntimeError(f"Unexpected outlined dual-positive nucleus counts: {outlined_nucleus_counts}")

    fig = plt.figure(figsize=(15.9, 6.35), facecolor="white")
    grid = fig.add_gridspec(
        3, 4, width_ratios=(1.0, 1.0, 1.0, 1.0),
        left=0.050, right=0.998, bottom=0.035, top=0.925,
        hspace=0.055, wspace=0.010,
    )
    headers = (
        ("NPY$^{GFP}$", "#00a000"),
        ("c-FOS", "#d00000"),
        ("Merge", "black"),
        ("Magnified merge", "black"),
    )
    axes: list[list[plt.Axes]] = []
    for row_index, item in enumerate(rendered):
        spec = item["spec"]
        row_axes: list[plt.Axes] = []
        for column_index in range(4):
            axis = fig.add_subplot(grid[row_index, column_index])
            axis.set_facecolor("black")
            axis.set_box_aspect(PANEL_A_TILE_BOX_ASPECT)
            axis.add_patch(Rectangle(
                (0.0, 0.0), 1.0, 1.0, transform=axis.transAxes,
                facecolor="black", edgecolor="none", linewidth=0.0,
                zorder=-100, clip_on=False,
            ))
            axis.set_axis_off()
            if row_index == 0:
                axis.set_title(
                    headers[column_index][0], color=headers[column_index][1],
                    fontsize=12.0, fontweight="bold", pad=4,
                )
            row_axes.append(axis)
        axes.append(row_axes)
        fig.text(
            0.008, 0.777 - row_index * 0.297,
            f"{spec['condition']}\n{spec['sample'].replace('_', ' ')}",
            ha="left", va="center", fontsize=10.5, fontweight="bold",
        )
        for column_index, marker_name in enumerate(("npy", "cfos")):
            row_axes[column_index].imshow(item["colors"][marker_name], interpolation="nearest")
            _add_scale_bar(
                row_axes[column_index], item["signals"][marker_name].shape,
                float(spec["um_per_px"]), length_um=100.0,
            )
        row_axes[2].imshow(item["merge"], interpolation="nearest")
        _add_scale_bar(
            row_axes[2], item["signals"]["npy"].shape,
            float(spec["um_per_px"]), length_um=100.0,
        )
        x0, y0, box_width, box_height = item["zoom_box"]
        row_axes[2].add_patch(Rectangle(
            (x0, y0), box_width, box_height, fill=False,
            edgecolor="white", linewidth=1.8, linestyle=(0, (3, 2)), zorder=13,
        ))
        row_axes[3].imshow(item["zoom_merge"], interpolation="nearest")
        _add_scale_bar(
            row_axes[3], item["zoom_npy"].shape,
            float(spec["um_per_px"]), length_um=50.0,
        )
        drawn_roi_count, drawn_path_count = _draw_npy_roi_contours(
            row_axes[3], item["npy_roi_layers"], item["zoom_merge"].shape[:2],
        )
        if drawn_roi_count != int(item["npy_roi_receipt"]["outlined_npy_roi_count"]):
            raise RuntimeError(
                "Panel A NPY ROI contour count changed for {}".format(spec["sample"])
            )
        if drawn_path_count < drawn_roi_count:
            raise RuntimeError(
                "Panel A lost an NPY ROI contour path for {}".format(spec["sample"])
            )
        item["npy_roi_receipt"]["drawn_contour_path_count"] = drawn_path_count
        # Keep the lower connector above the bottom-right scale-bar region.
        # The source remains the zoom rectangle; only its landing point is raised.
        for source_y, target_y in ((y0, 1.0), (y0 + box_height, 0.25)):
            fig.add_artist(ConnectionPatch(
                xyA=(x0 + box_width, source_y), coordsA="data", axesA=row_axes[2],
                xyB=(0.0, target_y), coordsB="axes fraction", axesB=row_axes[3],
                color="white", linewidth=1.2, linestyle=(0, (3, 2)), zorder=12,
            ))

    fig.canvas.draw()
    axis_boxes = [axis.get_position().bounds for row in axes for axis in row]
    tile_widths = np.asarray([box[2] for box in axis_boxes], dtype=float)
    tile_heights = np.asarray([box[3] for box in axis_boxes], dtype=float)
    if not (
        np.allclose(tile_widths, tile_widths[0], atol=1e-10, rtol=0.0)
        and np.allclose(tile_heights, tile_heights[0], atol=1e-10, rtol=0.0)
    ):
        raise RuntimeError(
            "Panel A microscopy axes are not exactly even-sized: "
            f"width range={tile_widths.min():.12g}-{tile_widths.max():.12g}, "
            f"height range={tile_heights.min():.12g}-{tile_heights.max():.12g}"
        )
    tile_geometry = {
        "grid_rows_columns": [3, 4],
        "column_width_ratios": [1.0, 1.0, 1.0, 1.0],
        "tile_box_aspect_height_over_width": PANEL_A_TILE_BOX_ASPECT,
        "tile_width_figure_fraction": float(tile_widths[0]),
        "tile_height_figure_fraction": float(tile_heights[0]),
        "all_twelve_axis_boxes_equal": True,
        "grid_bounds_normalized": {
            "left": 0.050, "right": 0.998, "bottom": 0.035, "top": 0.925,
        },
        "inter_column_wspace": 0.010,
        "image_aspect_policy": (
            "native pixel aspect preserved; black letterboxing within equal outer tiles where needed"
        ),
    }

    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, format="png", dpi=dpi, facecolor="white")
    fig.savefig(pdf_path, format="pdf", dpi=PDF_RASTER_DPI, facecolor="white")
    plt.close(fig)
    panel = read_rgb(png_path)

    provenance = {
        "schema_version": "figure3-panel-a-raw-microscopy-v5",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "panel_png": str(png_path.resolve()),
        "panel_pdf": str(pdf_path.resolve()),
        "source_policy": (
            "raw NPY and c-FOS TIFFs provide image pixels; native DAPI/c-FOS/NPY "
            "Cellpose masks provide analysis-defined nucleus calls and native NPY ROI contours; "
            "no PowerPoint or screenshot pixels"
        ),
        "display_policy": (
            "Per-field/per-channel 0.5th–99.8th percentile stretch followed by fixed gamma 0.85. "
            "Display intensities must not be compared across conditions."
        ),
        "colors": {"NPY": "green", "c-FOS": "red"},
        "zoom_policy": (
            "Fixed-size base source rectangle centered on the maximum Gaussian-smoothed NPY*c-FOS "
            "display co-intensity. Water uses its exact horizontal mirror, preserving y/width/height, "
            "to show the same viewer-left anatomical side as Sucrose and Allulose. View selection "
            "never defines a positive cell or NPY ROI contour."
        ),
        "npy_roi_contour_policy": (
            "Each native NPY ROI is assigned only when one native DAPI nucleus contains a strict >50% "
            "majority of its DAPI-overlapping pixels. All assigned ROIs whose nucleus is an analysis-defined "
            "c-FOS-positive/NPY-positive nucleus with transformed DAPI centroid inside the half-open inset "
            "are outlined. Ambiguous assignments fail closed; if several strict-majority ROIs map to one nucleus, "
            "the ROI with greatest native overlap is displayed, with lowest ROI ID as tie-break. Contours are "
            "derived on the complete nearest-neighbor registered display label frame, translated "
            "to inset coordinates, and clipped to the microscopy-image extent before axes letterboxing "
            "without adding a crop-edge closure. Display intensity is never used for assignment or selection. "
            "Each contour is rendered as one clean 1.56-pt black (#000000) stroke without a white "
            "backing or shadow."
        ),
        "npy_roi_contour_counts": {
            "inset_dual_positive_nuclei_by_condition": inset_dual_counts,
            "outlined_npy_rois_by_condition": outlined_roi_counts,
            "dual_positive_nuclei_with_outlined_roi_by_condition": outlined_nucleus_counts,
        },
        "tile_geometry": tile_geometry,
        "scale_bars": {
            "full_field_um": 100.0,
            "zoom_um": 50.0,
            "Water2_um_per_px": REPRESENTATIVE_FIELDS[0]["um_per_px"],
            "legacy_um_per_px": REPRESENTATIVE_FIELDS[1]["um_per_px"],
        },
        "representative_fields": [
            {
                "condition": item["spec"]["condition"],
                "sample": item["spec"]["sample"],
                "reference_match": item["spec"]["reference_match"],
                "output_shape_yx": list(item["spec"]["output_shape_yx"]),
                "display_to_raw": item["spec"]["display_to_raw"],
                "base_zoom_box_xywh_half_open": list(item["base_zoom_box"]),
                "zoom_box_xywh_half_open": list(item["zoom_box"]),
                "display_ranges": item["channel_ranges"],
                "majority_assigned_double_positive_npy_roi_contours": item["npy_roi_receipt"],
            }
            for item in rendered
        ],
        "sources": records,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return panel, source_paths + [provenance_path, pdf_path]


def read_rgb(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        return np.asarray(image.convert("RGB"))


def add_panel_letter(
    ax: plt.Axes,
    letter: str,
    *,
    x: float = -0.025,
    y: float = 1.015,
    ha: str = "right",
) -> Artist:
    return ax.text(
        x, y, letter, transform=ax.transAxes, ha=ha, va="top",
        fontsize=22, fontweight="bold", color="black", clip_on=False,
    )


def export_panel_crops(
    fig: plt.Figure,
    panel_axes: Mapping[str, Sequence[plt.Axes]],
    panel_texts: Mapping[str, Sequence[Artist]],
    output_dir: Path,
    figure_name: str,
    dpi: int,
) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    all_axes = list(fig.axes)
    all_text = list(fig.texts)
    boxes: dict[str, Bbox] = {}
    for letter, axes in panel_axes.items():
        pixel_boxes = []
        for axis in axes:
            box = axis.get_tightbbox(renderer)
            if box is not None and np.isfinite(box.extents).all():
                pixel_boxes.append(box)
        for artist in panel_texts.get(letter, ()):
            box = artist.get_window_extent(renderer)
            if box is not None and np.isfinite(box.extents).all():
                pixel_boxes.append(box)
        if not pixel_boxes:
            raise RuntimeError(f"Panel {letter} has no drawable bounds")
        inch = Bbox.union(pixel_boxes).transformed(fig.dpi_scale_trans.inverted())
        pad_left = 0.24 if letter in "DEFG" else 0.10
        boxes[letter] = Bbox.from_extents(inch.x0 - pad_left, inch.y0 - 0.10, inch.x1 + 0.10, inch.y1 + 0.12)

    output_dir.mkdir(parents=True, exist_ok=True)
    for letter, axes in panel_axes.items():
        target_axes = set(axes)
        target_text = set(panel_texts.get(letter, ()))
        axis_visibility = [(axis, axis.get_visible()) for axis in all_axes]
        text_visibility = [(artist, artist.get_visible()) for artist in all_text]
        try:
            for axis in all_axes:
                axis.set_visible(axis in target_axes)
            for artist in all_text:
                artist.set_visible(artist in target_text)
            stem = output_dir / f"{figure_name}_Panel_{letter}"
            fig.savefig(stem.with_suffix(".pdf"), format="pdf", dpi=PDF_RASTER_DPI, bbox_inches=boxes[letter], facecolor="white")
            fig.savefig(stem.with_suffix(".png"), format="png", dpi=dpi, bbox_inches=boxes[letter], facecolor="white")
        finally:
            for axis, visible in axis_visibility:
                axis.set_visible(visible)
            for artist, visible in text_visibility:
                artist.set_visible(visible)


def validate_quantitative_tick_clearance(
    fig: plt.Figure,
    quantitative_axes: Mapping[str, plt.Axes],
    spatial_axes: Mapping[str, plt.Axes],
    *,
    minimum_inches: float = 0.04,
) -> None:
    """Fail if the spatial row can paint over the D/E x tick labels."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    minimum_pixels = minimum_inches * fig.dpi
    for quantitative_letter, spatial_letter in (("D", "F"), ("E", "G")):
        tick_boxes = [
            label.get_window_extent(renderer)
            for label in quantitative_axes[quantitative_letter].get_xticklabels()
            if label.get_visible() and label.get_text()
        ]
        if not tick_boxes:
            raise RuntimeError(f"Panel {quantitative_letter} has no visible x tick labels")
        tick_bottom = min(box.y0 for box in tick_boxes)
        spatial_top = spatial_axes[spatial_letter].get_window_extent(renderer).y1
        clearance = tick_bottom - spatial_top
        if clearance < minimum_pixels:
            raise RuntimeError(
                f"Panel {quantitative_letter} x tick labels have only "
                f"{clearance / fig.dpi:.3f} in clearance above Panel {spatial_letter}; "
                f"at least {minimum_inches:.3f} in is required"
            )


def validate_spatial_panel_geometry(
    fig: plt.Figure,
    spatial_axes: Mapping[str, plt.Axes],
    spatial_images: Mapping[str, np.ndarray],
    *,
    minimum_width_inches: float = 7.70,
) -> dict[str, dict[str, float]]:
    """Fail unless F/G use the wider allocation at their native pixel aspect."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    geometry: dict[str, dict[str, float]] = {}
    for letter in "FG":
        axis = spatial_axes[letter]
        image = spatial_images[letter]
        box = axis.get_window_extent(renderer)
        width_inches = float(box.width / fig.dpi)
        height_inches = float(box.height / fig.dpi)
        source_aspect = float(image.shape[1] / image.shape[0])
        rendered_aspect = float(width_inches / height_inches)
        position = axis.get_position()
        if width_inches < minimum_width_inches:
            raise RuntimeError(
                f"Panel {letter} is only {width_inches:.3f} in wide; "
                f"at least {minimum_width_inches:.3f} in is required"
            )
        if not np.isclose(rendered_aspect, source_aspect, atol=1e-6, rtol=0.0):
            raise RuntimeError(
                f"Panel {letter} aspect changed from {source_aspect:.9f} "
                f"to {rendered_aspect:.9f}"
            )
        if position.y0 < SPATIAL_PANEL_MIN_BOTTOM_FRACTION:
            raise RuntimeError(
                f"Panel {letter} bottom {position.y0:.6f} is below the protected "
                f"{SPATIAL_PANEL_MIN_BOTTOM_FRACTION:.6f} figure fraction"
            )
        geometry[letter] = {
            "width_inches": width_inches,
            "height_inches": height_inches,
            "source_pixel_aspect_width_over_height": source_aspect,
            "rendered_aspect_width_over_height": rendered_aspect,
            "bottom_figure_fraction": float(position.y0),
            "top_figure_fraction": float(position.y1),
        }
    if not np.isclose(
        geometry["F"]["width_inches"], geometry["G"]["width_inches"],
        atol=1e-9, rtol=0.0,
    ):
        raise RuntimeError("Panels F and G do not have exactly matched displayed widths")
    return geometry


def write_text_and_provenance(
    output_dir: Path,
    analysis_dir: Path,
    inputs: Mapping[str, Path],
    values: pd.DataFrame,
    selected_stats: pd.DataFrame,
    supplemental_stats: pd.DataFrame,
    figure_name: str,
    spatial_layout_geometry: Mapping[str, Mapping[str, float]],
) -> None:
    legends_dir = output_dir / "legends"
    source_dir = output_dir / "source_data"
    provenance_dir = output_dir / "provenance"
    for directory in (legends_dir, source_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=True)
    identity_columns = [column for column in ("animal", "condition", "sample") if column in values]
    displayed_columns = [endpoint for _, endpoint, _, _, _ in ENDPOINTS]
    animal_plot_values_path = source_dir / f"{figure_name}_animal_plot_values.csv"
    values[identity_columns + displayed_columns].to_csv(animal_plot_values_path, index=False)
    statistics_used = source_dir / f"{figure_name}_anova_mwu_statistics_used.csv"
    selected_stats.to_csv(statistics_used, index=False)
    supplemental_endpoint = SUPPLEMENTAL_ENDPOINTS[0][0]
    supplemental_values = values[identity_columns + [supplemental_endpoint]].copy()
    supplemental_values[f"{supplemental_endpoint}_percent"] = 100.0 * supplemental_values[supplemental_endpoint]
    supplemental_values_path = source_dir / f"{figure_name}_supplemental_NPY_DAPI_animal_values.csv"
    supplemental_statistics_path = source_dir / f"{figure_name}_supplemental_NPY_DAPI_anova_mwu_statistics.csv"
    supplemental_note_path = source_dir / f"{figure_name}_supplemental_NPY_DAPI_README.txt"
    supplemental_values.to_csv(supplemental_values_path, index=False)
    supplemental_stats.to_csv(supplemental_statistics_path, index=False)
    supplemental_note_path.write_text(
        "Supplemental Figure 3 endpoint retained after removal from the main plate\n"
        "Endpoint: NPY-positive nuclei / DAPI nuclei (npy_over_dapi).\n"
        "The values are animal-level ratios; the percent column is 100 times the ratio.\n"
        "Statistics are the original ordinary one-way ANOVA plus raw exact two-sided MWU receipts.\n"
        "This endpoint is intentionally not displayed in panels D-G and has not been discarded.\n",
        encoding="utf-8",
    )
    for stale_path in (
        source_dir / f"{figure_name}_exact_statistics_used.csv",
    ):
        if stale_path.is_file():
            stale_path.unlink()

    endpoint_counts = {
        endpoint: {
            condition: int(np.isfinite(pd.to_numeric(
                values.loc[values["condition"].eq(condition), endpoint], errors="coerce"
            )).sum())
            for condition in CONDITIONS
        }
        for _, endpoint, _, _, _ in ENDPOINTS
    }
    supplemental_counts = {
        condition: int(np.isfinite(pd.to_numeric(
            values.loc[values["condition"].eq(condition), "npy_over_dapi"], errors="coerce"
        )).sum())
        for condition in CONDITIONS
    }
    spatial_stats = pd.read_csv(inputs["spatial_exact_permanova"])
    spatial_required = {
        "endpoint", "test", "n_animals", "feature_count", "feature_transform",
        "p_value_exact", "enumerated_labelings", "experimental_unit",
        "cells_as_independent_replicates", "bh_q_value_two_endpoint_family",
    }
    if not spatial_required.issubset(spatial_stats.columns):
        raise ValueError(
            "Spatial statistics table is missing required columns: "
            + str(sorted(spatial_required.difference(spatial_stats.columns)))
        )
    if set(spatial_stats["endpoint"].astype(str)) != {"cfos_occurrence", "double_occurrence"}:
        raise ValueError("Spatial statistics must contain exactly the declared F/G occurrence endpoints")
    spatial_stats = spatial_stats.set_index("endpoint", drop=False)
    for endpoint in ("cfos_occurrence", "double_occurrence"):
        row = spatial_stats.loc[endpoint]
        independent = str(row["cells_as_independent_replicates"]).strip().lower()
        if (
            int(row["n_animals"]) != 9 or int(row["feature_count"]) != 6
            or int(row["enumerated_labelings"]) != 1680
            or str(row["experimental_unit"]).strip().lower() != "biological animal"
            or independent not in {"false", "0"}
            or "permanova" not in str(row["test"]).lower()
        ):
            raise ValueError(f"Spatial endpoint {endpoint} violates the animal-level six-shell exact-test contract")
    dispersion_stats = pd.read_csv(inputs["spatial_exact_dispersion"])
    dispersion_required = {
        "endpoint", "test", "statistic", "dispersion_f", "p_value_exact",
        "n_animals", "feature_count", "enumerated_labelings", "experimental_unit",
        "cells_as_independent_replicates", "bh_q_value_two_endpoint_diagnostic_family",
    }
    if not dispersion_required.issubset(dispersion_stats.columns):
        raise ValueError(
            "Spatial dispersion table is missing required columns: "
            + str(sorted(dispersion_required.difference(dispersion_stats.columns)))
        )
    if set(dispersion_stats["endpoint"].astype(str)) != {"cfos_occurrence", "double_occurrence"}:
        raise ValueError("Spatial dispersion table must contain exactly the declared F/G endpoints")
    dispersion_stats = dispersion_stats.set_index("endpoint", drop=False)
    for endpoint in ("cfos_occurrence", "double_occurrence"):
        row = dispersion_stats.loc[endpoint]
        independent = str(row["cells_as_independent_replicates"]).strip().lower()
        if (
            int(row["n_animals"]) != 9 or int(row["feature_count"]) != 6
            or int(row["enumerated_labelings"]) != 1680
            or str(row["experimental_unit"]).strip().lower() != "biological animal"
            or independent not in {"false", "0"}
            or "dispersion" not in str(row["test"]).lower()
            or "anova" not in str(row["statistic"]).lower()
        ):
            raise ValueError(f"Spatial endpoint {endpoint} violates the exact dispersion-diagnostic contract")
    spatial_result_text = (
        "Exact 1,680-label animal PERMANOVA: panel F p = "
        f"{fmt_p(float(spatial_stats.loc['cfos_occurrence', 'p_value_exact']))}, BH q = "
        f"{fmt_p(float(spatial_stats.loc['cfos_occurrence', 'bh_q_value_two_endpoint_family']))}; "
        "panel G p = "
        f"{fmt_p(float(spatial_stats.loc['double_occurrence', 'p_value_exact']))}, BH q = "
        f"{fmt_p(float(spatial_stats.loc['double_occurrence', 'bh_q_value_two_endpoint_family']))}. "
        "Exact 1,680-label multivariate-dispersion diagnostic: panel F F = "
        f"{float(dispersion_stats.loc['cfos_occurrence', 'dispersion_f']):.2f}, p = "
        f"{fmt_p(float(dispersion_stats.loc['cfos_occurrence', 'p_value_exact']))}, BH q = "
        f"{fmt_p(float(dispersion_stats.loc['cfos_occurrence', 'bh_q_value_two_endpoint_diagnostic_family']))}; "
        "panel G F = "
        f"{float(dispersion_stats.loc['double_occurrence', 'dispersion_f']):.2f}, p = "
        f"{fmt_p(float(dispersion_stats.loc['double_occurrence', 'p_value_exact']))}, BH q = "
        f"{fmt_p(float(dispersion_stats.loc['double_occurrence', 'bh_q_value_two_endpoint_diagnostic_family']))}. "
        "The panel G dispersion result is heterogeneous, so its PERMANOVA cannot be interpreted as pure "
        "group-centroid separation; it may reflect centroid and/or dispersion differences."
    )
    spatial_source_copies = {
        "spatial_exact_permanova": source_dir / f"{figure_name}_spatial_exact_permanova.csv",
        "spatial_exact_permutation_receipts": source_dir / f"{figure_name}_spatial_exact_permutation_receipts.csv",
        "spatial_exact_dispersion": source_dir / f"{figure_name}_spatial_exact_dispersion.csv",
        "spatial_exact_dispersion_permutation_receipts": source_dir / f"{figure_name}_spatial_exact_dispersion_permutation_receipts.csv",
        "spatial_animal_features": source_dir / f"{figure_name}_spatial_animal_radial_occurrence_features.csv",
        "spatial_animal_counts": source_dir / f"{figure_name}_spatial_animal_counts.csv",
        "spatial_model_setting_audit": source_dir / f"{figure_name}_spatial_segmentation_model_setting_audit.csv",
    }
    for role, target in spatial_source_copies.items():
        shutil.copyfile(inputs[role], target)
    panel_a_sidecar = json.loads(
        inputs["panel_A_raw_reconstruction_provenance"].read_text(encoding="utf-8")
    )
    cartoon_sidecar = json.loads(inputs["cartoon_provenance"].read_text(encoding="utf-8"))
    if int(cartoon_sidecar.get("schema_version", -1)) != 9:
        raise RuntimeError("Matched Water/Allulose cartoon provenance schema must be v9")
    expected_cartoon_map = {
        "B": "Water marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
        "C": "Allulose marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
    }
    if cartoon_sidecar.get("panel_map") != expected_cartoon_map:
        raise RuntimeError(f"Unexpected matched-cartoon panel map: {cartoon_sidecar.get('panel_map')}")
    grammar = cartoon_sidecar.get("matched_panel_grammar", {})
    expected_view_titles = ["NPY-positive nuclei", "c-FOS-positive nuclei", "Combined status"]
    if (
        grammar.get("views_in_order") != [
            "NPY-positive DAPI nuclei", "c-FOS-positive DAPI nuclei", "combined/dual status",
        ]
        or grammar.get("view_title_text") != expected_view_titles
        or grammar.get("same_physical_canvas") is not True
        or grammar.get("same_colors") is not True
        or grammar.get("same_title_and_scale_bar_grammar") is not True
        or grammar.get("in_panel_method_footnotes") is not False
        or grammar.get("in_panel_status_counts") is not False
        or grammar.get("count_location") != "Figure 3 caption, panel legends, and provenance only"
        or grammar.get("panel_title_x_figure_fraction") != 0.5
        or grammar.get("panel_title_horizontal_alignment") != "center"
    ):
        raise RuntimeError(f"Matched-cartoon grammar changed: {grammar}")
    expected_source_lumen_hashes = {
        "Water": "1e462a1232813f3ed4b36f6893d12b114473f1ff3cf7c089bfccf363e2e70063",
        "Allulose": "c4521d0f26dc6c1eb160343ef6f91a1516d16c53f1e65f6b8336e8977490b973",
    }
    expected_display_counts = {
        "Water": {"npy_positive_nuclei": 205, "cfos_positive_nuclei": 113, "dual_positive_nuclei": 57},
        "Allulose": {"npy_positive_nuclei": 493, "cfos_positive_nuclei": 271, "dual_positive_nuclei": 131},
    }
    cartoon_samples = cartoon_sidecar.get("samples", {})
    if set(cartoon_samples) != set(expected_source_lumen_hashes):
        raise RuntimeError(f"Matched cartoons must contain Water and Allulose; got {sorted(cartoon_samples)}")
    display_counts = {}
    accepted_lumen_hashes = {}
    accepted_lumen_hil = {}
    for condition, expected_source_hash in expected_source_lumen_hashes.items():
        receipt = cartoon_samples[condition]
        transform = receipt.get("display_transform", {})
        lumen = receipt.get("ventricular_lumen_outline", {})
        if transform.get("output_shape_yx") != [866, 1374]:
            raise RuntimeError(f"{condition} cartoon common-frame shape changed")
        if not np.allclose(receipt.get("physical_size_xy_um_per_px", []), [0.664212669961432] * 2):
            raise RuntimeError(f"{condition} cartoon physical scale changed")
        if transform.get("label_interpolation") != "nearest neighbor only":
            raise RuntimeError(f"{condition} cartoon label registration is not nearest neighbor")
        source_lumen = lumen.get("source_final_after_completion_and_hole_fill", {})
        if source_lumen.get("binary_sha256") != expected_source_hash:
            raise RuntimeError(f"{condition} raw-DAPI ventricular source receipt changed")
        hil = lumen.get("human_in_loop_annotation", {})
        hil_geometry = hil.get("mask_geometry", {})
        selected_hash = lumen.get("selected_binary_sha256")
        if (
            hil.get("method")
            != "human-in-the-loop closed polygon in the exact registered common frame"
            or hil.get("analysis_role")
            != "display annotation only; never an analysis ROI or denominator"
            or not str(hil.get("reviewer", "")).strip()
            or not str(hil.get("session", "")).strip()
            or not str(hil.get("accepted_at_utc", "")).strip()
            or int(hil.get("polygon_vertex_count", 0)) < 3
            or len(str(hil.get("receipt_sha256", ""))) != 64
            or len(str(hil.get("mask_file_sha256", ""))) != 64
            or hil.get("source_binding", {}).get("raw_dapi_lumen_source_binary_sha256")
            != expected_source_hash
            or hil_geometry.get("binary_sha256") != selected_hash
            or hil_geometry.get("connected_components_8") != 1
            or hil_geometry.get("contour_path_counts") != {"external": 1, "tree": 1}
        ):
            raise RuntimeError(f"{condition} accepted HIL ventricle receipt is invalid: {hil}")
        if lumen.get("connected_components_8") != 1 or lumen.get("contour_path_counts") != {"external": 1, "tree": 1}:
            raise RuntimeError(f"{condition} ventricular boundary must have one connected outer path")
        accepted_lumen_hashes[condition] = selected_hash
        accepted_lumen_hil[condition] = hil
        tissue_outline = receipt.get("external_tissue_outline", {})
        if (
            tissue_outline.get("contour_path_counts") != {"external": 1, "tree": 1}
            or tissue_outline.get("internal_gray_contour_path_count") != 0
        ):
            raise RuntimeError(f"{condition} gray tissue outline contains an internal path")
        display_counts[condition] = receipt.get("analysis", {}).get("display_field_counts", {})
    water_completion = cartoon_samples["Water"]["ventricular_lumen_outline"].get("superior_completion", {})
    if (
        water_completion.get("applied") is not True
        or water_completion.get("centerline_top_xy") != [698, 0]
        or water_completion.get("fixed_apex_xy") != [708, 219]
        or water_completion.get("corridor_x_half_open") != [604, 769]
        or water_completion.get("centerline_x_sha256_i32_le")
        != "44372b1db1143fa5c1bf90278b663ba22e9ab2152698a7e8f0b5f8b683001c67"
        or water_completion.get("centerline_zero_clearance_pixel_count") != 8
        or water_completion.get("final_contour_path_counts") != {"external": 1, "tree": 1}
    ):
        raise RuntimeError(f"Water superior ventricular completion receipt changed: {water_completion}")
    if display_counts != expected_display_counts:
        raise RuntimeError(f"Unexpected Water/Allulose display-field counts: {display_counts}")

    if panel_a_sidecar.get("schema_version") != "figure3-panel-a-raw-microscopy-v5":
        raise RuntimeError("Panel A raw-microscopy provenance schema must be v5")
    expected_contour_counts = {
        "inset_dual_positive_nuclei_by_condition": {"Water": 24, "Sucrose": 7, "Allulose": 38},
        "outlined_npy_rois_by_condition": {"Water": 23, "Sucrose": 6, "Allulose": 36},
        "dual_positive_nuclei_with_outlined_roi_by_condition": {"Water": 23, "Sucrose": 6, "Allulose": 36},
    }
    contour_counts = panel_a_sidecar.get("npy_roi_contour_counts", {})
    if contour_counts != expected_contour_counts:
        raise RuntimeError(f"Unexpected Panel A NPY-ROI contour counts: {contour_counts}")
    panel_a_fields = {
        str(item["condition"]): item for item in panel_a_sidecar.get("representative_fields", [])
    }
    if set(panel_a_fields) != {"Water", "Sucrose", "Allulose"}:
        raise RuntimeError(f"Unexpected Panel A conditions: {sorted(panel_a_fields)}")
    contour_receipts = {
        condition: field["majority_assigned_double_positive_npy_roi_contours"]
        for condition, field in panel_a_fields.items()
    }
    strict_candidate_counts = {
        condition: int(receipt["strict_majority_inset_candidate_npy_roi_count"])
        for condition, receipt in contour_receipts.items()
    }
    if strict_candidate_counts != {"Water": 23, "Sucrose": 6, "Allulose": 37}:
        raise RuntimeError(f"Unexpected strict-majority inset candidates: {strict_candidate_counts}")
    for condition, receipt in contour_receipts.items():
        if receipt.get("schema_version") != "figure3-panel-a-npy-roi-contours-v1":
            raise RuntimeError(f"Unexpected {condition} contour receipt schema")
        if receipt.get("intensity_used_for_roi_assignment_or_display_selection") is not False:
            raise RuntimeError(f"{condition} contour receipt improperly uses display intensity")
        if int(receipt.get("outlined_npy_roi_count", -1)) != expected_contour_counts["outlined_npy_rois_by_condition"][condition]:
            raise RuntimeError(f"Unexpected {condition} outlined NPY ROI count")
        if int(receipt.get("outlined_dual_positive_nucleus_count", -1)) != expected_contour_counts["dual_positive_nuclei_with_outlined_roi_by_condition"][condition]:
            raise RuntimeError(f"Unexpected {condition} covered dual-positive nucleus count")
    water_field = panel_a_fields["Water"]
    if water_field.get("base_zoom_box_xywh_half_open") != [887, 516, 370, 303]:
        raise RuntimeError("Water2 base zoom box changed")
    if water_field.get("zoom_box_xywh_half_open") != [423, 516, 370, 303]:
        raise RuntimeError("Water2 homologous-side zoom box changed")
    contour_count_text = (
        "inset dual-positive nuclei Water2/Sucrose/Allulose n={}/{}/{}; "
        "outlined NPY ROIs and covered nuclei n={}/{}, {}/{}, and {}/{}"
    ).format(
        contour_counts["inset_dual_positive_nuclei_by_condition"]["Water"],
        contour_counts["inset_dual_positive_nuclei_by_condition"]["Sucrose"],
        contour_counts["inset_dual_positive_nuclei_by_condition"]["Allulose"],
        contour_counts["outlined_npy_rois_by_condition"]["Water"],
        contour_counts["dual_positive_nuclei_with_outlined_roi_by_condition"]["Water"],
        contour_counts["outlined_npy_rois_by_condition"]["Sucrose"],
        contour_counts["dual_positive_nuclei_with_outlined_roi_by_condition"]["Sucrose"],
        contour_counts["outlined_npy_rois_by_condition"]["Allulose"],
        contour_counts["dual_positive_nuclei_with_outlined_roi_by_condition"]["Allulose"],
    )
    shared = (
        "Bars are animal means and error bars are SD; circles are individual biological animals and crosses mark values "
        "outside within-condition Tukey 1.5-IQR fences. Marked values are retained. Omnibus values are ordinary one-way "
        "ANOVA F tests (scipy.stats.f_oneway) on animal-level ratios. Pairwise W-S, W-A, and S-A values are raw exact "
        "two-sided Mann-Whitney U p values from exhaustive independent-animal label enumeration, without multiplicity "
        "adjustment. Ratios were multiplied by 100 only for display; this positive scaling does not change ANOVA or rank-test "
        "p values. Panels D-E contain Water/Sucrose/Allulose n=3/3/3. Panels F-G are a Cellpose-derived exploratory "
        "spatial analysis whose six-shell animal vectors are tested by exhaustive animal-label PERMANOVA; cells and shells "
        "are feature-extraction units, not biological replicates. Each field is centered and covariance-whitened from its own "
        "DAPI-nucleus centroid cloud; shells are not registered to the ventricle or another anatomical landmark. The "
        "occurrence profiles jointly reflect prevalence and radial location and are not abundance-independent. PERMANOVA "
        "can reflect group-centroid and/or dispersion differences and must be interpreted with the exact dispersion audit. "
        "This exploratory legacy dataset spans acquisition batches and magnifications. WATER_NPY3 is an independent "
        "biological replicate with a provisional identifier and was acquired on Leica SP8 at 40×, whereas legacy fields "
        "came from other sessions, so condition is not separable from batch. Segmentation metadata also records WATER_NPY3 "
        "DAPI with dapi_Keyence4 rather than the modal dapi_channel setting, and E7_FR7-5/E8_FR6-1 NPY with dapi_channel "
        "rather than the modal cfos_channel3 setting."
    )
    caption = (
        "Figure 3. c-FOS activation and NPY-positive nuclei after Water, Sucrose, or Allulose.\n\n"
        "(A) Raw-TIFF reconstruction of representative Water2, Sucrose E7_FR7-5, and Allulose E8_FR6-4 fields: "
        "NPY-GFP is green and c-FOS is red, followed by the two-channel merge and a traceable magnified merge. "
        "The Water2 view uses the exact Image 1 source crop and orientation. Dashed rectangles/connectors identify each "
        "zoom source. The Water2 base zoom is mirrored horizontally to [x, y, width, height] = [423, 516, 370, 303] "
        "while y, width, and height remain unchanged, so its magnified field is on the same registered anatomical side as "
        "the Sucrose and Allulose fields. In each magnified merge, clean 1.56-pt black (#000000) contours trace native NPY Cellpose ROIs "
        "assigned to analysis-defined c-FOS-positive/NPY-positive DAPI nuclei whose nearest-neighbor transformed DAPI "
        "centroid lies inside the half-open zoom. Assignment requires one DAPI nucleus to contain a strict >50% majority "
        "of the NPY ROI's native DAPI-overlapping pixels; ties, no overlap, and no majority fail closed. When multiple "
        "strict-majority NPY ROIs map to one eligible nucleus, the ROI with the greatest native overlap is displayed, with "
        "lowest ROI ID as tie-break. Contours are derived on the complete registered ROI frame before inset translation and "
        "are clipped to the exact microscopy-image extent before black axes letterboxing, so no artificial crop-edge "
        "closure is introduced and no ROI stroke enters the letterbox. "
        f"Receipt counts are {contour_count_text}. "
        "Display intensity does not define a contour. Full-field and zoom scale bars are 100 and 50 µm, respectively. "
        "Every field/channel receives the same percentile/gamma display procedure but is stretched independently, so "
        "brightness must not be compared across conditions. No PowerPoint or screenshot pixels are used.\n"
        "(B) Water2 marker-positive DAPI nuclei in three registered views titled exactly NPY-positive nuclei, "
        "c-FOS-positive nuclei, and Combined status; display-field counts are 205, 113, and 57, respectively. "
        "(C) Allulose E8_FR6-4 marker-positive DAPI nuclei in the identical three-view grammar; display-field counts are "
        "493, 271, and 131. Counts and method details are intentionally absent from the artwork and retained here, in the "
        "panel legends, and in provenance. Both cartoons use fixed Panel-A-derived transforms into one 866×1374 display "
        "frame at 0.66421267 µm/px; label fields use nearest-neighbor registration, condition headings are centered, and "
        "the 100-µm scale-bar grammar is identical. Classification is performed on each complete native DAPI/c-FOS/NPY "
        "mask set before display registration using the established overlap rules. The gray tissue outline is drawn from "
        "the hole-filled tissue silhouette, yielding exactly one external path while the original tissue mask continues to "
        "control background fill. The base ventricular components are hash-bound raw-DAPI annotations. Water alone "
        "receives a display-only superior completion derived from the registered DAPI Cellpose-label clearance map: its "
        "original low-intensity component is extended from the unique apex (x=708,y=219) to y=0 by a deterministic "
        "dynamic-programming path maximizing cumulative Euclidean distance-transform values computed from the registered "
        "DAPI-label complement within x=[604,769), constrained to |dx|≤1 with leftmost "
        "ties; disks of radius max(1,floor(distance)-1) are unioned with the base component and holes are filled. The source "
        "Water boundary reaches y=0. In both B/C cartoons, the final ventricular lumen is a closed polygon accepted in the "
        "dedicated human-in-the-loop reviewer in this exact registered frame. Each binary mask is bound to its polygon "
        "vertices, reviewer/session/timestamp, raw DAPI source and registered DAPI-label field, and yields exactly one "
        "filled component and one outer contour. These boundaries are never "
        "analysis ROIs or denominators.\n"
        "(D) c-FOS-positive nuclei as a percentage of DAPI nuclei.\n"
        "(E) c-FOS-positive NPY nuclei as a percentage of NPY nuclei.\n"
        "(F) Cellpose-derived exploratory c-FOS-positive/DAPI field-centric radial occurrence profile across six "
        "equal-DAPI-density elliptical shells. The moderately narrowed animal radial-profile plot is on the left and the enlarged descriptive condition-mean heatmap is on the right.\n"
        "(G) Cellpose-derived exploratory c-FOS-positive NPY-positive/DAPI field-centric radial occurrence profile across "
        "the same six-shell construction. Its moderately narrowed animal radial-profile plot is on the left and its enlarged descriptive condition-mean heatmap is on the right. F and G heatmaps each use their "
        "own endpoint-specific scale from zero to that endpoint maximum and must not be compared by color across panels. "
        "Shell occurrence rates are arcsin-square-root transformed for exact multivariate testing. These profiles jointly "
        "reflect prevalence and radial location and are not abundance-independent; this is not a trained condition "
        "classifier. "
        f"{spatial_result_text}\n"
        "The formerly displayed NPY-positive/DAPI endpoint is intentionally omitted from the main plate but preserved, "
        f"with all animal values and ANOVA/exact-MWU receipts, in {supplemental_values_path.name} and "
        f"{supplemental_statistics_path.name}.\n"
        f"{shared}\n"
    )
    (legends_dir / f"{figure_name}_LEGEND.txt").write_text(caption, encoding="utf-8")
    (output_dir / f"{figure_name}_caption.txt").write_text(caption, encoding="utf-8")
    panel_descriptions = {
        "A": (
            "Equal-sized raw-TIFF NPY-GFP/c-FOS/merge/zoom tiles with calibrated scale bars; Water2 uses the "
            "homologous-side mirrored zoom [423,516,370,303]. Clean 1.56-pt black (#000000) native "
            "NPY-ROI contours, clipped to the microscopy-image extent without entering black letterboxing and without "
            "white backing or shadow, are assigned by "
            "strict-majority native DAPI overlap to analysis-defined dual-positive nuclei. Inset dual-positive counts are "
            "24/7/38; displayed ROI/covered-nucleus counts are 23/23, 6/6, and 36/36."
        ),
        "B": (
            "Water2 marker-positive DAPI nuclei in matched NPY-positive, c-FOS-positive, and combined-status views; "
            "display-field counts 205/113/57, centered condition heading, nearest-neighbor common-frame registration, "
            "100-µm bars, an external-only gray tissue outline, and a single-path black lumen boundary whose Water superior completion is derived from registered DAPI-label clearance."
        ),
        "C": (
            "Allulose E8_FR6-4 marker-positive DAPI nuclei in the identical grammar; display-field counts 493/271/131, "
            "centered condition heading, nearest-neighbor common-frame registration, 100-µm bars, and display-only "
            "external-only gray tissue and single-path black lumen-contour rules as panel B."
        ),
        "D": "Animal-level c-FOS-positive nuclei (% DAPI), n=3/3/3.",
        "E": "Animal-level c-FOS-positive NPY nuclei (% NPY), n=3/3/3.",
        "F": "Cellpose-derived exploratory c-FOS-positive/DAPI field-centric occurrence across six equal-DAPI-density elliptical shells with a moderately narrowed left-side animal radial-profile plot and an enlarged right-side descriptive condition-mean heatmap; animal-level exact PERMANOVA.",
        "G": "Cellpose-derived exploratory c-FOS-positive NPY-positive/DAPI field-centric occurrence across six equal-DAPI-density elliptical shells with a moderately narrowed left-side animal radial-profile plot and an enlarged right-side descriptive condition-mean heatmap; animal-level exact PERMANOVA.",
    }
    for letter, description in panel_descriptions.items():
        (legends_dir / f"{figure_name}_Panel_{letter}_LEGEND.txt").write_text(
            f"Figure 3, panel {letter}. {description}\n\n{shared}\n", encoding="utf-8"
        )

    input_records = []
    for role, path in inputs.items():
        input_records.append({
            "record_type": "input", "role": role, "path": str(path.resolve()),
            "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
        })
    output_paths: dict[str, Path] = {
        "figure_pdf": output_dir / f"{figure_name}.pdf",
        "figure_png": output_dir / f"{figure_name}.png",
        "figure_caption": output_dir / f"{figure_name}_caption.txt",
        "figure_legend": legends_dir / f"{figure_name}_LEGEND.txt",
        "displayed_animal_values": animal_plot_values_path,
        "displayed_anova_mwu_statistics": statistics_used,
        "supplemental_npy_dapi_animal_values": supplemental_values_path,
        "supplemental_npy_dapi_statistics": supplemental_statistics_path,
        "supplemental_npy_dapi_readme": supplemental_note_path,
        **{f"figure_{role}": path for role, path in spatial_source_copies.items()},
    }
    for letter in "ABCDEFG":
        output_paths[f"panel_{letter}_pdf"] = output_dir / "panels" / f"{figure_name}_Panel_{letter}.pdf"
        output_paths[f"panel_{letter}_png"] = output_dir / "panels" / f"{figure_name}_Panel_{letter}.png"
        output_paths[f"panel_{letter}_legend"] = legends_dir / f"{figure_name}_Panel_{letter}_LEGEND.txt"
    output_records = []
    for role, path in output_paths.items():
        require_file(path, f"generated Figure 3 {role}")
        output_records.append({
            "record_type": "output", "role": role, "path": str(path.resolve()),
            "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
        })
    provenance = {
        "schema_version": "figure3-cfos-npy-v6",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": str(Path(__file__).resolve()),
        "analysis_directory": str(analysis_dir.resolve()),
        "conditions": list(CONDITIONS),
        "condition_colors": COLORS,
        "panel_map": {
            "A": "raw NPY/c-FOS microscopy, merge, and traceable zoom",
            "B": "Water marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
            "C": "Allulose marker-positive DAPI nuclei: NPY, c-FOS, and combined/dual status",
            "D": "c-FOS-positive nuclei over DAPI nuclei",
            "E": "c-FOS-positive NPY nuclei over NPY nuclei",
            "F": "c-FOS-positive/DAPI field-centric six-shell radial occurrence",
            "G": "c-FOS-positive NPY-positive/DAPI field-centric six-shell radial occurrence",
        },
        "layout": {
            "nominal_figure_inches": [16.0, 17.1],
            "pdf_embedded_raster_dpi": PDF_RASTER_DPI,
            "outer_grid": "four bands: full-width equal-tile A; full-width matched side-by-side B/C; compact centered D/E; wide side-by-side F/G",
            "outer_bounds_normalized": {"left": 0.020, "right": 0.992, "bottom": 0.025, "top": 0.985},
            "outer_height_ratios": [6.20, 2.62, 3.10, 3.72],
            "outer_hspace": 0.055,
            "panel_A_width_inches": 15.55,
            "panels_B_C_width_inches_each": 7.60,
            "panels_D_E_slot_inches_each": [3.50, 3.10],
            "panels_F_G_width_inches_each": 7.72,
            "spatial_row_width_ratios": list(SPATIAL_ROW_WIDTH_RATIOS),
            "spatial_panel_geometry": spatial_layout_geometry,
            "spatial_native_aspect_preserved": True,
            "spatial_row_uses_unused_bottom_margin": True,
            "spatial_panel_label_placement": "inside the top-left source-image whitespace",
            "spatial_source_panel_layout": "moderately narrowed radial-profile plot left; enlarged condition-means heatmap right",
            "spatial_source_panel_figure_inches": list(SPATIAL_SOURCE_PANEL_FIGSIZE),
            "intentional_whitespace": "minimal F/G side gutters and center gap; compact D/E row preserves Figure4/5 plot proportions",
        },
        "display_scale": "100 * animal ratio in D/E; source ratios used for ANOVA and exact MWU",
        "endpoint_counts": endpoint_counts,
        "supplemental_removed_endpoint": {
            "endpoint": "npy_over_dapi",
            "reason": "removed from main plate at user request; values and receipts explicitly retained",
            "group_counts": supplemental_counts,
            "values_file": str(supplemental_values_path.resolve()),
            "statistics_file": str(supplemental_statistics_path.resolve()),
            "readme_file": str(supplemental_note_path.resolve()),
        },
        "panel_a": {
            "source": "raw TIFF channels plus native Cellpose masks listed and hash-bound in the Panel A sidecar",
            "render": "fixed audited field transforms; per-field/channel percentile stretch; calibrated scale bars",
            "colors": {"NPY": "green", "c-FOS": "red"},
            "water_representative": "WATER_NPY2; exact Image 1 integer crop and 90-degree counter-clockwise orientation",
            "water_homologous_side_zoom": {
                "base_xywh_half_open": [887, 516, 370, 303],
                "final_xywh_half_open": [423, 516, 370, 303],
                "rule": "exact horizontal mirror; y, width, and height unchanged",
            },
            "npy_roi_assignment_rule": (
                "one native DAPI nucleus must contain a strict >50% majority of a native NPY ROI's "
                "DAPI-overlapping pixels; ties, no-overlap, and no-majority cases fail closed"
            ),
            "contour_inclusion_rule": (
                "assigned nucleus is analysis-defined c-FOS-positive/NPY-positive and its nearest-neighbor "
                "transformed DAPI centroid lies in the half-open zoom"
            ),
            "one_displayed_roi_per_nucleus_rule": (
                "when multiple strict-majority NPY ROIs map to one eligible nucleus, display the ROI with "
                "greatest native overlap; lowest ROI ID breaks ties"
            ),
            "contour_geometry": (
                "derived on the complete registered ROI frame, translated to inset coordinates, and clipped "
                "to the exact microscopy-image extent before axes letterboxing without adding a crop-edge closure"
            ),
            "contour_render": {
                "stroke_color": "black",
                "stroke_color_hex": NPY_ROI_CONTOUR_COLOR,
                "linewidth_pt": NPY_ROI_CONTOUR_WIDTH_PT,
                "white_backing_or_shadow": False,
                "stroke_count_per_path": 1,
                "clip_target": "microscopy image extent before axes letterboxing",
                "black_letterbox_stroke_allowed": False,
            },
            "strict_majority_inset_candidate_npy_roi_counts_by_condition": strict_candidate_counts,
            "inset_dual_positive_nucleus_counts_by_condition": contour_counts["inset_dual_positive_nuclei_by_condition"],
            "outlined_npy_roi_counts_by_condition": contour_counts["outlined_npy_rois_by_condition"],
            "dual_positive_nuclei_with_outlined_roi_counts_by_condition": contour_counts["dual_positive_nuclei_with_outlined_roi_by_condition"],
            "roi_mapping_and_contour_receipts": str(inputs["panel_A_raw_reconstruction_provenance"].resolve()),
            "tile_geometry": panel_a_sidecar["tile_geometry"],
            "display_intensity_used_for_contour_assignment_or_selection": False,
            "no_powerpoint_or_screenshot_pixels": True,
        },
        "cartoons": {
            "panel_roles": {
                "B": "Water marker-positive DAPI nuclei",
                "C": "Allulose marker-positive DAPI nuclei",
            },
            "views_in_order": [
                "NPY-positive DAPI nuclei", "c-FOS-positive DAPI nuclei", "combined/dual status",
            ],
            "view_title_text": ["NPY-positive nuclei", "c-FOS-positive nuclei", "Combined status"],
            "common_frame_shape_yx": [866, 1374],
            "common_frame_um_per_px": 0.664212669961432,
            "classification_coordinates": "complete native DAPI/c-FOS/NPY label fields before display transform",
            "label_registration": "fixed Panel-A-derived homographies; nearest neighbor only",
            "same_physical_canvas_colors_titles_and_scale_bar_grammar": True,
            "panel_title_x_figure_fraction": 0.5,
            "panel_title_horizontal_alignment": "center",
            "external_tissue_boundaries": (
                "binary_fill_holes(display tissue silhouette) for contouring only; original tissue mask controls fill; "
                "exactly one external gray path and no internal gray paths per view"
            ),
            "ventricular_lumen_boundaries": (
                "hash-bound raw-DAPI source masks; Water adds registered DAPI-label-clearance superior completion; "
                "accepted HIL closed polygons are source-, reviewer-, timestamp- and mask-hash-bound; "
                "exactly one outer black path per view; display annotations only, never analysis ROIs or denominators"
            ),
            "water_superior_lumen_extension": {
                "base_apex_xy": [708, 219],
                "centerline_rule": (
                    "dynamic-programming path maximizing cumulative Euclidean distance-transform values computed from "
                    "the registered DAPI-label complement; fixed apex, free top endpoint, |dx|<=1, leftmost ties"
                ),
                "corridor_x_half_open": [604, 769],
                "centerline_top_xy": [698, 0],
                "centerline_path_sha256_i32_le": "44372b1db1143fa5c1bf90278b663ba22e9ab2152698a7e8f0b5f8b683001c67",
                "final_area_px": 35541,
                "final_bbox_yxyx_half_open": [0, 537, 606, 802],
                "final_binary_sha256": expected_source_lumen_hashes["Water"],
                "top_row_x_inclusive": [692, 704],
            },
            "accepted_ventricle_hil": accepted_lumen_hil,
            "ventricular_source_binary_sha256": expected_source_lumen_hashes,
            "in_panel_method_footnotes": False,
            "in_panel_status_counts": False,
            "count_location": "caption, panel legends, and cartoon provenance sidecar",
            "method_location": "caption, panel legends, and cartoon provenance sidecar",
            "source_receipt": str(inputs["cartoon_provenance"].resolve()),
            "ventricular_boundary_binary_sha256": accepted_lumen_hashes,
            "external_tissue_outline_binary_sha256": {
                condition: cartoon_samples[condition]["external_tissue_outline"]["outline_binary_sha256"]
                for condition in ("Water", "Allulose")
            },
            "display_field_counts": display_counts,
        },
        "statistics": (
            "Panels D/E: ordinary scipy.stats.f_oneway ANOVA on animal values; raw exact exhaustive two-sided "
            "Mann-Whitney U pairwise tests without multiplicity adjustment. Panels F/G use arcsin-square-root-transformed "
            "six-shell animal vectors with exhaustive 1,680-label PERMANOVA, exact multivariate-dispersion diagnostics, "
            "and separate two-endpoint BH q families. The PERMANOVA can reflect centroid and/or dispersion differences. "
            "This is a Cellpose-derived exploratory spatial analysis, not a trained condition CNN or pooled-cell inference."
        ),
        "spatial_analysis": {
            "panels": ["F", "G"],
            "shells": (
                "six equal-DAPI-density elliptical radial shells after centering and covariance whitening "
                "each animal field from its own DAPI-nucleus centroid cloud; not ventricle/anatomy registered"
            ),
            "feature_transform": "arcsin(sqrt(shell positive nuclei / shell DAPI nuclei))",
            "panel_inset": (
                "descriptive 3-condition by 6-shell heatmap of condition-mean occurrence percentages; "
                "F/G use separate endpoint-specific zero-to-maximum color scales"
            ),
            "inference": "exact exhaustive 1,680-allocation animal-label PERMANOVA; biological animal is the sole inferential unit",
            "profiles_are_abundance_independent": False,
            "trained_condition_classifier": False,
            "permanova_results": {
                endpoint: {
                    "pseudo_f": float(spatial_stats.loc[endpoint, "pseudo_f"]),
                    "p_value_exact": float(spatial_stats.loc[endpoint, "p_value_exact"]),
                    "bh_q_value_two_endpoint_family": float(spatial_stats.loc[endpoint, "bh_q_value_two_endpoint_family"]),
                }
                for endpoint in ("cfos_occurrence", "double_occurrence")
            },
            "dispersion_diagnostic": {
                "test": "exact animal-label ANOVA F on Euclidean distances to allocated group centroids",
                "permanova_centroid_only_interpretation_allowed": False,
                "results": {
                    endpoint: {
                        "dispersion_f": float(dispersion_stats.loc[endpoint, "dispersion_f"]),
                        "p_value_exact": float(dispersion_stats.loc[endpoint, "p_value_exact"]),
                        "bh_q_value_two_endpoint_diagnostic_family": float(
                            dispersion_stats.loc[endpoint, "bh_q_value_two_endpoint_diagnostic_family"]
                        ),
                    }
                    for endpoint in ("cfos_occurrence", "double_occurrence")
                },
            },
            "figure_source_data_copies": {role: str(path.resolve()) for role, path in spatial_source_copies.items()},
        },
        "caveats": [
            "exploratory legacy analysis",
            "acquisition/batch confounding",
            "WATER_NPY3 animal identifier provisional",
            "Panel A display stretches are independent and cannot support cross-condition intensity comparison",
            "Panels F/G occurrence profiles jointly reflect prevalence and radial location and are not abundance-independent",
            "PERMANOVA may reflect group-centroid and/or dispersion differences; interpret with the exact dispersion audit",
            "F/G field coordinates are DAPI-cloud centered/whitened and are not anatomy or ventricle aligned",
            "WATER_NPY3 DAPI model setting differs from modal DAPI; E7_FR7-5 and E8_FR6-1 NPY settings differ from modal NPY",
            "F/G heatmap colors use separate endpoint-specific scales and cannot be compared across panels",
        ],
        "inputs": input_records,
        "outputs": output_records,
    }
    provenance_path = provenance_dir / f"{figure_name}_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_rows = input_records + output_records + [{
        "record_type": "output", "role": "figure_provenance", "path": str(provenance_path.resolve()),
        "sha256": sha256_file(provenance_path), "size_bytes": provenance_path.stat().st_size,
    }]
    pd.DataFrame(manifest_rows).to_csv(provenance_dir / f"{figure_name}_manifest.csv", index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--spatial-dir", type=Path, default=DEFAULT_SPATIAL_DIR)
    parser.add_argument("--per-animal", type=Path, default=None)
    parser.add_argument(
        "--statistics", "--exact-statistics", dest="statistics", type=Path, default=None,
        help="ANOVA + exact-MWU statistics table (legacy alias retained for CLI compatibility).",
    )
    parser.add_argument("--panel-b", type=Path, default=None)
    parser.add_argument("--panel-c", type=Path, default=None)
    parser.add_argument("--panel-f", type=Path, default=None)
    parser.add_argument("--panel-g", type=Path, default=None)
    parser.add_argument("--figure-name", default=FIGURE_NAME)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--panel-a-dpi", type=int, default=600)
    parser.add_argument("--seed", type=int, default=31)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.dpi < 300 or args.panel_a_dpi < 300:
        raise ValueError("Publication master, panels, and Panel A must each be at least 300 dpi")
    analysis_dir = args.analysis_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    raw_root = args.raw_root.expanduser().resolve()
    spatial_dir = args.spatial_dir.expanduser().resolve()
    panels_input = analysis_dir / "panels"
    panels_input.mkdir(parents=True, exist_ok=True)
    if args.per_animal is not None:
        per_animal_path = require_file(args.per_animal, "animal-level values")
    else:
        candidates = (
            analysis_dir / "per_animal_cfos_npy.csv",
            analysis_dir / "per_animal_whole_field_cfos_npy.csv",
        )
        found = [path for path in candidates if path.is_file() and path.stat().st_size > 0]
        if len(found) != 1:
            raise FileNotFoundError(
                "Expected exactly one canonical animal-level table; found "
                + str([str(path) for path in found])
            )
        per_animal_path = found[0].resolve()
    statistics_path = require_file(
        args.statistics or analysis_dir / "anova_mwu_statistics.csv",
        "one-way ANOVA and exact-MWU statistics",
    )
    panel_b_path = require_file(
        args.panel_b or panels_input / "Figure3_Panel_B_Water_marker_positive_nuclei.png",
        "Panel B cartoon",
    )
    panel_c_path = require_file(
        args.panel_c or panels_input / "Figure3_Panel_C_Allulose_marker_positive_nuclei.png",
        "Panel C cartoon",
    )
    panel_f_path = require_file(
        args.panel_f or spatial_dir / "Figure3_Spatial_cFOS_occurrence.png",
        "Panel F c-FOS spatial-occurrence plot",
    )
    panel_g_path = require_file(
        args.panel_g or spatial_dir / "Figure3_Spatial_cFOS_NPY_occurrence.png",
        "Panel G c-FOS+NPY spatial-occurrence plot",
    )

    raw_values = pd.read_csv(per_animal_path)
    values = normalize_animal_values(raw_values)
    raw_stats = pd.read_csv(statistics_path)
    stats_by_endpoint: dict[
        str, tuple[float, float, dict[tuple[str, str], float]]
    ] = {}
    selected_stats: list[pd.DataFrame] = []
    for _, endpoint, _, _, _ in ENDPOINTS:
        anova_p, anova_f, pairs, selected = exact_statistics_for_endpoint(raw_stats, endpoint)
        stats_by_endpoint[endpoint] = anova_p, anova_f, pairs
        selected_stats.append(selected)
    supplemental_stats: list[pd.DataFrame] = []
    for endpoint, _, _ in SUPPLEMENTAL_ENDPOINTS:
        _, _, _, selected = exact_statistics_for_endpoint(raw_stats, endpoint)
        supplemental_stats.append(selected)

    panel_a_path = panels_input / "Figure3_Panel_A_raw_microscopy.png"
    panel_a_pdf = panels_input / "Figure3_Panel_A_raw_microscopy.pdf"
    panel_a_provenance = panels_input / "Figure3_Panel_A_raw_microscopy_provenance.json"
    panel_a, panel_a_inputs = build_raw_microscopy_panel(
        raw_root, panel_a_path, panel_a_pdf, panel_a_provenance, dpi=args.panel_a_dpi,
    )
    stale_panel_a = panels_input / "Figure3_Panel_A_reference_microscopy.png"
    if stale_panel_a.is_file():
        stale_panel_a.unlink()
    panel_b = read_rgb(panel_b_path)
    panel_c = read_rgb(panel_c_path)
    panel_f = read_rgb(panel_f_path)
    panel_g = read_rgb(panel_g_path)
    panel_a_aspect = panel_a.shape[1] / panel_a.shape[0]
    if not 2.35 <= panel_a_aspect <= 2.65:
        raise ValueError(f"Raw-reconstructed Panel A has unexpected compact 3x4 aspect ratio: {panel_a.shape}")
    for name, panel in (("B", panel_b), ("C", panel_c)):
        aspect = panel.shape[1] / panel.shape[0]
        if not 2.85 <= aspect <= 3.00:
            raise ValueError(
                f"Panel {name} has unexpected footnote-free three-view aspect ratio: {panel.shape}"
            )
    for name, panel in (("F", panel_f), ("G", panel_g)):
        aspect = panel.shape[1] / panel.shape[0]
        expected_aspect = SPATIAL_SOURCE_PANEL_FIGSIZE[0] / SPATIAL_SOURCE_PANEL_FIGSIZE[1]
        if not np.isclose(aspect, expected_aspect, rtol=0.0, atol=SPATIAL_SOURCE_ASPECT_ATOL):
            raise ValueError(
                f"Panel {name} has unexpected spatial-occurrence aspect ratio: "
                f"{panel.shape}; expected {expected_aspect:.6f}"
            )

    matplotlib.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans",
        "axes.linewidth": 1.6, "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(16.0, 17.1), facecolor="white")
    outer = fig.add_gridspec(
        4, 12,
        left=0.020, right=0.992, bottom=0.025, top=0.985,
        height_ratios=(6.20, 2.62, 3.10, 3.72), hspace=0.055,
    )
    panel_axes: dict[str, list[plt.Axes]] = {}
    panel_texts: dict[str, list[Artist]] = {}

    # Panel A spans the full content width. Its generator enforces equal 1:1:1:1
    # subpanel columns; imshow preserves the registered asset aspect without
    # set_position overrides.
    ax_a = fig.add_subplot(outer[0, :])
    ax_a.imshow(panel_a, interpolation="nearest")
    ax_a.set_axis_off()
    panel_axes["A"] = [ax_a]
    panel_texts["A"] = [add_panel_letter(ax_a, "A", x=0.018, y=1.008)]

    cartoon_row = outer[1, :].subgridspec(
        1, 5, width_ratios=(0.03, 7.60, 0.29, 7.60, 0.03), wspace=0.0,
    )
    for column_index, (letter, image) in zip((1, 3), (("B", panel_b), ("C", panel_c))):
        axis = fig.add_subplot(cartoon_row[0, column_index])
        axis.imshow(image, interpolation="nearest")
        axis.set_axis_off()
        panel_axes[letter] = [axis]
        # Panels B/C carry their letters inside the provenance-bound assets.
        panel_texts[letter] = []

    # Match Figure 4/5 plot proportions: the two barplots remain compact,
    # approximately 3.55 in wide with a 2.4-in plot axis plus the protected
    # tick-label gutter, rather than being stretched across the full plate.
    statistics_row = outer[2, :].subgridspec(
        1, 5, width_ratios=(3.75, 3.50, 1.05, 3.50, 3.75), wspace=0.0,
    )
    for column_index, (index, item) in zip((1, 3), enumerate(ENDPOINTS)):
        letter, endpoint, title, ylabel, _ = item
        # The final row is a protected gutter for the 24-degree x tick labels.
        # Without it, the later-created opaque F/G axes paint over the lower
        # portion of Water/Sucrose/Allulose in the master figure.
        endpoint_grid = statistics_row[0, column_index].subgridspec(
            3, 1, height_ratios=(0.40, 2.42, 0.28), hspace=0.02,
        )
        header_axis = fig.add_subplot(endpoint_grid[0, 0])
        header_axis.set_axis_off()
        header_axis.text(
            0.0, 0.95, letter, transform=header_axis.transAxes,
            fontsize=20, fontweight="bold", ha="left", va="top",
        )
        header_axis.text(
            0.5, 0.95, title, transform=header_axis.transAxes,
            fontsize=11.2, fontweight="bold", ha="center", va="top",
        )
        axis = fig.add_subplot(endpoint_grid[1, 0])
        anova_p, anova_f, pairs = stats_by_endpoint[endpoint]
        draw_quantitative_panel(
            axis, values, endpoint, ylabel, anova_p, anova_f, pairs,
            args.seed + 100 * index,
        )
        panel_axes[letter] = [header_axis, axis]
        panel_texts[letter] = []

    spatial_row = outer[3, :].subgridspec(
        1, 5, width_ratios=SPATIAL_ROW_WIDTH_RATIOS, wspace=0.0,
    )
    for column_index, (letter, image) in zip((1, 3), (("F", panel_f), ("G", panel_g))):
        axis = fig.add_subplot(spatial_row[0, column_index])
        axis.imshow(image, interpolation="nearest")
        # The GridSpec row is height-limited. Preserve its top edge (and D/E
        # tick-label clearance), then use the idle bottom margin for the exact
        # native-aspect height required by the wider slot.
        slot = spatial_row[0, column_index].get_position(fig)
        source_aspect = float(image.shape[1] / image.shape[0])
        native_height = (
            slot.width * fig.get_figwidth() / source_aspect / fig.get_figheight()
        )
        native_bottom = slot.y1 - native_height
        if native_bottom < SPATIAL_PANEL_MIN_BOTTOM_FRACTION:
            raise RuntimeError(
                f"Panel {letter} wider native-aspect box would exceed the protected bottom margin"
            )
        axis.set_position((slot.x0, native_bottom, slot.width, native_height))
        axis.set_axis_off()
        panel_axes[letter] = [axis]
        panel_texts[letter] = [add_panel_letter(
            axis, letter,
            x=SPATIAL_PANEL_LABEL_X, y=SPATIAL_PANEL_LABEL_Y, ha="left",
        )]

    if set(panel_axes) != set("ABCDEFG"):
        raise RuntimeError(f"Figure 3 panel registry must contain exactly A-G; found {sorted(panel_axes)}")
    validate_quantitative_tick_clearance(
        fig,
        {letter: panel_axes[letter][-1] for letter in "DE"},
        {letter: panel_axes[letter][0] for letter in "FG"},
    )
    spatial_layout_geometry = validate_spatial_panel_geometry(
        fig,
        {letter: panel_axes[letter][0] for letter in "FG"},
        {"F": panel_f, "G": panel_g},
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / args.figure_name
    # Preserve the intentional centered whitespace on the exact nominal plate.
    # Individual panel crops remain tight through export_panel_crops below.
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", dpi=PDF_RASTER_DPI, facecolor="white")
    # The master PNG was capped at 300 while the PDF rasters and every panel crop
    # took the full --dpi, so Figure 3's two deliverables disagreed. Both now
    # honour --dpi, which the run-book passes as 600.
    fig.savefig(
        stem.with_suffix(".png"), format="png", dpi=args.dpi,
        facecolor="white",
    )
    export_panel_crops(fig, panel_axes, panel_texts, panels_dir, args.figure_name, args.dpi)
    plt.close(fig)
    # Panel letters A-G are all live in this composition; no stale crop is retained.

    inputs: dict[str, Path] = {
        "animal_values": per_animal_path,
        "anova_mwu_statistics": statistics_path,
        "panel_A_raw_reconstruction": panel_a_path,
        "panel_A_raw_reconstruction_pdf": panel_a_pdf,
        "panel_A_raw_reconstruction_provenance": panel_a_provenance,
        "panel_B_Water_marker_positive_nuclei": panel_b_path,
        "panel_C_Allulose_marker_positive_nuclei": panel_c_path,
        "panel_F_cfos_spatial_occurrence": panel_f_path,
        "panel_G_cfos_npy_spatial_occurrence": panel_g_path,
        "panel_B_vector": require_file(
            panels_input / "Figure3_Panel_B_Water_marker_positive_nuclei.pdf", "Panel B vector cartoon",
        ),
        "panel_C_vector": require_file(
            panels_input / "Figure3_Panel_C_Allulose_marker_positive_nuclei.pdf", "Panel C vector cartoon",
        ),
        "cartoon_provenance": require_file(
            panels_input / "Figure3_cFOS_NPY_cartoon_provenance.json", "matched cartoon provenance",
        ),
        "cartoon_provenance_csv": require_file(
            panels_input / "Figure3_cFOS_NPY_cartoon_provenance.csv", "matched cartoon provenance table",
        ),
        "generator": Path(__file__).resolve(),
    }
    for index, source_path in enumerate(dict.fromkeys(path.resolve() for path in panel_a_inputs), start=1):
        inputs[f"panel_A_source_{index:02d}_{source_path.stem}"] = source_path
    for optional_role, filename in (
        ("segmentation_qc", "segmentation_qc.csv"),
        ("analysis_source_manifest", "source_manifest.csv"),
        ("statistical_test_receipts", "statistical_test_receipts.csv"),
    ):
        candidate = analysis_dir / filename
        if candidate.is_file():
            inputs[optional_role] = candidate
    for spatial_role, filename in (
        ("spatial_F_vector", "Figure3_Spatial_cFOS_occurrence.pdf"),
        ("spatial_G_vector", "Figure3_Spatial_cFOS_NPY_occurrence.pdf"),
        ("spatial_exact_permanova", "exact_spatial_permanova.csv"),
        ("spatial_exact_permutation_receipts", "exact_spatial_permutation_receipts.csv"),
        ("spatial_exact_dispersion", "exact_spatial_dispersion.csv"),
        ("spatial_exact_dispersion_permutation_receipts", "exact_spatial_dispersion_permutation_receipts.csv"),
        ("spatial_animal_features", "animal_radial_occurrence_features.csv"),
        ("spatial_animal_counts", "animal_spatial_counts.csv"),
        ("spatial_model_setting_audit", "segmentation_model_setting_audit.csv"),
        ("spatial_provenance", "spatial_distribution_provenance.json"),
        ("spatial_output_manifest", "output_manifest.csv"),
    ):
        inputs[spatial_role] = require_file(spatial_dir / filename, spatial_role.replace("_", " "))
    write_text_and_provenance(
        output_dir, analysis_dir, inputs, values,
        pd.concat(selected_stats, ignore_index=True),
        pd.concat(supplemental_stats, ignore_index=True), args.figure_name,
        spatial_layout_geometry,
    )
    print(f"Wrote {stem.with_suffix('.pdf')}")
    print(f"Wrote {stem.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
