#!/usr/bin/env python3
"""Render Figure S2 female single-bottle outcomes and live panel A-D crops."""

from __future__ import annotations

import argparse
import atexit
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# Stabilize PDF timestamps for byte-reproducible command-level reruns.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.35,
    "axes.spines.top": False,
    "axes.spines.right": False,
})
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
sys.path.insert(0, str(PAPER_ROOT / "scripts" / "shared"))
from spanish_panel_text import translate_figure_texts_to_spanish
ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "FigS2"
DEFAULT_ANALYSIS = ANALYSIS_ROOT / "results"
DEFAULT_OUTPUT = PAPER_ROOT / "FigS2"

FIGURE_NUMBER = "S2"
TARGET_SEX = "Female"
STEM_NAME = "Figure_S2_single_bottle_female"
CONDITIONS = ("Water", "Sucrose", "Allulose")
COLORS = {"Water": "#4C9BD6", "Sucrose": "#D95F5F", "Allulose": "#39A96B"}
OMNIBUS_TEST = "One_way_ANOVA"
EXACT_OMNIBUS_TEST = "exact cage-label permutation omnibus F"
PAIRWISE_TEST = "exact cage-label permutation mean difference"
CONSUMPTION_OUTCOME = "Day-6 cumulative cage consumption (mL/cage)"
WEIGHT_OUTCOME = "Day-6 cage-mean body-weight change (%)"
BOOTSTRAP_REPLICATES = 20_000
P_VALUE_GLOBAL_FONTSIZE_PT = 10.1
P_VALUE_PAIRWISE_FONTSIZE_PT = 9.2
P_VALUE_FONTWEIGHT = "bold"
FIGURE_SIZE_INCHES = (10.8, 10.2)
OUTPUT_SUBDIRS = ("panels", "legends", "source_data", "provenance")
GENERATED_TREE_MARKER = ".paper_generated_tree.json"
OUTPUT_PRODUCER = "Paper/scripts/FigS2/02_make_figure_s2_female.py"
PANEL_MAP = {
    (TARGET_SEX, "consumption_trajectory"): "A",
    (TARGET_SEX, "consumption_endpoint"): "B",
    (TARGET_SEX, "weight_trajectory"): "C",
    (TARGET_SEX, "weight_endpoint"): "D",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20260816)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_inventory(root: Path) -> tuple[list[str], dict[str, str]]:
    """Return the complete, symlink-free directory and file inventory."""
    directories: list[str] = []
    files: dict[str, str] = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise RuntimeError(f"Generated output must not contain symlinks: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_dir():
            directories.append(relative)
        elif path.is_file() and relative != GENERATED_TREE_MARKER:
            files[relative] = sha256(path)
        elif not path.is_file():
            raise RuntimeError(f"Unsupported output entry: {path}")
    return directories, files


def validate_replacement_target(final_output: Path) -> None:
    """Fail closed unless an existing nonempty tree is wholly script-owned."""
    if final_output.is_symlink():
        raise RuntimeError(f"Refusing symlink output directory: {final_output}")
    if not final_output.exists():
        return
    if not final_output.is_dir():
        raise RuntimeError(f"Output path is not a directory: {final_output}")
    if not any(final_output.iterdir()):
        return
    marker_path = final_output / GENERATED_TREE_MARKER
    if not marker_path.is_file() or marker_path.is_symlink():
        raise RuntimeError(
            f"Refusing to replace nonempty unowned directory {final_output}: it "
            f"holds no trusted {GENERATED_TREE_MARKER}, so this script cannot "
            f"tell that it produced what is in there and will not delete it. "
            f"Point --output-dir at a new or empty directory, or remove that one "
            f"yourself if you know it is disposable."
        )
    try:
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Invalid generated-tree marker: {marker_path}") from exc
    if marker.get("schema") != "paper-generated-tree-v1" or marker.get("producer") != OUTPUT_PRODUCER:
        raise RuntimeError(f"Output marker does not authorize {OUTPUT_PRODUCER}: {marker_path}")
    directories, files = tree_inventory(final_output)
    if directories != marker.get("directories") or files != marker.get("files_sha256"):
        raise RuntimeError(
            f"Refusing to replace modified or untracked output tree: "
            f"{final_output}. It carries this script's marker, but its contents no "
            f"longer match what the marker recorded, so something else edited or "
            f"added to it. Point --output-dir somewhere new, or remove that tree "
            f"yourself if the changes were yours and are disposable."
        )


def write_generated_tree_marker(staging: Path) -> None:
    directories, files = tree_inventory(staging)
    marker = {
        "schema": "paper-generated-tree-v1",
        "producer": OUTPUT_PRODUCER,
        "directories": directories,
        "files_sha256": files,
    }
    (staging / GENERATED_TREE_MARKER).write_text(
        json.dumps(marker, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def prepare_staging_output(final_output: Path) -> Path:
    """Create a complete sibling staging tree after validating target ownership."""
    final_output.parent.mkdir(parents=True, exist_ok=True)
    validate_replacement_target(final_output)
    staging = Path(tempfile.mkdtemp(prefix=f".{final_output.name}.stage-", dir=final_output.parent))
    for name in OUTPUT_SUBDIRS:
        (staging / name).mkdir()
    atexit.register(shutil.rmtree, staging, ignore_errors=True)
    return staging


def publish_staging_output(staging: Path, final_output: Path) -> None:
    """Publish a complete owned tree atomically, restoring the old tree on error."""
    write_generated_tree_marker(staging)
    validate_replacement_target(final_output)
    backup = staging.with_name(staging.name + ".previous")
    if backup.exists():
        raise RuntimeError(f"Refusing pre-existing publication backup: {backup}")
    had_previous = final_output.exists()
    try:
        if had_previous:
            os.replace(final_output, backup)
        os.replace(staging, final_output)
    except Exception:
        if backup.exists() and not final_output.exists():
            os.replace(backup, final_output)
        raise
    else:
        if backup.exists():
            shutil.rmtree(backup)

def write_local_readme(output_dir: Path) -> None:
    text = f"""FIGURE {FIGURE_NUMBER} - {TARGET_SEX.upper()} SINGLE-BOTTLE EXPERIMENT

The combined publication PDF and PNG are at the parent Fig{FIGURE_NUMBER} directory's root.
Exact live A-D crops are in panels/, legends and the caption are in legends/, plotted
values/statistics/confidence intervals/source hashes/audit CSVs are in
source_data/, and the deterministic manifest, scientific audit, and this
README are in provenance/. Every publication artifact is a real file.

Panels A-D show the {TARGET_SEX.lower()} consumption trajectory/endpoint and
body-weight trajectory/endpoint. Behavior-specific exact cage-label
permutation tests are retained because this is a small cage-randomized design;
pairwise tests are Holm-adjusted within sex and endpoint. No treatment-by-sex
interaction is claimed.

Reproduce from the Paper directory:

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig{FIGURE_NUMBER}/01_analyze_single_bottle_{TARGET_SEX.lower()}.py
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig{FIGURE_NUMBER}/02_make_figure_{FIGURE_NUMBER.lower()}_{TARGET_SEX.lower()}.py

The renderer validates a complete sibling staging tree and replaces only an empty or
previously marker-owned tree. It refuses symlinks, modified files, and untracked
content, so the recovered live Fig{FIGURE_NUMBER} directory must not be used as
--output-dir; render to a fresh directory and promote reviewed leaf artifacts.
"""
    (output_dir / "README.txt").write_text(text, encoding="utf-8")


def p_text(value: float) -> str:
    if not np.isfinite(value):
        return "p = NA"
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.3f}"


def condition_counts(data: pd.DataFrame, unit: str) -> dict[str, int]:
    return (
        data[[unit, "Treatment"]].drop_duplicates()
        .groupby("Treatment")[unit].nunique().to_dict()
    )


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.13, 1.08, letter, transform=ax.transAxes, fontsize=19,
        fontweight="bold", ha="left", va="top", clip_on=False,
    )


def trajectory_summary(
    data: pd.DataFrame,
    value_column: str,
    sex: str,
    condition: str,
    panel: str,
    seed: int,
) -> pd.DataFrame:
    if data.duplicated(["E", "study_day"]).any():
        raise ValueError(f"Panel {panel} has duplicate cage/day rows")
    pivot = (
        data.pivot(index="E", columns="study_day", values=value_column)
        .sort_index(axis=0).sort_index(axis=1)
    )
    if pivot.empty or pivot.isna().any().any() or tuple(pivot.columns) != (1, 3, 6):
        raise ValueError(f"Panel {panel}, {condition}: incomplete cage trajectory")
    values = pivot.to_numpy(float)
    means = values.mean(axis=0)
    if values.shape[0] < 2:
        low = np.full(values.shape[1], np.nan)
        high = np.full(values.shape[1], np.nan)
        estimable = False
        method = "not estimable: fewer than 2 cages"
    else:
        rng = np.random.default_rng(seed)
        draw_indices = rng.integers(
            0, values.shape[0], size=(BOOTSTRAP_REPLICATES, values.shape[0])
        )
        bootstrap_means = values[draw_indices, :].mean(axis=1)
        low, high = np.percentile(bootstrap_means, [2.5, 97.5], axis=0)
        estimable = True
        method = "pointwise percentile cluster bootstrap of whole cage trajectories"
    return pd.DataFrame({
        "panel": panel,
        "Sex": sex,
        "Treatment": condition,
        "study_day": pivot.columns.to_numpy(int),
        "n_cages": values.shape[0],
        "cage_balanced_mean": means,
        "ci95_low": low,
        "ci95_high": high,
        "ci_estimable": estimable,
        "ci_method": method,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES if estimable else 0,
        "bootstrap_seed": seed,
    })


def style_trajectory(
    axes: list[plt.Axes],
    anchor: plt.Axes,
    sex: str,
    title: str,
    key: str,
    ylabel: str,
    letter: str,
) -> None:
    anchor.set_axis_off()
    anchor.text(
        -0.16, 1.20, letter, transform=anchor.transAxes, fontsize=19,
        fontweight="bold", ha="left", va="top", clip_on=False,
    )
    anchor.text(
        0.5, 1.20, f"{sex} · {title}", transform=anchor.transAxes,
        fontsize=11.0, fontweight="bold", ha="center", va="top", clip_on=False,
    )
    anchor.text(
        0.5, -0.235, key, transform=anchor.transAxes, fontsize=6.85,
        color="#454545", ha="center", va="top", clip_on=False,
    )
    anchor.text(
        0.5, -0.14, "Study day", transform=anchor.transAxes, fontsize=10.0,
        ha="center", va="top", clip_on=False,
    )
    for index, axis in enumerate(axes):
        axis.set_xticks([1, 3, 6])
        axis.grid(axis="y", color="#dddddd", linewidth=0.65)
        axis.tick_params(axis="both", labelsize=7.9)
        if index:
            axis.tick_params(axis="y", labelleft=False)
            axis.spines["left"].set_color("#b5b5b5")
            axis.spines["left"].set_linewidth(0.75)
    axes[0].set_ylabel(ylabel, fontsize=9.1)


def draw_consumption_trajectory(
    axes: list[plt.Axes],
    anchor: plt.Axes,
    data: pd.DataFrame,
    sex: str,
    letter: str,
    seed: int,
) -> pd.DataFrame:
    value_column = "cumulative_cage_volume_removed_ml"
    subset_sex = data.loc[data["Sex"].eq(sex)].copy()
    counts = condition_counts(subset_sex, "E")
    summaries: list[pd.DataFrame] = []
    for index, (condition, axis) in enumerate(zip(CONDITIONS, axes)):
        subset = subset_sex.loc[subset_sex["Treatment"].eq(condition)].copy()
        for _, cage in subset.groupby("E", sort=True):
            cage = cage.sort_values("study_day")
            axis.plot(
                cage["study_day"], cage[value_column], color=COLORS[condition],
                alpha=0.38, linewidth=1.0, marker="o", markersize=3.0,
                markerfacecolor="white", markeredgewidth=0.7, zorder=2,
            )
        summary = trajectory_summary(
            subset, value_column, sex, condition, letter, seed + 100 * index
        )
        summaries.append(summary)
        x = summary["study_day"].to_numpy(float)
        mean = summary["cage_balanced_mean"].to_numpy(float)
        if bool(summary["ci_estimable"].iloc[0]):
            axis.fill_between(
                x,
                summary["ci95_low"].to_numpy(float),
                summary["ci95_high"].to_numpy(float),
                color=COLORS[condition], alpha=0.20, linewidth=0, zorder=1,
            )
        else:
            axis.text(
                0.5, 0.93, "95% CI not estimable\n(one cage)",
                transform=axis.transAxes, ha="center", va="top", fontsize=6.9,
                color="#555555", fontweight="bold",
            )
        axis.plot(
            x, mean, color=COLORS[condition], linewidth=2.5, marker="o",
            markersize=5.0, markeredgecolor="white", markeredgewidth=0.75, zorder=5,
        )
        cage_count = counts.get(condition, 0)
        axis.set_title(
            f"{condition}\n$n$={cage_count} cage{'s' if cage_count != 1 else ''}",
            color=COLORS[condition], fontsize=8.7, fontweight="bold", pad=4,
        )
    style_trajectory(
        axes, anchor, sex, "Cumulative consumption",
        "thin: cages   |   thick/band: equal-cage mean/95% CI",
        "Consumption from Day 1 (mL/cage)", letter,
    )
    return pd.concat(summaries, ignore_index=True)


def draw_weight_trajectory(
    axes: list[plt.Axes],
    anchor: plt.Axes,
    mouse_data: pd.DataFrame,
    cage_day_data: pd.DataFrame,
    sex: str,
    letter: str,
    seed: int,
) -> pd.DataFrame:
    value_column = "cage_mean_body_weight_change_pct"
    mice_sex = mouse_data.loc[mouse_data["Sex"].eq(sex)].copy()
    cages_sex = cage_day_data.loc[cage_day_data["Sex"].eq(sex)].copy()
    mouse_counts = condition_counts(mice_sex, "subject_id")
    cage_counts = condition_counts(cages_sex, "E")
    summaries: list[pd.DataFrame] = []
    for index, (condition, axis) in enumerate(zip(CONDITIONS, axes)):
        mouse_subset = mice_sex.loc[mice_sex["Treatment"].eq(condition)].copy()
        cage_subset = cages_sex.loc[cages_sex["Treatment"].eq(condition)].copy()
        animals = list(mouse_subset.groupby("subject_id", sort=True))
        offsets = np.linspace(-0.08, 0.08, max(len(animals), 1))
        for offset, (_, animal) in zip(offsets, animals):
            animal = animal.sort_values("study_day")
            axis.scatter(
                animal["study_day"].to_numpy(float) + offset,
                animal["body_weight_change_pct"], s=15,
                facecolor=COLORS[condition], edgecolor="none", alpha=0.25, zorder=2,
            )
        for _, cage in cage_subset.groupby("E", sort=True):
            cage = cage.sort_values("study_day")
            axis.plot(
                cage["study_day"], cage[value_column], color=COLORS[condition],
                alpha=0.38, linewidth=1.0, marker="o", markersize=3.0,
                markerfacecolor="white", markeredgewidth=0.7, zorder=3,
            )
        summary = trajectory_summary(
            cage_subset, value_column, sex, condition, letter, seed + 100 * index
        )
        summaries.append(summary)
        x = summary["study_day"].to_numpy(float)
        mean = summary["cage_balanced_mean"].to_numpy(float)
        if bool(summary["ci_estimable"].iloc[0]):
            axis.fill_between(
                x,
                summary["ci95_low"].to_numpy(float),
                summary["ci95_high"].to_numpy(float),
                color=COLORS[condition], alpha=0.20, linewidth=0, zorder=1,
            )
        else:
            axis.text(
                0.5, 0.93, "95% CI not estimable\n(one cage)",
                transform=axis.transAxes, ha="center", va="top", fontsize=6.9,
                color="#555555", fontweight="bold",
            )
        axis.plot(
            x, mean, color=COLORS[condition], linewidth=2.5, marker="o",
            markersize=5.0, markeredgecolor="white", markeredgewidth=0.75, zorder=5,
        )
        axis.axhline(0, color="#777777", linewidth=0.8, linestyle="--", zorder=0)
        axis.set_title(
            f"{condition}\n{mouse_counts.get(condition, 0)} mice / "
            f"{cage_counts.get(condition, 0)} cages",
            color=COLORS[condition], fontsize=8.0, fontweight="bold", pad=4,
        )
    finite = pd.to_numeric(cages_sex[value_column], errors="coerce").dropna().to_numpy(float)
    if finite.size:
        span = max(float(np.max(finite) - np.min(finite)), 1.0)
        axes[0].set_ylim(float(np.min(finite)) - 0.08 * span, float(np.max(finite)) + 0.15 * span)
    style_trajectory(
        axes, anchor, sex, "Cage-balanced body-weight change",
        "pale: mice   |   thin: cage-day   |   thick/band: equal-cage mean/95% CI",
        "Change from Day 1 (%)", letter,
    )
    return pd.concat(summaries, ignore_index=True)


def draw_endpoint(
    axis: plt.Axes,
    data: pd.DataFrame,
    value_column: str,
    statistics: pd.DataFrame,
    sex: str,
    outcome: str,
    title: str,
    ylabel: str,
    letter: str,
    seed: int,
) -> None:
    subset = data.loc[data["Sex"].eq(sex)].copy()
    if subset["E"].duplicated().any():
        raise ValueError(f"Panel {letter} endpoint must contain one value per cage")
    stats = statistics.loc[
        statistics["Sex"].eq(sex) & statistics["outcome"].eq(outcome)
    ].copy()
    descriptions = {
        str(row["comparison"]): row
        for _, row in stats.loc[stats["test"].eq("descriptive cage summary")].iterrows()
    }
    rng = np.random.default_rng(seed)
    for index, condition in enumerate(CONDITIONS):
        values = pd.to_numeric(
            subset.loc[subset["Treatment"].eq(condition), value_column], errors="coerce"
        ).dropna().to_numpy(float)
        jitter = rng.uniform(-0.10, 0.10, values.size)
        axis.scatter(
            np.full(values.size, index) + jitter, values, s=39, marker="o",
            facecolor=COLORS[condition], edgecolor="black", linewidth=0.62,
            alpha=0.93, zorder=4,
        )
        row = descriptions[condition]
        mean = float(row["estimate"])
        low = float(row["ci95_low"])
        high = float(row["ci95_high"])
        if np.isfinite(low) and np.isfinite(high):
            axis.errorbar(
                index, mean, yerr=[[mean - low], [high - mean]], fmt="D", markersize=6.1,
                color="black", markerfacecolor="white", markeredgewidth=1.25,
                linewidth=1.5, capsize=3.5, zorder=6,
            )
        else:
            axis.plot(
                index, mean, marker="D", markersize=6.1, color="black",
                markerfacecolor="white", markeredgewidth=1.25, zorder=6,
            )
            axis.annotate(
                "95% CI\nnot estimable", xy=(index, mean), xytext=(10, 8),
                textcoords="offset points", fontsize=6.7, color="#4d4d4d",
                fontweight="bold", ha="left", va="bottom", zorder=8,
            )

    omnibus = stats.loc[stats["test"].eq(OMNIBUS_TEST)].iloc[0]
    pairs = stats.loc[stats["test"].eq(PAIRWISE_TEST)]
    abbreviations = {"Water": "W", "Sucrose": "S", "Allulose": "A"}
    pair_lines: list[str] = []
    for row in pairs.itertuples(index=False):
        left, right = str(row.comparison).split("|")
        pair_lines.append(
            f"{abbreviations[left]}-{abbreviations[right]} Holm exact "
            f"{p_text(float(row.p_adjusted_holm))}"
        )
    values = pd.to_numeric(subset[value_column], errors="coerce").dropna().to_numpy(float)
    data_min = float(np.min(values))
    data_max = float(np.max(values))
    span = max(data_max - data_min, max(abs(data_max), 1.0) * 0.12)
    axis.set_ylim(data_min - 0.10 * span, data_max + 0.56 * span)
    box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 0.7}
    axis.text(
        0.03, 0.97, f"One-way ANOVA {p_text(float(omnibus.p_raw))}",
        transform=axis.transAxes, ha="left", va="top",
        fontsize=P_VALUE_GLOBAL_FONTSIZE_PT, fontweight=P_VALUE_FONTWEIGHT,
        bbox=box, zorder=12,
    )
    axis.text(
        0.03, 0.855, "\n".join(pair_lines),
        transform=axis.transAxes, ha="left", va="top",
        fontsize=P_VALUE_PAIRWISE_FONTSIZE_PT, fontweight=P_VALUE_FONTWEIGHT,
        linespacing=1.14, bbox=box, zorder=11,
    )
    counts = condition_counts(subset, "E")
    axis.set_title(f"{sex} · {title}", fontsize=10.6, fontweight="bold", pad=7)
    axis.set_ylabel(ylabel, fontsize=9.5)
    axis.set_xticks(np.arange(len(CONDITIONS)))
    axis.set_xticklabels([
        f"{condition}\n(n={counts.get(condition, 0)} "
        f"cage{'s' if counts.get(condition, 0) != 1 else ''})"
        for condition in CONDITIONS
    ], fontsize=8.2)
    axis.grid(axis="y", color="#dddddd", linewidth=0.7)
    panel_label(axis, letter)


def make_panel_axes(
    fig: plt.Figure, slot: object
) -> tuple[plt.Axes, list[plt.Axes]]:
    anchor = fig.add_subplot(slot)
    grid = slot.subgridspec(1, 3, wspace=0.13)
    axes: list[plt.Axes] = []
    for index in range(3):
        axes.append(fig.add_subplot(grid[0, index], sharey=axes[0] if axes else None))
    return anchor, axes


def export_panel_crops(
    fig: plt.Figure,
    panel_axes: dict[str, plt.Axes | list[plt.Axes]],
    output_dir: Path,
    dpi: int,
    suffix: str = "",
    crop_boxes: dict[str, Bbox] | None = None,
) -> dict[str, Bbox]:
    """Export only target axes from the combined live figure; never redraw panels."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    all_axes = list(fig.axes)
    all_figure_texts = list(fig.texts)
    result: dict[str, Bbox] = {}
    for letter, axes in panel_axes.items():
        axis_list = axes if isinstance(axes, list) else [axes]
        if crop_boxes is None:
            boxes = [axis.get_tightbbox(renderer) for axis in axis_list]
            boxes.extend(
                text.get_window_extent(renderer)
                for axis in axis_list for text in axis.texts if text.get_visible()
            )
            finite_boxes = [
                box for box in boxes if box is not None and np.isfinite(box.extents).all()
            ]
            pixel_box = Bbox.union(finite_boxes)
            inch_box = pixel_box.transformed(fig.dpi_scale_trans.inverted()).padded(0.055)
        else:
            inch_box = crop_boxes[letter]
        result[letter] = inch_box.frozen()
        stem = output_dir / f"{STEM_NAME}_Panel_{letter}{suffix}"
        axis_visibility = [(axis, axis.get_visible()) for axis in all_axes]
        text_visibility = [(artist, artist.get_visible()) for artist in all_figure_texts]
        try:
            for axis in all_axes:
                axis.set_visible(axis in axis_list)
            for artist in all_figure_texts:
                artist.set_visible(False)
            fig.savefig(
                stem.with_suffix(".pdf"), bbox_inches=inch_box, facecolor="white",
                metadata={"Title": f"Figure {FIGURE_NUMBER} live panel {letter}{suffix}"},
            )
            fig.savefig(
                stem.with_suffix(".png"), dpi=dpi, bbox_inches=inch_box, facecolor="white",
            )
        finally:
            for axis, visible in axis_visibility:
                axis.set_visible(visible)
            for artist, visible in text_visibility:
                artist.set_visible(visible)
    return result


def write_spanish_panel_legends(output_dir: Path) -> None:
    panels = {
        "A": "Trayectorias femeninas de consumo por jaula desde el día 1 hasta el día 6; las líneas gruesas y bandas resumen jaulas con el mismo peso.",
        "B": "Consumo por jaula femenina en el día 6; se muestra ANOVA de una vía y comparaciones exactas con ajuste de Holm.",
        "C": "Cambio de peso corporal femenino desde basal. Los puntos de ratón son descriptivos y la inferencia permanece al nivel de jaula.",
        "D": "Cambio medio de peso corporal por jaula femenina en el día 6; se muestra ANOVA de una vía y comparaciones exactas con ajuste de Holm.",
    }
    for letter, text in panels.items():
        (output_dir / f"{STEM_NAME}_Panel_{letter}_spanish_LEGEND.txt").write_text(
            f"Figura {FIGURE_NUMBER}, panel {letter}. {text}\n", encoding="utf-8"
        )
    sex_es = "masculino" if TARGET_SEX == "Male" else "femenino"
    cage_counts = (
        "Agua n=2, Sacarosa n=3 y Alulosa n=2"
        if TARGET_SEX == "Male"
        else "Agua n=1, Sacarosa n=2 y Alulosa n=2"
    )
    master = (
        f"Figura {FIGURE_NUMBER}. Experimento de una sola botella, grupo {sex_es} "
        "(análisis exploratorio).\n\n"
        + "\n".join(f"({letter}) {panels[letter]}" for letter in "ABCD")
        + f"\n\nTodas las jaulas de enero fueron verificadas como de un solo sexo; "
          f"los tamaños por tratamiento son {cage_counts}. La jaula es la unidad "
          "experimental y de inferencia; los puntos de ratones en los paneles de peso "
          "son descriptivos. Los valores globales mostrados usan ANOVA de una vía, las "
          "comparaciones exactas por pares se ajustan por Holm y el ómnibus exacto por "
          "etiquetas de jaula se conserva como sensibilidad. No se estimó una interacción "
          "tratamiento por sexo y las diferencias de significación entre análisis por "
          "sexo no deben interpretarse como una diferencia sexual.\n"
    )
    (output_dir / f"{STEM_NAME}_LEGEND_spanish.txt").write_text(
        master, encoding="utf-8"
    )


def make_plotted_values(
    consumption: pd.DataFrame,
    weight_mice: pd.DataFrame,
    weight_cage_day: pd.DataFrame,
    weight_endpoint: pd.DataFrame,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for sex in (TARGET_SEX,):
        c = consumption.loc[consumption["Sex"].eq(sex)].copy()
        c_endpoint = c.loc[c["study_day"].eq(6)].copy()
        w_mouse = weight_mice.loc[weight_mice["Sex"].eq(sex)].copy()
        w_cage = weight_cage_day.loc[weight_cage_day["Sex"].eq(sex)].copy()
        w_endpoint = weight_endpoint.loc[weight_endpoint["Sex"].eq(sex)].copy()
        a_letter = PANEL_MAP[(sex, "consumption_trajectory")]
        b_letter = PANEL_MAP[(sex, "consumption_endpoint")]
        c_letter = PANEL_MAP[(sex, "weight_trajectory")]
        d_letter = PANEL_MAP[(sex, "weight_endpoint")]
        frames.extend([
            c.assign(
                panel=a_letter, value_role="cage_trajectory", experimental_unit="cage",
                unit_id=c["E"], value_name="cumulative_cage_volume_removed_ml",
                plotted_value=c["cumulative_cage_volume_removed_ml"], subject_id=pd.NA,
            ),
            c_endpoint.assign(
                panel=b_letter, value_role="cage_endpoint", experimental_unit="cage",
                unit_id=c_endpoint["E"], value_name="cumulative_cage_volume_removed_ml",
                plotted_value=c_endpoint["cumulative_cage_volume_removed_ml"], subject_id=pd.NA,
            ),
            w_mouse.assign(
                panel=c_letter, value_role="mouse_trajectory_descriptive",
                experimental_unit="mouse_descriptive_only", unit_id=w_mouse["subject_id"],
                value_name="body_weight_change_pct",
                plotted_value=w_mouse["body_weight_change_pct"],
            ),
            w_cage.assign(
                panel=c_letter, value_role="cage_day_trajectory", experimental_unit="cage",
                unit_id=w_cage["E"], value_name="cage_mean_body_weight_change_pct",
                plotted_value=w_cage["cage_mean_body_weight_change_pct"], subject_id=pd.NA,
            ),
            w_endpoint.assign(
                panel=d_letter, value_role="cage_endpoint", experimental_unit="cage",
                unit_id=w_endpoint["E"], value_name="cage_mean_body_weight_change_pct",
                plotted_value=w_endpoint["cage_mean_body_weight_change_pct"], subject_id=pd.NA,
            ),
        ])
    columns = [
        "panel", "value_role", "experimental_unit", "unit_id", "E", "subject_id",
        "Treatment", "Sex", "study_day", "value_name", "plotted_value",
    ]
    normalized = []
    for frame in frames:
        for column in columns:
            if column not in frame.columns:
                frame[column] = pd.NA
        normalized.append(frame[columns])
    return pd.concat(normalized, ignore_index=True)


def write_legends(output_dir: Path, statistics: pd.DataFrame) -> None:
    def global_p(sex: str, outcome: str) -> float:
        return float(statistics.loc[
            statistics["Sex"].eq(sex)
            & statistics["outcome"].eq(outcome)
            & statistics["test"].eq(OMNIBUS_TEST),
            "p_raw",
        ].iloc[0])

    panels: dict[str, str] = {}
    for sex in (TARGET_SEX,):
        consumption_letter = PANEL_MAP[(sex, "consumption_trajectory")]
        consumption_endpoint_letter = PANEL_MAP[(sex, "consumption_endpoint")]
        weight_letter = PANEL_MAP[(sex, "weight_trajectory")]
        weight_endpoint_letter = PANEL_MAP[(sex, "weight_endpoint")]
        panels[consumption_letter] = (
            f"{sex} cage consumption from Day 1 through Day 6. "
            "Thin lines are complete cage trajectories; thick lines are equal-cage means. "
            "Bands are pointwise descriptive 95% intervals from whole-cage trajectory "
            f"resampling ({BOOTSTRAP_REPLICATES:,} draws). "
            + ("Water has one cage, so its 95% interval is not estimable." if sex == "Female" else "")
        )
        panels[consumption_endpoint_letter] = (
            f"{sex} Day-6 cage-consumption endpoint. Dots are cages; diamonds are cage means. "
            f"One-way ANOVA {p_text(global_p(sex, CONSUMPTION_OUTCOME))}; "
            "the three pairwise exact p values are Holm-adjusted within this sex and outcome. "
            + ("The Water mean has no interval because it is one cage." if sex == "Female" else "")
        )
        panels[weight_letter] = (
            f"{sex} body-weight change from each mouse's Day-1 baseline. Pale points are mice, "
            "thin lines are cage-day means, and thick lines are equal-cage means; inferential "
            "weight remains at cage level. Bands use whole-cage trajectory resampling. "
            + ("Water has one cage, so its 95% interval is not estimable." if sex == "Female" else "")
        )
        panels[weight_endpoint_letter] = (
            f"{sex} Day-6 cage-mean body-weight endpoint. Each dot is one cage. "
            f"One-way ANOVA {p_text(global_p(sex, WEIGHT_OUTCOME))}; "
            "the three pairwise exact p values are Holm-adjusted within this sex and outcome. "
            + ("The Water mean has no interval because it is one cage." if sex == "Female" else "")
        )

    master = (
        f"Figure {FIGURE_NUMBER}. {TARGET_SEX} single bottle experiment (exploratory).\n\n"
        + "\n".join(f"({letter}) {panels[letter]}" for letter in "ABCD")
        + "\n\nAll January cages were verified as single-sex in both source tables before the "
          f"cage-level bottle data were assigned to sex. {TARGET_SEX} cage N is "
          + ("Water 1, Sucrose 2, Allulose 2. " if TARGET_SEX == "Female" else
             "Water 2, Sucrose 3, Allulose 2. ")
          + "Displayed omnibus p values use ordinary one-way ANOVA; exact cage-label omnibus results are retained as sensitivity analyses. Pairwise exact p values are low-resolution. Mouse rows in body-weight "
          "panels are descriptive. No treatment-by-sex interaction was estimated or tested; "
          "differences in within-sex significance must not be interpreted as a sex difference. "
          "Treatment randomization was not documented, so the analyses are exploratory.\n"
    )
    (output_dir / f"{STEM_NAME}_LEGEND.txt").write_text(master, encoding="utf-8")
    (output_dir / f"{STEM_NAME}_caption.txt").write_text(master, encoding="utf-8")
    for letter in "ABCD":
        (output_dir / f"{STEM_NAME}_Panel_{letter}_LEGEND.txt").write_text(
            f"Figure {FIGURE_NUMBER} panel {letter}. {panels[letter]}\n", encoding="utf-8"
        )


def main() -> int:
    args = parse_args()
    args.analysis_dir = args.analysis_dir.resolve()
    final_output = Path(os.path.abspath(os.path.expanduser(args.output_dir)))
    staging_output = prepare_staging_output(final_output)
    panel_output = staging_output / "panels"
    legend_output = staging_output / "legends"
    source_output = staging_output / "source_data"
    provenance_output = staging_output / "provenance"
    paths = {
        "consumption_values": args.analysis_dir / "single_bottle_consumption_primary_values.csv",
        "consumption_endpoints": args.analysis_dir / "single_bottle_consumption_cage_endpoint_values.csv",
        "weight_mouse_values": args.analysis_dir / "single_bottle_body_weight_mouse_values.csv",
        "weight_cage_day_values": args.analysis_dir / "single_bottle_body_weight_cage_day_values.csv",
        "weight_cage_endpoints": args.analysis_dir / "single_bottle_body_weight_cage_endpoint_values.csv",
        "statistics": args.analysis_dir / "single_bottle_endpoint_statistics.csv",
        "cage_sex_audit": args.analysis_dir / "cage_sex_assignment_audit.csv",
        "analysis_sources": args.analysis_dir / "analysis_sources.csv",
        "analysis_manifest": args.analysis_dir / "run_manifest.json",
        "value_correspondence": args.analysis_dir / "figure1_value_correspondence.json",
        "statistical_audit": args.analysis_dir / f"FIGURE_{FIGURE_NUMBER}_STATISTICAL_AUDIT.md",
        "analysis_script": PAPER_ROOT / "scripts" / f"Fig{FIGURE_NUMBER}" / f"01_analyze_single_bottle_{TARGET_SEX.lower()}.py",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    consumption = pd.read_csv(paths["consumption_values"])
    consumption_endpoint = pd.read_csv(paths["consumption_endpoints"])
    weight_mice = pd.read_csv(paths["weight_mouse_values"])
    weight_cage_day = pd.read_csv(paths["weight_cage_day_values"])
    weight_endpoint = pd.read_csv(paths["weight_cage_endpoints"])
    statistics = pd.read_csv(paths["statistics"])
    cage_audit = pd.read_csv(paths["cage_sex_audit"])
    if not cage_audit["eligible_for_sex_stratified_cage_analysis"].all():
        raise ValueError("Renderer refused: not every January cage passed the single-sex gate")
    expected_counts = (
        {"Water": 1, "Sucrose": 2, "Allulose": 2}
        if TARGET_SEX == "Female" else {"Water": 2, "Sucrose": 3, "Allulose": 2}
    )
    observed_sexes = set(consumption["Sex"].dropna().astype(str))
    if observed_sexes != {TARGET_SEX}:
        raise ValueError(f"Expected only {TARGET_SEX} values, found {observed_sexes}")
    observed = condition_counts(consumption, "E")
    if observed != expected_counts:
        raise ValueError(f"Unexpected {TARGET_SEX} cage counts: {observed}")
    if TARGET_SEX == "Female":
        female_water = statistics.loc[
            statistics["Sex"].eq("Female")
            & statistics["outcome"].eq(CONSUMPTION_OUTCOME)
            & statistics["test"].eq("descriptive cage summary")
            & statistics["comparison"].eq("Water")
        ].iloc[0]
        if pd.notna(female_water["ci95_low"]) or pd.notna(female_water["ci95_high"]):
            raise ValueError("Female Water one-cage confidence interval must be not estimable")

    fig = plt.figure(figsize=FIGURE_SIZE_INCHES, facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.095, right=0.985, bottom=0.085, top=0.875,
        wspace=0.28, hspace=0.54,
    )
    axes_by_panel: dict[str, plt.Axes | list[plt.Axes]] = {}
    intervals: list[pd.DataFrame] = []

    for sex_index, sex in enumerate((TARGET_SEX,)):
        row_base = 0
        consumption_letter = PANEL_MAP[(sex, "consumption_trajectory")]
        endpoint_letter = PANEL_MAP[(sex, "consumption_endpoint")]
        weight_letter = PANEL_MAP[(sex, "weight_trajectory")]
        weight_endpoint_letter = PANEL_MAP[(sex, "weight_endpoint")]

        consumption_anchor, consumption_axes = make_panel_axes(fig, outer[row_base, 0])
        consumption_endpoint_axis = fig.add_subplot(outer[row_base, 1])
        weight_anchor, weight_axes = make_panel_axes(fig, outer[row_base + 1, 0])
        weight_endpoint_axis = fig.add_subplot(outer[row_base + 1, 1])

        intervals.append(draw_consumption_trajectory(
            consumption_axes, consumption_anchor, consumption, sex,
            consumption_letter, args.seed + sex_index * 10_000 + 1_000,
        ))
        draw_endpoint(
            consumption_endpoint_axis, consumption_endpoint,
            "cumulative_cage_volume_removed_ml", statistics, sex, CONSUMPTION_OUTCOME,
            "Day-6 consumption endpoint", "Consumption (mL/cage)",
            endpoint_letter, args.seed + sex_index * 10_000 + 10,
        )
        intervals.append(draw_weight_trajectory(
            weight_axes, weight_anchor, weight_mice, weight_cage_day, sex,
            weight_letter, args.seed + sex_index * 10_000 + 2_000,
        ))
        draw_endpoint(
            weight_endpoint_axis, weight_endpoint,
            "cage_mean_body_weight_change_pct", statistics, sex, WEIGHT_OUTCOME,
            "Day-6 body-weight endpoint", "Cage-mean change from baseline (%)",
            weight_endpoint_letter, args.seed + sex_index * 10_000 + 20,
        )
        axes_by_panel[consumption_letter] = [consumption_anchor, *consumption_axes]
        axes_by_panel[endpoint_letter] = consumption_endpoint_axis
        axes_by_panel[weight_letter] = [weight_anchor, *weight_axes]
        axes_by_panel[weight_endpoint_letter] = weight_endpoint_axis

    fig.suptitle(
        f"Single bottle experiment · {TARGET_SEX.lower()} cages",
        fontsize=15.2, fontweight="bold", y=0.978,
    )
    fig.text(
        0.5, 0.951,
        "Exploratory within-sex cage analyses · no treatment-by-sex interaction test",
        ha="center", va="top", fontsize=10.0, color="#4d4d4d", fontweight="bold",
    )
    stem = staging_output / STEM_NAME
    fig.savefig(
        stem.with_suffix(".pdf"), facecolor="white",
        metadata={"Title": f"Single bottle experiment: {TARGET_SEX.lower()} cages"},
    )
    fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, facecolor="white")
    english_crop_boxes = export_panel_crops(fig, axes_by_panel, panel_output, args.dpi)
    translation_receipt = translate_figure_texts_to_spanish(fig)
    export_panel_crops(
        fig, axes_by_panel, panel_output, args.dpi,
        suffix="_spanish", crop_boxes=english_crop_boxes,
    )
    plt.close(fig)

    interval_values = pd.concat(intervals, ignore_index=True)
    interval_path = source_output / f"{STEM_NAME}_trajectory_confidence_intervals.csv"
    interval_values.to_csv(interval_path, index=False)
    if TARGET_SEX == "Female":
        female_water_intervals = interval_values.loc[
            interval_values["Treatment"].eq("Water")
        ]
        if female_water_intervals["ci_estimable"].any():
            raise RuntimeError("Female Water trajectory CI was incorrectly estimated")

    plotted_values = make_plotted_values(
        consumption, weight_mice, weight_cage_day, weight_endpoint
    )
    plotted_values.to_csv(source_output / f"{STEM_NAME}_plotted_values.csv", index=False)
    stats_for_figure = statistics.copy()
    stats_for_figure["panel"] = stats_for_figure.apply(
        lambda row: PANEL_MAP[(
            str(row["Sex"]),
            "consumption_endpoint" if str(row["outcome"]) == CONSUMPTION_OUTCOME
            else "weight_endpoint",
        )],
        axis=1,
    )
    stats_for_figure.to_csv(source_output / f"{STEM_NAME}_statistics.csv", index=False)
    cage_audit.to_csv(source_output / f"{STEM_NAME}_cage_sex_assignment_audit.csv", index=False)
    shutil.copyfile(
        paths["statistical_audit"],
        provenance_output / f"FIGURE_{FIGURE_NUMBER}_AUDIT.md",
    )
    shutil.copyfile(
        paths["value_correspondence"],
        provenance_output / f"{STEM_NAME}_value_correspondence.json",
    )
    write_legends(legend_output, statistics)
    write_spanish_panel_legends(legend_output)
    translation_receipt_path = provenance_output / f"{STEM_NAME}_spanish_translation_receipt.json"
    translation_receipt_path.write_text(
        json.dumps(translation_receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_local_readme(provenance_output)

    with (source_output / f"{STEM_NAME}_sources.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(("role", "path", "sha256"))
        for role, path in paths.items():
            writer.writerow((role, str(path), sha256(path)))
        writer.writerow(("renderer", str(Path(__file__).resolve()), sha256(Path(__file__).resolve())))

    expected_graphics = [
        f"{STEM_NAME}.pdf",
        f"{STEM_NAME}.png",
        *[
            f"panels/{STEM_NAME}_Panel_{letter}.{extension}"
            for letter in "ABCD" for extension in ("pdf", "png")
        ],
        *[
            f"panels/{STEM_NAME}_Panel_{letter}_spanish.{extension}"
            for letter in "ABCD" for extension in ("pdf", "png")
        ],
    ]
    missing = [name for name in expected_graphics if not (staging_output / name).is_file()]
    if missing:
        raise RuntimeError(f"Missing rendered graphics: {missing}")
    manifest = {
        "figure": STEM_NAME,
        "analysis_status": "exploratory_low_resolution",
        "analysis_directory": str(args.analysis_dir),
        "output_directory": str(final_output),
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha256(Path(__file__).resolve()),
        "seed": args.seed,
        "dpi": args.dpi,
        "figure_size_inches": list(FIGURE_SIZE_INCHES),
        "figure_orientation": "portrait; two rows by two macro columns",
        "subpanels": list("ABCD"),
        "formats": ["pdf", "png"],
        "language_contract": {
            "combined_figure": "English only",
            "standalone_panels": ["English", "Spanish"],
            "spanish_translation_receipt": str(translation_receipt_path.relative_to(staging_output)),
            "spanish_translation_receipt_sha256": sha256(translation_receipt_path),
            "crop_geometry": "Spanish and English pairs use identical fixed crop boxes",
        },
        "panel_export_method": "exact crops of live axes in the combined rendered figure",
        "panel_assignment": {
            "A-D": f"{TARGET_SEX}: consumption trajectory/endpoint; body-weight trajectory/endpoint",
        },
        "experimental_unit": "cage",
        "sex": TARGET_SEX,
        "cage_counts": expected_counts,
        "statistics": {
            "endpoint": "Day 6",
            "omnibus": "ordinary equal-variance one-way ANOVA within sex/outcome",
            "omnibus_sensitivity": "exact exhaustive cage-label allocation of the same F statistic",
            "pairwise": "exact exhaustive cage-label permutation; Holm within sex/outcome",
            "interaction_test": "not performed",
            "interpretation": "exploratory; exact p values have low resolution",
        },
        "trajectory_confidence_intervals": {
            "level": 0.95,
            "type": "pointwise descriptive percentile interval",
            "resampling_unit": "whole cage trajectory",
            "same_cage_draw_across_days": True,
            "bootstrap_replicates_when_estimable": BOOTSTRAP_REPLICATES,
            "group_limit": (
                "Female Water not estimable (one cage); no band drawn"
                if TARGET_SEX == "Female" else "all groups estimable"
            ),
            "values_file": str(interval_path.relative_to(staging_output)),
            "values_sha256": sha256(interval_path),
        },
        "p_value_annotation_style": {
            "layout": "vertical multiline list; text itself is not rotated",
            "omnibus_font_size_pt": P_VALUE_GLOBAL_FONTSIZE_PT,
            "pairwise_font_size_pt": P_VALUE_PAIRWISE_FONTSIZE_PT,
            "fontweight": P_VALUE_FONTWEIGHT,
        },
        "graphics": expected_graphics,
        "graphics_sha256": {
            name: sha256(staging_output / name) for name in expected_graphics
        },
    }
    (provenance_output / f"{STEM_NAME}_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    publish_staging_output(staging_output, final_output)
    print(f"[FIGURE] {final_output / stem.with_suffix('.pdf').name}")
    print(f"[FIGURE] {final_output / stem.with_suffix('.png').name}")
    print("[PANELS] live A-D crops exported as PDF and PNG")
    print(
        f"[P-VALUE FONT] omnibus={P_VALUE_GLOBAL_FONTSIZE_PT} pt; "
        f"pairwise={P_VALUE_PAIRWISE_FONTSIZE_PT} pt; bold; vertical stack"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
