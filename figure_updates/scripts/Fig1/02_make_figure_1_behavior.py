#!/usr/bin/env python3
"""Render Figure 1 and exact PDF/PNG crops of its live A-D axes."""

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

# Stabilize PDF CreationDate so command-level reruns are byte-reproducible.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.4,
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
if PAPER_ROOT is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
PROJECT_ROOT = PAPER_ROOT.parent
sys.path.insert(0, str(PAPER_ROOT / "scripts" / "shared"))
from spanish_panel_text import translate_figure_texts_to_spanish
FIGURE_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig1"
DEFAULT_ANALYSIS = FIGURE_ANALYSIS_ROOT / "results"
DEFAULT_OUTPUT = PAPER_ROOT / "Fig1"
CONDITIONS = ("Water", "Sucrose", "Allulose")
COLORS = {"Water": "#4C9BD6", "Sucrose": "#D95F5F", "Allulose": "#39A96B"}
OMNIBUS_TEST = "One_way_ANOVA"
EXACT_OMNIBUS_TEST = "exact cage-label permutation omnibus F"
P_VALUE_GLOBAL_FONTSIZE_PT = 10.1
P_VALUE_PAIRWISE_FONTSIZE_PT = 9.2
P_VALUE_FONTWEIGHT = "bold"
BOOTSTRAP_REPLICATES = 20_000
FIGURE_SIZE_INCHES = (11.2, 12.2)
OUTPUT_SUBDIRS = ("panels", "legends", "source_data", "provenance")
GENERATED_TREE_MARKER = ".paper_generated_tree.json"
OUTPUT_PRODUCER = "Paper/scripts/Fig1/02_make_figure_1_behavior.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20260814)
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
    text = """FIGURE 1 - SINGLE BOTTLE EXPERIMENT

The combined publication PDF and PNG are at the parent Fig1 directory's root. Exact live
panel crops are in panels/, legends and the caption are in legends/, plotted
values/statistics/confidence intervals/source hashes are in source_data/, and
the deterministic manifest, statistical audit, and this README are in
provenance/. Every publication artifact is a real file; symlinks are forbidden.

Reproduce from the Paper directory:

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig1/01_analyze_behavior_experiment1.py
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig1/02_make_figure_1_behavior.py

The renderer builds a complete sibling staging tree and replaces only an empty or
previously marker-owned tree. It refuses symlinks, modified files, and untracked
content, so the recovered live Fig1 directory must not be used as --output-dir;
render to a fresh directory and promote reviewed leaf artifacts.
"""
    (output_dir / "README.txt").write_text(text, encoding="utf-8")


def p_text(value: float) -> str:
    if not np.isfinite(value):
        return "p = NA"
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.3f}"


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.13, 1.08, letter, transform=ax.transAxes, fontsize=20,
        fontweight="bold", ha="left", va="top", clip_on=False,
    )


def condition_counts(data: pd.DataFrame, unit: str) -> dict[str, int]:
    return (
        data[[unit, "Treatment"]].drop_duplicates()
        .groupby("Treatment")[unit].nunique().to_dict()
    )


def descriptive_lookup(stats: pd.DataFrame) -> dict[str, pd.Series]:
    subset = stats.loc[stats["test"].eq("descriptive cage summary")]
    return {str(row["comparison"]): row for _, row in subset.iterrows()}


def whole_cage_trajectory_bootstrap(
    data: pd.DataFrame,
    value_column: str,
    seed: int,
    panel: str,
    condition: str,
) -> pd.DataFrame:
    """Pointwise percentile CIs after resampling complete cage trajectories.

    A bootstrap replicate draws cage IDs with replacement once and uses that
    same draw at every study day. Thus longitudinal dependence is retained and
    every cage has equal weight. These intervals are descriptive because each
    treatment contains only three to five cages.
    """
    required = {"E", "study_day", value_column}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Panel {panel} bootstrap data are missing {sorted(missing)}")
    if data.duplicated(["E", "study_day"]).any():
        raise ValueError(f"Panel {panel} needs one value per cage and study day")
    pivot = (
        data.pivot(index="E", columns="study_day", values=value_column)
        .sort_index(axis=0).sort_index(axis=1)
    )
    if pivot.empty or pivot.isna().any().any():
        raise ValueError(f"Panel {panel} has an incomplete cage trajectory for {condition}")
    values = pivot.to_numpy(float)
    rng = np.random.default_rng(seed)
    resampled_indices = rng.integers(0, values.shape[0], size=(BOOTSTRAP_REPLICATES, values.shape[0]))
    bootstrap_means = values[resampled_indices, :].mean(axis=1)
    low, high = np.percentile(bootstrap_means, [2.5, 97.5], axis=0)
    return pd.DataFrame({
        "panel": panel,
        "Treatment": condition,
        "study_day": pivot.columns.to_numpy(float),
        "n_cages": values.shape[0],
        "cage_balanced_mean": values.mean(axis=0),
        "ci95_low": low,
        "ci95_high": high,
        "ci_method": "pointwise percentile cluster bootstrap of whole cage trajectories",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": seed,
    })


def style_trajectory_facets(
    axes: list[plt.Axes],
    anchor: plt.Axes,
    title: str,
    ylabel: str,
    letter: str,
) -> None:
    anchor.set_axis_off()
    anchor.text(
        -0.16, 1.21, letter, transform=anchor.transAxes, fontsize=20,
        fontweight="bold", ha="left", va="top", clip_on=False,
    )
    anchor.text(
        0.5, 1.21, title, transform=anchor.transAxes, fontsize=11.2,
        fontweight="bold", ha="center", va="top", clip_on=False,
    )
    anchor.text(
        0.5, -0.145, "Study day", transform=anchor.transAxes, fontsize=10.5,
        ha="center", va="top", clip_on=False,
    )
    for index, ax in enumerate(axes):
        ax.set_xticks([1, 3, 6])
        ax.grid(axis="y", color="#dddddd", linewidth=0.65)
        ax.tick_params(axis="both", labelsize=8.2)
        if index:
            ax.tick_params(axis="y", labelleft=False)
            ax.spines["left"].set_color("#b5b5b5")
            ax.spines["left"].set_linewidth(0.75)
    axes[0].set_ylabel(ylabel, fontsize=9.5)


def draw_consumption_trajectory(
    axes: list[plt.Axes], anchor: plt.Axes, data: pd.DataFrame, seed: int,
) -> pd.DataFrame:
    value_column = "cumulative_cage_volume_removed_ml"
    counts = condition_counts(data, "E")
    summaries: list[pd.DataFrame] = []
    for condition_index, (condition, ax) in enumerate(zip(CONDITIONS, axes)):
        subset = data.loc[data["Treatment"].eq(condition)].copy()
        for _, cage in subset.groupby("E", sort=True):
            cage = cage.sort_values("study_day")
            ax.plot(
                cage["study_day"], cage[value_column], color=COLORS[condition],
                alpha=0.34, linewidth=0.95, marker="o", markersize=3.2,
                markerfacecolor="white", markeredgewidth=0.7, zorder=2,
            )
        summary = whole_cage_trajectory_bootstrap(
            subset, value_column, seed + 100 * condition_index,
            panel="A", condition=condition,
        )
        summaries.append(summary)
        x = summary["study_day"].to_numpy(float)
        mean = summary["cage_balanced_mean"].to_numpy(float)
        ax.fill_between(
            x, summary["ci95_low"], summary["ci95_high"],
            color=COLORS[condition], alpha=0.20, linewidth=0, zorder=1,
        )
        ax.plot(
            x, mean, color=COLORS[condition], linewidth=2.7, marker="o",
            markersize=5.4, markeredgecolor="white", markeredgewidth=0.75, zorder=5,
        )
        ax.set_title(
            f"{condition}\n$n$={counts.get(condition, 0)} cages",
            color=COLORS[condition], fontsize=9.1, fontweight="bold", pad=5,
        )
    style_trajectory_facets(
        axes, anchor, "Cumulative consumption",
        "Consumption from Day 1 (mL/cage)", "A",
    )
    return pd.concat(summaries, ignore_index=True)


def draw_cage_endpoint(
    ax: plt.Axes,
    data: pd.DataFrame,
    value_column: str,
    stats: pd.DataFrame,
    title: str,
    ylabel: str,
    letter: str,
    seed: int,
) -> None:
    if data["E"].duplicated().any():
        raise ValueError(f"Panel {letter} endpoint must have one value per cage")
    rng = np.random.default_rng(seed)
    descriptions = descriptive_lookup(stats)
    for index, condition in enumerate(CONDITIONS):
        subset = data.loc[data["Treatment"].eq(condition)]
        values = pd.to_numeric(subset[value_column], errors="coerce").dropna().to_numpy(float)
        jitter = rng.uniform(-0.11, 0.11, values.size)
        ax.scatter(
            np.full(values.size, index) + jitter, values, s=42, marker="o",
            facecolor=COLORS[condition], edgecolor="black", linewidth=0.65, alpha=0.92, zorder=4,
        )
        row = descriptions[condition]
        mean = float(row["estimate"])
        low, high = float(row["ci95_low"]), float(row["ci95_high"])
        ax.errorbar(
            index, mean, yerr=[[mean - low], [high - mean]], fmt="D", markersize=6.5,
            color="black", markerfacecolor="white", markeredgewidth=1.3,
            linewidth=1.6, capsize=4, zorder=6,
        )

    omnibus = stats.loc[stats["test"].eq(OMNIBUS_TEST)].iloc[0]
    pairs = stats.loc[
        stats["comparison"].isin(["Water|Sucrose", "Water|Allulose", "Sucrose|Allulose"])
    ]
    abbreviations = {"Water": "W", "Sucrose": "S", "Allulose": "A"}
    pair_text: list[str] = []
    for row in pairs.itertuples(index=False):
        left, right = str(row.comparison).split("|")
        pair_text.append(
            f"{abbreviations[left]}-{abbreviations[right]} exact, Holm-adjusted "
            f"{p_text(float(row.p_adjusted_holm))}"
        )
    finite_values = pd.to_numeric(data[value_column], errors="coerce").dropna().to_numpy(float)
    if finite_values.size:
        data_min = float(np.min(finite_values))
        data_max = float(np.max(finite_values))
        data_span = max(data_max - data_min, max(abs(data_max), 1.0) * 0.10)
        ax.set_ylim(data_min - 0.08 * data_span, data_max + 0.42 * data_span)
    annotation_box = {
        "facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.8,
    }
    ax.text(
        0.03, 0.97, f"One-way ANOVA {p_text(float(omnibus.p_raw))}",
        transform=ax.transAxes, ha="left", va="top", fontsize=P_VALUE_GLOBAL_FONTSIZE_PT,
        fontweight=P_VALUE_FONTWEIGHT, bbox=annotation_box, zorder=12,
    )
    ax.text(
        0.03, 0.855, "\n".join(pair_text),
        transform=ax.transAxes, ha="left", va="top", fontsize=P_VALUE_PAIRWISE_FONTSIZE_PT,
        fontweight=P_VALUE_FONTWEIGHT, linespacing=1.15, bbox=annotation_box, zorder=11,
    )
    counts = condition_counts(data, "E")
    ax.set_title(title, fontweight="bold")
    ax.set_ylabel(ylabel)
    ax.set_xticks(np.arange(len(CONDITIONS)))
    ax.set_xticklabels([f"{condition}\n(n={counts.get(condition, 0)} cages)" for condition in CONDITIONS])
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    panel_label(ax, letter)


def draw_weight_trajectory(
    axes: list[plt.Axes], anchor: plt.Axes, data: pd.DataFrame, seed: int,
) -> pd.DataFrame:
    value_column = "body_weight_change_pct"
    mouse_counts = condition_counts(data, "subject_id")
    cage_counts = condition_counts(data, "E")
    summaries: list[pd.DataFrame] = []
    for condition_index, (condition, ax) in enumerate(zip(CONDITIONS, axes)):
        subset = data.loc[data["Treatment"].eq(condition)].copy()
        animals = list(subset.groupby("subject_id", sort=True))
        offsets = np.linspace(-0.085, 0.085, max(len(animals), 1))
        for offset, (_, animal) in zip(offsets, animals):
            animal = animal.sort_values("study_day")
            ax.scatter(
                animal["study_day"].to_numpy(float) + offset, animal[value_column], s=17,
                marker="o", facecolor=COLORS[condition], edgecolor="none",
                alpha=0.25, zorder=2,
            )
        cage_day = subset.groupby(["E", "study_day"], as_index=False)[value_column].mean()
        for _, cage in cage_day.groupby("E", sort=True):
            cage = cage.sort_values("study_day")
            ax.plot(
                cage["study_day"], cage[value_column], color=COLORS[condition],
                alpha=0.35, linewidth=0.95, marker="o", markersize=3.2,
                markerfacecolor="white", markeredgewidth=0.7, zorder=3,
            )
        summary = whole_cage_trajectory_bootstrap(
            cage_day, value_column, seed + 100 * condition_index,
            panel="C", condition=condition,
        )
        summaries.append(summary)
        x = summary["study_day"].to_numpy(float)
        mean = summary["cage_balanced_mean"].to_numpy(float)
        ax.fill_between(
            x, summary["ci95_low"], summary["ci95_high"],
            color=COLORS[condition], alpha=0.20, linewidth=0, zorder=1,
        )
        ax.plot(
            x, mean, color=COLORS[condition], linewidth=2.7, marker="o",
            markersize=5.4, markeredgecolor="white", markeredgewidth=0.75, zorder=5,
        )
        ax.axhline(0, color="#777777", linewidth=0.85, linestyle="--", zorder=0)
        ax.set_title(
            f"{condition}\n{mouse_counts.get(condition, 0)} mice / "
            f"{cage_counts.get(condition, 0)} cages",
            color=COLORS[condition], fontsize=8.6, fontweight="bold", pad=5,
        )
    finite_values = pd.to_numeric(data[value_column], errors="coerce").dropna().to_numpy(float)
    if finite_values.size:
        data_min = float(np.min(finite_values))
        data_max = float(np.max(finite_values))
        data_span = max(data_max - data_min, 1.0)
        axes[0].set_ylim(data_min - 0.06 * data_span, data_max + 0.13 * data_span)
    style_trajectory_facets(
        axes, anchor, "Cage-balanced body-weight change",
        "Change from Day 1 (%)", "C",
    )
    return pd.concat(summaries, ignore_index=True)


def export_panel_crops(
    fig: plt.Figure,
    axes: dict[str, plt.Axes | list[plt.Axes]],
    output_dir: Path,
    dpi: int,
    suffix: str = "",
    crop_boxes: dict[str, Bbox] | None = None,
) -> dict[str, Bbox]:
    """Crop only target axes from the combined live figure; panels are not redrawn."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    all_axes = list(fig.axes)
    all_figure_texts = list(fig.texts)
    result: dict[str, Bbox] = {}
    for letter, panel_axes in axes.items():
        axis_list = panel_axes if isinstance(panel_axes, list) else [panel_axes]
        if crop_boxes is None:
            boxes = [axis.get_tightbbox(renderer) for axis in axis_list]
            boxes.extend(
                text.get_window_extent(renderer)
                for axis in axis_list for text in axis.texts if text.get_visible()
            )
            pixel_box = Bbox.union([
                box for box in boxes if box is not None and np.isfinite(box.extents).all()
            ])
            inch_box = pixel_box.transformed(fig.dpi_scale_trans.inverted()).padded(0.06)
        else:
            inch_box = crop_boxes[letter]
        result[letter] = inch_box.frozen()
        stem = output_dir / f"Figure_1_behavior_experiment_1_Panel_{letter}{suffix}"
        axis_visibility = [(axis, axis.get_visible()) for axis in all_axes]
        text_visibility = [(artist, artist.get_visible()) for artist in all_figure_texts]
        try:
            for axis in all_axes:
                axis.set_visible(axis in axis_list)
            for artist in all_figure_texts:
                artist.set_visible(False)
            fig.savefig(
                stem.with_suffix(".pdf"), bbox_inches=inch_box, facecolor="white",
                metadata={"Title": f"Figure 1 live panel {letter}{suffix}"},
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
        "A": "Trayectorias de consumo por jaula desde el día 1 hasta el día 6; las líneas gruesas y las bandas resumen jaulas con el mismo peso.",
        "B": "Consumo por jaula en el día 6. Los puntos son jaulas; se muestra ANOVA de una vía y comparaciones exactas con ajuste de Holm.",
        "C": "Cambio de peso corporal desde basal. Los puntos pálidos son ratones; la inferencia y los resúmenes permanecen equilibrados por jaula.",
        "D": "Cambio medio de peso corporal por jaula en el día 6. Los puntos son jaulas; se muestra ANOVA de una vía y comparaciones exactas con ajuste de Holm.",
    }
    for letter, text in panels.items():
        path = output_dir / f"Figure_1_behavior_experiment_1_Panel_{letter}_spanish_LEGEND.txt"
        path.write_text(f"Figura 1, panel {letter}. {text}\n", encoding="utf-8")
    master = (
        "Figura 1. Experimento de una sola botella: análisis exploratorio de la cohorte de enero.\n\n"
        + "\n".join(f"({letter}) {panels[letter]}" for letter in "ABCD")
        + "\n\nLa unidad experimental y de inferencia es la jaula. El análisis principal "
          "usa un valor del día 6 por jaula; las observaciones individuales de ratones en "
          "los paneles de peso son descriptivas. Las bandas de A y C son intervalos "
          "descriptivos del 95 % obtenidos por remuestreo de trayectorias completas de "
          "jaulas, conservando la dependencia longitudinal y ponderando cada jaula por "
          "igual. Los valores globales corresponden a ANOVA de una vía y las tres "
          "comparaciones exactas por pares se ajustan por Holm dentro de cada variable. "
          "La aleatorización del tratamiento no fue documentada, los tamaños muestrales "
          "son pequeños y los resultados deben interpretarse como exploratorios.\n"
    )
    (output_dir / "Figure_1_behavior_experiment_1_LEGEND_spanish.txt").write_text(
        master, encoding="utf-8"
    )


def make_figure_values(
    consumption: pd.DataFrame,
    weights: pd.DataFrame,
    weight_cages: pd.DataFrame,
) -> pd.DataFrame:
    a = consumption.assign(
        panel="A", value_role="cage_trajectory", experimental_unit="cage",
        unit_id=consumption["E"], value_name="cumulative_cage_volume_removed_ml",
        plotted_value=consumption["cumulative_cage_volume_removed_ml"], subject_id=pd.NA,
    )
    b = consumption.loc[consumption["study_day"].eq(6)].copy().assign(
        panel="B", value_role="cage_endpoint", experimental_unit="cage",
        unit_id=lambda frame: frame["E"], value_name="cumulative_cage_volume_removed_ml",
        plotted_value=lambda frame: frame["cumulative_cage_volume_removed_ml"], subject_id=pd.NA,
    )
    c_mouse = weights.assign(
        panel="C", value_role="mouse_trajectory_descriptive", experimental_unit="mouse_descriptive_only",
        unit_id=weights["subject_id"], value_name="body_weight_change_pct",
        plotted_value=weights["body_weight_change_pct"],
    )
    c_cage = (
        weights.groupby(["E", "Treatment", "study_day"], as_index=False)["body_weight_change_pct"]
        .mean()
        .assign(
            panel="C", value_role="cage_day_trajectory", experimental_unit="cage",
            unit_id=lambda frame: frame["E"], value_name="body_weight_change_pct",
            plotted_value=lambda frame: frame["body_weight_change_pct"],
            subject_id=pd.NA, Sex=pd.NA, Genotype=pd.NA,
        )
    )
    d = weight_cages.assign(
        panel="D", value_role="cage_endpoint", experimental_unit="cage",
        unit_id=weight_cages["E"], value_name="cage_mean_body_weight_change_pct",
        plotted_value=weight_cages["cage_mean_body_weight_change_pct"],
        subject_id=pd.NA, Sex=pd.NA, Genotype=pd.NA,
    )
    columns = [
        "panel", "value_role", "experimental_unit", "unit_id", "E", "subject_id",
        "Treatment", "Sex", "Genotype", "study_day", "value_name", "plotted_value",
    ]
    return pd.concat(
        [a[columns], b[columns], c_mouse[columns], c_cage[columns], d[columns]],
        ignore_index=True,
    )


def write_legends(output_dir: Path, consumption_stats: pd.DataFrame, weight_stats: pd.DataFrame) -> None:
    intake_global = float(consumption_stats.loc[consumption_stats["test"].eq(OMNIBUS_TEST), "p_raw"].iloc[0])
    weight_global = float(weight_stats.loc[weight_stats["test"].eq(OMNIBUS_TEST), "p_raw"].iloc[0])
    panels = {
        "A": (
            "January E1-E12 cage consumption from Day 1 through Day 6. "
            "Treatments are split into facets. Visual key: thin: cages; thick/band: equal-cage mean/95% CI. "
            "Thin hollow-marker lines are cages; thick lines are equal-cage means and shaded bands are "
            "pointwise descriptive 95% percentile intervals from "
            f"{BOOTSTRAP_REPLICATES:,} whole-cage trajectory bootstrap resamples. "
            "Consumption was calculated from cage-bottle volume change and may include "
            "unmeasured leakage or evaporation."
        ),
        "B": (
            f"Day-6 cage-consumption endpoint. Dots are cages; diamonds and intervals are bootstrap cage means "
            f"and descriptive 95% intervals. One-way ANOVA {p_text(intake_global)}; the three "
            "pairwise exact p values are Holm-adjusted. The exact cage-label omnibus is retained in source data as a sensitivity analysis."
        ),
        "C": (
            "January mouse body-weight changes from each mouse's Day-1 baseline, split into treatment "
            "facets. Visual key: pale: mice; thin: cage-day; thick/band: equal-cage mean/95% CI. Pale points "
            "are individual mice, thin hollow-marker lines are cage-day means, and thick lines are equal-cage "
            "treatment means after first averaging mice within cage. Shaded bands are "
            f"pointwise descriptive 95% percentile intervals from {BOOTSTRAP_REPLICATES:,} whole-cage "
            "trajectory bootstrap resamples. Individual mouse observations are descriptive."
        ),
        "D": (
            f"Day-6 cage-mean body-weight percent change. Each dot is one cage, regardless of its mouse "
            f"count. One-way ANOVA {p_text(weight_global)}; the three pairwise exact p values "
            "are Holm-adjusted. The exact cage-label omnibus is retained in source data as a sensitivity analysis."
        ),
    }
    master = (
        "Figure 1. Single bottle experiment: exploratory analysis of the January cohort.\n\n"
        + "\n".join(f"({letter}) {text}" for letter, text in panels.items())
        + "\n\nPrimary inference uses one Day-6 value per cage. August E13/E14 are a separate "
          "Water-only cohort with bottle-weight but no bottle-volume records and are not pooled. Nonpositive "
          "weight rows are excluded; E10/FR6-2 (source label Control) is included as Water. FR5-4 was transferred "
          "from E2 to E10 for dehydration after the day-6 measurement; its primary measurements retain the "
          "original Allulose assignment and subsequent observations are flagged separately. E9 is retained with genotype Unknown. Pairwise families are "
          "Holm-adjusted within outcome. Treatment randomization was not documented, sample sizes are small, "
          "and results should be presented as exploratory. Panel A/C confidence bands are descriptive: each "
          "replicate resamples entire cages with replacement and applies the same sampled cage identities at "
          "all three days, preserving within-cage longitudinal dependence; cages are equally weighted.\n"
    )
    (output_dir / "Figure_1_behavior_experiment_1_LEGEND.txt").write_text(master, encoding="utf-8")
    (output_dir / "Figure_1_behavior_experiment_1_caption.txt").write_text(master, encoding="utf-8")
    for letter, text in panels.items():
        (output_dir / f"Figure_1_behavior_experiment_1_Panel_{letter}_LEGEND.txt").write_text(
            f"Figure 1 panel {letter}. {text}\n", encoding="utf-8"
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
        "consumption_values": args.analysis_dir / "consumption_primary_values.csv",
        "consumption_statistics": args.analysis_dir / "consumption_primary_statistics.csv",
        "body_weight_mouse_values": args.analysis_dir / "body_weight_primary_values.csv",
        "body_weight_cage_values": args.analysis_dir / "body_weight_cage_endpoint_values.csv",
        "body_weight_statistics": args.analysis_dir / "body_weight_primary_statistics.csv",
        "cage_design_qc": args.analysis_dir / "cage_design_and_size_qc.csv",
        "cohort_inventory": args.analysis_dir / "cohort_inventory.csv",
        "analysis_sources": args.analysis_dir / "analysis_sources.csv",
        "analysis_manifest": args.analysis_dir / "run_manifest.json",
        "statistical_audit": args.analysis_dir / "FIGURE1_STATISTICAL_AUDIT.md",
        "analysis_script": PAPER_ROOT / "scripts" / "Fig1" / "01_analyze_behavior_experiment1.py",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    consumption = pd.read_csv(paths["consumption_values"])
    consumption_stats = pd.read_csv(paths["consumption_statistics"])
    weights = pd.read_csv(paths["body_weight_mouse_values"])
    weight_cages = pd.read_csv(paths["body_weight_cage_values"])
    weight_stats = pd.read_csv(paths["body_weight_statistics"])
    if set(consumption["cohort"]) != {"January_E1-E12"}:
        raise ValueError("Consumption figure values must contain only the January cohort")
    if weights["E"].isin(["E13", "E14"]).any() or weight_cages["E"].isin(["E13", "E14"]).any():
        raise ValueError("August E13/E14 must not enter Figure 1")
    if not consumption_stats["experimental_unit"].eq("cage").all():
        raise ValueError("Consumption statistics are not cage-level")
    if not weight_stats["experimental_unit"].eq("cage").all():
        raise ValueError("Body-weight statistics are not cage-level")

    consumption_endpoint = consumption.loc[consumption["study_day"].eq(6)].copy()
    fig = plt.figure(figsize=FIGURE_SIZE_INCHES, facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.10, right=0.985, bottom=0.075, top=0.86,
        wspace=0.28, hspace=0.55,
    )
    panel_a_anchor = fig.add_subplot(outer[0, 0])
    panel_a_grid = outer[0, 0].subgridspec(1, 3, wspace=0.13)
    panel_a_axes: list[plt.Axes] = []
    for index in range(3):
        panel_a_axes.append(
            fig.add_subplot(
                panel_a_grid[0, index],
                sharey=panel_a_axes[0] if panel_a_axes else None,
            )
        )
    panel_b_axis = fig.add_subplot(outer[0, 1])
    panel_c_anchor = fig.add_subplot(outer[1, 0])
    panel_c_grid = outer[1, 0].subgridspec(1, 3, wspace=0.13)
    panel_c_axes: list[plt.Axes] = []
    for index in range(3):
        panel_c_axes.append(
            fig.add_subplot(
                panel_c_grid[0, index],
                sharey=panel_c_axes[0] if panel_c_axes else None,
            )
        )
    panel_d_axis = fig.add_subplot(outer[1, 1])

    consumption_intervals = draw_consumption_trajectory(
        panel_a_axes, panel_a_anchor, consumption, args.seed + 10_000,
    )
    draw_cage_endpoint(
        panel_b_axis, consumption_endpoint, "cumulative_cage_volume_removed_ml", consumption_stats,
        "Primary Day-6 consumption endpoint", "Consumption (mL/cage)", "B", args.seed,
    )
    weight_intervals = draw_weight_trajectory(
        panel_c_axes, panel_c_anchor, weights, args.seed + 20_000,
    )
    draw_cage_endpoint(
        panel_d_axis, weight_cages, "cage_mean_body_weight_change_pct", weight_stats,
        "Primary Day-6 body-weight endpoint", "Cage-mean change from baseline (%)", "D", args.seed + 1,
    )
    fig.suptitle(
        "Single bottle experiment",
        fontsize=15.5, fontweight="bold", y=0.975,
    )

    stem = staging_output / "Figure_1_behavior_experiment_1"
    fig.savefig(
        stem.with_suffix(".pdf"), facecolor="white",
        metadata={"Title": "Single bottle experiment"},
    )
    fig.savefig(stem.with_suffix(".png"), dpi=args.dpi, facecolor="white")
    panel_axes = {
        "A": [panel_a_anchor, *panel_a_axes],
        "B": panel_b_axis,
        "C": [panel_c_anchor, *panel_c_axes],
        "D": panel_d_axis,
    }
    english_crop_boxes = export_panel_crops(
        fig,
        panel_axes,
        panel_output,
        args.dpi,
    )
    translation_receipt = translate_figure_texts_to_spanish(fig)
    export_panel_crops(
        fig, panel_axes, panel_output, args.dpi,
        suffix="_spanish", crop_boxes=english_crop_boxes,
    )
    plt.close(fig)

    write_legends(legend_output, consumption_stats, weight_stats)
    write_spanish_panel_legends(legend_output)
    figure_values = make_figure_values(consumption, weights, weight_cages)
    figure_values.to_csv(source_output / "Figure_1_behavior_experiment_1_values.csv", index=False)
    trajectory_intervals = pd.concat(
        [consumption_intervals, weight_intervals], ignore_index=True,
    )
    intervals_path = source_output / "Figure_1_behavior_experiment_1_trajectory_confidence_intervals.csv"
    trajectory_intervals.to_csv(intervals_path, index=False)
    figure_stats = pd.concat([
        consumption_stats.assign(panel="B"),
        weight_stats.assign(panel="D"),
    ], ignore_index=True, sort=False)
    figure_stats.to_csv(source_output / "Figure_1_behavior_experiment_1_statistics.csv", index=False)

    with (source_output / "Figure_1_behavior_experiment_1_sources.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(("role", "path", "sha256"))
        for role, path in paths.items():
            writer.writerow((role, str(path), sha256(path)))
        writer.writerow(("renderer", str(Path(__file__).resolve()), sha256(Path(__file__).resolve())))

    shutil.copyfile(paths["statistical_audit"], provenance_output / "FIGURE1_STATISTICAL_AUDIT.md")
    translation_receipt_path = provenance_output / "Figure_1_behavior_experiment_1_spanish_translation_receipt.json"
    translation_receipt_path.write_text(
        json.dumps(translation_receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_local_readme(provenance_output)

    expected_graphics = [
        "Figure_1_behavior_experiment_1.pdf",
        "Figure_1_behavior_experiment_1.png",
        *[
            f"panels/Figure_1_behavior_experiment_1_Panel_{letter}.{extension}"
            for letter in "ABCD" for extension in ("pdf", "png")
        ],
        *[
            f"panels/Figure_1_behavior_experiment_1_Panel_{letter}_spanish.{extension}"
            for letter in "ABCD" for extension in ("pdf", "png")
        ],
    ]
    missing_graphics = [name for name in expected_graphics if not (staging_output / name).is_file()]
    if missing_graphics:
        raise RuntimeError(f"Missing rendered graphics: {missing_graphics}")
    manifest = {
        "figure": "Figure_1_behavior_experiment_1",
        "analysis_status": "exploratory",
        "primary_cohort": "January_E1-E12",
        "excluded_separate_cohort": "August_E13-E14",
        "seed": args.seed,
        "dpi": args.dpi,
        "analysis_directory": str(args.analysis_dir),
        "output_directory": str(final_output),
        "subpanels": ["A", "B", "C", "D"],
        "formats": ["pdf", "png"],
        "language_contract": {
            "combined_figure": "English only",
            "standalone_panels": ["English", "Spanish"],
            "spanish_translation_receipt": str(translation_receipt_path.relative_to(staging_output)),
            "spanish_translation_receipt_sha256": sha256(translation_receipt_path),
            "crop_geometry": "Spanish and English pairs use identical fixed crop boxes",
        },
        "panel_export_method": "exact crops of live axes in the combined rendered figure",
        "figure_size_inches": list(FIGURE_SIZE_INCHES),
        "figure_orientation": "portrait",
        "primary_experimental_unit": "cage for all treatment summaries, intervals, and inferential tests",
        "panel_C_mouse_rows": "descriptive; treatment summaries are cage-balanced",
        "trajectory_confidence_intervals": {
            "panels": ["A", "C"],
            "level": 0.95,
            "type": "pointwise descriptive percentile interval",
            "resampling_unit": "whole cage trajectory",
            "longitudinal_rule": "same resampled cage identities at Days 1, 3, and 6",
            "cage_weighting": "equal",
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "small_sample_limitation": "only 3-5 cages per treatment; intervals are descriptive",
            "values_file": str(intervals_path.relative_to(staging_output)),
            "values_sha256": sha256(intervals_path),
        },
        "p_value_annotation_style": {
            "layout": "vertical multiline list; text itself is not rotated",
            "omnibus_font_size_pt": P_VALUE_GLOBAL_FONTSIZE_PT,
            "pairwise_font_size_pt": P_VALUE_PAIRWISE_FONTSIZE_PT,
            "fontweight": P_VALUE_FONTWEIGHT,
        },
        "statistics": {
            "displayed_omnibus": "ordinary equal-variance one-way ANOVA on independent cage endpoints",
            "omnibus_sensitivity": "exact exhaustive cage-label allocation of the same F statistic",
            "pairwise": "exact cage-label mean-difference tests; Holm within outcome",
        },
        "legend_only_visual_keys": {
            "A": "thin: cages; thick/band: equal-cage mean/95% CI",
            "C": "pale: mice; thin: cage-day; thick/band: equal-cage mean/95% CI",
        },
        "in_panel_explanatory_annotations": [],
        "renderer_sha256": sha256(Path(__file__).resolve()),
        "graphics": expected_graphics,
        "graphics_sha256": {
            name: sha256(staging_output / name) for name in expected_graphics
        },
    }
    (provenance_output / "Figure_1_behavior_experiment_1_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    publish_staging_output(staging_output, final_output)
    print(f"[FIGURE] {final_output / stem.with_suffix('.pdf').name}")
    print(f"[FIGURE] {final_output / stem.with_suffix('.png').name}")
    print("[PANELS] live A-D crops exported as PDF and PNG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
