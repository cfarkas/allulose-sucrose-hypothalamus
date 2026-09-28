#!/usr/bin/env python3
"""Build Figure S4: cage-aware Experiment 2 behavior preference."""

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
from typing import Any, Sequence

# Matplotlib honors SOURCE_DATE_EPOCH for PDF metadata.  A fixed value makes
# repeated publication renders byte-reproducible instead of embedding run time.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.4,
    "axes.titlesize": 12.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "legend.fontsize": 8.5,
})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.transforms import Bbox


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
sys.path.insert(0, str(PAPER_ROOT / "scripts" / "shared"))
from spanish_panel_text import translate_figure_texts_to_spanish
PROJECT_ROOT = PAPER_ROOT.parent
FIGURE_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "FigS4"
DEFAULT_ANALYSIS = FIGURE_ANALYSIS_ROOT / "results"
DEFAULT_OUTPUT = PAPER_ROOT / "FigS4"
def _first_existing(*candidates: Path) -> Path:
    """First candidate that exists, else the first, so defaults name real files.

    The analyses/FigS3 location is a private working copy that a downloaded tree
    does not contain; the workbook itself ships as FigS4/raw_data.
    """
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


SHIPPED_WORKBOOK = PAPER_ROOT / "FigS4" / "raw_data" / "UCSC_MICE.xlsx"
DEFAULT_WORKBOOK = _first_existing(
    FIGURE_ANALYSIS_ROOT / "raw" / "UCSC_MICE.xlsx", SHIPPED_WORKBOOK,
)
FIGURE_NUMBER = "S3"
TREATMENTS = ("Water", "Allulose", "Sucrose")
COLORS = {"Water": "#56B4E9", "Allulose": "#009E73", "Sucrose": "#D55E00"}
ABBR = {"Water": "W", "Allulose": "A", "Sucrose": "S"}
MARKERS = {
    "J1": "o", "J2": "s", "J3": "^", "J4": "D",
    "J5": "P", "J6": "<", "J7": "v", "J8": "h",
}
PANELS = {
    "A": "solution",
    "B": "water",
    "C": "solution_fraction",
    "D": "food",
    "E": "body_weight",
}
ENDPOINT_STYLE = {
    "solution": ("Solution-side bottle", "Consumption (g/cage/day)"),
    "water": ("Water-side bottle", "Consumption (g/cage/day)"),
    "solution_fraction": ("Two-bottle solution-side share", "Solution / total fluid"),
    "food": ("Food consumption", "Food (g/cage/day)"),
}
STEM_NAME = "Figure_S4_behavior_preference"
OUTPUT_SUBDIRS = ("panels", "legends", "source_data", "provenance")
GENERATED_TREE_MARKER = ".paper_generated_tree.json"
OUTPUT_PRODUCER = "Paper/scripts/FigS4/02_make_figure_s4_behavior_preference.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--formats", default="pdf,png")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
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
    """Fail closed unless an existing nonempty tree is wholly script-owned.

    Publication is an atomic directory replacement, so without this the default
    --output-dir, which is the live FigS4 directory, would take the numbered
    scripts, raw_data and legacy down with it. Its three sibling behaviour
    renderers have always had this check; this one did not.
    """
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
    """Create a complete sibling staging tree without touching the live figure."""
    final_output.parent.mkdir(parents=True, exist_ok=True)
    validate_replacement_target(final_output)
    staging = Path(tempfile.mkdtemp(prefix=f".{final_output.name}.stage-", dir=final_output.parent))
    for name in OUTPUT_SUBDIRS:
        (staging / name).mkdir()
    atexit.register(shutil.rmtree, staging, ignore_errors=True)
    return staging


def publish_staging_output(staging: Path, final_output: Path) -> None:
    """Atomically expose a validated tree, restoring the previous tree on error."""
    links = [path for path in staging.rglob("*") if path.is_symlink()]
    if links:
        raise RuntimeError(f"Publication output contains forbidden symlinks: {links}")
    write_generated_tree_marker(staging)
    validate_replacement_target(final_output)
    backup = staging.with_name(staging.name + ".previous")
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
    text = 'FIGURE S4 — TWO-BOTTLE CONSUMPTION, FOOD REMOVAL AND BODY WEIGHT\n\nFigure S4 presents a two-bottle assay with eight cages and two mice per cage:\n2 Water, 3 Sucrose and 3 Allulose cages. Control cages receive two water bottles;\nthe other cages receive water and the assigned sugar solution.\n\nPanels\n  A  Consumption from the bottle on the assigned-solution side.\n  B  Water-bottle consumption.\n  C  Assigned-solution fraction of total fluid consumption.\n  D  Food removal.\n  E  Percentage body-weight change.\n\nBottle mass changes and food removal are reported in g/cage/day. Summaries\nexclude initialization records and use seven daily observations per cage.\nThe solution fraction in controls describes the two water-bottle positions.\nWeights are expressed relative to the pre-exposure baseline. Cage is the\nindependent unit. Displayed endpoints use ordinary one-way ANOVA and raw exact\ntwo-sided Mann–Whitney comparisons; those pairwise values are unadjusted.\nSmall cage counts and unmeasured fluid losses constrain interpretation.\n\nFiles\n  Figure_S4_behavior_preference.pdf/.png  Complete English figure.\n  panels/       Individual A–E panels in English and Spanish, PNG and PDF.\n  legends/      Figure and panel captions in both languages.\n  source_data/  Plotted values, endpoint tests and source-audit tables.\n  raw_data/UCSC_MICE.xlsx  Input workbook, supplied through Zenodo.\n  provenance/   Input and output hashes and validation records.\nAll publication PNGs are exported at 600 dpi.\n\nReproduce from Paper/ with the paper environment and fresh output paths:\n  python FigS4/01_analyze_behavior_experiment2.py \\\n    --workbook FigS4/raw_data/UCSC_MICE.xlsx --output-dir /tmp/FigS4_analysis\n  python FigS4/02_make_figure_s4_behavior_preference.py \\\n    --workbook FigS4/raw_data/UCSC_MICE.xlsx \\\n    --analysis-dir /tmp/FigS4_analysis --output-dir /tmp/FigS4_render --dpi 600\n\nRAW_DATA_MANIFEST.csv identifies the packaged workbook and its scientific\nsource. The analysis reads its cells, formulas and recorded measurements.\n'
    (output_dir / "README.txt").write_text(text, encoding="utf-8")


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.13, 1.08, letter, transform=ax.transAxes, fontsize=19,
        fontweight="bold", ha="right", va="top", clip_on=False,
    )


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", width=1.1, length=4)
    ax.grid(axis="y", color="#d8dde0", linewidth=0.65, alpha=0.7, zorder=0)


def condition_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color=COLORS[treatment], marker="o", markersize=5.5,
               linewidth=2.0, label=treatment)
        for treatment in TREATMENTS
    ]


def compact_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.0001:
        return f"{value:.1e}"
    return f"{value:.4f}".rstrip("0").rstrip(".")


def pvalue_text(stats: pd.DataFrame, endpoint: str) -> str:
    subset = stats.loc[stats["endpoint"].eq(endpoint)]
    global_row = subset.loc[subset["scope"].eq("global")].iloc[0]
    pairs = subset.loc[subset["scope"].eq("pairwise")].set_index("contrast")
    return "\n".join((
        f"One-way ANOVA p = {compact_p(float(global_row['p_value_raw']))}",
        f"MWU W-S p = {compact_p(float(pairs.loc['Water - Sucrose', 'p_value_raw']))}",
        f"MWU W-A p = {compact_p(float(pairs.loc['Water - Allulose', 'p_value_raw']))}",
        f"MWU S-A p = {compact_p(float(pairs.loc['Sucrose - Allulose', 'p_value_raw']))}",
    ))


def tukey_outlier_mask(values: pd.Series) -> np.ndarray:
    array = values.astype(float).to_numpy()
    q1, q3 = np.quantile(array, (0.25, 0.75))
    iqr = q3 - q1
    return (array < q1 - 1.5 * iqr) | (array > q3 + 1.5 * iqr)


def draw_summary_points(ax: plt.Axes, values: pd.DataFrame, stats: pd.DataFrame, endpoint: str) -> None:
    clean_axes(ax)
    all_values = values["value"].astype(float).to_numpy()
    data_min, data_max = float(np.min(all_values)), float(np.max(all_values))
    span = max(data_max - data_min, 0.1 if endpoint == "solution_fraction" else 1.0)
    lower = data_min - 0.12 * span
    upper = data_max + 1.15 * span
    if endpoint == "solution_fraction":
        lower = min(0.38, lower)
        upper = max(1.12, upper)
        ax.axhline(0.5, color="#5e676a", linestyle="--", linewidth=1.0, zorder=1)
    ax.set_ylim(lower, upper)
    for x, treatment in enumerate(TREATMENTS):
        subset = values.loc[values["treatment"].eq(treatment)].sort_values("cage")
        offsets = np.linspace(-0.11, 0.11, len(subset)) if len(subset) > 1 else np.array([0.0])
        flagged = tukey_outlier_mask(subset["value"])
        for offset, is_flagged, row in zip(offsets, flagged, subset.itertuples(index=False)):
            if is_flagged:
                ax.scatter(
                    x + offset, row.value, s=54, marker="x", color="black",
                    linewidths=1.55, zorder=6,
                )
            else:
                ax.scatter(
                    x + offset, row.value, s=42, marker=MARKERS[row.cage],
                    facecolor=COLORS[treatment], edgecolor="black", linewidth=0.55, zorder=4,
                )
        mean = float(subset["value"].mean())
        sd = float(subset["value"].std(ddof=1)) if len(subset) > 1 else 0.0
        ax.errorbar(
            x, mean, yerr=sd, fmt="none", ecolor="black", elinewidth=1.3,
            capsize=4, capthick=1.3, zorder=3,
        )
        ax.plot([x - 0.22, x + 0.22], [mean, mean], color="black", linewidth=2.2, zorder=5)
    ax.set_xlim(-0.50, 2.50)
    ax.set_xticks(range(3))
    ax.set_xticklabels([ABBR[item] for item in TREATMENTS], fontweight="bold")
    summary_title = "Cage summaries, mean +/- SD"
    if endpoint == "solution_fraction":
        # Invisible discriminator for the one Spanish title that needs a shorter fit.
        summary_title += " "
    ax.set_title(summary_title, fontsize=9.0, fontweight="bold", pad=3)
    ax.text(
        0.03, 0.97, pvalue_text(stats, endpoint), transform=ax.transAxes,
        ha="left", va="top", fontsize=8.0, fontweight="bold", linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 1.0},
        zorder=8,
    )


def draw_endpoint(
    fig: plt.Figure,
    spec: Any,
    letter: str,
    endpoint: str,
    daily: pd.DataFrame,
    summary: pd.DataFrame,
    stats: pd.DataFrame,
) -> list[plt.Axes]:
    inner = spec.subgridspec(1, 2, width_ratios=(3.15, 1.55), wspace=0.34)
    ax = fig.add_subplot(inner[0, 0])
    sax = fig.add_subplot(inner[0, 1])
    clean_axes(ax)
    panel_label(ax, letter)
    title, ylabel = ENDPOINT_STYLE[endpoint]
    ax.set_title(title, loc="left", pad=7)
    subset = daily.loc[
        daily["endpoint"].eq(endpoint)
        & daily["analysis_period"].eq("primary")
        & np.isfinite(pd.to_numeric(daily["value"], errors="coerce"))
    ].copy()
    for cage, cage_data in subset.groupby("cage", observed=True):
        cage_data = cage_data.sort_values("study_day")
        treatment = str(cage_data["treatment"].iloc[0])
        ax.plot(
            cage_data["study_day"], cage_data["value"],
            color=COLORS[treatment], marker=MARKERS[cage], markersize=4.5,
            linewidth=1.25, alpha=0.78, label=cage,
        )
    if endpoint == "solution_fraction":
        ax.axhline(0.5, color="#5e676a", linestyle="--", linewidth=1.0, zorder=1)
        ax.set_ylim(0.35, 1.0)
    else:
        ax.set_ylim(bottom=0)
    ax.set_xlim(0.75, 7.25)
    ax.set_xticks(range(1, 8))
    ax.set_xlabel("Study day")
    ax.set_ylabel(ylabel)
    ax.legend(handles=condition_handles(), loc="best", frameon=False, ncol=1, handlelength=1.5)
    endpoint_summary = summary.loc[summary["endpoint"].eq(endpoint)].copy()
    draw_summary_points(sax, endpoint_summary, stats, endpoint)
    return [ax, sax]


def draw_weights(
    fig: plt.Figure,
    spec: Any,
    weights: pd.DataFrame,
    cage_summary: pd.DataFrame,
    stats: pd.DataFrame,
) -> list[plt.Axes]:
    inner = spec.subgridspec(1, 2, width_ratios=(3.15, 1.55), wspace=0.34)
    ax = fig.add_subplot(inner[0, 0])
    sax = fig.add_subplot(inner[0, 1])
    clean_axes(ax)
    clean_axes(sax)
    panel_label(ax, "E")
    ax.set_title("Body-weight change", loc="left", pad=7)
    primary = weights.loc[
        weights["analysis_period"].eq("primary")
        & np.isfinite(pd.to_numeric(weights["weight_change_g"], errors="coerce"))
    ].copy()
    cage_daily = primary.groupby(
        ["cage", "treatment", "study_day"], as_index=False, observed=True
    )["weight_change_g"].mean()
    for cage, subset in cage_daily.groupby("cage", observed=True):
        subset = subset.sort_values("study_day")
        treatment = str(subset["treatment"].iloc[0])
        ax.plot(
            subset["study_day"], subset["weight_change_g"], color=COLORS[treatment],
            marker=MARKERS[cage], markersize=4.7, markeredgewidth=0.8,
            linewidth=1.65, alpha=0.78,
        )
    ax.axhline(0, color="#444b4e", linewidth=1.25, linestyle="--")
    ax.set_xlim(-0.2, 7.2)
    ax.set_xticks(range(0, 8))
    ax.set_xlabel("Study day")
    ax.set_ylabel("Change from baseline (g)")
    ax.legend(handles=condition_handles(), loc="best", frameon=False)

    draw_summary_points(sax, cage_summary, stats, "weight_change")
    sax.axhline(0, color="#444b4e", linewidth=0.9, linestyle="--")
    return [ax, sax]


def export_panel_crops(
    fig: plt.Figure,
    panel_axes: dict[str, list[plt.Axes]],
    output_dir: Path,
    dpi: int,
    suffix: str = "",
    crop_boxes: dict[str, Bbox] | None = None,
) -> dict[str, Bbox]:
    """Save exact live crops of the combined figure; no panel is redrawn."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    result: dict[str, Bbox] = {}
    for letter, slug in PANELS.items():
        if crop_boxes is None:
            boxes = []
            for axis in panel_axes[letter]:
                box = axis.get_tightbbox(renderer)
                if box is not None and np.isfinite(box.extents).all():
                    boxes.append(box)
            if not boxes:
                raise RuntimeError(f"No drawable axes for panel {letter}")
            box = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted())
            pad = 0.07
            box = Bbox.from_extents(box.x0 - pad, box.y0 - pad, box.x1 + pad, box.y1 + pad)
        else:
            box = crop_boxes[letter]
        result[letter] = box.frozen()
        stem = output_dir / f"Figure_S4_behavior_preference_Panel_{letter}_{slug}{suffix}"
        fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches=box, facecolor="white")
        fig.savefig(stem.with_suffix(".png"), format="png", dpi=dpi, bbox_inches=box, facecolor="white")
    return result


def write_spanish_panel_legends(output_dir: Path) -> None:
    texts = {
        "A_solution": "Consumo de la botella de solución por jaula y promedio de siete días por jaula; las barras muestran media ± DE.",
        "B_water": "Consumo de la botella de agua por jaula y promedio de siete días por jaula; las barras muestran media ± DE.",
        "C_solution_fraction": "Fracción diaria asignada a la botella de solución. En controles Agua corresponde a la asignación al lado de la botella 2; la línea discontinua marca 0,5.",
        "D_food": "Consumo compartido de alimento por jaula y promedio de siete días por jaula; la jaula es la unidad experimental.",
        "E_body_weight": "Cambio de peso corporal desde basal, promediado entre los dos ratones por jaula; la inferencia usa una media final-menos-basal por jaula.",
    }
    for key, text in texts.items():
        letter = key.split("_", 1)[0]
        (output_dir / f"{STEM_NAME}_Panel_{letter}_spanish_LEGEND.txt").write_text(
            f"Figura {FIGURE_NUMBER}, panel {letter}. {text} Se muestra ANOVA de una vía; las comparaciones por pares son pruebas exactas de Mann–Whitney sin ajuste de multiplicidad.\n",
            encoding="utf-8",
        )
    master = (
        "Figura S4. Experimento 2 de conducta con dos botellas y peso corporal a nivel de jaula.\n\n"
        + "\n".join(f"({key[0]}) {texts[key]}" for key in texts)
        + "\n\nLa jaula es la unidad experimental en todos los análisis (Agua n=2, "
          "Alulosa n=3 y Sacarosa n=3). Los valores globales son valores p de ANOVA "
          "ordinario de una vía. Las comparaciones Agua–Sacarosa, Agua–Alulosa y "
          "Sacarosa–Alulosa usan valores p exactos bilaterales de Mann–Whitney sin "
          "ajuste por multiplicidad. Los puntos son jaulas, las barras muestran media "
          "± DE y las cruces señalan valores fuera de las cercas de Tukey de 1,5 RIC, "
          "que se conservaron en todos los análisis. Las entradas del día 0 son masas "
          "de inicialización y no consumo; los días 1–7 forman la ventana principal. "
          "En controles Agua, la fracción de la botella 2 mide asignación lateral, no "
          "preferencia por edulcorante. Los tamaños de muestra son pequeños y no se "
          "formula una conclusión confirmatoria.\n"
    )
    (output_dir / f"{STEM_NAME}_LEGEND_spanish.txt").write_text(
        master, encoding="utf-8"
    )


def write_legends(
    output_dir: Path,
    stats: pd.DataFrame,
    analysis_dir: Path,
    workbook: Path,
) -> None:
    panel_text = {
        "A_solution": (
            "(A) Solution-side bottle consumption. For Allulose and Sucrose this is the sweetener bottle; for Water "
            "controls it is the second water bottle. Lines are cages. The inset shows each cage's seven-day mean and "
            "group mean +/- SD."
        ),
        "B_water": (
            "(B) Water-side bottle consumption. Lines are cages; inset points are independent cage summaries with "
            "group mean +/- SD."
        ),
        "C_solution_fraction": (
            "(C) Daily solution-side share, solution/(solution+water). In Water controls this is bottle-2 share and "
            "therefore measures side allocation rather than sweetener preference. The dashed line marks 0.5."
        ),
        "D_food": (
            "(D) Shared cage food consumption. Lines are cages; inset points are seven-day cage means with group mean +/- SD."
        ),
        "E_body_weight": (
            "(E) Body-weight change from baseline, averaged across the two mice in each cage at each time point. Lines "
            "and inset points are cages; inference uses one final-minus-baseline mean per cage."
        ),
    }
    shared = (
        "Cage is the experimental unit throughout (Water n=2, Allulose n=3, Sucrose n=3). Global values are ordinary "
        "one-way ANOVA p values. Pairwise W-S, W-A, and S-A "
        "values are raw exact two-sided Mann-Whitney U p values; no pairwise multiplicity adjustment is applied. "
        "With these cage counts, ANOVA assumptions are weakly assessable and no confirmatory claim is made. Dots are cages; "
        "bars are mean +/- SD. Crosses mark within-condition Tukey 1.5-IQR outliers, which remain in all analyses. "
        "Day 0 intake entries are initialization masses, not consumption, and were excluded; days 1-7 form the "
        "primary window. "
        "Hoja 3 versus Hoja 4 reconciliation found refill-adjacent solution discrepancies and one food discrepancy; "
        "Hoja 4 values were retained as the intended daily sheet, and a source-discrepancy-robust sensitivity excludes "
        "implicated dates. Sex/genotype adjustment was not possible from Hoja 4."
    )
    stats_lines = ["Primary exploratory statistics:"]
    for endpoint in ("solution", "water", "solution_fraction", "food", "weight_change"):
        subset = stats.loc[stats["endpoint"].eq(endpoint)]
        global_row = subset.loc[subset["scope"].eq("global")].iloc[0]
        stats_lines.append(
            f"  {endpoint}: one-way ANOVA raw p={float(global_row['p_value_raw']):.6g} "
            f"(F={float(global_row['statistic']):.6g}, df={int(global_row['df_between'])},"
            f"{int(global_row['df_within'])})."
        )
        for row in subset.loc[subset["scope"].eq("pairwise")].itertuples(index=False):
            stats_lines.append(
                f"    {row.contrast}: U={row.statistic:.6g}, exact raw p={row.p_value_raw:.6g}."
            )
    master = "\n\n".join([
        "Figure S4. Experiment 2 two-bottle behavior and cage-level body weight.",
        "\n\n".join(panel_text.values()),
        shared,
        "\n".join(stats_lines),
        f"Source workbook: {workbook.resolve()}",
        f"Analysis directory: {analysis_dir.resolve()}",
    ]) + "\n"
    stem = "Figure_S4_behavior_preference"
    (output_dir / f"{stem}_LEGEND.txt").write_text(master, encoding="utf-8")
    (output_dir / f"{stem}_caption.txt").write_text(master, encoding="utf-8")
    for slug, text in panel_text.items():
        (output_dir / f"{stem}_panel_{slug}_LEGEND.txt").write_text(
            "Figure S4. Experiment 2 two-bottle behavior and cage-level body weight.\n\n"
            + text + "\n\n" + shared + "\n",
            encoding="utf-8",
        )


def main() -> int:
    args = parse_args()
    args.analysis_dir = args.analysis_dir.resolve()
    args.workbook = args.workbook.resolve()
    final_output = Path(os.path.abspath(os.path.expanduser(args.output_dir)))
    staging_output = prepare_staging_output(final_output)
    panel_output = staging_output / "panels"
    legend_output = staging_output / "legends"
    source_output = staging_output / "source_data"
    provenance_output = staging_output / "provenance"
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    if set(formats) != {"pdf", "png"}:
        raise ValueError("Publication output requires exactly PDF and PNG formats")

    daily_path = args.analysis_dir / "daily_endpoints_long.csv"
    summary_path = args.analysis_dir / "cage_endpoint_summaries.csv"
    stats_path = args.analysis_dir / "one_way_anova_mann_whitney_statistics.csv"
    weights_path = args.analysis_dir / "mouse_weight_long.csv"
    mouse_summary_path = args.analysis_dir / "mouse_weight_summary_primary.csv"
    weight_cage_path = args.analysis_dir / "weight_cage_cluster_summary.csv"
    validation_path = args.analysis_dir / "hoja3_hoja4_intake_validation.csv"
    audit_path = args.analysis_dir / "legacy_analysis_audit.csv"
    required = [daily_path, summary_path, stats_path, weights_path, mouse_summary_path, weight_cage_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing analysis inputs: " + ", ".join(missing))

    daily = pd.read_csv(daily_path)
    summary = pd.read_csv(summary_path)
    stats_all = pd.read_csv(stats_path)
    weights = pd.read_csv(weights_path)
    mouse_summary = pd.read_csv(mouse_summary_path)
    weight_cage = pd.read_csv(weight_cage_path)
    primary_summary = summary.loc[summary["analysis_set"].eq("primary")].copy()
    primary_stats = stats_all.loc[
        stats_all["analysis_set"].isin(("primary", "primary_weight_cluster_sensitivity"))
    ].copy()
    if primary_summary["cage"].nunique() != 8 or mouse_summary["mouse_id"].nunique() != 16:
        raise ValueError("Unexpected design counts in analysis outputs")

    expected_endpoints = {"solution", "water", "solution_fraction", "food", "weight_change"}
    if set(primary_stats["endpoint"].unique()).intersection(expected_endpoints) != expected_endpoints:
        raise ValueError("Missing one-way ANOVA/exact MWU rows for one or more displayed endpoints")
    if any(
        str(column).lower().startswith("p_")
        and any(token in str(column).lower() for token in ("adjust", "correct"))
        for column in primary_stats.columns
    ):
        raise ValueError("Figure S4 forbids adjusted pairwise columns in the current statistics schema")

    fig = plt.figure(figsize=(14.2, 11.2), facecolor="white")
    outer = fig.add_gridspec(
        3, 2, left=0.065, right=0.985, bottom=0.065, top=0.965,
        wspace=0.20, hspace=0.37, height_ratios=(1.0, 1.0, 1.0),
    )
    panel_axes: dict[str, list[plt.Axes]] = {}
    panel_axes["A"] = draw_endpoint(fig, outer[0, 0], "A", "solution", daily, primary_summary, primary_stats)
    panel_axes["B"] = draw_endpoint(fig, outer[0, 1], "B", "water", daily, primary_summary, primary_stats)
    panel_axes["C"] = draw_endpoint(fig, outer[1, 0], "C", "solution_fraction", daily, primary_summary, primary_stats)
    panel_axes["D"] = draw_endpoint(fig, outer[1, 1], "D", "food", daily, primary_summary, primary_stats)
    panel_axes["E"] = draw_weights(fig, outer[2, :], weights, weight_cage, primary_stats)

    stem = staging_output / STEM_NAME
    for extension in formats:
        fig.savefig(
            stem.with_suffix(f".{extension}"), dpi=args.dpi,
            facecolor="white", bbox_inches="tight",
        )
    english_crop_boxes = export_panel_crops(fig, panel_axes, panel_output, args.dpi)
    translation_receipt = translate_figure_texts_to_spanish(
        fig,
        extra={
            "Solution-side bottle": "Botella del lado de la solución",
            "Water-side bottle": "Botella del lado del agua",
            "Consumption (g/cage/day)": "Consumo (g/jaula/día)",
            "Two-bottle solution-side share": "Fracción de solución (dos botellas)",
            "Solution / total fluid": "Solución / líquido total",
            "Food consumption": "Consumo de alimento",
            "Food (g/cage/day)": "Alimento (g/jaula/día)",
            "Group mean ± SD": "Media del grupo ± DE",
            "Cage summaries, mean +/- SD ": "Media por jaula ± DE",
            "Cage summaries, mean +/- SD": "Media por jaula ± DE",
        },
    )
    export_panel_crops(
        fig, panel_axes, panel_output, args.dpi,
        suffix="_spanish", crop_boxes=english_crop_boxes,
    )
    plt.close(fig)

    primary_stats = primary_stats.loc[primary_stats["endpoint"].isin(expected_endpoints)].copy()
    primary_stats["multiplicity_adjustment"] = "none; raw p values displayed"
    primary_stats.to_csv(source_output / f"{STEM_NAME}_statistics.csv", index=False)
    plotted_cages = pd.concat([
        primary_summary.loc[primary_summary["endpoint"].isin(
            ("solution", "water", "solution_fraction", "food")
        )],
        weight_cage,
    ], ignore_index=True, sort=False)
    plotted_cages["outlier_tukey_1_5_iqr"] = False
    for (_, treatment), index in plotted_cages.groupby(
        ["endpoint", "treatment"], observed=True
    ).groups.items():
        plotted_cages.loc[index, "outlier_tukey_1_5_iqr"] = tukey_outlier_mask(
            plotted_cages.loc[index, "value"]
        )
    plotted_cages.to_csv(source_output / f"{STEM_NAME}_plotted_cage_values.csv", index=False)
    cage_weight_daily = weights.loc[weights["analysis_period"].eq("primary")].groupby(
        ["cage", "treatment", "date", "study_day"], as_index=False, observed=True
    )["weight_change_g"].mean()
    cage_weight_daily.to_csv(
        source_output / f"{STEM_NAME}_plotted_cage_weight_trajectories.csv", index=False
    )
    write_legends(legend_output, primary_stats, args.analysis_dir, args.workbook)
    write_spanish_panel_legends(legend_output)
    translation_receipt_path = provenance_output / f"{STEM_NAME}_spanish_translation_receipt.json"
    translation_receipt_path.write_text(
        json.dumps(translation_receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    sources = [
        ("workbook", args.workbook),
        ("analysis_script", PAPER_ROOT / "scripts" / "FigS4" / "01_analyze_behavior_experiment2.py"),
        ("figure_script", Path(__file__).resolve()),
        ("daily_endpoints", daily_path),
        ("cage_summaries", summary_path),
        ("anova_mwu_statistics", stats_path),
        ("mouse_weights", weights_path),
        ("mouse_weight_summary", mouse_summary_path),
        ("weight_cage_summary", weight_cage_path),
        ("Hoja3_Hoja4_validation", validation_path),
        ("legacy_audit", audit_path),
    ]
    with (source_output / f"{STEM_NAME}_sources.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(("role", "path", "exists", "bytes", "sha256"))
        for role, path in sources:
            writer.writerow((
                role, str(path.resolve()), path.is_file(),
                path.stat().st_size if path.is_file() else "",
                sha256(path) if path.is_file() else "",
            ))

    if validation_path.is_file():
        shutil.copyfile(
            validation_path,
            source_output / f"{STEM_NAME}_hoja3_hoja4_intake_validation.csv",
        )
    if audit_path.is_file():
        shutil.copyfile(
            audit_path,
            source_output / f"{STEM_NAME}_legacy_analysis_audit.csv",
        )
    write_local_readme(provenance_output)

    expected_graphics = [
        f"{STEM_NAME}.pdf",
        f"{STEM_NAME}.png",
        *[
            f"panels/{STEM_NAME}_Panel_{letter}_{slug}.{extension}"
            for letter, slug in PANELS.items()
            for extension in ("pdf", "png")
        ],
        *[
            f"panels/{STEM_NAME}_Panel_{letter}_{slug}_spanish.{extension}"
            for letter, slug in PANELS.items()
            for extension in ("pdf", "png")
        ],
    ]
    missing_graphics = [
        name for name in expected_graphics if not (staging_output / name).is_file()
    ]
    if missing_graphics:
        raise RuntimeError(f"Missing rendered graphics: {missing_graphics}")
    manifest = {
        "figure": STEM_NAME,
        "analysis_status": "exploratory cage-level one-way ANOVA with exact pairwise MWU",
        "analysis_directory": str(args.analysis_dir),
        "output_directory": str(final_output),
        "workbook": str(args.workbook),
        "workbook_sha256": sha256(args.workbook),
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha256(Path(__file__).resolve()),
        "dpi": args.dpi,
        "subpanels": list(PANELS),
        "formats": ["pdf", "png"],
        "language_contract": {
            "combined_figure": "English only",
            "standalone_panels": ["English", "Spanish"],
            "spanish_translation_receipt": str(translation_receipt_path.relative_to(staging_output)),
            "spanish_translation_receipt_sha256": sha256(translation_receipt_path),
            "crop_geometry": "Spanish and English pairs use identical fixed crop boxes",
        },
        "panel_export_method": "exact crops of live axes in the combined rendered figure",
        "experimental_units": {
            "intake": "cage",
            "body_weight": "cage mean of two mouse changes",
        },
        "statistics": {
            "global": "ordinary equal-variance one-way ANOVA with cage as experimental unit",
            "pairwise": "raw exact exhaustive two-sided Mann-Whitney U; W-S, W-A, S-A",
            "multiplicity_adjustment": "none",
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
    print(f"Wrote Figure S4 and five live-cropped panels to {final_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
