#!/usr/bin/env python3
"""Exact spatial-occurrence analysis of ARC/ME c-FOS and c-FOS+POMC nuclei.

This is the Figure 3 spatial analysis applied to Figure 4's reconstructed
ARC/ME sections. Within each reconstructed section the DAPI centroids define
six covariance-normalized radial shells of equal DAPI density. c-FOS/DAPI and
c-FOS+POMC/DAPI occurrence rates are calculated in every shell, transformed by
arcsin(sqrt(p)), and compared between conditions by an exhaustive animal-label
PERMANOVA over every allocation that preserves the observed group sizes.

Two things differ from Figure 3, both forced by the data rather than chosen.

Figure 3 has one registered field per animal, so its shells come straight from
that field. Figure 4 has two to five reconstructed sections per animal, each in
its own registered coordinate frame, so pooling their centroids would mix
frames. Shells are therefore built inside each section, and an animal's shell
counts are the sums over its own sections. The six-shell vector is still one
vector per animal.

Figure 3 is balanced 3/3/3, so its enumeration is the 1,680 allocations of nine
animals. Figure 4 is unbalanced, and the two endpoints do not even share a
denominator: animals whose POMC channel is incomplete cannot contribute a
c-FOS+POMC rate and are dropped from that endpoint only. Each endpoint is
therefore enumerated exhaustively over its own animal set and group sizes.

Cells and shells are feature-extraction units only. The biological animal is
the experimental unit for every p-value reported here, and both endpoints stay
exploratory: the same DAPI centroids define every shell and denominator.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

PAPER_ROOT = next(
    path for path in Path(__file__).resolve().parents
    if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                                 and (path / "README.txt").is_file()))
)
SHARED_SCRIPTS = PAPER_ROOT / "scripts" / "shared"
if str(SHARED_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SHARED_SCRIPTS))
from spanish_panel_text import translate_figure_texts_to_spanish

DEFAULT_ANALYSIS = PAPER_ROOT / "analyses" / "Fig4" / "results" / "human_final_run"
PER_CELL_NAME = "per_cell_arc_me_inferential_measurements.csv"

CONDITIONS = ("Water", "Sucrose", "Allulose")
CONDITION_COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
MEAN_LINE_COLORS = {"Water": "#2f91b8", "Sucrose": "#c51b1f", "Allulose": "#159447"}
N_RADIAL_SHELLS = 6
ANALYSIS_VERSION = "fig4_spatial_distribution_v1.1_2026-08-27"
PLOT_WIDTH_RATIOS = (2.30, 1.55)
PLOT_WSPACE = 0.26

OUTPUT_NAMES = (
    "animal_spatial_counts.csv",
    "animal_radial_occurrence_features.csv",
    "section_radial_shell_audit.csv",
    "exact_spatial_permanova.csv",
    "exact_spatial_permutation_receipts.csv",
    "Figure4_Spatial_cFOS_occurrence.png",
    "Figure4_Spatial_cFOS_occurrence.pdf",
    "Figure4_Spatial_cFOS_occurrence_spanish.png",
    "Figure4_Spatial_cFOS_occurrence_spanish.pdf",
    "Figure4_Spatial_cFOS_POMC_occurrence.png",
    "Figure4_Spatial_cFOS_POMC_occurrence.pdf",
    "Figure4_Spatial_cFOS_POMC_occurrence_spanish.png",
    "Figure4_Spatial_cFOS_POMC_occurrence_spanish.pdf",
    "spatial_spanish_translation_receipt.json",
    "spatial_distribution_provenance.json",
    "output_manifest.csv",
)


@dataclass(frozen=True)
class EndpointSpec:
    endpoint: str
    label: str
    plot_title: str
    y_label: str
    output_stem: str
    requires_pomc_channel: bool


ENDPOINTS = (
    EndpointSpec(
        "cfos_occurrence",
        "c-FOS-positive nuclei / DAPI nuclei spatial occurrence in ARC and ME",
        "c-FOS spatial occurrence",
        "c-FOS⁺ nuclei / DAPI nuclei (%)",
        "Figure4_Spatial_cFOS_occurrence",
        False,
    ),
    EndpointSpec(
        "double_occurrence",
        "c-FOS and POMC double-positive nuclei / DAPI nuclei spatial occurrence in ARC and ME",
        "c-FOS⁺POMC⁺ spatial occurrence",
        "c-FOS⁺POMC⁺ nuclei / DAPI nuclei (%)",
        "Figure4_Spatial_cFOS_POMC_occurrence",
        True,
    ),
)

SPANISH_PANEL_TEXT = {
    "c-FOS spatial occurrence": "Distribución espacial de c-FOS",
    "c-FOS⁺POMC⁺ spatial occurrence": "Distribución espacial de c-FOS⁺POMC⁺",
    "c-FOS⁺ nuclei / DAPI nuclei (%)": "Núcleos c-FOS⁺ / núcleos DAPI (%)",
    "c-FOS⁺POMC⁺ nuclei / DAPI nuclei (%)": "Núcleos c-FOS⁺POMC⁺ / núcleos DAPI (%)",
    "Equal-DAPI-density radial shell": "Capa radial de igual densidad DAPI",
    "Exact PERMANOVA p": "PERMANOVA exacta p",
    "Condition means": "Medias por condición",
    "radial shell": "capa radial",
    "inner": "interior",
    "outer": "exterior",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS,
                        help=f"Directory holding {PER_CELL_NAME} from 02_analyze_pomc_cfos.py.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Fresh output directory. Refuses to overwrite without --force.")
    parser.add_argument("--dpi", type=int, default=600, help="Panel raster resolution.")
    parser.add_argument("--force", action="store_true",
                        help="Replace an existing output directory this script generated.")
    return parser


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def covariance_normalized_radial_shells(
    dapi_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Six equal-DAPI-density shells, identical in method to Figure 3's."""
    centered = dapi_points - np.mean(dapi_points, axis=0, keepdims=True)
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


def load_per_cell(analysis_dir: Path) -> pd.DataFrame:
    path = analysis_dir / PER_CELL_NAME
    if not path.is_file():
        raise SystemExit(f"Per-cell table not found: {path}")
    frame = pd.read_csv(path, low_memory=False)
    for column in ("animal_id", "cond", "section_index", "centroid_y", "centroid_x",
                   "region", "is_cfos", "is_pomc", "pomc_channel_complete",
                   "included_in_arc_me_analysis"):
        if column not in frame.columns:
            raise SystemExit(f"{PER_CELL_NAME} is missing required column {column}")
    frame = frame.loc[frame["included_in_arc_me_analysis"].astype(bool)]
    frame = frame.loc[frame["cond"].isin(CONDITIONS)].copy()
    if frame.empty:
        raise SystemExit("No ARC/ME cells in Water, Sucrose or Allulose")
    for column in ("is_cfos", "is_pomc", "pomc_channel_complete"):
        frame[column] = frame[column].astype(bool)
    frame["animal_id"] = frame["animal_id"].astype(str)
    return frame


def extract_features(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-section shells summed into one six-shell vector per animal."""
    audit_rows: list[dict[str, Any]] = []
    shell_counts: dict[tuple[str, int], dict[str, int]] = {}
    animal_condition: dict[str, str] = {}
    animal_pomc: dict[str, bool] = {}

    for (animal_id, section_index), section in frame.groupby(
        ["animal_id", "section_index"], sort=True
    ):
        condition = str(section["cond"].iloc[0])
        animal_condition.setdefault(animal_id, condition)
        if animal_condition[animal_id] != condition:
            raise RuntimeError(f"{animal_id} appears under two conditions")
        animal_pomc[animal_id] = bool(
            animal_pomc.get(animal_id, True) and section["pomc_channel_complete"].all()
        )
        points = section[["centroid_y", "centroid_x"]].to_numpy(dtype=np.float64)
        if points.shape[0] < N_RADIAL_SHELLS:
            raise RuntimeError(
                f"{animal_id} section {section_index} has fewer nuclei than shells"
            )
        shell_index, edges, eigenvalues = covariance_normalized_radial_shells(points)
        is_cfos = section["is_cfos"].to_numpy(dtype=bool)
        is_double = is_cfos & section["is_pomc"].to_numpy(dtype=bool)
        denominators: list[int] = []
        for shell_zero in range(N_RADIAL_SHELLS):
            selected = shell_index == shell_zero
            denominator = int(np.count_nonzero(selected))
            if denominator == 0:
                raise RuntimeError(
                    f"Empty radial shell {shell_zero + 1}: {animal_id}/{section_index}"
                )
            denominators.append(denominator)
            bucket = shell_counts.setdefault(
                (animal_id, shell_zero + 1),
                {"dapi": 0, "cfos_occurrence": 0, "double_occurrence": 0},
            )
            bucket["dapi"] += denominator
            bucket["cfos_occurrence"] += int(np.count_nonzero(is_cfos & selected))
            bucket["double_occurrence"] += int(np.count_nonzero(is_double & selected))
        if max(denominators) - min(denominators) > 1:
            raise RuntimeError(
                f"DAPI quantile shells are not equal-density for "
                f"{animal_id}/{section_index}: {denominators}"
            )
        audit_rows.append({
            "animal_id": animal_id,
            "condition": condition,
            "section_index": int(section_index),
            "section_dapi_nuclei": int(points.shape[0]),
            "shell_min_dapi_nuclei": min(denominators),
            "shell_max_dapi_nuclei": max(denominators),
            "whitening_eigenvalue_small": float(np.min(eigenvalues)),
            "whitening_eigenvalue_large": float(np.max(eigenvalues)),
            "innermost_whitened_edge": float(edges[0]),
            "outermost_whitened_edge": float(edges[-1]),
            "shells_built_in": "this reconstructed section's own registered frame",
        })

    counts_rows: list[dict[str, Any]] = []
    feature_rows: list[dict[str, Any]] = []
    for animal_id in sorted(animal_condition):
        condition = animal_condition[animal_id]
        totals = {"dapi": 0, "cfos_occurrence": 0, "double_occurrence": 0}
        for shell in range(1, N_RADIAL_SHELLS + 1):
            bucket = shell_counts[(animal_id, shell)]
            for key in totals:
                totals[key] += bucket[key]
            for endpoint in ENDPOINTS:
                positive = bucket[endpoint.endpoint]
                occurrence = positive / bucket["dapi"]
                feature_rows.append({
                    "animal_id": animal_id,
                    "condition": condition,
                    "endpoint": endpoint.endpoint,
                    "endpoint_label": endpoint.label,
                    "radial_shell": shell,
                    "shell_dapi_nuclei": bucket["dapi"],
                    "shell_positive_nuclei": positive,
                    "occurrence": occurrence,
                    "occurrence_percent": 100.0 * occurrence,
                    "arcsin_sqrt_occurrence": float(np.arcsin(np.sqrt(occurrence))),
                    "pomc_channel_complete": animal_pomc[animal_id],
                    "inference_role": "one component of one animal-level six-shell vector",
                })
        counts_rows.append({
            "animal_id": animal_id,
            "condition": condition,
            "sections": int(frame.loc[frame["animal_id"].eq(animal_id), "section_index"].nunique()),
            "dapi_nuclei": totals["dapi"],
            "cfos_positive_nuclei": totals["cfos_occurrence"],
            "cfos_pomc_positive_nuclei": totals["double_occurrence"],
            "pomc_channel_complete": animal_pomc[animal_id],
            "experimental_unit": "biological animal",
        })

    counts = pd.DataFrame(counts_rows).sort_values(
        ["condition", "animal_id"], key=lambda s: s.map(
            {c: i for i, c in enumerate(CONDITIONS)}) if s.name == "condition" else s
    ).reset_index(drop=True)
    features = pd.DataFrame(feature_rows).sort_values(
        ["endpoint", "animal_id", "radial_shell"]).reset_index(drop=True)
    audit = pd.DataFrame(audit_rows).sort_values(
        ["animal_id", "section_index"]).reset_index(drop=True)
    return counts, features, audit


def feature_matrix(features: pd.DataFrame, animal_ids: Sequence[str], endpoint: str) -> np.ndarray:
    selected = features.loc[features["endpoint"].eq(endpoint)]
    pivot = selected.pivot(index="animal_id", columns="radial_shell",
                           values="arcsin_sqrt_occurrence").reindex(list(animal_ids))
    if list(pivot.columns) != list(range(1, N_RADIAL_SHELLS + 1)) or pivot.isna().any().any():
        raise RuntimeError(f"Incomplete six-shell feature matrix for {endpoint}")
    matrix = pivot.to_numpy(dtype=np.float64)
    if matrix.shape != (len(animal_ids), N_RADIAL_SHELLS) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid feature matrix for {endpoint}: {matrix.shape}")
    return matrix


def pseudo_f_and_r_squared(matrix: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    grand_mean = np.mean(matrix, axis=0)
    total_sum_squares = float(np.sum((matrix - grand_mean) ** 2))
    between_sum_squares = 0.0
    for group_index in range(len(CONDITIONS)):
        group = matrix[labels == group_index]
        if group.shape[0] == 0:
            raise ValueError("A PERMANOVA allocation contains an empty condition")
        between_sum_squares += float(
            group.shape[0] * np.sum((np.mean(group, axis=0) - grand_mean) ** 2)
        )
    within_sum_squares = total_sum_squares - between_sum_squares
    df_between = len(CONDITIONS) - 1
    df_within = matrix.shape[0] - len(CONDITIONS)
    pseudo_f = (math.inf if within_sum_squares <= 0
                else (between_sum_squares / df_between) / (within_sum_squares / df_within))
    r_squared = between_sum_squares / total_sum_squares if total_sum_squares > 0 else math.nan
    return float(pseudo_f), float(r_squared)


def enumerate_condition_allocations(group_sizes: Sequence[int]) -> Iterable[np.ndarray]:
    """Every labelling that preserves the observed group sizes, in order."""
    total = int(sum(group_sizes))
    indices = tuple(range(total))
    for first in itertools.combinations(indices, group_sizes[0]):
        remaining = tuple(index for index in indices if index not in first)
        for second in itertools.combinations(remaining, group_sizes[1]):
            labels = np.full(total, 2, dtype=np.int8)
            labels[list(first)] = 0
            labels[list(second)] = 1
            yield labels


def expected_allocation_count(group_sizes: Sequence[int]) -> int:
    total = int(sum(group_sizes))
    return (math.comb(total, group_sizes[0])
            * math.comb(total - group_sizes[0], group_sizes[1]))


def allocation_text(labels: np.ndarray, animal_ids: Sequence[str], group_index: int) -> str:
    return ";".join(a for a, l in zip(animal_ids, labels) if int(l) == group_index)


def benjamini_hochberg(raw_p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(raw_p_values, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or np.any((values < 0) | (values > 1)):
        raise ValueError("BH adjustment requires a non-empty one-dimensional p-value vector")
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = ranked * values.size / np.arange(1, values.size + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def exact_permanova(counts: pd.DataFrame,
                    features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    condition_to_index = {condition: index for index, condition in enumerate(CONDITIONS)}
    summary_rows: list[dict[str, Any]] = []
    receipt_rows: list[dict[str, Any]] = []

    for endpoint in ENDPOINTS:
        eligible = counts
        if endpoint.requires_pomc_channel:
            eligible = counts.loc[counts["pomc_channel_complete"]]
        animal_ids = eligible["animal_id"].astype(str).tolist()
        observed_labels = eligible["condition"].map(condition_to_index).to_numpy(dtype=np.int8)
        group_sizes = [int(np.count_nonzero(observed_labels == index))
                       for index in range(len(CONDITIONS))]
        if min(group_sizes) < 2:
            raise RuntimeError(
                f"{endpoint.endpoint}: every condition needs at least two animals, got {group_sizes}"
            )
        expected = expected_allocation_count(group_sizes)
        matrix = feature_matrix(features, animal_ids, endpoint.endpoint)
        observed_f, observed_r_squared = pseudo_f_and_r_squared(matrix, observed_labels)

        extreme = 0
        enumerated = 0
        observed_seen = 0
        for labeling_index, labels in enumerate(
            enumerate_condition_allocations(group_sizes), start=1
        ):
            permuted_f, permuted_r_squared = pseudo_f_and_r_squared(matrix, labels)
            enumerated += 1
            is_observed = bool(np.array_equal(labels, observed_labels))
            observed_seen += int(is_observed)
            as_extreme = bool(permuted_f >= observed_f - 1e-12)
            extreme += int(as_extreme)
            # The full receipt table for both endpoints is ~228,000 rows, which is
            # provenance rather than a result. Keep the observed allocation and the
            # ones that drive the p-value, which is what a reader checks.
            if is_observed or as_extreme:
                receipt_rows.append({
                    "endpoint": endpoint.endpoint,
                    "endpoint_label": endpoint.label,
                    "labeling_index": labeling_index,
                    "water_animals": allocation_text(labels, animal_ids, 0),
                    "sucrose_animals": allocation_text(labels, animal_ids, 1),
                    "allulose_animals": allocation_text(labels, animal_ids, 2),
                    "pseudo_f": permuted_f,
                    "r_squared": permuted_r_squared,
                    "is_observed_allocation": is_observed,
                    "is_as_or_more_extreme": as_extreme,
                    "permutation_unit": "biological animal",
                })
        if enumerated != expected:
            raise RuntimeError(
                f"{endpoint.endpoint}: expected {expected} allocations, enumerated {enumerated}"
            )
        if observed_seen != 1:
            raise RuntimeError(
                f"Observed allocation appears {observed_seen} times for {endpoint.endpoint}"
            )
        summary_rows.append({
            "endpoint": endpoint.endpoint,
            "endpoint_label": endpoint.label,
            "test": "exact animal-label PERMANOVA",
            "comparison": "Water|Sucrose|Allulose",
            "n_animals": len(animal_ids),
            "n_water": group_sizes[0],
            "n_sucrose": group_sizes[1],
            "n_allulose": group_sizes[2],
            "animals_excluded_incomplete_pomc": int(len(counts) - len(eligible)),
            "feature_count": N_RADIAL_SHELLS,
            "feature_transform": "arcsin(sqrt(shell positive nuclei / shell DAPI nuclei))",
            "shell_construction": "per reconstructed section, summed over an animal's sections",
            "distance": "Euclidean on six-shell animal vectors",
            "pseudo_f": observed_f,
            "r_squared": observed_r_squared,
            "p_value_exact": float(extreme / enumerated),
            "extreme_labelings": extreme,
            "enumerated_labelings": enumerated,
            "permutation_exhaustive": True,
            "experimental_unit": "biological animal",
            "cells_as_independent_replicates": False,
            "analysis_tier": "exploratory",
        })

    summary = pd.DataFrame(summary_rows)
    summary["bh_q_value_two_endpoint_family"] = benjamini_hochberg(
        summary["p_value_exact"].to_numpy(dtype=float)
    )
    summary["multiple_testing_family"] = (
        "two declared ARC/ME spatial occurrence profiles: c-FOS/DAPI and c-FOS+POMC/DAPI"
    )
    return summary, pd.DataFrame(receipt_rows)


def p_text(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}"
    if value < 0.1:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return f"{value:.3f}".rstrip("0").rstrip(".")


def configure_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.0, "axes.titlesize": 10.0,
        "axes.labelsize": 8.6, "xtick.labelsize": 7.4, "ytick.labelsize": 7.4,
        "legend.fontsize": 7.2, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.8, "savefig.dpi": "figure",
    })


def save_panel(features: pd.DataFrame, summary: pd.DataFrame,
               endpoint: EndpointSpec, output_dir: Path, dpi: int) -> dict[str, str]:
    configure_plot_style()
    selected = features.loc[features["endpoint"].eq(endpoint.endpoint)].copy()
    if endpoint.requires_pomc_channel:
        selected = selected.loc[selected["pomc_channel_complete"]]
    statistic = summary.loc[summary["endpoint"].eq(endpoint.endpoint)].iloc[0]

    fig = plt.figure(figsize=(6.15, 3.12), facecolor="white")
    # Give the 6-column condition heatmap materially more horizontal room. The
    # total panel aspect stays 6.15:3.12 so I/J still fit the publication grid
    # without stretching; only the internal allocation changes.
    grid = fig.add_gridspec(
        1, 2, width_ratios=PLOT_WIDTH_RATIOS, wspace=PLOT_WSPACE,
    )
    axis = fig.add_subplot(grid[0, 0])
    heat_axis = fig.add_subplot(grid[0, 1])
    # A little more bottom room than Figure 3 needs: these y-axis labels are
    # longer, so the legend row has to clear the two-line shell tick labels.
    fig.subplots_adjust(left=0.128, right=0.970, bottom=0.255, top=0.79)
    x_values = np.arange(1, N_RADIAL_SHELLS + 1, dtype=float)

    for condition in CONDITIONS:
        condition_rows = selected.loc[selected["condition"].eq(condition)]
        animal_profiles: list[np.ndarray] = []
        for animal_id in condition_rows["animal_id"].drop_duplicates().tolist():
            profile = (condition_rows.loc[condition_rows["animal_id"].eq(animal_id)]
                       .sort_values("radial_shell")["occurrence_percent"].to_numpy(dtype=float))
            if profile.size != N_RADIAL_SHELLS:
                raise RuntimeError(f"Incomplete plot profile for {animal_id}/{endpoint.endpoint}")
            animal_profiles.append(profile)
            axis.plot(x_values, profile, color=MEAN_LINE_COLORS[condition], alpha=0.20,
                      linewidth=0.8, marker="o", markersize=2.5, markeredgewidth=0.0, zorder=2)
        mean_profile = np.mean(np.stack(animal_profiles), axis=0)
        axis.plot(x_values, mean_profile, color=MEAN_LINE_COLORS[condition], linewidth=2.25,
                  marker="o", markersize=4.5, markerfacecolor=CONDITION_COLORS[condition],
                  markeredgecolor="#172126", markeredgewidth=0.55, zorder=4)

    ymax = float(selected["occurrence_percent"].max())
    axis.set_ylim(0.0, max(1.0, ymax * 1.20))
    axis.set_xlim(0.72, N_RADIAL_SHELLS + 0.28)
    axis.set_xticks(x_values)
    axis.set_xticklabels(["1\ninner", "2", "3", "4", "5", "6\nouter"])
    axis.set_xlabel("Equal-DAPI-density radial shell", labelpad=3.0)
    axis.set_ylabel(endpoint.y_label, labelpad=4.0)
    axis.set_title(endpoint.plot_title, loc="left", fontweight="bold", pad=8.0)
    axis.text(0.98, 0.955,
              (f"Exact PERMANOVA p = {p_text(float(statistic['p_value_exact']))}"
               f"  ·  BH q = {p_text(float(statistic['bh_q_value_two_endpoint_family']))}"),
              transform=axis.transAxes, ha="right", va="top", fontsize=7.3, fontweight="bold",
              color="#26353b",
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.7},
              zorder=8)
    axis.grid(axis="y", color="#d9dee1", linewidth=0.6, alpha=0.9)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(width=0.75, length=3.0, color="#30393d")

    legend_handles = [
        Line2D([0], [0], color=MEAN_LINE_COLORS[condition], marker="o",
               markerfacecolor=CONDITION_COLORS[condition], markeredgecolor="#172126",
               markeredgewidth=0.45, linewidth=2.0, markersize=4.0, label=condition)
        for condition in CONDITIONS
    ]
    axis.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
                ncol=3, frameon=False, columnspacing=1.3, handlelength=1.6)

    condition_means = np.vstack([
        selected.loc[selected["condition"].eq(condition)]
        .groupby("radial_shell", sort=True)["occurrence_percent"].mean()
        .reindex(range(1, N_RADIAL_SHELLS + 1)).to_numpy(dtype=float)
        for condition in CONDITIONS
    ])
    if condition_means.shape != (len(CONDITIONS), N_RADIAL_SHELLS):
        raise RuntimeError(f"Invalid descriptive heatmap for {endpoint.endpoint}")
    heat_image = heat_axis.imshow(condition_means, cmap="inferno", vmin=0.0,
                                  vmax=max(float(np.max(condition_means)), 1e-12),
                                  aspect="auto", interpolation="nearest")
    heat_axis.set_title("Condition means", fontsize=8.7, fontweight="bold", pad=7.0)
    heat_axis.set_xticks(np.arange(N_RADIAL_SHELLS))
    heat_axis.set_xticklabels(range(1, N_RADIAL_SHELLS + 1), fontsize=6.8)
    heat_axis.set_xlabel("radial shell", fontsize=7.2, labelpad=2.0)
    heat_axis.set_yticks(np.arange(len(CONDITIONS)))
    heat_axis.set_yticklabels(("W", "S", "A"), fontsize=7.4, fontweight="bold")
    for tick, condition in zip(heat_axis.get_yticklabels(), CONDITIONS):
        tick.set_color(MEAN_LINE_COLORS[condition])
    for row_index in range(len(CONDITIONS)):
        for column_index in range(N_RADIAL_SHELLS):
            value = float(condition_means[row_index, column_index])
            red, green, blue, _a = heat_image.cmap(heat_image.norm(value))
            luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
            heat_axis.text(column_index, row_index, f"{value:.1f}", ha="center", va="center",
                           fontsize=5.8, fontweight="bold",
                           color="#172126" if luminance >= 0.56 else "white")
    heat_axis.set_xticks(np.arange(-0.5, N_RADIAL_SHELLS, 1), minor=True)
    heat_axis.set_yticks(np.arange(-0.5, len(CONDITIONS), 1), minor=True)
    heat_axis.grid(which="minor", color="white", linewidth=0.55, alpha=0.72)
    heat_axis.tick_params(which="minor", bottom=False, left=False)
    heat_axis.tick_params(which="major", length=0)
    for spine in heat_axis.spines.values():
        spine.set_color("#45565d")
        spine.set_linewidth(0.65)
    colorbar = fig.colorbar(heat_image, ax=heat_axis, fraction=0.085, pad=0.055)
    colorbar.set_label("%", fontsize=7.0, labelpad=1.5)
    colorbar.ax.tick_params(labelsize=6.2, width=0.6, length=2.0)
    colorbar.outline.set_linewidth(0.55)

    fig.savefig(output_dir / f"{endpoint.output_stem}.png", dpi=dpi, facecolor="white",
                metadata={"Software": ANALYSIS_VERSION, "Title": endpoint.plot_title})
    fig.savefig(output_dir / f"{endpoint.output_stem}.pdf", dpi=dpi, facecolor="white",
                metadata={"Title": endpoint.plot_title,
                          "Author": "Figure 4 reproducible analysis",
                          "Creator": ANALYSIS_VERSION,
                          "CreationDate": None, "ModDate": None})
    translation = translate_figure_texts_to_spanish(fig, extra=SPANISH_PANEL_TEXT)
    spanish_title = SPANISH_PANEL_TEXT.get(endpoint.plot_title, endpoint.plot_title)
    fig.savefig(output_dir / f"{endpoint.output_stem}_spanish.png", dpi=dpi, facecolor="white",
                metadata={"Software": ANALYSIS_VERSION, "Title": spanish_title})
    fig.savefig(output_dir / f"{endpoint.output_stem}_spanish.pdf", dpi=dpi, facecolor="white",
                metadata={"Title": spanish_title,
                          "Author": "Figure 4 reproducible analysis",
                          "Creator": ANALYSIS_VERSION,
                          "CreationDate": None, "ModDate": None})
    plt.close(fig)
    return translation


def main() -> int:
    args = build_parser().parse_args()
    if int(args.dpi) != 600:
        raise SystemExit(f"Figure 4 spatial panels are fixed at 600 dpi, got {args.dpi}")
    analysis_dir = args.analysis_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists():
        if not args.force:
            raise SystemExit(f"Output directory already exists: {output_dir} (use --force)")
        if not (output_dir / "output_manifest.csv").is_file():
            raise SystemExit(f"Refusing to --force a directory this script did not write: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    frame = load_per_cell(analysis_dir)
    counts, features, audit = extract_features(frame)
    summary, receipts = exact_permanova(counts, features)

    counts.to_csv(output_dir / "animal_spatial_counts.csv", index=False)
    features.to_csv(output_dir / "animal_radial_occurrence_features.csv", index=False)
    audit.to_csv(output_dir / "section_radial_shell_audit.csv", index=False)
    summary.to_csv(output_dir / "exact_spatial_permanova.csv", index=False)
    receipts.to_csv(output_dir / "exact_spatial_permutation_receipts.csv", index=False)
    translation_receipt: dict[str, dict[str, str]] = {}
    for endpoint in ENDPOINTS:
        translation_receipt[endpoint.output_stem] = save_panel(
            features, summary, endpoint, output_dir, args.dpi
        )
    translation_receipt_path = output_dir / "spatial_spanish_translation_receipt.json"
    translation_receipt_path.write_text(
        json.dumps(translation_receipt, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    per_cell_path = analysis_dir / PER_CELL_NAME
    provenance = {
        "schema": "fig4_spatial_distribution_provenance_v1",
        "analysis_version": ANALYSIS_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "method": (
            "Figure 3's spatial occurrence analysis applied to Figure 4's reconstructed "
            "ARC/ME sections: six covariance-normalized equal-DAPI-density radial shells, "
            "arcsin(sqrt(p)) shell occurrence, exhaustive animal-label PERMANOVA, "
            "Benjamini-Hochberg across the two endpoints"
        ),
        "shell_construction": (
            "built inside each reconstructed section, because each section is its own "
            "registered coordinate frame; an animal's shell counts are the sums over its "
            "own sections, giving one six-shell vector per animal"
        ),
        "input": {
            "per_cell_table": str(per_cell_path),
            "per_cell_sha256": sha256_file(per_cell_path),
            "cells_analyzed": int(len(frame)),
            "sections_analyzed": int(len(audit)),
            "animals_analyzed": int(len(counts)),
        },
        "experimental_unit": "biological animal",
        "analysis_tier": "exploratory",
        "raster_dpi": int(args.dpi),
        "master_language": "en",
        "isolated_spatial_panels_bilingual": True,
        "isolated_spatial_panel_locales": ["en", "es"],
        "plot_width_ratios_line_to_heatmap": list(PLOT_WIDTH_RATIOS),
        "plot_horizontal_spacing": PLOT_WSPACE,
        "plot_right_margin": 0.970,
        "spanish_translation_receipt": str(translation_receipt_path),
        "endpoints": [
            {
                "endpoint": row["endpoint"],
                "n_animals": int(row["n_animals"]),
                "group_sizes": [int(row["n_water"]), int(row["n_sucrose"]), int(row["n_allulose"])],
                "enumerated_labelings": int(row["enumerated_labelings"]),
                "p_value_exact": float(row["p_value_exact"]),
                "bh_q_value_two_endpoint_family": float(row["bh_q_value_two_endpoint_family"]),
                "pseudo_f": float(row["pseudo_f"]),
                "r_squared": float(row["r_squared"]),
            }
            for _, row in summary.iterrows()
        ],
    }
    (output_dir / "spatial_distribution_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    manifest = pd.DataFrame([
        {"name": name, "bytes": (output_dir / name).stat().st_size,
         "sha256": sha256_file(output_dir / name)}
        for name in OUTPUT_NAMES if name != "output_manifest.csv"
    ])
    manifest.to_csv(output_dir / "output_manifest.csv", index=False)

    missing = [name for name in OUTPUT_NAMES if not (output_dir / name).is_file()]
    if missing:
        raise SystemExit(f"Expected outputs are missing: {missing}")
    for _, row in summary.iterrows():
        print(f"[SPATIAL] {row['endpoint']}: n={row['n_animals']} "
              f"({row['n_water']}/{row['n_sucrose']}/{row['n_allulose']}), "
              f"{row['enumerated_labelings']} labelings, "
              f"pseudo-F={row['pseudo_f']:.4f}, p={row['p_value_exact']:.6f}, "
              f"BH q={row['bh_q_value_two_endpoint_family']:.6f}")
    print(f"[PASS] Figure 4 spatial occurrence analysis: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
