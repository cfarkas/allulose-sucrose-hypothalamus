#!/usr/bin/env python3
"""Audit and analyze the January single-bottle behavior experiment.

The primary drinking outcome is cage consumption from Day 1, calculated from
the change in bottle volume.
The primary body-weight outcome is the Day-6 cage mean of mouse-level percent
changes. Repeated observations are retained for trajectories but are not
treated as independent replicates. August E13/E14 records and all exclusions
remain visible in the QC outputs.
"""

from __future__ import annotations

import argparse
import atexit
import hashlib
import itertools
import json
import os
import platform
import shutil
import tempfile
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import scipy


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
PROJECT_ROOT = PAPER_ROOT.parent
FIGURE_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig1"
DEFAULT_DATA = FIGURE_ANALYSIS_ROOT / "raw"
DEFAULT_OUT = FIGURE_ANALYSIS_ROOT / "results"
GENERATED_TREE_MARKER = ".paper_generated_tree.json"
OUTPUT_PRODUCER = "Paper/scripts/Fig1/01_analyze_behavior_experiment1.py"
CONDITIONS = ("Water", "Sucrose", "Allulose")
JANUARY_CAGES = tuple(f"E{number}" for number in range(1, 13))
AUGUST_CAGES = ("E13", "E14")
PRIMARY_DAYS = (1, 3, 6)
DATE_DAY_MAP = {
    "2025-01-13": 1,
    "2025-01-15": 3,
    "2025-01-17": 6,
    "2025-01-20": 9,
    "2025-01-22": 11,
    "2025-08-18": 1,
    "2025-08-21": 3,
    "2025-08-24": 6,
}

# Optional sensitivity metadata copied from the historical plot_trends.py.
# It is absent from the raw CSV and conflicts with observed weight-record mice.
NOMINAL_CAGE_SIZES = {
    "E1": 2, "E2": 3, "E3": 3, "E4": 2, "E5": 3, "E6": 1,
    "E7": 2, "E8": 3, "E9": 2, "E10": 3, "E11": 2, "E12": 2,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consumption", type=Path, default=DEFAULT_DATA / "consumption.csv")
    parser.add_argument("--weights", type=Path, default=DEFAULT_DATA / "weights.csv")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--bootstrap", type=int, default=20_000)
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


def clean_text(series: pd.Series, unknown: str = "Unknown") -> pd.Series:
    cleaned = series.astype("string").str.strip()
    cleaned = cleaned.mask(cleaned.str.lower().isin(["", "nan", "na", "none"]))
    return cleaned.fillna(unknown)


def normalize_treatment(series: pd.Series) -> pd.Series:
    cleaned = clean_text(series)
    replacements = {
        "water": "Water", "agua": "Water", "control": "Control",
        "sucrose": "Sucrose", "sacarosa": "Sucrose", "sucrosa": "Sucrose",
        "allulose": "Allulose", "alulosa": "Allulose",
    }
    return cleaned.str.lower().map(replacements).fillna(cleaned)


def assign_cohort(cage: pd.Series, date: pd.Series) -> pd.Series:
    date_text = date.dt.strftime("%Y-%m-%d")
    january = cage.isin(JANUARY_CAGES) & date_text.str.startswith("2025-01", na=False)
    august = cage.isin(AUGUST_CAGES) & date_text.str.startswith("2025-08", na=False)
    return pd.Series(
        np.select(
            [january, august],
            ["January_E1-E12", "August_E13-E14"],
            default="Unassigned_or_mismatched",
        ),
        index=cage.index,
        dtype="string",
    )


def natural_cage_sort(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["_cage_number"] = pd.to_numeric(
        result["E"].astype("string").str.extract(r"(\d+)", expand=False), errors="coerce"
    )
    return result.sort_values(["_cage_number", "E"], kind="mergesort").drop(columns="_cage_number")


def joined_reasons(frame: pd.DataFrame, rules: Sequence[tuple[pd.Series, str]]) -> pd.Series:
    result: list[str] = []
    for index in frame.index:
        labels = [label for mask, label in rules if bool(mask.loc[index])]
        result.append(";".join(labels) if labels else "included_primary")
    return pd.Series(result, index=frame.index, dtype="string")


def holm_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if finite.size == 0:
        return adjusted
    order = finite[np.argsort(values[finite], kind="mergesort")]
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (order.size - rank) * values[index]))
        adjusted[index] = running
    return adjusted


def anova_f(values: np.ndarray, labels: np.ndarray) -> float:
    groups = [values[labels == label] for label in np.unique(labels)]
    overall = float(np.mean(values))
    ss_between = sum(group.size * (float(np.mean(group)) - overall) ** 2 for group in groups)
    ss_within = sum(float(np.sum((group - np.mean(group)) ** 2)) for group in groups)
    df_between = len(groups) - 1
    df_within = values.size - len(groups)
    if df_within <= 0 or ss_within <= 0:
        return float("inf") if ss_between > 0 else 0.0
    return (ss_between / df_between) / (ss_within / df_within)


def exact_global_permutation(values_by_group: dict[str, np.ndarray]) -> tuple[float, int, float]:
    """Enumerate all cage-label allocations preserving the three group sizes."""
    arrays = [np.asarray(values_by_group[group], dtype=float) for group in CONDITIONS]
    if any(array.size == 0 for array in arrays):
        raise ValueError("All three treatment groups are required for exact permutation")
    values = np.concatenate(arrays)
    observed_labels = np.concatenate([[group] * array.size for group, array in zip(CONDITIONS, arrays)])
    observed = anova_f(values, observed_labels)
    indices = tuple(range(values.size))
    exceed = 0
    total = 0
    for first in itertools.combinations(indices, arrays[0].size):
        first_set = set(first)
        remainder = tuple(index for index in indices if index not in first_set)
        for second in itertools.combinations(remainder, arrays[1].size):
            labels = np.full(values.size, CONDITIONS[2], dtype=object)
            labels[list(first)] = CONDITIONS[0]
            labels[list(second)] = CONDITIONS[1]
            exceed += int(anova_f(values, labels) >= observed - 1e-12)
            total += 1
    return exceed / total, total, observed


def ordinary_one_way_anova(
    values_by_group: dict[str, np.ndarray],
) -> tuple[float, float, int, int]:
    """Ordinary equal-variance one-way ANOVA on independent cage endpoints."""
    arrays = [np.asarray(values_by_group[group], dtype=float) for group in CONDITIONS]
    if any(array.size < 1 or not np.isfinite(array).all() for array in arrays):
        raise ValueError("One-way ANOVA requires finite values in all three treatment groups")
    values = np.concatenate(arrays)
    labels = np.concatenate([[group] * array.size for group, array in zip(CONDITIONS, arrays)])
    statistic = anova_f(values, labels)
    df_between = len(arrays) - 1
    df_within = values.size - len(arrays)
    if df_within <= 0:
        raise ValueError("One-way ANOVA has no residual degrees of freedom")
    p_value = float(scipy.stats.f.sf(statistic, df_between, df_within))
    return statistic, p_value, df_between, df_within


def exact_pairwise_permutation(left: np.ndarray, right: np.ndarray) -> tuple[float, int, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    pooled = np.concatenate([left, right])
    observed = float(np.mean(left) - np.mean(right))
    exceed = 0
    total = 0
    for left_indices in itertools.combinations(tuple(range(pooled.size)), left.size):
        mask = np.zeros(pooled.size, dtype=bool)
        mask[list(left_indices)] = True
        difference = float(np.mean(pooled[mask]) - np.mean(pooled[~mask]))
        exceed += int(abs(difference) >= abs(observed) - 1e-12)
        total += 1
    return exceed / total, total, observed


def bootstrap_mean_ci(
    values: Iterable[float], rng: np.random.Generator, repetitions: int
) -> tuple[float, float]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    draws = rng.choice(array, size=(repetitions, array.size), replace=True).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def exact_endpoint_statistics(
    endpoint: pd.DataFrame,
    value_column: str,
    outcome: str,
    role: str,
    rng: np.random.Generator,
    bootstrap: int,
) -> pd.DataFrame:
    if endpoint["E"].duplicated().any():
        raise ValueError(f"Endpoint {outcome} must contain one row per cage")
    values_by_group = {
        condition: endpoint.loc[endpoint["Treatment"].eq(condition), value_column].to_numpy(float)
        for condition in CONDITIONS
    }
    statistic, anova_p, df_between, df_within = ordinary_one_way_anova(values_by_group)
    exact_p, permutations, exact_statistic = exact_global_permutation(values_by_group)
    if not np.isclose(statistic, exact_statistic, rtol=0.0, atol=1e-12):
        raise RuntimeError("Parametric and permutation omnibus F statistics disagree")
    rows: list[dict[str, object]] = [{
        "outcome": outcome,
        "cohort": "January_E1-E12",
        "endpoint_day": 6,
        "experimental_unit": "cage",
        "analysis_role": role,
        "test": "One_way_ANOVA",
        "comparison": "Water|Sucrose|Allulose",
        "estimate": statistic,
        "estimate_type": "F statistic",
        "p_raw": anova_p,
        "p_adjusted_holm": np.nan,
        "n_left": endpoint.shape[0],
        "n_right": np.nan,
        "permutations": np.nan,
        "df_between": df_between,
        "df_within": df_within,
        "p_value_method": "ordinary equal-variance one-way ANOVA F distribution",
    }, {
        "outcome": outcome,
        "cohort": "January_E1-E12",
        "endpoint_day": 6,
        "experimental_unit": "cage",
        "analysis_role": role,
        "test": "exact cage-label permutation omnibus F",
        "comparison": "Water|Sucrose|Allulose",
        "estimate": exact_statistic,
        "estimate_type": "F statistic",
        "p_raw": exact_p,
        "p_adjusted_holm": np.nan,
        "n_left": endpoint.shape[0],
        "n_right": np.nan,
        "permutations": permutations,
        "df_between": df_between,
        "df_within": df_within,
        "p_value_method": "exact exhaustive cage-label allocation of the ANOVA F statistic",
    }]
    pair_rows: list[dict[str, object]] = []
    for left, right in itertools.combinations(CONDITIONS, 2):
        p_value, pair_permutations, difference = exact_pairwise_permutation(
            values_by_group[left], values_by_group[right]
        )
        pair_rows.append({
            "outcome": outcome,
            "cohort": "January_E1-E12",
            "endpoint_day": 6,
            "experimental_unit": "cage",
            "analysis_role": role,
            "test": "exact cage-label permutation mean difference",
            "comparison": f"{left}|{right}",
            "estimate": difference,
            "estimate_type": f"mean({left})-mean({right})",
            "p_raw": p_value,
            "n_left": values_by_group[left].size,
            "n_right": values_by_group[right].size,
            "permutations": pair_permutations,
        })
    for row, adjusted in zip(pair_rows, holm_adjust([row["p_raw"] for row in pair_rows])):
        row["p_adjusted_holm"] = adjusted
    rows.extend(pair_rows)
    for condition in CONDITIONS:
        values = values_by_group[condition]
        low, high = bootstrap_mean_ci(values, rng, bootstrap)
        rows.append({
            "outcome": outcome,
            "cohort": "January_E1-E12",
            "endpoint_day": 6,
            "experimental_unit": "cage",
            "analysis_role": role,
            "test": "descriptive cage summary",
            "comparison": condition,
            "estimate": float(np.mean(values)),
            "estimate_type": "cage mean",
            "sd": float(np.std(values, ddof=1)),
            "ci95_low": low,
            "ci95_high": high,
            "p_raw": np.nan,
            "p_adjusted_holm": np.nan,
            "n_left": values.size,
            "n_right": np.nan,
            "permutations": np.nan,
        })
    return pd.DataFrame(rows)


def prepare_consumption(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path)
    required = {"E", "Treatment", "Sex", "Genotype", "Date", "MeasurementType", "Value"}
    missing = sorted(required.difference(raw.columns))
    if missing:
        raise ValueError(f"Consumption CSV is missing columns: {missing}")
    raw.insert(0, "raw_row_id", np.arange(1, len(raw) + 1))
    raw["E"] = clean_text(raw["E"])
    raw["Treatment"] = normalize_treatment(raw["Treatment"])
    raw["Sex"] = clean_text(raw["Sex"]).str.title()
    raw["Genotype"] = clean_text(raw["Genotype"])
    raw["MeasurementType"] = clean_text(raw["MeasurementType"])
    raw["Date"] = pd.to_datetime(raw["Date"], errors="coerce")
    raw["Value"] = pd.to_numeric(raw["Value"], errors="coerce")
    raw["cohort"] = assign_cohort(raw["E"], raw["Date"])
    raw["study_day"] = raw["Date"].dt.strftime("%Y-%m-%d").map(DATE_DAY_MAP)
    raw["is_measured_bottle_volume"] = raw["MeasurementType"].eq("Bottle_Volume")
    raw["measured_zero_value"] = raw["Value"].eq(0)
    raw["legacy_zero_imputation_risk"] = (
        raw["cohort"].eq("August_E13-E14") & ~raw["is_measured_bottle_volume"]
    )
    raw["primary_consumption_row"] = (
        raw["cohort"].eq("January_E1-E12")
        & raw["is_measured_bottle_volume"]
        & raw["study_day"].isin(PRIMARY_DAYS)
    )
    raw["row_disposition"] = joined_reasons(raw, [
        (raw["cohort"].eq("August_E13-E14"), "separate_August_cohort_not_pooled"),
        (raw["cohort"].eq("Unassigned_or_mismatched"), "unassigned_or_mismatched_cohort"),
        (~raw["is_measured_bottle_volume"], "not_Bottle_Volume"),
        (~raw["study_day"].isin(PRIMARY_DAYS), "outside_Day1_3_6_window"),
    ])

    volume = raw.loc[
        raw["cohort"].eq("January_E1-E12") & raw["is_measured_bottle_volume"]
    ].copy()
    if volume.duplicated(["E", "Date"]).any():
        raise ValueError("Duplicate January cage/date Bottle_Volume rows")
    volume = volume.sort_values(["E", "Date"], kind="mergesort")
    volume["previous_volume_ml"] = volume.groupby("E")["Value"].shift()
    volume["increase_from_previous_ml"] = volume["Value"] - volume["previous_volume_ml"]
    volume["possible_refill_reset_or_measurement_increase"] = volume["increase_from_previous_ml"].gt(2.0)
    baseline = volume.loc[volume["study_day"].eq(1), ["E", "Value"]].rename(
        columns={"Value": "baseline_volume_ml"}
    )
    if baseline["E"].duplicated().any() or set(baseline["E"]) != set(JANUARY_CAGES):
        raise ValueError("Each January cage E1-E12 must have one Day-1 volume baseline")
    volume = volume.merge(baseline, on="E", how="left", validate="many_to_one")
    volume["cumulative_cage_volume_removed_ml"] = volume["baseline_volume_ml"] - volume["Value"]
    volume["nominal_cage_size_legacy"] = volume["E"].map(NOMINAL_CAGE_SIZES)
    volume["nominal_per_mouse_removed_ml_sensitivity"] = (
        volume["cumulative_cage_volume_removed_ml"] / volume["nominal_cage_size_legacy"]
    )
    volume["primary_window"] = volume["study_day"].isin(PRIMARY_DAYS)
    primary = volume.loc[volume["primary_window"]].copy()
    completeness = primary.groupby("E")["study_day"].nunique()
    if set(completeness.index) != set(JANUARY_CAGES) or not completeness.eq(len(PRIMARY_DAYS)).all():
        raise ValueError(f"Incomplete primary consumption series: {completeness.to_dict()}")
    if primary.shape[0] != len(JANUARY_CAGES) * len(PRIMARY_DAYS):
        raise ValueError("Primary consumption must have exactly one row per cage and day")
    return raw, volume, primary


def prepare_weights(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path)
    required = {"E", "Treatment", "Sex", "Genotype", "AnimalID", "Date", "BodyWeight"}
    missing = sorted(required.difference(raw.columns))
    if missing:
        raise ValueError(f"Weights CSV is missing columns: {missing}")
    raw.insert(0, "raw_row_id", np.arange(1, len(raw) + 1))
    raw["E"] = clean_text(raw["E"])
    raw["Treatment"] = normalize_treatment(raw["Treatment"])
    raw["Sex"] = clean_text(raw["Sex"]).str.title()
    raw["Genotype"] = clean_text(raw["Genotype"])
    raw["AnimalID"] = clean_text(raw["AnimalID"])
    raw["subject_id"] = raw["E"] + "|" + raw["AnimalID"]
    raw["Date"] = pd.to_datetime(raw["Date"], errors="coerce")
    # Author clarification, 2026-09-08: FR5-4 moved E2 -> E10 for dehydration
    # on January 17 AFTER the day-6 measurement. Keep source labels and records.
    raw["Treatment_source"] = raw["Treatment"]
    raw["animal_id_reconciled"] = raw["subject_id"]
    moved = raw["AnimalID"].eq("FR5-4") & raw["E"].isin(["E2", "E10"])
    raw.loc[moved, "animal_id_reconciled"] = "E2|FR5-4"
    control = raw["E"].eq("E10") & raw["AnimalID"].eq("FR6-2") & raw["Treatment"].eq("Control")
    raw.loc[control, "Treatment"] = "Water"
    after_transfer = moved & raw["E"].eq("E10") & raw["Date"].gt(pd.Timestamp("2025-01-17"))
    raw.loc[after_transfer, "Treatment"] = "Water"
    raw["author_clarification"] = ""
    raw.loc[control, "author_clarification"] = "Control label reconciled to Water in control cage E10"
    raw.loc[moved, "author_clarification"] = "Same animal; humane transfer for dehydration after day-6 measurement on 2025-01-17"
    raw["post_humane_transfer"] = after_transfer
    raw["BodyWeight"] = pd.to_numeric(raw["BodyWeight"], errors="coerce")
    raw["cohort"] = assign_cohort(raw["E"], raw["Date"])
    raw["study_day"] = raw["Date"].dt.strftime("%Y-%m-%d").map(DATE_DAY_MAP)
    raw["invalid_nonpositive_weight"] = raw["BodyWeight"].isna() | raw["BodyWeight"].le(0)
    raw["BodyWeight_clean"] = raw["BodyWeight"].mask(raw["invalid_nonpositive_weight"])
    raw["primary_treatment"] = raw["Treatment"].isin(CONDITIONS)
    raw["primary_body_weight_row"] = (
        raw["cohort"].eq("January_E1-E12")
        & raw["primary_treatment"]
        & raw["study_day"].isin(PRIMARY_DAYS)
        & ~raw["invalid_nonpositive_weight"]
    )
    raw["row_disposition"] = joined_reasons(raw, [
        (raw["cohort"].eq("August_E13-E14"), "separate_August_cohort_not_pooled"),
        (raw["cohort"].eq("Unassigned_or_mismatched"), "unassigned_or_mismatched_cohort"),
        (~raw["primary_treatment"], "Control_or_nonprimary_treatment"),
        (~raw["study_day"].isin(PRIMARY_DAYS), "outside_Day1_3_6_window"),
        (raw["invalid_nonpositive_weight"], "nonpositive_or_missing_weight"),
    ])
    if raw.duplicated(["subject_id", "Date", "Treatment"]).any():
        raise ValueError("Duplicate subject/date/treatment body-weight rows")

    primary = raw.loc[raw["primary_body_weight_row"]].copy()
    baseline = primary.loc[primary["study_day"].eq(1), ["subject_id", "BodyWeight_clean"]].rename(
        columns={"BodyWeight_clean": "baseline_weight_g"}
    )
    if baseline["subject_id"].duplicated().any():
        raise ValueError("Duplicate Day-1 body weights for an E|AnimalID subject")
    primary = primary.merge(baseline, on="subject_id", how="left", validate="many_to_one")
    primary["body_weight_change_g"] = primary["BodyWeight_clean"] - primary["baseline_weight_g"]
    primary["body_weight_change_pct"] = 100.0 * primary["body_weight_change_g"] / primary["baseline_weight_g"]
    completeness = primary.groupby("subject_id")["study_day"].nunique()
    if not completeness.eq(len(PRIMARY_DAYS)).all():
        raise ValueError(f"Incomplete January primary body-weight series: {completeness.to_dict()}")
    e9 = primary.loc[primary["E"].eq("E9")]
    if e9.empty or not e9["Genotype"].eq("Unknown").all():
        raise ValueError("E9 must be retained with genotype Unknown")
    exclusions = raw.loc[~raw["primary_body_weight_row"]].copy()
    return raw, primary, exclusions


def make_weight_cage_endpoint(primary: pd.DataFrame) -> pd.DataFrame:
    endpoint = primary.loc[primary["study_day"].eq(6)].copy()
    treatment_counts = endpoint.groupby("E")["Treatment"].nunique()
    if not treatment_counts.eq(1).all():
        raise ValueError("Primary weight endpoint has more than one treatment per cage")
    cage = endpoint.groupby(["E", "Treatment"], as_index=False).agg(
        cage_mean_body_weight_change_pct=("body_weight_change_pct", "mean"),
        cage_sd_body_weight_change_pct=("body_weight_change_pct", "std"),
        cage_mean_baseline_weight_g=("baseline_weight_g", "mean"),
        n_mice=("subject_id", "nunique"),
        sexes=("Sex", lambda values: ";".join(sorted(set(map(str, values))))),
        genotypes=("Genotype", lambda values: ";".join(sorted(set(map(str, values))))),
    )
    cage["cohort"] = "January_E1-E12"
    cage["study_day"] = 6
    cage["experimental_unit"] = "cage"
    return natural_cage_sort(cage)


def make_cage_design(consumption_raw: pd.DataFrame, weights_raw: pd.DataFrame) -> pd.DataFrame:
    cages = sorted(set(consumption_raw["E"]).union(weights_raw["E"]), key=lambda value: int(value[1:]))
    rows: list[dict[str, object]] = []
    for cage in cages:
        c = consumption_raw.loc[consumption_raw["E"].eq(cage)]
        w = weights_raw.loc[weights_raw["E"].eq(cage)]
        january_w = w.loc[w["cohort"].eq("January_E1-E12")]
        primary_w = january_w.loc[january_w["primary_body_weight_row"]]
        observed_all = january_w["subject_id"].nunique()
        observed_primary = primary_w["subject_id"].nunique()
        nominal = NOMINAL_CAGE_SIZES.get(cage, np.nan)
        rows.append({
            "E": cage,
            "cohort_consumption": ";".join(sorted(set(map(str, c["cohort"])))) if len(c) else "none",
            "cohort_weights": ";".join(sorted(set(map(str, w["cohort"])))) if len(w) else "none",
            "consumption_treatments": ";".join(sorted(set(map(str, c["Treatment"])))) if len(c) else "none",
            "weight_treatments": ";".join(sorted(set(map(str, w["Treatment"])))) if len(w) else "none",
            "consumption_genotypes": ";".join(sorted(set(map(str, c["Genotype"])))) if len(c) else "none",
            "weight_genotypes": ";".join(sorted(set(map(str, w["Genotype"])))) if len(w) else "none",
            "has_measured_bottle_volume": bool(c["is_measured_bottle_volume"].any()) if len(c) else False,
            "has_bottle_weight": bool(c["MeasurementType"].eq("Bottle_Weight").any()) if len(c) else False,
            "nominal_cage_size_legacy": nominal,
            "observed_january_weight_subjects_all_labels": observed_all,
            "observed_january_primary_condition_subjects": observed_primary,
            "nominal_matches_observed_all_labels": bool(nominal == observed_all) if pd.notna(nominal) else False,
            "included_consumption_primary": cage in JANUARY_CAGES,
            "included_body_weight_primary": bool(primary_w.shape[0]),
            "design_note": (
                "mixed mouse-level treatment labels within cage" if cage == "E10" else
                "cross-file genotype metadata disagreement; weight genotype retained as Unknown" if cage == "E9" else ""
            ),
        })
    return pd.DataFrame(rows)


def make_cohort_inventory(consumption_raw: pd.DataFrame, weights_raw: pd.DataFrame) -> pd.DataFrame:
    c = consumption_raw.assign(source_table="consumption", record_type=consumption_raw["MeasurementType"])
    c_summary = c.groupby(
        ["source_table", "record_type", "cohort"], dropna=False, as_index=False
    ).agg(
        n_rows=("E", "size"),
        n_cages=("E", "nunique"),
        first_date=("Date", "min"),
        last_date=("Date", "max"),
        n_measured_zeros=("Value", lambda values: int(pd.to_numeric(values, errors="coerce").eq(0).sum())),
    )
    c_summary["n_subjects"] = 0

    w = weights_raw.assign(source_table="weights", record_type="BodyWeight")
    w_summary = w.groupby(
        ["source_table", "record_type", "cohort"], dropna=False, as_index=False
    ).agg(
        n_rows=("E", "size"),
        n_cages=("E", "nunique"),
        n_subjects=("subject_id", "nunique"),
        first_date=("Date", "min"),
        last_date=("Date", "max"),
        n_measured_zeros=("BodyWeight", lambda values: int(pd.to_numeric(values, errors="coerce").eq(0).sum())),
    )
    columns = [
        "source_table", "record_type", "cohort", "n_rows", "n_cages", "n_subjects",
        "first_date", "last_date", "n_measured_zeros",
    ]
    return pd.concat([c_summary[columns], w_summary[columns]], ignore_index=True)


def write_audit(
    outdir: Path,
    consumption_raw: pd.DataFrame,
    consumption_primary: pd.DataFrame,
    consumption_stats: pd.DataFrame,
    weights_raw: pd.DataFrame,
    weights_primary: pd.DataFrame,
    weight_endpoint: pd.DataFrame,
    weight_stats: pd.DataFrame,
    cages: pd.DataFrame,
) -> None:
    cage_counts = consumption_primary.groupby("Treatment")["E"].nunique().to_dict()
    mouse_counts = weights_primary.groupby("Treatment")["subject_id"].nunique().to_dict()
    weight_cage_counts = weight_endpoint.groupby("Treatment")["E"].nunique().to_dict()
    intake_p = float(consumption_stats.loc[
        consumption_stats["test"].eq("One_way_ANOVA"), "p_raw"
    ].iloc[0])
    weight_p = float(weight_stats.loc[
        weight_stats["test"].eq("One_way_ANOVA"), "p_raw"
    ].iloc[0])
    intake_exact_p = float(consumption_stats.loc[
        consumption_stats["test"].eq("exact cage-label permutation omnibus F"), "p_raw"
    ].iloc[0])
    weight_exact_p = float(weight_stats.loc[
        weight_stats["test"].eq("exact cage-label permutation omnibus F"), "p_raw"
    ].iloc[0])
    mismatch = cages.loc[
        cages["included_consumption_primary"] & ~cages["nominal_matches_observed_all_labels"], "E"
    ].tolist()
    august_c = consumption_raw.loc[consumption_raw["cohort"].eq("August_E13-E14")]
    august_w = weights_raw.loc[weights_raw["cohort"].eq("August_E13-E14")]
    invalid_weights = int(weights_raw["invalid_nonpositive_weight"].sum())
    primary_increase = int(consumption_primary["possible_refill_reset_or_measurement_increase"].sum())
    e9_rows = int(weights_primary["E"].eq("E9").sum())
    text = f"""# Figure 1 scientific and statistical audit

## Decision

Figure 1 is an exploratory January-cohort analysis. It does not reproduce the legacy pooled-cohort MANOVA, row-level OLS, PCA, LDA, or HC3 mouse-level inference. Repeated rows are used for trajectories; each inferential endpoint contributes exactly one value per cage.

## Cohorts and source integrity

- January cohort: E1-E12, dated January 13-22, 2025. Primary window is Day 1, 3, and 6 (January 13, 15, and 17).
- August add-on: E13/E14, dated August 18, 21, and 24, 2025. It contains {len(august_c)} consumption-table rows, all Bottle_Weight, and {len(august_w)} mouse-weight rows. It contains zero measured Bottle_Volume rows and is not pooled with January.
- The raw consumption CSV contains {int(consumption_raw['measured_zero_value'].sum())} measured zero values. The apparent August zero-volume/intake values in the legacy combined output were join/fill artifacts: absent Bottle_Volume was converted to zero. This pipeline never imputes absent bottle volume.
- E9 is retained: {e9_rows} primary rows for one mouse, with genotype explicitly represented as Unknown. The consumption file separately labels its volume record NPY, so the cross-file genotype conflict is preserved rather than guessed away.
- Author clarification (8 September 2026): FR5-4 was transferred from allulose cage E2 to water-control cage E10 for dehydration on 17 January, after the day-6 measurement. Both cage records identify the same animal. The primary window retains its pre-transfer allulose observations; later observations are flagged as post-transfer. E10/FR6-2, originally labelled Control, is reconciled to Water. The original CSV and source labels are retained.
- Cages containing multiple mouse-level treatment labels are flagged in cage_design_and_size_qc.csv. For the body-weight endpoint, only valid rows matching the cage bottle treatment enter the cage-balanced primary summary.
- There are {invalid_weights} nonpositive/missing weight rows. They are flagged as invalid and never interpreted as 0-g mice. None is needed to complete the January Day 1/3/6 primary series.

## Consumption (panels A-B)

- Experimental unit: cage.
- Primary quantity: measured Day-1 bottle volume minus the current measured bottle volume, in mL per cage. It is reported as consumption from Day 1; without leakage/evaporation controls it must not be interpreted as individually ingested volume.
- Primary group sizes: Water {cage_counts.get('Water', 0)}, Sucrose {cage_counts.get('Sucrose', 0)}, Allulose {cage_counts.get('Allulose', 0)} cages.
- Day-6 omnibus inference: ordinary equal-variance one-way ANOVA on one value per cage, p={intake_p:.6g}. Exact enumeration of all 27,720 cage-label allocations is retained as a sensitivity analysis, p={intake_exact_p:.6g}; three exact pairwise tests use Holm correction.
- One primary-window reading (E5, Day 6) increases by more than 2 mL from the prior reading ({primary_increase} flagged row). It may reflect measurement variation or an undocumented intervention and is retained transparently; the endpoint is still baseline minus Day-6 volume.
- The optional per-mouse sensitivity divides each cage value once by historical nominal occupancy. It never duplicates cage values to mouse rows and is not the primary figure outcome.

## Body weight (panels C-D)

- Panel C displays {weights_primary['subject_id'].nunique()} mice descriptively, keyed by their primary-window cage and AnimalID, with a separate reconciled identity linking the transferred animal. Thick summaries are cage-balanced means, not mouse-level inferential estimates.
- Day-6 inference first averages mouse percent changes within each cage, then compares Water {weight_cage_counts.get('Water', 0)}, Sucrose {weight_cage_counts.get('Sucrose', 0)}, and Allulose {weight_cage_counts.get('Allulose', 0)} cages.
- Ordinary one-way ANOVA omnibus p={weight_p:.6g}; exact cage-label sensitivity p={weight_exact_p:.6g}. Three pairwise exact tests use Holm correction. No HC3 mouse-level model is used.
- Descriptive mouse counts are Water {mouse_counts.get('Water', 0)}, Sucrose {mouse_counts.get('Sucrose', 0)}, Allulose {mouse_counts.get('Allulose', 0)}. They do not replace cage N for inference.

## Design limitations and interpretation

1. Historical nominal cage sizes disagree with observed January weight-record subject counts for {', '.join(mismatch)}. Consequently, cage-total consumption is primary; nominal per-mouse results are sensitivity only.
2. Treatment allocation/randomization documentation was not found. Exact label-permutation p values therefore rely on exchangeability and are exploratory, not proof from a documented randomized design.
3. There are only 3-5 cages per treatment. Bootstrap intervals are descriptive, and exact p values are discrete.
4. Sex and genotype are sparse, imbalanced, and partly confounded with cage/treatment. No treatment-by-sex/genotype claim is supported here.
5. Repeated Day 1/3/6 observations are not entered into ordinary independent-row tests. Inference uses one prespecified Day-6 endpoint per cage.
6. Pairwise multiplicity is controlled by Holm across the three treatment contrasts separately for each outcome.
7. No IQR deletion is used. All valid prespecified-window observations are retained.
"""
    (outdir / "FIGURE1_STATISTICAL_AUDIT.md").write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.consumption = args.consumption.resolve()
    args.weights = args.weights.resolve()
    final_outdir = Path(os.path.abspath(os.path.expanduser(args.outdir)))
    final_outdir.parent.mkdir(parents=True, exist_ok=True)
    validate_replacement_target(final_outdir)
    args.outdir = Path(tempfile.mkdtemp(
        prefix=f".{final_outdir.name}.stage-", dir=final_outdir.parent
    ))
    atexit.register(shutil.rmtree, args.outdir, ignore_errors=True)
    if args.bootstrap < 1_000:
        raise ValueError("Use at least 1,000 bootstrap repetitions")
    for source in (args.consumption, args.weights):
        if not source.is_file():
            raise FileNotFoundError(source)
    rng = np.random.default_rng(args.seed)

    consumption_raw, consumption_all, consumption_primary = prepare_consumption(args.consumption)
    weights_raw, weights_primary, weight_exclusions = prepare_weights(args.weights)
    cage_design = make_cage_design(consumption_raw, weights_raw)
    cohort_inventory = make_cohort_inventory(consumption_raw, weights_raw)

    consumption_endpoint = consumption_primary.loc[consumption_primary["study_day"].eq(6)].copy()
    consumption_stats = exact_endpoint_statistics(
        consumption_endpoint,
        "cumulative_cage_volume_removed_ml",
        "Day-6 cumulative cage consumption (mL/cage)",
        "primary",
        rng,
        args.bootstrap,
    )
    consumption_sensitivity_stats = exact_endpoint_statistics(
        consumption_endpoint,
        "nominal_per_mouse_removed_ml_sensitivity",
        "Day-6 nominal per-mouse consumption (mL; sensitivity only)",
        "sensitivity_only_legacy_nominal_occupancy",
        rng,
        args.bootstrap,
    )
    weight_endpoint = make_weight_cage_endpoint(weights_primary)
    weight_stats = exact_endpoint_statistics(
        weight_endpoint,
        "cage_mean_body_weight_change_pct",
        "Day-6 cage-mean body-weight change (%)",
        "primary",
        rng,
        args.bootstrap,
    )

    consumption_raw.to_csv(args.outdir / "consumption_raw_normalized.csv", index=False)
    consumption_all.to_csv(args.outdir / "consumption_all_readings_qc.csv", index=False)
    consumption_primary.to_csv(args.outdir / "consumption_primary_values.csv", index=False)
    consumption_stats.to_csv(args.outdir / "consumption_primary_statistics.csv", index=False)
    consumption_primary.to_csv(
        args.outdir / "consumption_nominal_per_mouse_sensitivity_values.csv", index=False
    )
    consumption_sensitivity_stats.to_csv(
        args.outdir / "consumption_nominal_per_mouse_sensitivity_statistics.csv", index=False
    )
    weights_raw.to_csv(args.outdir / "body_weight_all_readings_qc.csv", index=False)
    weights_primary.to_csv(args.outdir / "body_weight_primary_values.csv", index=False)
    weight_exclusions.to_csv(args.outdir / "body_weight_exclusions_and_invalid_values.csv", index=False)
    weight_endpoint.to_csv(args.outdir / "body_weight_cage_endpoint_values.csv", index=False)
    weight_stats.to_csv(args.outdir / "body_weight_primary_statistics.csv", index=False)
    cage_design.to_csv(args.outdir / "cage_design_and_size_qc.csv", index=False)
    cohort_inventory.to_csv(args.outdir / "cohort_inventory.csv", index=False)

    factor_balance = (
        weights_primary[["subject_id", "E", "Treatment", "Sex", "Genotype"]]
        .drop_duplicates()
        .groupby(["Treatment", "Sex", "Genotype"], dropna=False)["subject_id"]
        .nunique().rename("n_animals").reset_index()
    )
    factor_balance.to_csv(args.outdir / "body_weight_factor_balance.csv", index=False)
    write_audit(
        args.outdir,
        consumption_raw,
        consumption_primary,
        consumption_stats,
        weights_raw,
        weights_primary,
        weight_endpoint,
        weight_stats,
        cage_design,
    )

    sources = pd.DataFrame([
        {"role": "raw_consumption", "path": str(args.consumption), "sha256": sha256(args.consumption)},
        {"role": "raw_weights", "path": str(args.weights), "sha256": sha256(args.weights)},
        {"role": "analysis_script", "path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
    ])
    sources.to_csv(args.outdir / "analysis_sources.csv", index=False)
    stale_hc3 = args.outdir / "body_weight_HC3_model_coefficients.csv"
    if stale_hc3.exists():
        stale_hc3.unlink()

    output_names = sorted(
        path.name for path in args.outdir.iterdir()
        if path.is_file() and path.name != "run_manifest.json"
    )
    manifest = {
        "analysis": "behavior experiment 1: January cohort, cage-level primary inference",
        "analysis_status": "exploratory",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "consumption_source": str(args.consumption),
        "consumption_sha256": sha256(args.consumption),
        "weights_source": str(args.weights),
        "weights_sha256": sha256(args.weights),
        "seed": args.seed,
        "bootstrap_repetitions": args.bootstrap,
        "primary_cohort": "January_E1-E12",
        "excluded_separate_cohort": "August_E13-E14",
        "primary_days": list(PRIMARY_DAYS),
        "conditions": list(CONDITIONS),
        "primary_consumption_unit": "cage",
        "primary_consumption_outcome": "measured cage consumption from Day 1",
        "primary_weight_unit": "cage",
        "primary_weight_outcome": "cage mean of mouse Day-6 percent changes",
        "pairwise_multiplicity": "Holm within each outcome",
        "nominal_per_mouse_analysis": "sensitivity only; historical occupancy conflicts with weight records",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "outputs_present_before_manifest": output_names,
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    publish_staging_output(args.outdir, final_outdir)
    print(f"[ANALYSIS] {final_outdir}")
    print(f"[AUDIT] {final_outdir / 'FIGURE1_STATISTICAL_AUDIT.md'}")
    print(f"[PRIMARY N] consumption cages={consumption_primary['E'].nunique()}; weight cages={weight_endpoint['E'].nunique()}; weight mice={weights_primary['animal_id_reconciled'].nunique()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
