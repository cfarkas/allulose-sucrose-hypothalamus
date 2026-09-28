#!/usr/bin/env python3
"""Analyze Figure S1: male single-bottle outcomes.

This is an exploratory, sex-stratified companion to Figure 1.  Cage consumption
is calculated from bottle-volume change and assigned to sex only after verifying
that every January cage is single-sex.  Body-weight trajectories retain mouse
rows for display, but all treatment summaries and tests use cages. Day-6
omnibus p values use ordinary one-way ANOVA; exhaustive cage-label omnibus
results are retained as sensitivity analyses and pairwise tests remain exact.
This script does not test a treatment-by-sex interaction.
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
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "FigS1"
DEFAULT_RAW = ANALYSIS_ROOT / "raw"
DEFAULT_OUT = ANALYSIS_ROOT / "results"
DEFAULT_FIGURE1_REFERENCE = PAPER_ROOT / "analyses" / "Fig1" / "results"
GENERATED_TREE_MARKER = ".paper_generated_tree.json"
OUTPUT_PRODUCER = "Paper/scripts/FigS1/01_analyze_single_bottle_male.py"
FIGURE_NUMBER = "S1"
TARGET_SEX = "Male"

CONDITIONS = ("Water", "Sucrose", "Allulose")
SEXES = ("Female", "Male")
JANUARY_CAGES = tuple(f"E{number}" for number in range(1, 13))
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
EXPECTED_CAGE_COUNTS = {
    ("Female", "Water"): 1,
    ("Female", "Sucrose"): 2,
    ("Female", "Allulose"): 2,
    ("Male", "Water"): 2,
    ("Male", "Sucrose"): 3,
    ("Male", "Allulose"): 2,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consumption", type=Path, default=DEFAULT_RAW / "consumption.csv")
    parser.add_argument("--weights", type=Path, default=DEFAULT_RAW / "weights.csv")
    parser.add_argument(
        "--figure1-reference-dir", type=Path,
        default=DEFAULT_FIGURE1_REFERENCE,
        help="Figure 1 analysis used for the sex-split value correspondence audit.",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=20260816)
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


def cohort_for(cage: pd.Series, date: pd.Series) -> pd.Series:
    date_text = date.dt.strftime("%Y-%m-%d")
    january = cage.isin(JANUARY_CAGES) & date_text.str.startswith("2025-01", na=False)
    august = cage.isin(["E13", "E14"]) & date_text.str.startswith("2025-08", na=False)
    return pd.Series(
        np.select(
            [january, august],
            ["January_E1-E12", "August_E13-E14"],
            default="Unassigned_or_mismatched",
        ),
        index=cage.index,
        dtype="string",
    )


def natural_cage(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["_cage_number"] = pd.to_numeric(
        result["E"].astype("string").str.extract(r"(\d+)", expand=False), errors="coerce"
    )
    sort_columns = ["_cage_number"]
    for column in ("study_day", "subject_id"):
        if column in result.columns:
            sort_columns.append(column)
    return result.sort_values(sort_columns, kind="mergesort").drop(columns="_cage_number")


def holm_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
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


def exact_global(values_by_group: dict[str, np.ndarray]) -> tuple[float, int, float]:
    arrays = [np.asarray(values_by_group[condition], dtype=float) for condition in CONDITIONS]
    if any(array.size == 0 for array in arrays):
        raise ValueError("All three treatments are required within each sex")
    values = np.concatenate(arrays)
    observed_labels = np.concatenate(
        [[condition] * array.size for condition, array in zip(CONDITIONS, arrays)]
    )
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


def exact_pair(left: np.ndarray, right: np.ndarray) -> tuple[float, int, float]:
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
) -> tuple[float, float, str]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if array.size < 2:
        return float("nan"), float("nan"), "not_estimable_fewer_than_2_cages"
    draws = rng.choice(array, size=(repetitions, array.size), replace=True).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high), "percentile_cage_bootstrap"


def endpoint_statistics(
    endpoint: pd.DataFrame,
    value_column: str,
    outcome: str,
    sex: str,
    rng: np.random.Generator,
    repetitions: int,
) -> pd.DataFrame:
    subset = endpoint.loc[endpoint["Sex"].eq(sex)].copy()
    if subset["E"].duplicated().any():
        raise ValueError(f"{sex} {outcome}: endpoint must contain one row per cage")
    values_by_group = {
        condition: subset.loc[subset["Treatment"].eq(condition), value_column].to_numpy(float)
        for condition in CONDITIONS
    }
    statistic, anova_p, df_between, df_within = ordinary_one_way_anova(values_by_group)
    exact_p, global_permutations, exact_statistic = exact_global(values_by_group)
    if not np.isclose(statistic, exact_statistic, rtol=0.0, atol=1e-12):
        raise RuntimeError("Parametric and permutation omnibus F statistics disagree")
    common = {
        "Sex": sex,
        "outcome": outcome,
        "cohort": "January_E1-E12",
        "endpoint_day": 6,
        "experimental_unit": "cage",
        "analysis_role": "exploratory_sex_stratified_not_interaction",
    }
    rows: list[dict[str, object]] = [{
        **common,
        "test": "One_way_ANOVA",
        "comparison": "Water|Sucrose|Allulose",
        "estimate": statistic,
        "estimate_type": "F statistic",
        "p_raw": anova_p,
        "p_adjusted_holm": np.nan,
        "n_left": subset.shape[0],
        "n_right": np.nan,
        "permutations": np.nan,
        "df_between": df_between,
        "df_within": df_within,
        "p_value_method": "ordinary equal-variance one-way ANOVA F distribution",
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "ci_status": "not_applicable",
    }, {
        **common,
        "test": "exact cage-label permutation omnibus F",
        "comparison": "Water|Sucrose|Allulose",
        "estimate": exact_statistic,
        "estimate_type": "F statistic",
        "p_raw": exact_p,
        "p_adjusted_holm": np.nan,
        "n_left": subset.shape[0],
        "n_right": np.nan,
        "permutations": global_permutations,
        "df_between": df_between,
        "df_within": df_within,
        "p_value_method": "exact exhaustive cage-label allocation of the ANOVA F statistic",
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "ci_status": "not_applicable",
    }]
    pair_rows: list[dict[str, object]] = []
    for left, right in itertools.combinations(CONDITIONS, 2):
        p_value, pair_permutations, difference = exact_pair(
            values_by_group[left], values_by_group[right]
        )
        pair_rows.append({
            **common,
            "test": "exact cage-label permutation mean difference",
            "comparison": f"{left}|{right}",
            "estimate": difference,
            "estimate_type": f"mean({left})-mean({right})",
            "p_raw": p_value,
            "n_left": values_by_group[left].size,
            "n_right": values_by_group[right].size,
            "permutations": pair_permutations,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
            "ci_status": "not_computed_for_contrast",
        })
    adjusted = holm_adjust([float(row["p_raw"]) for row in pair_rows])
    for row, value in zip(pair_rows, adjusted):
        row["p_adjusted_holm"] = value
    rows.extend(pair_rows)
    for condition in CONDITIONS:
        values = values_by_group[condition]
        low, high, status = bootstrap_mean_ci(values, rng, repetitions)
        rows.append({
            **common,
            "test": "descriptive cage summary",
            "comparison": condition,
            "estimate": float(np.mean(values)),
            "estimate_type": "cage mean",
            "sd": float(np.std(values, ddof=1)) if values.size >= 2 else np.nan,
            "p_raw": np.nan,
            "p_adjusted_holm": np.nan,
            "n_left": values.size,
            "n_right": np.nan,
            "permutations": np.nan,
            "ci95_low": low,
            "ci95_high": high,
            "ci_status": status,
        })
    return pd.DataFrame(rows)


def prepare_consumption(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    raw["cohort"] = cohort_for(raw["E"], raw["Date"])
    raw["study_day"] = raw["Date"].dt.strftime("%Y-%m-%d").map(DATE_DAY_MAP)
    raw["is_bottle_volume"] = raw["MeasurementType"].eq("Bottle_Volume")

    volume = raw.loc[
        raw["cohort"].eq("January_E1-E12") & raw["is_bottle_volume"]
    ].copy()
    if volume.duplicated(["E", "Date"]).any():
        raise ValueError("Duplicate January cage/date Bottle_Volume rows")
    cage_design = volume.groupby("E", as_index=False).agg(
        n_sexes=("Sex", "nunique"),
        n_treatments=("Treatment", "nunique"),
        Sex=("Sex", "first"),
        Treatment=("Treatment", "first"),
    )
    if set(cage_design["E"]) != set(JANUARY_CAGES):
        raise ValueError("January bottle-volume data must contain exactly E1-E12")
    if not cage_design["n_sexes"].eq(1).all():
        raise ValueError("Cage consumption cannot be sex-stratified: a January cage is mixed-sex")
    if not cage_design["n_treatments"].eq(1).all():
        raise ValueError("A January cage has multiple bottle treatments")
    if not set(cage_design["Sex"]).issubset(set(SEXES)):
        raise ValueError("January cage sex must be Female or Male")

    volume = volume.sort_values(["E", "Date"], kind="mergesort")
    volume["previous_volume_ml"] = volume.groupby("E")["Value"].shift()
    volume["increase_from_previous_ml"] = volume["Value"] - volume["previous_volume_ml"]
    volume["possible_refill_reset_or_measurement_increase"] = volume[
        "increase_from_previous_ml"
    ].gt(2.0)
    baseline = volume.loc[volume["study_day"].eq(1), ["E", "Value"]].rename(
        columns={"Value": "baseline_volume_ml"}
    )
    if baseline["E"].duplicated().any() or set(baseline["E"]) != set(JANUARY_CAGES):
        raise ValueError("Every January cage must have exactly one Day-1 bottle-volume baseline")
    volume = volume.merge(baseline, on="E", how="left", validate="many_to_one")
    volume["cumulative_cage_volume_removed_ml"] = volume["baseline_volume_ml"] - volume["Value"]
    primary = volume.loc[volume["study_day"].isin(PRIMARY_DAYS)].copy()
    if primary.duplicated(["E", "study_day"]).any():
        raise ValueError("Consumption primary window contains duplicate cage/day rows")
    completeness = primary.groupby("E")["study_day"].nunique()
    if not completeness.eq(len(PRIMARY_DAYS)).all() or primary.shape[0] != 36:
        raise ValueError(f"Incomplete January consumption trajectories: {completeness.to_dict()}")
    return raw, natural_cage(primary)


def prepare_weights(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    raw["cohort"] = cohort_for(raw["E"], raw["Date"])
    raw["study_day"] = raw["Date"].dt.strftime("%Y-%m-%d").map(DATE_DAY_MAP)
    raw["invalid_nonpositive_weight"] = raw["BodyWeight"].isna() | raw["BodyWeight"].le(0)
    raw["primary_row"] = (
        raw["cohort"].eq("January_E1-E12")
        & raw["Treatment"].isin(CONDITIONS)
        & raw["study_day"].isin(PRIMARY_DAYS)
        & ~raw["invalid_nonpositive_weight"]
    )
    if raw.duplicated(["subject_id", "Date", "Treatment"]).any():
        raise ValueError("Duplicate subject/date/treatment body-weight rows")
    primary = raw.loc[raw["primary_row"]].copy()
    subject_design = primary.groupby("subject_id", as_index=False).agg(
        n_sexes=("Sex", "nunique"),
        n_cages=("E", "nunique"),
        n_treatments=("Treatment", "nunique"),
    )
    if not subject_design[["n_sexes", "n_cages", "n_treatments"]].eq(1).all().all():
        raise ValueError("A primary subject changes cage, sex, or treatment")
    baseline = primary.loc[primary["study_day"].eq(1), ["subject_id", "BodyWeight"]].rename(
        columns={"BodyWeight": "baseline_weight_g"}
    )
    if baseline["subject_id"].duplicated().any():
        raise ValueError("Duplicate Day-1 body-weight baseline")
    primary = primary.merge(baseline, on="subject_id", how="left", validate="many_to_one")
    primary["body_weight_change_g"] = primary["BodyWeight"] - primary["baseline_weight_g"]
    primary["body_weight_change_pct"] = (
        100.0 * primary["body_weight_change_g"] / primary["baseline_weight_g"]
    )
    completeness = primary.groupby("subject_id")["study_day"].nunique()
    if not completeness.eq(len(PRIMARY_DAYS)).all():
        raise ValueError(f"Incomplete primary mouse trajectories: {completeness.to_dict()}")

    cage_day = primary.groupby(["E", "Treatment", "Sex", "study_day"], as_index=False).agg(
        cage_mean_body_weight_change_pct=("body_weight_change_pct", "mean"),
        n_mice=("subject_id", "nunique"),
    )
    if cage_day.duplicated(["E", "study_day"]).any():
        raise ValueError("Weight cage-day table is not one row per cage/day")
    endpoint = cage_day.loc[cage_day["study_day"].eq(6)].copy()
    return raw, natural_cage(primary), natural_cage(cage_day), natural_cage(endpoint)


def make_cage_sex_audit(
    consumption: pd.DataFrame, weights: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cage in JANUARY_CAGES:
        c = consumption.loc[consumption["E"].eq(cage)]
        w = weights.loc[weights["E"].eq(cage)]
        c_sexes = sorted(set(c["Sex"].astype(str)))
        w_sexes = sorted(set(w["Sex"].astype(str)))
        c_treatments = sorted(set(c["Treatment"].astype(str)))
        w_treatments = sorted(set(w["Treatment"].astype(str)))
        rows.append({
            "E": cage,
            "consumption_sexes": ";".join(c_sexes),
            "weight_primary_sexes": ";".join(w_sexes),
            "single_sex_consumption": len(c_sexes) == 1,
            "single_sex_weight_primary": len(w_sexes) == 1,
            "sex_agrees_across_sources": c_sexes == w_sexes,
            "consumption_treatments": ";".join(c_treatments),
            "weight_primary_treatments": ";".join(w_treatments),
            "treatment_agrees_across_sources": c_treatments == w_treatments,
            "Sex": c_sexes[0] if len(c_sexes) == 1 else "Ambiguous",
            "Treatment": c_treatments[0] if len(c_treatments) == 1 else "Ambiguous",
            "n_primary_weight_mice": w["subject_id"].nunique(),
            "eligible_for_sex_stratified_cage_analysis": (
                len(c_sexes) == len(w_sexes) == 1
                and c_sexes == w_sexes
                and len(c_treatments) == len(w_treatments) == 1
                and c_treatments == w_treatments
            ),
        })
    audit = pd.DataFrame(rows)
    if not audit["eligible_for_sex_stratified_cage_analysis"].all():
        failures = audit.loc[~audit["eligible_for_sex_stratified_cage_analysis"], "E"].tolist()
        raise ValueError(f"Sex-stratification audit failed for cages: {failures}")
    observed = audit.groupby(["Sex", "Treatment"])["E"].nunique().to_dict()
    if observed != EXPECTED_CAGE_COUNTS:
        raise ValueError(f"Unexpected cage counts by sex/treatment: {observed}")
    return audit


def write_audit(
    outdir: Path,
    cage_audit: pd.DataFrame,
    consumption: pd.DataFrame,
    weights: pd.DataFrame,
    statistics: pd.DataFrame,
) -> None:
    lines = [
        f"# Figure {FIGURE_NUMBER} scientific and statistical audit",
        "",
        "## Decision",
        "",
        f"Figure {FIGURE_NUMBER} is an exploratory sex-stratified companion to Figure 1. It repeats the "
        f"January single-bottle outcomes within {TARGET_SEX.lower()} cages. It is not a "
        "treatment-by-sex interaction analysis and must not be interpreted as evidence that a "
        "treatment effect differs between sexes.",
        "",
        "## Eligibility of cage-level consumption for sex stratification",
        "",
        f"- All {cage_audit.shape[0]} included {TARGET_SEX.lower()} cages contain one recorded sex in the "
        "bottle-volume file and one recorded sex among primary weight mice.",
        "- Sex and treatment agree across the consumption and primary body-weight sources for "
        "every January cage. Therefore cage consumption can be labeled by cage sex without "
        "assigning a shared bottle measurement to individual mice.",
        "- Consumption is calculated from cage-bottle volume change and can include "
        "unmeasured leakage or evaporation; it is not individual intake.",
        "",
        "## Experimental unit and sample sizes",
        "",
        "- The experimental unit is the cage for consumption, cage-balanced body-weight "
        "summaries, confidence intervals, and all p values.",
        f"- {TARGET_SEX} cage N: " + ", ".join(
            f"{condition} {EXPECTED_CAGE_COUNTS[(TARGET_SEX, condition)]}"
            for condition in CONDITIONS
        ) + ".",
        f"- Mouse rows ({weights['subject_id'].nunique()} mice total) appear only as pale "
        "descriptive points in body-weight trajectories; mice do not replace cages for inference.",
        "",
        "## Statistics",
        "",
        "- The prespecified endpoint is Day 6. The displayed omnibus p value is an ordinary "
        "equal-variance one-way ANOVA on one value per cage. Exact cage-label enumeration of "
        "the same F statistic is retained as a sensitivity analysis.",
        "- Pairwise tests exhaustively enumerate treatment labels for the two groups being "
        "compared. Holm adjustment is applied to the three pairwise p values separately within "
        "each sex and outcome.",
        "- Pointwise trajectory intervals resample complete cage trajectories so the same sampled "
        "cage identities are used at Days 1, 3, and 6.",
        *( ["- Female Water has one cage. Its trajectory and endpoint 95% confidence interval are "
            "explicitly not estimable; no degenerate one-cage bootstrap interval is reported."]
           if TARGET_SEX == "Female" else [] ),
        "",
        "## Interpretation limits",
        "",
        "1. Treatment randomization documentation was not found, so label exchangeability is an "
        "assumption and all p values are exploratory.",
        "2. Within-sex cage counts are extremely small. Exact p values have low resolution and "
        "bootstrap intervals with two or three cages are descriptive and unstable.",
        "3. No treatment-by-sex interaction is estimated or tested. Comparing significance in "
        "females with significance in males would be invalid.",
        "4. Sex is assigned at cage level only because every January cage passed the single-sex "
        "cross-source audit recorded in cage_sex_assignment_audit.csv.",
        "5. Repeated observations are used for trajectories, not treated as independent rows in "
        "the inferential tests.",
        "6. All valid prespecified-window observations are retained; no IQR deletion is used.",
        "",
        "## Recorded omnibus results",
        "",
    ]
    omnibus = statistics.loc[statistics["test"].eq("One_way_ANOVA")]
    for row in omnibus.itertuples(index=False):
        exact_row = statistics.loc[
            statistics["Sex"].eq(row.Sex)
            & statistics["outcome"].eq(row.outcome)
            & statistics["test"].eq("exact cage-label permutation omnibus F")
        ].iloc[0]
        lines.append(
            f"- {row.Sex}, {row.outcome}: one-way ANOVA p={float(row.p_raw):.12g}; "
            f"exact cage-label sensitivity p={float(exact_row.p_raw):.12g} "
            f"({int(exact_row.permutations)} allocations)."
        )
    lines.append("")
    (outdir / f"FIGURE_{FIGURE_NUMBER}_STATISTICAL_AUDIT.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def validate_figure1_correspondence(
    consumption: pd.DataFrame, weights: pd.DataFrame, outdir: Path,
    reference_root: Path,
) -> None:
    """Prove that this sex split preserves Figure 1's plotted raw values."""
    consumption_reference = reference_root / "consumption_primary_values.csv"
    weight_reference = reference_root / "body_weight_primary_values.csv"
    report: dict[str, object] = {
        "sex": TARGET_SEX,
        "status": "not_checked_reference_missing",
        "references": [str(consumption_reference), str(weight_reference)],
    }
    if consumption_reference.is_file() and weight_reference.is_file():
        consumption_columns = [
            "E", "Treatment", "Sex", "study_day", "cumulative_cage_volume_removed_ml"
        ]
        weight_columns = [
            "E", "subject_id", "Treatment", "Sex", "study_day", "body_weight_change_pct"
        ]
        expected_consumption = pd.read_csv(consumption_reference)
        expected_consumption = expected_consumption.loc[
            expected_consumption["Sex"].eq(TARGET_SEX), consumption_columns
        ].sort_values(consumption_columns[:-1], kind="mergesort").reset_index(drop=True)
        observed_consumption = consumption[consumption_columns].sort_values(
            consumption_columns[:-1], kind="mergesort"
        ).reset_index(drop=True)
        expected_weights = pd.read_csv(weight_reference)
        expected_weights = expected_weights.loc[
            expected_weights["Sex"].eq(TARGET_SEX), weight_columns
        ].sort_values(weight_columns[:-1], kind="mergesort").reset_index(drop=True)
        observed_weights = weights[weight_columns].sort_values(
            weight_columns[:-1], kind="mergesort"
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(
            observed_consumption, expected_consumption, check_dtype=False,
            check_exact=False, rtol=0.0, atol=1e-12,
        )
        pd.testing.assert_frame_equal(
            observed_weights, expected_weights, check_dtype=False,
            check_exact=False, rtol=0.0, atol=1e-12,
        )
        report.update({
            "status": "keys_exact_numeric_match_atol_1e-12",
            "consumption_rows": int(observed_consumption.shape[0]),
            "body_weight_rows": int(observed_weights.shape[0]),
            "maximum_absolute_difference": {
                "cumulative_cage_volume_removed_ml": float(np.max(np.abs(
                    observed_consumption["cumulative_cage_volume_removed_ml"].to_numpy(float)
                    - expected_consumption["cumulative_cage_volume_removed_ml"].to_numpy(float)
                ))),
                "body_weight_change_pct": float(np.max(np.abs(
                    observed_weights["body_weight_change_pct"].to_numpy(float)
                    - expected_weights["body_weight_change_pct"].to_numpy(float)
                ))),
            },
            "reference_sha256": {
                "consumption_primary_values": sha256(consumption_reference),
                "body_weight_primary_values": sha256(weight_reference),
            },
        })
    (outdir / "figure1_value_correspondence.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    args = parse_args()
    args.consumption = args.consumption.resolve()
    args.weights = args.weights.resolve()
    args.figure1_reference_dir = args.figure1_reference_dir.resolve()
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

    consumption_raw, consumption_all = prepare_consumption(args.consumption)
    weights_raw, weights_all, weight_cage_day_all, weight_endpoint_all = prepare_weights(args.weights)
    cage_audit_all = make_cage_sex_audit(consumption_all, weights_all)
    consumption = consumption_all.loc[consumption_all["Sex"].eq(TARGET_SEX)].copy()
    weights = weights_all.loc[weights_all["Sex"].eq(TARGET_SEX)].copy()
    weight_cage_day = weight_cage_day_all.loc[
        weight_cage_day_all["Sex"].eq(TARGET_SEX)
    ].copy()
    weight_endpoint = weight_endpoint_all.loc[
        weight_endpoint_all["Sex"].eq(TARGET_SEX)
    ].copy()
    cage_audit = cage_audit_all.loc[cage_audit_all["Sex"].eq(TARGET_SEX)].copy()
    consumption_endpoint = consumption.loc[consumption["study_day"].eq(6)].copy()

    stats_frames: list[pd.DataFrame] = []
    stats_frames.append(endpoint_statistics(
        consumption_endpoint,
        "cumulative_cage_volume_removed_ml",
        "Day-6 cumulative cage consumption (mL/cage)",
        TARGET_SEX,
        np.random.default_rng(args.seed + 100),
        args.bootstrap,
    ))
    stats_frames.append(endpoint_statistics(
        weight_endpoint,
        "cage_mean_body_weight_change_pct",
        "Day-6 cage-mean body-weight change (%)",
        TARGET_SEX,
        np.random.default_rng(args.seed + 200),
        args.bootstrap,
    ))
    statistics = pd.concat(stats_frames, ignore_index=True, sort=False)

    consumption.to_csv(args.outdir / "single_bottle_consumption_primary_values.csv", index=False)
    consumption_endpoint.to_csv(
        args.outdir / "single_bottle_consumption_cage_endpoint_values.csv", index=False
    )
    weights.to_csv(args.outdir / "single_bottle_body_weight_mouse_values.csv", index=False)
    weight_cage_day.to_csv(
        args.outdir / "single_bottle_body_weight_cage_day_values.csv", index=False
    )
    weight_endpoint.to_csv(
        args.outdir / "single_bottle_body_weight_cage_endpoint_values.csv", index=False
    )
    statistics.to_csv(args.outdir / "single_bottle_endpoint_statistics.csv", index=False)
    cage_audit.to_csv(args.outdir / "cage_sex_assignment_audit.csv", index=False)

    excluded_weights = weights_raw.loc[
        weights_raw["Sex"].eq(TARGET_SEX) & ~weights_raw["primary_row"]
    ].copy()
    excluded_weights.to_csv(args.outdir / "excluded_weight_rows_qc.csv", index=False)
    raw_inventory = pd.DataFrame([
        {
            "role": "raw_consumption",
            "path": str(args.consumption),
            "sha256": sha256(args.consumption),
            "rows": pd.read_csv(args.consumption).shape[0],
        },
        {
            "role": "raw_weights",
            "path": str(args.weights),
            "sha256": sha256(args.weights),
            "rows": pd.read_csv(args.weights).shape[0],
        },
        {
            "role": "analysis_script",
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
            "rows": np.nan,
        },
    ])
    raw_inventory.to_csv(args.outdir / "analysis_sources.csv", index=False)
    write_audit(args.outdir, cage_audit, consumption, weights, statistics)
    validate_figure1_correspondence(
        consumption, weights, args.outdir, args.figure1_reference_dir
    )

    result_files = sorted(
        path.name for path in args.outdir.iterdir()
        if path.is_file() and path.name != "run_manifest.json"
    )
    manifest = {
        "analysis": f"Figure {FIGURE_NUMBER}: {TARGET_SEX.lower()} single-bottle outcomes",
        "analysis_status": "exploratory_low_resolution",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "raw_sources": {
            "consumption": {"path": str(args.consumption), "sha256": sha256(args.consumption)},
            "weights": {"path": str(args.weights), "sha256": sha256(args.weights)},
        },
        "seed": args.seed,
        "bootstrap_repetitions": args.bootstrap,
        "primary_cohort": "January_E1-E12",
        "primary_days": list(PRIMARY_DAYS),
        "sex": TARGET_SEX,
        "conditions": list(CONDITIONS),
        "cage_counts": {
            condition: EXPECTED_CAGE_COUNTS[(TARGET_SEX, condition)]
            for condition in CONDITIONS
        },
        "single_sex_gate": "all 12 January cages passed cross-source audit",
        "experimental_unit": "cage",
        "endpoint": "Day 6",
        "omnibus": "ordinary equal-variance one-way ANOVA within sex/outcome",
        "omnibus_sensitivity": "exact exhaustive cage-label allocation of the same F statistic",
        "pairwise": "exact exhaustive cage-label permutation; Holm within sex/outcome",
        "interaction_test": "not performed; no treatment-by-sex claim",
        "confidence_interval_limit": (
            "Female Water not estimable because n=1 cage"
            if TARGET_SEX == "Female" else "all treatment groups have at least 2 cages"
        ),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "outputs_present_before_manifest": result_files,
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    publish_staging_output(args.outdir, final_outdir)
    print(f"[ANALYSIS] {final_outdir}")
    print(
        f"[CAGE N] {TARGET_SEX} W/S/A="
        + "/".join(str(EXPECTED_CAGE_COUNTS[(TARGET_SEX, c)]) for c in CONDITIONS)
    )
    if TARGET_SEX == "Female":
        print("[LIMIT] Female Water 95% CI not estimable (one cage)")
    print("[INTERACTION] not tested")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
