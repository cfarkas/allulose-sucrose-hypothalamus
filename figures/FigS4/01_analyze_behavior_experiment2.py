#!/usr/bin/env python3
"""Audit and analyse Experiment 2 two-bottle behaviour from ``UCSC_MICE.xlsx``.

The publication analysis deliberately uses cages for every shared intake outcome.
Body weights are first reduced to one mean change per cage because treatment was
delivered to cages.  Every displayed endpoint uses the same inferential family:
an ordinary one-way ANOVA followed by raw exact two-sided Mann-Whitney U
comparisons (W-S, W-A, and S-A).  Inference is exploratory because there are
just 2 Water, 3 Allulose, and 3 Sucrose cages.

``Hoja 4`` is the intended daily-intake sheet.  Its first value in each intake
row is an initial bottle/food mass, not consumption, and is therefore retained
for provenance but excluded from intake summaries.  The primary window ends on
2025-10-31, reproducing the explicit late-date exclusion in Commands_2026.txt;
all-observed and refill-robust sensitivity outputs are also written.
"""

from __future__ import annotations

import argparse
import atexit
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import re
import shutil
import tempfile
import unicodedata
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import f_oneway


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
PROJECT_ROOT = PAPER_ROOT.parent
FIGURE_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "FigS4"
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
DEFAULT_LEGACY_ROOT = FIGURE_ANALYSIS_ROOT / "raw" / "reference_analysis_snapshot"
DEFAULT_OUTPUT = FIGURE_ANALYSIS_ROOT / "results"
TREATMENT_ORDER = ("Water", "Allulose", "Sucrose")
EXPECTED_CAGE_TREATMENT = {
    "J1": "Allulose", "J2": "Sucrose", "J3": "Allulose", "J4": "Sucrose",
    "J5": "Allulose", "J6": "Water", "J7": "Water", "J8": "Sucrose",
}
ENDPOINT_ORDER = ("solution", "water", "total_fluid", "solution_fraction", "food")
ENDPOINT_UNITS = {
    "solution": "g/cage/day",
    "water": "g/cage/day",
    "total_fluid": "g/cage/day",
    "solution_fraction": "fraction",
    "food": "g/cage/day",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--sheet", default="Hoja 4")
    parser.add_argument("--validation-sheet", default="Hoja 3")
    parser.add_argument("--primary-end-date", default="2025-10-31")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--legacy-root", type=Path, default=DEFAULT_LEGACY_ROOT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def simplify(value: Any) -> str:
    if value is None:
        return ""
    text = unicodedata.normalize("NFKD", str(value))
    text = "".join(char for char in text if not unicodedata.combining(char))
    return re.sub(r"\s+", " ", text.strip().lower())


def normalize_treatment(value: Any) -> str:
    text = simplify(value)
    if "alulosa" in text or "allulose" in text:
        return "Allulose"
    if "sucrosa" in text or "sucrose" in text or "sacarosa" in text:
        return "Sucrose"
    if "agua" in text or "water" in text:
        return "Water"
    return "Unknown"


def endpoint_from_label(value: Any) -> str | None:
    text = simplify(value)
    if "solucion" in text or "solution" in text:
        return "solution"
    if re.search(r"\bagua\b", text) or "water" in text:
        return "water"
    if "comida" in text or "food" in text:
        return "food"
    return None


def numeric(value: Any) -> float:
    if isinstance(value, bool) or value is None:
        return float("nan")
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value) if np.isfinite(value) else float("nan")
    text = simplify(value)
    if text in {"", "-", "--", "n/m", "nm", "na", "n/a", "nan"}:
        return float("nan")
    try:
        return float(text.replace(",", "."))
    except ValueError:
        return float("nan")


def as_date(value: Any) -> date | None:
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    return None


def find_date_header(ws: openpyxl.worksheet.worksheet.Worksheet) -> tuple[int, list[tuple[int, date]]]:
    best_row = -1
    best: list[tuple[int, date]] = []
    for row in range(1, ws.max_row + 1):
        candidates = []
        for col in range(1, ws.max_column + 1):
            parsed = as_date(ws.cell(row, col).value)
            if parsed is not None:
                candidates.append((col, parsed))
        if len(candidates) > len(best):
            best_row, best = row, candidates
    if best_row < 0 or len(best) < 2:
        raise ValueError(f"Could not locate a date header in sheet {ws.title!r}")
    return best_row, sorted(best)


def cage_starts(ws: openpyxl.worksheet.worksheet.Worksheet, after_row: int) -> list[int]:
    starts = []
    for row in range(after_row + 1, ws.max_row + 1):
        if re.fullmatch(r"J\d+", str(ws.cell(row, 1).value or "").strip().upper()):
            starts.append(row)
    return starts


def source_cells(ws: openpyxl.worksheet.worksheet.Worksheet) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in ws.iter_rows():
        for cell in row:
            if cell.value is None:
                continue
            value = cell.value
            if isinstance(value, (datetime, date)):
                rendered = value.isoformat()
            else:
                rendered = str(value)
            rows.append({
                "sheet": ws.title,
                "coordinate": cell.coordinate,
                "row": cell.row,
                "column": cell.column,
                "value": rendered,
                "data_type": cell.data_type,
            })
    return rows


def parse_hoja4(
    ws: openpyxl.worksheet.worksheet.Worksheet,
    primary_end: date,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    header_row, date_columns = find_date_header(ws)
    starts = cage_starts(ws, header_row)
    if not starts:
        raise ValueError("No J# cage blocks were found")
    dates = [item[1] for item in date_columns]
    baseline_date = min(dates)
    intake_rows: list[dict[str, Any]] = []
    weight_rows: list[dict[str, Any]] = []
    cages: list[dict[str, Any]] = []
    qc: list[dict[str, Any]] = []

    for index, start in enumerate(starts):
        stop = starts[index + 1] - 1 if index + 1 < len(starts) else ws.max_row
        cage = str(ws.cell(start, 1).value).strip().upper()
        treatment_values = [normalize_treatment(ws.cell(row, 2).value) for row in range(start, stop + 1)]
        treatment = next((item for item in treatment_values if item != "Unknown"), "Unknown")
        expected = EXPECTED_CAGE_TREATMENT.get(cage)
        if expected and treatment != expected:
            raise ValueError(f"{cage}: workbook treatment {treatment!r} does not match expected {expected!r}")

        animal_source_rows: list[int] = []
        measure_source_rows: dict[str, int] = {}
        for row in range(start, stop + 1):
            endpoint = endpoint_from_label(ws.cell(row, 3).value)
            if endpoint:
                measure_source_rows[endpoint] = row
                continue
            animal_id = str(ws.cell(row, 3).value or "").strip()
            plausible = [numeric(ws.cell(row, col).value) for col, _ in date_columns]
            if animal_id and sum(np.isfinite(value) and 5 <= value <= 80 for value in plausible) >= 2:
                animal_source_rows.append(row)
        missing_endpoints = sorted({"water", "solution", "food"}.difference(measure_source_rows))
        if missing_endpoints or len(animal_source_rows) != 2:
            raise ValueError(
                f"{cage}: expected two mouse rows and water/solution/food rows; "
                f"found {len(animal_source_rows)} mice and missing {missing_endpoints}"
            )
        euthanasia = next(
            (as_date(ws.cell(row, 15).value) for row in animal_source_rows if as_date(ws.cell(row, 15).value)),
            None,
        )
        cages.append({
            "cage": cage,
            "treatment": treatment,
            "n_mice": len(animal_source_rows),
            "baseline_date": baseline_date.isoformat(),
            "primary_end_date": primary_end.isoformat(),
            "euthanasia_date": euthanasia.isoformat() if euthanasia else "",
            "source_rows": f"{start}:{stop}",
        })

        for endpoint, row in measure_source_rows.items():
            label = str(ws.cell(row, 3).value or "").strip()
            values = [numeric(ws.cell(row, col).value) for col, _ in date_columns]
            finite_after_baseline = [value for value in values[1:] if np.isfinite(value)]
            if not (np.isfinite(values[0]) and values[0] > 80 and finite_after_baseline and np.nanmedian(finite_after_baseline) < 80):
                raise ValueError(f"{cage}/{endpoint}: could not verify first value as initialization mass")
            for (col, day), value in zip(date_columns, values):
                study_day = (day - baseline_date).days
                if study_day == 0:
                    period = "baseline_initialization"
                    included = False
                    status = "initial_mass_not_consumption"
                elif day <= primary_end:
                    period = "primary"
                    included = bool(np.isfinite(value) and value >= 0)
                    status = "included" if included else "missing"
                else:
                    period = "post_primary"
                    included = False
                    status = "available_sensitivity_only" if np.isfinite(value) and value >= 0 else "missing"
                record = {
                    "cage": cage,
                    "treatment": treatment,
                    "endpoint": endpoint,
                    "source_label": label,
                    "source_sheet": ws.title,
                    "source_cell": ws.cell(row, col).coordinate,
                    "date": day.isoformat(),
                    "study_day": study_day,
                    "analysis_period": period,
                    "recorded_value": value,
                    "unit": "g/cage/day" if study_day else "g initial mass",
                    "primary_included": included,
                    "qc_status": status,
                }
                intake_rows.append(record)
                if status == "missing":
                    qc.append({
                        "severity": "warning",
                        "scope": "intake",
                        "cage": cage,
                        "subject": "",
                        "endpoint": endpoint,
                        "date": day.isoformat(),
                        "source_cell": record["source_cell"],
                        "issue": "missing recorded daily value",
                        "analysis_action": "endpoint/date excluded; other endpoints retained",
                    })

        for row in animal_source_rows:
            animal_id = str(ws.cell(row, 3).value).strip().upper().replace("_", "-")
            color = str(ws.cell(row, 4).value or "").strip().title()
            observed = []
            for col, day in date_columns:
                value = numeric(ws.cell(row, col).value)
                if np.isfinite(value) and value > 0:
                    observed.append((day, value))
            if not observed:
                raise ValueError(f"{cage}/{animal_id}: no body weights")
            animal_baseline_date, baseline_weight = observed[0]
            for col, day in date_columns:
                value = numeric(ws.cell(row, col).value)
                finite = bool(np.isfinite(value) and value > 0)
                period = "primary" if day <= primary_end else "post_primary"
                weight_rows.append({
                    "mouse_id": animal_id,
                    "cage": cage,
                    "treatment": treatment,
                    "color_id": color,
                    "source_sheet": ws.title,
                    "source_cell": ws.cell(row, col).coordinate,
                    "date": day.isoformat(),
                    "study_day": (day - animal_baseline_date).days,
                    "analysis_period": period,
                    "weight_g": value if finite else np.nan,
                    "baseline_weight_g": baseline_weight,
                    "weight_change_g": value - baseline_weight if finite else np.nan,
                    "primary_included": bool(period == "primary" and finite),
                    "qc_status": "included" if finite else "missing",
                })
                if not finite:
                    qc.append({
                        "severity": "warning",
                        "scope": "body_weight",
                        "cage": cage,
                        "subject": animal_id,
                        "endpoint": "weight",
                        "date": day.isoformat(),
                        "source_cell": ws.cell(row, col).coordinate,
                        "issue": "missing body weight",
                        "analysis_action": "date excluded for this mouse only",
                    })

    cage_frame = pd.DataFrame(cages).sort_values("cage").reset_index(drop=True)
    if cage_frame["cage"].tolist() != sorted(EXPECTED_CAGE_TREATMENT):
        raise ValueError(f"Unexpected cages: {cage_frame['cage'].tolist()}")
    return (
        pd.DataFrame(intake_rows),
        pd.DataFrame(weight_rows),
        cage_frame,
        pd.DataFrame(qc),
    )


def add_derived_endpoints(intake: pd.DataFrame) -> pd.DataFrame:
    base = intake.loc[intake["endpoint"].isin(["solution", "water", "food"])].copy()
    output = base.rename(columns={"recorded_value": "value"})
    output["endpoint_interpretation"] = output["endpoint"].map({
        "solution": "sweetener-side bottle; second water bottle in Water controls",
        "water": "water-side bottle",
        "food": "shared cage food",
    })
    wide = base.pivot_table(
        index=["cage", "treatment", "date", "study_day", "analysis_period"],
        columns="endpoint",
        values="recorded_value",
        aggfunc="first",
        dropna=True,
        observed=True,
    ).reset_index()
    derived: list[dict[str, Any]] = []
    for row in wide.itertuples(index=False):
        solution = float(row.solution) if np.isfinite(row.solution) else np.nan
        water = float(row.water) if np.isfinite(row.water) else np.nan
        if row.study_day == 0:
            # Initial bottle masses are retained in the recorded-intake table,
            # but their sum/ratio must not masquerade as daily behavior.
            total = fraction = np.nan
        elif np.isfinite(solution) and np.isfinite(water):
            total = solution + water
            fraction = solution / total if total > 0 else np.nan
        else:
            total = fraction = np.nan
        for endpoint, value, interpretation in (
            ("total_fluid", total, "sum of both cage bottles"),
            (
                "solution_fraction", fraction,
                "solution-side share; for Water controls this is bottle-2 share, not sweetener preference",
            ),
        ):
            derived.append({
                "cage": row.cage,
                "treatment": row.treatment,
                "endpoint": endpoint,
                "source_label": "derived",
                "source_sheet": "Hoja 4",
                "source_cell": "derived from solution+water",
                "date": row.date,
                "study_day": row.study_day,
                "analysis_period": row.analysis_period,
                "value": value,
                "unit": ENDPOINT_UNITS[endpoint],
                "primary_included": bool(row.analysis_period == "primary" and np.isfinite(value)),
                "qc_status": (
                    "initialization_not_analyzed" if row.study_day == 0
                    else "included" if np.isfinite(value)
                    else "missing_component"
                ),
                "endpoint_interpretation": interpretation,
            })
    output["unit"] = output["endpoint"].map(ENDPOINT_UNITS)
    keep = [
        "cage", "treatment", "endpoint", "source_label", "source_sheet", "source_cell",
        "date", "study_day", "analysis_period", "value", "unit", "primary_included",
        "qc_status", "endpoint_interpretation",
    ]
    return pd.concat([output[keep], pd.DataFrame(derived)[keep]], ignore_index=True).sort_values(
        ["endpoint", "cage", "date"]
    ).reset_index(drop=True)


def complete_dates(frame: pd.DataFrame, expected_cages: int) -> list[str]:
    return sorted(
        date_value
        for date_value, subset in frame.groupby("date", observed=True)
        if subset.loc[np.isfinite(pd.to_numeric(subset["value"], errors="coerce")), "cage"].nunique() == expected_cages
    )


def summarize_cages(
    daily: pd.DataFrame,
    cages: pd.DataFrame,
    period: str,
    excluded_dates: Iterable[str] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expected_cages = cages["cage"].nunique()
    exclusions = set(excluded_dates)
    summaries: list[dict[str, Any]] = []
    date_rows: list[dict[str, Any]] = []
    for endpoint in ENDPOINT_ORDER:
        subset = daily.loc[daily["endpoint"].eq(endpoint)].copy()
        if period == "primary":
            subset = subset.loc[subset["analysis_period"].eq("primary")]
        elif period == "all_observed":
            subset = subset.loc[subset["study_day"].gt(0)]
        else:
            raise ValueError(period)
        subset = subset.loc[~subset["date"].isin(exclusions)]
        dates = complete_dates(subset, expected_cages)
        date_rows.extend({
            "analysis_set": period,
            "endpoint": endpoint,
            "date": item,
            "included_common_complete_date": True,
        } for item in dates)
        used = subset.loc[subset["date"].isin(dates) & np.isfinite(subset["value"])].copy()
        for row in cages.itertuples(index=False):
            values = used.loc[used["cage"].eq(row.cage), "value"].astype(float)
            summaries.append({
                "analysis_set": period,
                "endpoint": endpoint,
                "unit": ENDPOINT_UNITS[endpoint],
                "cage": row.cage,
                "treatment": row.treatment,
                "value": float(values.mean()) if len(values) else np.nan,
                "sd_across_days": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
                "n_common_days": int(len(values)),
                "dates_used": ";".join(dates),
                "experimental_unit": "cage",
            })
    return pd.DataFrame(summaries), pd.DataFrame(date_rows)


def group_descriptives(cage_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (analysis_set, endpoint, unit, treatment), subset in cage_summary.groupby(
        ["analysis_set", "endpoint", "unit", "treatment"], observed=True
    ):
        values = subset["value"].dropna().astype(float)
        rows.append({
            "analysis_set": analysis_set,
            "endpoint": endpoint,
            "unit": unit,
            "treatment": treatment,
            "n_cages": int(len(values)),
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
            "median": float(values.median()),
            "minimum": float(values.min()),
            "maximum": float(values.max()),
        })
    return pd.DataFrame(rows)


def average_ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks with deterministic tie handling."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
        start = stop
    return ranks


def one_way_anova(values: np.ndarray, labels: Sequence[str]) -> tuple[float, float, int, int]:
    """Ordinary equal-variance one-way ANOVA on independent cage summaries."""
    values = np.asarray(values, dtype=float)
    labels_array = np.asarray(labels, dtype=object)
    groups = [values[labels_array == treatment] for treatment in TREATMENT_ORDER]
    if any(len(group) < 2 for group in groups):
        raise ValueError("One-way ANOVA requires at least two finite cages per treatment")
    result = f_oneway(*groups)
    return float(result.statistic), float(result.pvalue), len(groups) - 1, len(values) - len(groups)


def exact_pairwise(values_a: np.ndarray, values_b: np.ndarray) -> tuple[float, float, int, int]:
    """Exact two-sided Mann-Whitney U by exhaustive group-label permutation."""
    pooled = np.concatenate([values_a, values_b])
    n_a = len(values_a)
    n_b = len(values_b)
    ranks = average_ranks(pooled)
    observed_u = float(ranks[:n_a].sum() - n_a * (n_a + 1) / 2.0)
    center = n_a * n_b / 2.0
    observed_distance = abs(observed_u - center)
    statistics: list[float] = []
    for indices_a in itertools.combinations(range(len(pooled)), n_a):
        u_value = float(ranks[list(indices_a)].sum() - n_a * (n_a + 1) / 2.0)
        statistics.append(u_value)
    extreme = sum(abs(item - center) >= observed_distance - 1e-12 for item in statistics)
    return observed_u, extreme / len(statistics), extreme, len(statistics)


def permutation_statistics(
    cage_summary: pd.DataFrame,
    analysis_label: str,
    endpoints: Sequence[str] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    endpoint_list = tuple(endpoints) if endpoints is not None else ENDPOINT_ORDER
    for endpoint in endpoint_list:
        subset = cage_summary.loc[cage_summary["endpoint"].eq(endpoint)].dropna(subset=["value"]).sort_values("cage")
        if subset.empty:
            raise ValueError(f"No cage-summary values for endpoint {endpoint!r} in {analysis_label!r}")
        unit = str(subset["unit"].iloc[0]) if "unit" in subset else ENDPOINT_UNITS.get(endpoint, "")
        values = subset["value"].to_numpy(float)
        labels = subset["treatment"].astype(str).tolist()
        statistic, pvalue, df_between, df_within = one_way_anova(values, labels)
        rows.append({
            "analysis_set": analysis_label,
            "endpoint": endpoint,
            "unit": unit,
            "scope": "global",
            "contrast": "Water vs Allulose vs Sucrose",
            "test": "One_way_ANOVA",
            "statistic": statistic,
            "statistic_definition": "ordinary equal-variance one-way ANOVA F on cage summaries",
            "p_value_raw": pvalue,
            "df_between": df_between,
            "df_within": df_within,
            "extreme_permutations": np.nan,
            "total_permutations": np.nan,
            "n_cages_group_1": int((subset["treatment"] == "Water").sum()),
            "n_cages_group_2": int((subset["treatment"] != "Water").sum()),
            "method": "ordinary parametric one-way ANOVA with cage as the experimental unit",
            "interpretation": "exploratory; small and unbalanced 2/3/3 cage design; assumptions are weakly assessable",
        })
        for first, second in (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose")):
            first_values = subset.loc[subset["treatment"].eq(first), "value"].to_numpy(float)
            second_values = subset.loc[subset["treatment"].eq(second), "value"].to_numpy(float)
            u_value, pair_p, pair_extreme, pair_total = exact_pairwise(first_values, second_values)
            rows.append({
                "analysis_set": analysis_label,
                "endpoint": endpoint,
                "unit": unit,
                "scope": "pairwise",
                "contrast": f"{first} - {second}",
                "test": "Mann_Whitney_U_two_sided",
                "statistic": u_value,
                "statistic_definition": f"Mann-Whitney U for {first} (group 1)",
                "p_value_raw": pair_p,
                "df_between": np.nan,
                "df_within": np.nan,
                "extreme_permutations": pair_extreme,
                "total_permutations": pair_total,
                "n_cages_group_1": len(first_values),
                "n_cages_group_2": len(second_values),
                "method": "exact exhaustive two-sided Mann-Whitney U label permutation of cage summaries",
                "interpretation": "exploratory; discrete and underpowered",
            })
    return pd.DataFrame(rows)


def mouse_summaries(weights: pd.DataFrame, primary_end: date) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = weights.loc[weights["primary_included"]].copy()
    primary["date_parsed"] = pd.to_datetime(primary["date"])
    rows = []
    for mouse_id, subset in primary.groupby("mouse_id", observed=True):
        subset = subset.sort_values("date_parsed")
        first, last = subset.iloc[0], subset.iloc[-1]
        rows.append({
            "mouse_id": mouse_id,
            "cage": first["cage"],
            "treatment": first["treatment"],
            "baseline_date": first["date"],
            "primary_final_date": last["date"],
            "baseline_weight_g": first["weight_g"],
            "primary_final_weight_g": last["weight_g"],
            "primary_weight_change_g": last["weight_g"] - first["weight_g"],
            "primary_weight_change_percent": 100.0 * (last["weight_g"] - first["weight_g"]) / first["weight_g"],
            "n_weight_measurements": len(subset),
            "measurement_unit": "individual mouse",
            "treatment_assignment_unit": "cage",
        })
    mouse = pd.DataFrame(rows).sort_values(["treatment", "cage", "mouse_id"])
    descriptive_rows = []
    for treatment, subset in mouse.groupby("treatment", observed=True):
        values = subset["primary_weight_change_g"].astype(float)
        descriptive_rows.append({
            "treatment": treatment,
            "n_mice": len(values),
            "n_cages": subset["cage"].nunique(),
            "mean_change_g": values.mean(),
            "sd_change_g": values.std(ddof=1),
            "median_change_g": values.median(),
            "minimum_change_g": values.min(),
            "maximum_change_g": values.max(),
            "inference": "descriptive only; mice share cage-delivered treatment",
        })
    cage = mouse.groupby(["cage", "treatment"], as_index=False, observed=True).agg(
        value=("primary_weight_change_g", "mean"), n_mice=("mouse_id", "nunique")
    )
    cage["analysis_set"] = "primary_weight_cluster_sensitivity"
    cage["endpoint"] = "weight_change"
    cage["unit"] = "g"
    cage["experimental_unit"] = "cage mean of two mouse changes"
    return mouse, pd.DataFrame(descriptive_rows), cage


def parse_plus_components(value: Any) -> list[float]:
    if not isinstance(value, str) or "+" not in value:
        return []
    return [float(item.replace(",", ".")) for item in re.findall(r"\d+(?:[\.,]\d+)?", value)]


def parse_validation_sheet(
    ws: openpyxl.worksheet.worksheet.Worksheet,
    hoja4_daily: pd.DataFrame,
) -> pd.DataFrame:
    """Compare Hoja 4 against candidate consecutive differences in Hoja 3.

    Refill cells such as ``35,5 + 134,8`` are intrinsically ambiguous.  We
    expose both a pre-refill and next-day candidate; discrepancies are audit
    flags, not automatic data replacement.
    """
    header_row, date_columns = find_date_header(ws)
    starts = cage_starts(ws, header_row)
    output: list[dict[str, Any]] = []
    h4_lookup = hoja4_daily.set_index(["cage", "endpoint", "date"])["recorded_value"].to_dict()
    for index, start in enumerate(starts):
        stop = starts[index + 1] - 1 if index + 1 < len(starts) else ws.max_row
        cage = str(ws.cell(start, 1).value).strip().upper()
        endpoint_rows = {
            endpoint_from_label(ws.cell(row, 3).value): row
            for row in range(start, stop + 1)
            if endpoint_from_label(ws.cell(row, 3).value)
        }
        for endpoint, row in endpoint_rows.items():
            raw = [ws.cell(row, col).value for col, _ in date_columns]
            for position in range(1, len(date_columns)):
                col, current_date = date_columns[position]
                previous = raw[position - 1]
                current = raw[position]
                previous_number = numeric(previous)
                current_number = numeric(current)
                previous_plus = parse_plus_components(previous)
                current_plus = parse_plus_components(current)
                candidate = np.nan
                rule = "not_comparable"
                if current_plus and np.isfinite(previous_number):
                    candidate = previous_number - current_plus[0]
                    rule = "pre_refill_remainder: previous - first component"
                elif previous_plus and len(previous_plus) >= 2 and np.isfinite(current_number):
                    candidate = previous_plus[-1] - current_number
                    rule = "post_refill_candidate: last component - current"
                elif np.isfinite(previous_number) and np.isfinite(current_number) and current_number <= previous_number:
                    candidate = previous_number - current_number
                    rule = "simple consecutive difference"
                elif np.isfinite(previous_number) and np.isfinite(current_number) and current_number > previous_number:
                    rule = "unannotated refill/increase; not comparable"
                h4_value = h4_lookup.get((cage, endpoint, current_date.isoformat()), np.nan)
                if np.isfinite(candidate) and np.isfinite(h4_value):
                    difference = float(h4_value - candidate)
                    status = "matches_candidate" if abs(difference) <= 0.11 else "discrepancy_or_refill_ambiguity"
                elif np.isfinite(h4_value):
                    difference = np.nan
                    status = "Hoja4_value_not_independently_checkable"
                elif np.isfinite(candidate):
                    difference = np.nan
                    status = "Hoja4_missing_despite_candidate"
                else:
                    difference = np.nan
                    status = "not_comparable"
                output.append({
                    "cage": cage,
                    "endpoint": endpoint,
                    "date": current_date.isoformat(),
                    "hoja3_previous_cell": ws.cell(row, date_columns[position - 1][0]).coordinate,
                    "hoja3_current_cell": ws.cell(row, col).coordinate,
                    "hoja3_previous_raw": str(previous),
                    "hoja3_current_raw": str(current),
                    "candidate_daily_value_g": candidate,
                    "candidate_rule": rule,
                    "hoja4_recorded_value_g": h4_value,
                    "hoja4_minus_candidate_g": difference,
                    "audit_status": status,
                    "analysis_action": "Hoja 4 retained as intended daily sheet; flagged results remain exploratory",
                })
    return pd.DataFrame(output)


def legacy_audit(legacy_root: Path) -> pd.DataFrame:
    script = legacy_root / "digest_ucsc_mice.py"
    outputs = legacy_root / "UCSC_MICE_sheet4_analysis"
    return pd.DataFrame([
        {
            "issue_id": "L01",
            "severity": "major",
            "legacy_behavior": "Initial bottle/food masses were converted to zero and kept as Day 1 response values.",
            "evidence": "infer_daily_values() sets the first >80 value to 0; command output reports eight dates including baseline.",
            "impact": "Creates an artificial zero-consumption day and dominates Day/interaction tests.",
            "current_resolution": "Initialization masses are retained in tidy provenance but excluded from every intake summary/statistic.",
        },
        {
            "issue_id": "L02",
            "severity": "major",
            "legacy_behavior": "Each cage-date intake was copied to both mice and OLS/MANOVA used the duplicated animal-date rows.",
            "evidence": "parse_ucsc_sheet() attaches crecord intake to every animal; ANOVA/MANOVA consume the unified frame.",
            "impact": "Pseudoreplication and anticonservative p-values (for example 10^-57 to 10^-68).",
            "current_resolution": "Shared solution/water/food outcomes are reduced to one mean per cage before exploratory inference.",
        },
        {
            "issue_id": "L03",
            "severity": "major",
            "legacy_behavior": "Repeated animal dates were treated as independent by OLS/MANOVA.",
            "evidence": "Formulas use C(Day)*C(Treatment) without animal/cage correlation for the OLS/MANOVA pages.",
            "impact": "Nominal residual degrees of freedom are much larger than the number of independent treatment units.",
            "current_resolution": "Primary tests use one time-averaged value per cage; trajectories are descriptive.",
        },
        {
            "issue_id": "L04",
            "severity": "major",
            "legacy_behavior": "The mixed model used categorical Day-by-Treatment effects with only eight cage clusters.",
            "evidence": "64 cage-date rows, 8 cages, and a 24-parameter fixed-effect structure; asymptotic LRTs use up to 21 df.",
            "impact": "Random-effect and chi-square approximations are unreliable; tiny p-values are not confirmatory.",
            "current_resolution": "No asymptotic mixed model is claimed; all 560 global cage-label allocations are enumerated exactly.",
        },
        {
            "issue_id": "L05",
            "severity": "moderate",
            "legacy_behavior": "The command manually excluded 2025-11-01 and 2025-11-02 without a recorded scientific rationale.",
            "evidence": "Commands_2026.txt and user_excluded_dates.csv list both dates.",
            "impact": "Outcome-sensitive date selection cannot be ruled out.",
            "current_resolution": "Primary output reproduces the command-defined 7-day window; all-observed sensitivity tables retain endpoint-complete late dates.",
        },
        {
            "issue_id": "L06",
            "severity": "major",
            "legacy_behavior": "The primary drinking variable compared the solution bottle in sweetener cages with one water bottle in Water cages.",
            "evidence": "choose_intake(..., drink_source='treatment') switches source by treatment.",
            "impact": "The outcome is not the same physical construct across groups and ignores the second bottle.",
            "current_resolution": "Solution-side, water-side, total fluid, and two-bottle fraction are reported separately; Water fraction is explicitly bottle-2 share.",
        },
        {
            "issue_id": "L07",
            "severity": "major",
            "legacy_behavior": "Daily values around solution-bottle refill notation were accepted without a reconciliation audit.",
            "evidence": "Hoja 3 contains plus-form refill cells; several Hoja 4 values disagree with simple pre/post-refill candidates.",
            "impact": "Some Sucrose daily solution values cannot be independently reconstructed from remaining-mass entries.",
            "current_resolution": "Every Hoja 3/Hoja 4 comparison is exported; refill-adjacent dates receive a separate robust sensitivity analysis.",
        },
        {
            "issue_id": "L08",
            "severity": "moderate",
            "legacy_behavior": "Many raw pairwise tests were unadjusted.",
            "evidence": "Legacy CSV explicitly records MultipleTestingCorrection=none.",
            "impact": "Reported raw p-values overstate evidence across endpoints/contrasts.",
            "current_resolution": "Raw exact two-sided MWU p-values are reported explicitly for W-S, W-A, and S-A; no adjusted pairwise claim is made.",
        },
        {
            "issue_id": "L09",
            "severity": "moderate",
            "legacy_behavior": "Mouse body weights were modeled without acknowledging cage-shared treatment delivery.",
            "evidence": "Animal random intercepts were used, but treatment was assigned through cage bottles.",
            "impact": "Individual measurements do not create independent treatment assignments.",
            "current_resolution": "Mouse changes are averaged within cage before the same one-way ANOVA/exact MWU inference used for intake endpoints.",
        },
        {
            "issue_id": "L10",
            "severity": "moderate",
            "legacy_behavior": "Sex and genotype were recorded as unknown in the behavior model.",
            "evidence": "Legacy unified_data.csv has Sex=unknown and Genotype=unknown.",
            "impact": "Balance/confounding cannot be evaluated from Hoja 4 alone.",
            "current_resolution": "No adjusted causal treatment claim is made; missing design covariates are a stated limitation.",
        },
    ]).assign(
        legacy_script=str(script.resolve()),
        legacy_outputs=str(outputs.resolve()),
    )


def write_manifest(output_dir: Path, files: Sequence[tuple[str, Path]]) -> None:
    rows = []
    for role, path in files:
        rows.append({
            "role": role,
            "path": str(path.resolve()),
            "exists": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else np.nan,
            "sha256": sha256(path) if path.is_file() else "",
        })
    pd.DataFrame(rows).to_csv(output_dir / "source_manifest.csv", index=False)


def prepare_staging_output(final_output: Path) -> Path:
    """Create an empty sibling results tree, leaving the live tree untouched."""
    final_output.parent.mkdir(parents=True, exist_ok=True)
    if final_output.is_symlink():
        raise RuntimeError(f"Refusing symlink output directory: {final_output}")
    if final_output.exists() and not final_output.is_dir():
        raise RuntimeError(f"Output path is not a directory: {final_output}")
    staging = Path(tempfile.mkdtemp(prefix=f".{final_output.name}.stage-", dir=final_output.parent))
    atexit.register(shutil.rmtree, staging, ignore_errors=True)
    return staging


def publish_staging_output(staging: Path, final_output: Path) -> None:
    """Atomically replace a complete results tree, restoring on failure."""
    links = [path for path in staging.rglob("*") if path.is_symlink()]
    if links:
        raise RuntimeError(f"Analysis output contains forbidden symlinks: {links}")
    backup = staging.with_name(staging.name + ".previous")
    try:
        if final_output.exists():
            os.replace(final_output, backup)
        os.replace(staging, final_output)
    except Exception:
        if backup.exists() and not final_output.exists():
            os.replace(backup, final_output)
        raise
    else:
        if backup.exists():
            shutil.rmtree(backup)


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.4f}"


def main() -> int:
    args = parse_args()
    args.workbook = args.workbook.resolve()
    final_output = args.output_dir.resolve()
    args.legacy_root = args.legacy_root.resolve()
    primary_end = date.fromisoformat(args.primary_end_date)
    if not args.workbook.is_file():
        raise FileNotFoundError(args.workbook)
    args.output_dir = prepare_staging_output(final_output)

    workbook = openpyxl.load_workbook(args.workbook, data_only=False, read_only=False)
    if args.sheet not in workbook.sheetnames or args.validation_sheet not in workbook.sheetnames:
        raise ValueError(f"Required sheets not found; available sheets: {workbook.sheetnames}")
    hoja4 = workbook[args.sheet]
    hoja3 = workbook[args.validation_sheet]
    intake, weights, cages, qc = parse_hoja4(hoja4, primary_end)
    daily = add_derived_endpoints(intake)
    validation = parse_validation_sheet(hoja3, intake)

    discrepancy = validation.loc[validation["audit_status"].eq("discrepancy_or_refill_ambiguity")].copy()
    for row in discrepancy.itertuples(index=False):
        qc = pd.concat([qc, pd.DataFrame([{
            "severity": "warning",
            "scope": "Hoja3_vs_Hoja4",
            "cage": row.cage,
            "subject": "",
            "endpoint": row.endpoint,
            "date": row.date,
            "source_cell": row.hoja3_current_cell,
            "issue": f"Hoja 4 differs from candidate consecutive/refill difference by {row.hoja4_minus_candidate_g:.3g} g",
            "analysis_action": "retain intended Hoja 4 value; flag and run refill-robust sensitivity",
        }])], ignore_index=True)

    primary_summary, primary_dates = summarize_cages(daily, cages, "primary")
    observed_summary, observed_dates = summarize_cages(daily, cages, "all_observed")
    # Exclude every primary date implicated in a Hoja3/Hoja4 discrepancy so all cages
    # retain the same dates within each endpoint comparison.
    refill_dates = sorted(set(discrepancy.loc[discrepancy["date"].le(primary_end.isoformat()), "date"]))
    robust_summary, robust_dates = summarize_cages(daily, cages, "primary", excluded_dates=refill_dates)
    robust_summary["analysis_set"] = "primary_refill_robust"
    robust_dates["analysis_set"] = "primary_refill_robust"

    primary_stats = permutation_statistics(primary_summary, "primary")
    observed_stats = permutation_statistics(observed_summary, "all_observed")
    robust_stats = permutation_statistics(robust_summary, "primary_refill_robust")
    all_stats = pd.concat([primary_stats, observed_stats, robust_stats], ignore_index=True)

    mouse_summary, weight_descriptive, weight_cage = mouse_summaries(weights, primary_end)
    weight_stats = permutation_statistics(
        weight_cage,
        "primary_weight_cluster_sensitivity",
        endpoints=("weight_change",),
    )
    weight_stats["interpretation"] = "exploratory primary inference on one mean body-weight change per cage"
    complete_stats = pd.concat([all_stats, weight_stats], ignore_index=True)

    audit = legacy_audit(args.legacy_root)
    cells = pd.DataFrame(source_cells(hoja4) + source_cells(hoja3))
    summaries = pd.concat([primary_summary, observed_summary, robust_summary], ignore_index=True)
    included_dates = pd.concat([primary_dates, observed_dates, robust_dates], ignore_index=True)
    descriptives = group_descriptives(summaries)

    intake.to_csv(args.output_dir / "hoja4_recorded_intake_long.csv", index=False)
    daily.to_csv(args.output_dir / "daily_endpoints_long.csv", index=False)
    weights.to_csv(args.output_dir / "mouse_weight_long.csv", index=False)
    cages.to_csv(args.output_dir / "cage_design.csv", index=False)
    summaries.to_csv(args.output_dir / "cage_endpoint_summaries.csv", index=False)
    included_dates.to_csv(args.output_dir / "endpoint_complete_dates.csv", index=False)
    descriptives.to_csv(args.output_dir / "group_descriptive_statistics.csv", index=False)
    complete_stats.to_csv(
        args.output_dir / "one_way_anova_mann_whitney_statistics.csv", index=False
    )
    mouse_summary.to_csv(args.output_dir / "mouse_weight_summary_primary.csv", index=False)
    weight_descriptive.to_csv(args.output_dir / "mouse_weight_group_descriptives.csv", index=False)
    weight_cage.to_csv(args.output_dir / "weight_cage_cluster_summary.csv", index=False)
    weight_stats.to_csv(
        args.output_dir / "weight_cage_anova_mann_whitney.csv", index=False
    )
    validation.to_csv(args.output_dir / "hoja3_hoja4_intake_validation.csv", index=False)
    qc.sort_values(["severity", "scope", "cage", "date"]).to_csv(args.output_dir / "qc_issues.csv", index=False)
    audit.to_csv(args.output_dir / "legacy_analysis_audit.csv", index=False)
    cells.to_csv(args.output_dir / "source_cells_long.csv", index=False)

    primary_global = primary_stats.loc[primary_stats["scope"].eq("global")]
    report_lines = [
        "Experiment 2 behavior: cage-aware audit and exploratory analysis",
        "",
        f"Workbook: {args.workbook}",
        f"Workbook SHA256: {sha256(args.workbook)}",
        f"Source sheet: {args.sheet}; validation sheet: {args.validation_sheet}",
        f"Primary intake window: day 1 through {primary_end.isoformat()} (7 daily observations/cage).",
        "The baseline date is an initialization mass and is not consumption.",
        f"Design: {cages['cage'].nunique()} cages / {mouse_summary['mouse_id'].nunique()} mice; "
        + ", ".join(
            f"{trt}={int((cages['treatment'] == trt).sum())} cages"
            for trt in TREATMENT_ORDER
        ),
        "",
        "Primary cage-level one-way ANOVA tests (exploratory):",
    ]
    for row in primary_global.itertuples(index=False):
        report_lines.append(
            f"  {row.endpoint}: F({int(row.df_between)},{int(row.df_within)})={row.statistic:.4g}, "
            f"p={format_p(row.p_value_raw)}"
        )
    weight_global = weight_stats.loc[weight_stats["scope"].eq("global")].iloc[0]
    report_lines.append(
        f"  weight_change: F({int(weight_global.df_between)},{int(weight_global.df_within)})="
        f"{weight_global.statistic:.4g}, p={format_p(weight_global.p_value_raw)}"
    )
    report_lines.extend([
        "",
        "Scientific interpretation boundary:",
        "  - Intake is a shared cage outcome. Daily rows and the two mice in a cage are not independent replicates.",
        "  - Water has only two cages; Allulose and Sucrose have three each. ANOVA assumptions are weakly assessable and inference is exploratory.",
        "  - Body-weight inference uses one mean final-minus-baseline change per cage; mouse trajectories are descriptive.",
        "  - Hoja 4 does not contain usable sex/genotype covariates, so adjusted causal claims are not supported.",
        "  - Hoja 3 refill arithmetic is ambiguous for several solution-bottle dates; see validation and sensitivity outputs.",
        "  - Water-control 'solution fraction' is the second-water-bottle share, not sweetener preference.",
        "",
        "Legacy analysis audit:",
    ])
    for row in audit.itertuples(index=False):
        report_lines.append(f"  {row.issue_id} [{row.severity}]: {row.impact}")
    (args.output_dir / "ANALYSIS_REPORT.txt").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    metadata = {
        "analysis": "Experiment 2 behavior preference",
        "created_by": str(Path(__file__).resolve()),
        "workbook": str(args.workbook),
        "workbook_sha256": sha256(args.workbook),
        "sheet": args.sheet,
        "validation_sheet": args.validation_sheet,
        "primary_end_date": primary_end.isoformat(),
        "primary_window_provenance": "Commands_2026.txt explicitly excluded 2025-11-01 and 2025-11-02; rationale not recorded",
        "cage_counts": cages.groupby("treatment")["cage"].nunique().to_dict(),
        "mouse_counts": mouse_summary.groupby("treatment")["mouse_id"].nunique().to_dict(),
        "refill_robust_excluded_dates": refill_dates,
        "statistics": (
            "ordinary one-way ANOVA globally; raw exact two-sided Mann-Whitney U "
            "for W-S, W-A, and S-A; one summary per cage"
        ),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "openpyxl": openpyxl.__version__,
    }
    (args.output_dir / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    legacy_outputs = args.legacy_root / "UCSC_MICE_sheet4_analysis"
    manifest_files = [
        ("current_workbook", args.workbook),
        ("analysis_script", Path(__file__).resolve()),
        ("legacy_script", args.legacy_root / "digest_ucsc_mice.py"),
        ("legacy_results", legacy_outputs / "UCSC_MICE_Hoja_4_v9_analysis_results.txt"),
        ("legacy_cage_intakes", legacy_outputs / "UCSC_MICE_Hoja_4_v9_cage_intakes.csv"),
        ("legacy_anova", legacy_outputs / "UCSC_MICE_Hoja_4_v9_anova_pvalues.csv"),
        ("legacy_mixed_lrt", legacy_outputs / "UCSC_MICE_Hoja_4_v9_mixed_effects_lrt_pvalues.csv"),
    ]
    write_manifest(args.output_dir, manifest_files)
    required_outputs = (
        "daily_endpoints_long.csv",
        "cage_endpoint_summaries.csv",
        "one_way_anova_mann_whitney_statistics.csv",
        "mouse_weight_long.csv",
        "mouse_weight_summary_primary.csv",
        "weight_cage_cluster_summary.csv",
        "ANALYSIS_REPORT.txt",
        "analysis_metadata.json",
        "source_manifest.csv",
    )
    missing_outputs = [name for name in required_outputs if not (args.output_dir / name).is_file()]
    if missing_outputs:
        raise RuntimeError(f"Incomplete staged analysis: {missing_outputs}")
    publish_staging_output(args.output_dir, final_output)
    print(f"Wrote Experiment 2 analysis to {final_output}")
    print(f"Primary design: {len(cages)} cages, {len(mouse_summary)} mice, 7 post-initialization days")
    print(f"Hoja 3/4 discrepancies flagged: {len(discrepancy)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
