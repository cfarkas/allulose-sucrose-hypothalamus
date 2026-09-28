#!/usr/bin/env python3
"""Audit and re-analyse the pooled four-region HIL-annotated global c-FOS data.

The upstream combined file accidentally includes a pseudoanimal whose literal
animal ID is ``horizontal``.  This script never edits the upstream results.  It
removes only that exact ID, retaining the real horizontal-section animals from
the separately quantified run (for example E7_FR7-5H and E9_FR8-1M).

Primary endpoint and inference
------------------------------
Animal-level c-FOS+/DAPI ratio, recomputed after collapsing the raw source rows
to one animal/condition/region value, for ARC, ME, VMN and OTHERS. All finite
values are retained. ARC and ME are the prespecified anatomical endpoints; VMN
is exploratory anatomy and OTHERS is unpainted residual tissue/QC, not an ROI.

The inferential machinery is intentionally identical to Figure 4 panels D-E:
within each region, the global Water/Sucrose/Allulose comparison uses ordinary
equal-variance one-way ANOVA. Pairwise comparisons use a two-sided
rank-sum/Mann-Whitney U statistic with
exhaustive independent-label permutation. The three pairwise p values are
reported unadjusted, as requested; global and sensitivity four-region Holm
values are reported separately and never substituted for the visible raw p.
Effect sizes are Cliff's delta in the displayed second-minus-first direction.
Group summaries remain arithmetic mean +/- sample SD, matching Figure 4.

The area-normalised density is emitted only as a descriptive sensitivity/QC
endpoint.  It is not tested because annotated area and acquisition batch vary
substantially and are partly condition-linked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import f_oneway, rankdata


HERE = Path(__file__).resolve()
PAPER = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
ROOT = PAPER.parent
FIGURE_ANALYSIS_ROOT = PAPER / "analyses" / "Fig2"
DEFAULT_INPUT = FIGURE_ANALYSIS_ROOT / "raw" / "combined_raw_per_animal_human_region_summary.csv"
SHIPPED_INPUT = (
    PAPER / "Fig2" / "raw_data" / "Apotome" / "combined_human_region_plots"
    / "combined_raw_per_animal_human_region_summary.csv"
)
DEFAULT_OUT = FIGURE_ANALYSIS_ROOT / "results"
CONDITIONS = ("Water", "Sucrose", "Allulose")
REGIONS = ("ARC", "ME", "VMN", "OTHERS")
ACQUISITION_FAMILIES = ("June", "August", "September")
REGION_ROLES = {
    "ARC": "primary_anatomical_endpoint",
    "ME": "primary_anatomical_endpoint",
    "VMN": "exploratory_anatomical_endpoint",
    "OTHERS": "unpainted_residual_tissue_qc_not_anatomical_roi",
}
EXPECTED_N = {
    "ARC": {"Water": 4, "Sucrose": 6, "Allulose": 6},
    "ME": {"Water": 4, "Sucrose": 4, "Allulose": 6},
    "VMN": {"Water": 4, "Sucrose": 6, "Allulose": 6},
    "OTHERS": {"Water": 4, "Sucrose": 6, "Allulose": 6},
}
EXPECTED_GLOBAL_P = {
    "ARC": 0.13263744932024776,
    "ME": 0.46649542038267056,
    "VMN": 0.5208878895377919,
    "OTHERS": 0.16343912014165743,
}
EXPECTED_PAIRWISE = {
    "ARC": {
        "Water|Sucrose": (0.6095238095238096, -0.25),
        "Water|Allulose": (0.6095238095238096, 0.25),
        "Sucrose|Allulose": (0.04112554112554113, 0.7222222222222222),
    },
    "ME": {
        "Water|Sucrose": (0.34285714285714286, -0.5),
        "Water|Allulose": (0.9142857142857143, 0.08333333333333333),
        "Sucrose|Allulose": (0.47619047619047616, 0.3333333333333333),
    },
    "VMN": {
        "Water|Sucrose": (0.7619047619047619, -0.16666666666666666),
        "Water|Allulose": (0.7619047619047619, -0.16666666666666666),
        "Sucrose|Allulose": (0.5887445887445888, 0.2222222222222222),
    },
    "OTHERS": {
        "Water|Sucrose": (0.3523809523809524, -0.4166666666666667),
        "Water|Allulose": (0.7619047619047619, 0.16666666666666666),
        "Sucrose|Allulose": (0.09307359307359307, 0.6111111111111112),
    },
}
EXPECTED_IQR_FLAGS = {
    ("VMN", "Sucrose", "E3_FR2-1"),
    ("VMN", "Allulose", "E8_FR6-4H"),
}
CANONICAL_INPUT_SHA256 = "d80bfe4dd5b259f1d309dcf59f54cb0effd62d7911f72de769cfc2b9b8170ebf"
CONDITION_OVERRIDES = {"E3_FR2-1": "Sucrose"}
DATASETS_ROOT = FIGURE_ANALYSIS_ROOT / "raw" / "datasets"
HIL_RUNS_ROOT = FIGURE_ANALYSIS_ROOT / "hil" / "runs"
# The same June 2025 evidence files, as this tree actually publishes them. The
# analyses/Fig2 locations are a private working copy that a downloaded tree does
# not contain, and the condition-provenance audit reads these files
# unconditionally, so naming only the working copy made Figure 2 unbuildable
# anywhere but the machine that first produced it. raw_data/SOURCE_ROOT_MAPPING
# .csv records the two roots as the same rsync-verified source.
FIG2_RAW_DATA_ROOT = PAPER / "Fig2" / "raw_data" / "Apotome"
SHIPPED_DATASETS_ROOT = FIG2_RAW_DATA_ROOT / "24_06_2025" / "analysis"
SHIPPED_HIL_RUN_ROOT = SHIPPED_DATASETS_ROOT / "human_regions_v9"


def _first_existing(*candidates: Path) -> Path:
    """First candidate that exists, else the first, so defaults name real files."""
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


CONDITION_EVIDENCE = {
    "behavior_weights": _first_existing(
        PAPER / "analyses" / "Fig1" / "raw" / "weights.csv",
        PAPER / "Fig1" / "raw_data" / "weights.csv",
    ),
    "binning_metadata": _first_existing(
        DATASETS_ROOT / "june_2025_primary" / "bins_heuristic_stats_all.csv",
        SHIPPED_DATASETS_ROOT / "bins_heuristic_stats_all.csv",
    ),
    "assay_samplesheet": _first_existing(
        DATASETS_ROOT / "june_2025_primary" / "apotome_samplesheet_template.csv",
        SHIPPED_DATASETS_ROOT / "apotome_samplesheet_template.csv",
    ),
    "image_discovery": _first_existing(
        HIL_RUNS_ROOT / "june_2025_primary" / "discovered_images.csv",
        SHIPPED_HIL_RUN_ROOT / "discovered_images.csv",
    ),
    "annotation_status": _first_existing(
        HIL_RUNS_ROOT / "june_2025_primary" / "annotation_status.csv",
        SHIPPED_HIL_RUN_ROOT / "annotation_status.csv",
    ),
}
def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def holm_adjust(p_values: list[float]) -> list[float]:
    """Figure 4's Holm step-down family-wise-error adjustment."""

    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(p))
    if finite.size:
        order = finite[np.argsort(p[finite], kind="mergesort")]
        running = 0.0
        for rank, index in enumerate(order):
            running = max(running, min(1.0, (order.size - rank) * p[index]))
            adjusted[index] = running
    return adjusted.tolist()


def acquisition_family(source_dataset: object, source_csv: object) -> tuple[str, str]:
    """Map audited source paths to the documented acquisition month.

    The mapping deliberately uses acquisition-directory evidence, never animal
    names or outcomes. Multiple source rows for one animal are permitted only
    when every source belongs to the same acquisition family.
    """

    evidence = f"{source_dataset};{source_csv}".casefold()
    matches: list[tuple[str, str]] = []
    if "29.09.2025" in evidence or "september_2025_water" in evidence:
        matches.append(("September", "29.09.2025 or september_2025_water source"))
    if (
        "aug_2025_npy" in evidence
        or "selected_human_regions" in evidence
        or "selected_horizontal_human_regions" in evidence
        or "august_2025_" in evidence
    ):
        matches.append(("August", "Aug_2025_NPY or selected_human_regions source"))
    if (
        "24_06_2025" in evidence
        or "2025.06.24" in evidence
        or "2025.06.24_human_regions" in evidence
        or "june_2025_" in evidence
    ):
        matches.append(("June", "24_06_2025 or 2025.06.24 source directory/dataset"))
    matches = list(dict.fromkeys(matches))
    if len(matches) != 1:
        raise ValueError(
            "Each animal must map to exactly one documented acquisition family; "
            f"found {matches} for source evidence {evidence!r}"
        )
    return matches[0]


def canonical_condition(value: object) -> str:
    key = str(value).strip().casefold()
    aliases = {
        "water": "Water",
        "agua": "Water",
        "sucrose": "Sucrose",
        "sacarosa": "Sucrose",
        "allulose": "Allulose",
        "alulosa": "Allulose",
    }
    return aliases.get(key, str(value).strip())


def audit_and_prepare(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = {
        "animal",
        "cond",
        "region",
        "dapi_nuclei",
        "cfos_cells",
        "area_um2",
        "source_dataset",
        "source_csv",
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    audit = raw.copy()
    audit.insert(0, "raw_row_id", np.arange(len(audit), dtype=int))
    audit["animal"] = audit["animal"].astype(str).str.strip()
    audit["cond_original"] = audit["cond"].astype(str).str.strip()
    audit["cond_canonical_before_override"] = audit["cond"].map(canonical_condition)
    audit["condition_override"] = audit["animal"].map(CONDITION_OVERRIDES)
    audit["condition_override_applied"] = audit["condition_override"].notna()
    audit["cond_canonical"] = audit["condition_override"].fillna(
        audit["cond_canonical_before_override"]
    )
    audit["region_canonical"] = audit["region"].astype(str).str.strip().str.upper()
    for col in ("dapi_nuclei", "cfos_cells", "area_um2", "area_px", "n_images", "cfos_over_dapi"):
        if col in audit.columns:
            audit[col] = pd.to_numeric(audit[col], errors="coerce")

    audit["is_exact_pseudoanimal_horizontal"] = audit["animal"].str.casefold().eq("horizontal")
    audit["condition_in_scope"] = audit["cond_canonical"].isin(CONDITIONS)
    audit["region_in_scope"] = audit["region_canonical"].isin(REGIONS)
    audit["included_for_animal_region_collapse"] = (
        ~audit["is_exact_pseudoanimal_horizontal"]
        & audit["condition_in_scope"]
        & audit["region_in_scope"]
    )

    def raw_reason(row: pd.Series) -> str:
        reasons: list[str] = []
        if row["is_exact_pseudoanimal_horizontal"]:
            reasons.append("exclude_exact_pseudoanimal_horizontal")
        if not row["condition_in_scope"]:
            reasons.append("exclude_condition_out_of_scope")
        if not row["region_in_scope"]:
            reasons.append("exclude_region_out_of_scope")
        return ";".join(reasons) if reasons else "include_for_animal_condition_region_collapse"

    audit["decision_reason"] = audit.apply(raw_reason, axis=1)
    audit["raw_ratio_recomputed"] = np.where(
        audit["dapi_nuclei"] > 0,
        audit["cfos_cells"] / audit["dapi_nuclei"],
        np.nan,
    )
    audit["raw_ratio_abs_difference"] = (
        audit["cfos_over_dapi"] - audit["raw_ratio_recomputed"]
        if "cfos_over_dapi" in audit.columns
        else np.nan
    )
    audit["raw_ratio_abs_difference"] = pd.to_numeric(
        audit["raw_ratio_abs_difference"], errors="coerce"
    ).abs()

    collapse_input = audit.loc[audit["included_for_animal_region_collapse"]].copy()

    def sum_finite(series: pd.Series) -> float:
        return float(pd.to_numeric(series, errors="coerce").sum(min_count=1))

    def join_unique(series: pd.Series) -> str:
        return ";".join(sorted({str(value) for value in series.dropna() if str(value)}))

    aggregations: dict[str, tuple[str, object]] = {
        "dapi_nuclei": ("dapi_nuclei", sum_finite),
        "cfos_cells": ("cfos_cells", sum_finite),
        "area_um2": ("area_um2", sum_finite),
        "n_images": ("n_images", sum_finite),
        "sample": ("sample", join_unique),
        "source_dataset": ("source_dataset", join_unique),
        "source_csv": ("source_csv", join_unique),
        "source_row_count": ("raw_row_id", "size"),
        "max_raw_reported_ratio_abs_difference": ("raw_ratio_abs_difference", "max"),
        "condition_override_applied": ("condition_override_applied", "any"),
        "source_condition_labels": ("cond_original", join_unique),
    }
    if "area_px" in collapse_input.columns:
        aggregations["area_px"] = ("area_px", sum_finite)
    collapsed = (
        collapse_input.groupby(
            ["animal", "cond_canonical", "region_canonical"],
            sort=True,
            observed=True,
            as_index=False,
        )
        .agg(**aggregations)
    )
    collapsed["valid_counts_after_collapse"] = (
        np.isfinite(collapsed["dapi_nuclei"])
        & np.isfinite(collapsed["cfos_cells"])
        & (collapsed["dapi_nuclei"] > 0)
        & (collapsed["cfos_cells"] >= 0)
    )
    collapsed["collapsed_decision_reason"] = np.where(
        collapsed["valid_counts_after_collapse"],
        "included_all_finite_animal_region_value",
        "exclude_missing_or_nonpositive_DAPI_denominator_after_collapse",
    )

    values = collapsed.loc[collapsed["valid_counts_after_collapse"]].copy()
    if values.duplicated(["animal", "cond_canonical", "region_canonical"]).any():
        raise RuntimeError("Animal/condition/region collapse did not produce unique rows")
    values["condition"] = pd.Categorical(values["cond_canonical"], CONDITIONS, ordered=True)
    values["region"] = pd.Categorical(values["region_canonical"], REGIONS, ordered=True)
    values["region_role"] = values["region_canonical"].map(REGION_ROLES)
    values["cfos_over_dapi_recomputed"] = values["cfos_cells"] / values["dapi_nuclei"]
    values["cfos_percent_dapi"] = 100.0 * values["cfos_over_dapi_recomputed"]
    values["cfos_density_per_100k_um2"] = np.where(
        values["area_um2"] > 0,
        100_000.0 * values["cfos_cells"] / values["area_um2"],
        np.nan,
    )
    values["upstream_ratio_abs_difference"] = values["max_raw_reported_ratio_abs_difference"]

    keep = [
        "animal", "condition", "region", "region_role", "dapi_nuclei", "cfos_cells",
        "area_um2", "area_px", "cfos_over_dapi_recomputed", "cfos_percent_dapi",
        "cfos_density_per_100k_um2", "upstream_ratio_abs_difference", "n_images",
        "source_row_count", "sample", "source_dataset", "source_csv",
        "condition_override_applied", "source_condition_labels",
    ]
    keep = [col for col in keep if col in values.columns]
    values = values[keep].sort_values(["region", "condition", "animal"]).reset_index(drop=True)

    excluded_raw = audit.loc[~audit["included_for_animal_region_collapse"]].copy()
    excluded_raw.insert(0, "audit_level", "raw_source_row")
    excluded_collapsed = collapsed.loc[~collapsed["valid_counts_after_collapse"]].copy()
    excluded_collapsed.insert(0, "audit_level", "collapsed_animal_condition_region")
    excluded = pd.concat([excluded_raw, excluded_collapsed], ignore_index=True, sort=False)
    return values, audit, excluded


def descriptives(values: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for region in REGIONS:
        for condition in CONDITIONS:
            x = values.loc[
                (values["region"].astype(str) == region)
                & (values["condition"].astype(str) == condition),
                endpoint,
            ].dropna().to_numpy(float)
            rows.append(
                {
                    "endpoint": endpoint,
                    "region": region,
                    "region_role": REGION_ROLES[region],
                    "condition": condition,
                    "n_animals": len(x),
                    "mean": np.mean(x) if len(x) else np.nan,
                    "sd": np.std(x, ddof=1) if len(x) > 1 else np.nan,
                    "median": np.median(x) if len(x) else np.nan,
                    "q1": np.quantile(x, 0.25) if len(x) else np.nan,
                    "q3": np.quantile(x, 0.75) if len(x) else np.nan,
                    "min": np.min(x) if len(x) else np.nan,
                    "max": np.max(x) if len(x) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def build_outlier_audit(values: pd.DataFrame) -> pd.DataFrame:
    """Flag, but never exclude, values beyond within-group 1.5 x IQR fences.

    The rule is applied independently within each treatment x region cell and
    is descriptive only. A flag is not evidence of a measurement error. The
    validated HIL counts and denominators remain the primary inferential data.
    """

    endpoint = "cfos_over_dapi_recomputed"
    rows: list[dict[str, object]] = []
    for region in REGIONS:
        for condition in CONDITIONS:
            group = values.loc[
                (values["region"].astype(str) == region)
                & (values["condition"].astype(str) == condition)
            ].copy()
            x = group[endpoint].to_numpy(float)
            if len(x) == 0:
                continue
            q1 = float(np.quantile(x, 0.25))
            q3 = float(np.quantile(x, 0.75))
            iqr = q3 - q1
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            for record in group.itertuples(index=False):
                value = float(getattr(record, endpoint))
                flagged = bool(value < lower or value > upper)
                rows.append(
                    {
                        "animal": record.animal,
                        "condition": str(record.condition),
                        "region": str(record.region),
                        "region_role": record.region_role,
                        "endpoint": endpoint,
                        "value": value,
                        "dapi_nuclei": float(record.dapi_nuclei),
                        "cfos_cells": float(record.cfos_cells),
                        "group_n": len(x),
                        "q1": q1,
                        "q3": q3,
                        "iqr": iqr,
                        "lower_fence_1_5_iqr": lower,
                        "upper_fence_1_5_iqr": upper,
                        "descriptive_iqr_flag": flagged,
                        "inferential_inclusion": True,
                        "plot_marker": "x" if flagged else "circle",
                        "decision_reason": (
                            "plausible validated HIL value; retained; no outcome-dependent exclusion"
                            if flagged
                            else "within 1.5 x IQR fences; retained"
                        ),
                        "source_dataset": record.source_dataset,
                        "source_csv": record.source_csv,
                    }
                )
    return pd.DataFrame(rows).sort_values(
        ["region", "condition", "animal"]
    ).reset_index(drop=True)


def cliffs_delta(second: np.ndarray, first: np.ndarray) -> float:
    """Cliff's delta for the displayed direction second minus first."""

    left = np.asarray(second, dtype=float)
    right = np.asarray(first, dtype=float)
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if not left.size or not right.size:
        return float("nan")
    differences = left[:, None] - right[None, :]
    return float(
        (np.count_nonzero(differences > 0) - np.count_nonzero(differences < 0))
        / differences.size
    )


def one_way_anova(groups: list[np.ndarray]) -> tuple[float, float, int, int]:
    """Ordinary equal-variance one-way ANOVA on independent animal values."""
    usable = [np.asarray(group, dtype=float) for group in groups]
    usable = [
        group[np.isfinite(group)]
        for group in usable
        if np.count_nonzero(np.isfinite(group)) >= 2
    ]
    if len(usable) < 2:
        return float("nan"), float("nan"), 0, 0
    result = f_oneway(*usable)
    return (
        float(result.statistic), float(result.pvalue), len(usable) - 1,
        int(sum(group.size for group in usable) - len(usable)),
    )


def exact_permutation_mwu(
    first: np.ndarray, second: np.ndarray, max_partitions: int = 1_000_000,
) -> float:
    """Exact two-sided independent-label permutation p for rank-sum/U."""

    left = np.asarray(first, dtype=float)
    right = np.asarray(second, dtype=float)
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if left.size < 2 or right.size < 2:
        return float("nan")
    pooled = np.concatenate((left, right))
    n_first = int(left.size)
    total_partitions = math.comb(int(pooled.size), n_first)
    if total_partitions > max_partitions:
        raise RuntimeError(
            f"Exact rank-sum enumeration requires {total_partitions} partitions; "
            f"cap is {max_partitions}"
        )
    ranks = rankdata(pooled, method="average")
    expected_rank_sum = n_first * (pooled.size + 1.0) / 2.0
    observed = abs(float(np.sum(ranks[:n_first])) - expected_rank_sum)
    extreme = sum(
        abs(float(np.sum(ranks[list(selection)])) - expected_rank_sum)
        >= observed - 1e-12
        for selection in combinations(range(int(pooled.size)), n_first)
    )
    return float(extreme / total_partitions)


def inferential_statistics(values: pd.DataFrame) -> pd.DataFrame:
    """Animal-level one-way ANOVA and exact-MWU statistics."""

    endpoint = "cfos_over_dapi_recomputed"
    pair_order = (
        ("Water", "Sucrose"),
        ("Water", "Allulose"),
        ("Sucrose", "Allulose"),
    )
    rows: list[dict[str, object]] = []
    for region in REGIONS:
        groups = {
            condition: values.loc[
                (values["region"].astype(str) == region)
                & (values["condition"].astype(str) == condition),
                endpoint,
            ].dropna().to_numpy(float)
            for condition in CONDITIONS
        }
        if any(groups[condition].size < 2 for condition in CONDITIONS):
            raise ValueError(f"{region} requires at least two animals in every group")
        for condition in CONDITIONS:
            group = groups[condition]
            rows.append(
                {
                    "endpoint": endpoint,
                    "region": region,
                    "test": "descriptive",
                    "comparison": condition,
                    "n": int(group.size),
                    "mean": float(np.mean(group)),
                    "sd": float(np.std(group, ddof=1)),
                    "p_value_raw": np.nan,
                    "effect_size": np.nan,
                    "p_value_method": "",
                    "n_left": np.nan,
                    "n_right": np.nan,
                    "effect_size_definition": "",
                    "permutation_exhaustive": np.nan,
                    "permutations_evaluated": np.nan,
                    "permutations_total": np.nan,
                }
            )
        global_statistic, global_p, df_between, df_within = one_way_anova(
            [groups[condition] for condition in CONDITIONS]
        )
        rows.append(
            {
                "endpoint": endpoint,
                "region": region,
                "test": "One_way_ANOVA",
                "comparison": "Water|Sucrose|Allulose",
                "n": int(sum(groups[condition].size for condition in CONDITIONS)),
                "mean": np.nan,
                "sd": np.nan,
                "p_value_raw": global_p,
                "statistic": global_statistic,
                "df_between": df_between,
                "df_within": df_within,
                "effect_size": np.nan,
                "p_value_method": "ordinary equal-variance one-way ANOVA F",
                "n_left": np.nan,
                "n_right": np.nan,
                "effect_size_definition": "",
                "permutation_exhaustive": np.nan,
                "permutations_evaluated": np.nan,
                "permutations_total": np.nan,
            }
        )
        pair_raw = {
            pair: exact_permutation_mwu(groups[pair[0]], groups[pair[1]])
            for pair in pair_order
        }
        for first, second in pair_order:
            rows.append(
                {
                    "endpoint": endpoint,
                    "region": region,
                    "test": "Mann_Whitney_U_two_sided",
                    "comparison": f"{first}|{second}",
                    "n": int(groups[first].size + groups[second].size),
                    "mean": np.nan,
                    "sd": np.nan,
                    "p_value_raw": pair_raw[(first, second)],
                    "effect_size": cliffs_delta(groups[second], groups[first]),
                    "p_value_method": (
                        "exhaustive independent-label permutation of rank-sum/U"
                    ),
                    "n_left": int(groups[first].size),
                    "n_right": int(groups[second].size),
                    "effect_size_definition": (
                        f"Cliff's delta ({second} minus {first})"
                    ),
                    "permutation_exhaustive": True,
                    "permutations_evaluated": int(
                        math.comb(groups[first].size + groups[second].size, groups[first].size)
                    ),
                    "permutations_total": int(
                        math.comb(groups[first].size + groups[second].size, groups[first].size)
                    ),
                }
            )
    columns = [
        "endpoint", "region", "test", "comparison", "n", "mean", "sd",
        "p_value_raw", "statistic", "df_between", "df_within", "effect_size",
        "p_value_method", "n_left", "n_right", "effect_size_definition",
        "permutation_exhaustive", "permutations_evaluated", "permutations_total",
    ]
    frame = pd.DataFrame(rows, columns=columns)
    frame["p_value_holm_across_four_regions"] = np.nan
    global_mask = frame["test"].eq("One_way_ANOVA")
    frame.loc[global_mask, "p_value_holm_across_four_regions"] = holm_adjust(
        frame.loc[global_mask, "p_value_raw"].astype(float).tolist()
    )
    return frame


def build_acquisition_strata(
    values: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Attach the documented June/August/September acquisition family."""

    strata = values.copy()
    assignments = strata.apply(
        lambda row: acquisition_family(row["source_dataset"], row["source_csv"]),
        axis=1,
        result_type="expand",
    )
    assignments.columns = ["acquisition_family", "acquisition_mapping_rule"]
    strata = pd.concat([strata, assignments], axis=1)
    strata["acquisition_family"] = pd.Categorical(
        strata["acquisition_family"], ACQUISITION_FAMILIES, ordered=True,
    )
    per_animal = strata[
        ["animal", "condition", "acquisition_family", "acquisition_mapping_rule"]
    ].drop_duplicates()
    if per_animal["animal"].duplicated().any():
        raise RuntimeError("An animal maps to multiple condition/acquisition assignments")

    mapping_columns = [
        "animal", "condition", "region", "acquisition_family",
        "acquisition_mapping_rule", "dapi_nuclei", "cfos_cells",
        "cfos_over_dapi_recomputed", "source_dataset", "source_csv",
    ]
    mapping = strata[mapping_columns].copy().sort_values(
        ["region", "acquisition_family", "condition", "animal"]
    ).reset_index(drop=True)

    count_rows: list[dict[str, object]] = []
    for region in REGIONS:
        for family in ACQUISITION_FAMILIES:
            group = strata.loc[
                strata["region"].astype(str).eq(region)
                & strata["acquisition_family"].astype(str).eq(family)
            ]
            present = int(group["condition"].astype(str).nunique())
            informative = present >= 2
            for condition in CONDITIONS:
                members = group.loc[group["condition"].astype(str).eq(condition)]
                count_rows.append({
                    "region": region,
                    "acquisition_family": family,
                    "condition": condition,
                    "n_animals": int(len(members)),
                    "n_conditions_present_in_stratum": present,
                    "informative_for_condition_test": informative,
                    "statistic_contribution_if_single_condition": (
                        "not_applicable" if informative else "zero_no_information"
                    ),
                    "animals": ";".join(sorted(members["animal"].astype(str))),
                })
    counts = pd.DataFrame(count_rows)
    return strata, mapping, counts


def multinomial_labelings(counts: tuple[int, ...]) -> list[np.ndarray]:
    """Enumerate unique label vectors with fixed counts in stable order."""

    n_total = int(sum(counts))
    indices = tuple(range(n_total))
    output: list[np.ndarray] = []

    def recurse(
        condition_index: int,
        remaining_indices: tuple[int, ...],
        labels: np.ndarray,
    ) -> None:
        if condition_index == len(counts) - 1:
            labels[list(remaining_indices)] = condition_index
            output.append(labels.copy())
            return
        for selected in combinations(remaining_indices, counts[condition_index]):
            updated = labels.copy()
            updated[list(selected)] = condition_index
            selected_set = set(selected)
            remainder = tuple(i for i in remaining_indices if i not in selected_set)
            recurse(condition_index + 1, remainder, updated)

    recurse(0, indices, np.full(n_total, -1, dtype=np.int8))
    expected = 1
    remaining = n_total
    for count in counts[:-1]:
        expected *= math.comb(remaining, count)
        remaining -= count
    if len(output) != expected or any(np.any(labels < 0) for labels in output):
        raise RuntimeError(
            f"Fixed-count label enumeration failed: {len(output)} != {expected}"
        )
    return output


def _xlog_ratio(observed: float, expected: float) -> float:
    if observed == 0:
        return 0.0
    if expected <= 0:
        return float("inf")
    return float(observed * math.log(observed / expected))


def within_stratum_binomial_deviance(
    cfos: np.ndarray, dapi: np.ndarray, labels: np.ndarray,
) -> float:
    """Condition-rate versus stratum-only binomial deviance for one stratum."""

    present = [label for label in range(len(CONDITIONS)) if np.any(labels == label)]
    if len(present) < 2:
        return 0.0
    cfos = np.asarray(cfos, dtype=float)
    dapi = np.asarray(dapi, dtype=float)
    if (
        np.any(cfos < 0)
        or np.any(dapi <= 0)
        or np.any(cfos > dapi)
        or not np.allclose(cfos, np.round(cfos), rtol=0, atol=1e-9)
        or not np.allclose(dapi, np.round(dapi), rtol=0, atol=1e-9)
    ):
        raise ValueError("Count-aware sensitivity requires valid integer subset counts")
    pooled_rate = float(np.sum(cfos) / np.sum(dapi))
    if pooled_rate <= 0 or pooled_rate >= 1:
        return 0.0
    deviance = 0.0
    for label in present:
        selected = labels == label
        successes = float(np.sum(cfos[selected]))
        trials = float(np.sum(dapi[selected]))
        failures = trials - successes
        deviance += 2.0 * (
            _xlog_ratio(successes, trials * pooled_rate)
            + _xlog_ratio(failures, trials * (1.0 - pooled_rate))
        )
    return float(max(0.0, deviance))


def acquisition_stratified_sensitivities(
    strata: pd.DataFrame,
) -> pd.DataFrame:
    """Exhaustive whole-animal permutations constrained within acquisition family.

    The condition-rate versus stratum-only binomial deviance uses the actual
    c-FOS/DAPI counts. Its p values arise only from whole-animal label
    permutation; no cell/nucleus is treated as an independent experimental unit.
    """

    rows: list[dict[str, object]] = []
    for region in REGIONS:
        region_values = strata.loc[strata["region"].astype(str).eq(region)].copy()
        stratum_records: list[dict[str, object]] = []
        for family in ACQUISITION_FAMILIES:
            group = region_values.loc[
                region_values["acquisition_family"].astype(str).eq(family)
            ].sort_values("animal").reset_index(drop=True)
            if group.empty:
                continue
            counts = tuple(
                int(group["condition"].astype(str).eq(condition).sum())
                for condition in CONDITIONS
            )
            labelings = multinomial_labelings(counts)
            observed_labels = np.asarray(
                [CONDITIONS.index(str(value)) for value in group["condition"]],
                dtype=np.int8,
            )
            cfos = group["cfos_cells"].to_numpy(float)
            dapi = group["dapi_nuclei"].to_numpy(float)
            deviance_statistics = np.asarray(
                [
                    within_stratum_binomial_deviance(cfos, dapi, labels)
                    for labels in labelings
                ],
                dtype=float,
            )
            stratum_records.append({
                "family": family,
                "counts": counts,
                "labelings": len(labelings),
                "informative": sum(count > 0 for count in counts) >= 2,
                "deviance_observed": within_stratum_binomial_deviance(
                    cfos, dapi, observed_labels,
                ),
                "deviance_statistics": deviance_statistics,
            })
        if not stratum_records:
            raise RuntimeError(f"No acquisition strata available for {region}")

        total = int(math.prod(record["labelings"] for record in stratum_records))
        labeling_receipt = ";".join(
            f"{record['family']}:{record['labelings']}"
            for record in stratum_records
        )
        informative = ";".join(
            record["family"] for record in stratum_records if record["informative"]
        )
        no_information = ";".join(
            record["family"] for record in stratum_records if not record["informative"]
        )
        count_receipt = ";".join(
            f"{record['family']}:" + "/".join(
                f"{CONDITIONS[index][0]}={count}"
                for index, count in enumerate(record["counts"])
            )
            for record in stratum_records
        )
        for family_key, test, observed_key, statistics_key, formula in (
            (
                "acquisition_stratified_count",
                "sum_within_acquisition_binomial_deviance",
                "deviance_observed",
                "deviance_statistics",
                "T_count = sum_s 2*sum_g[y_g log(y_g/(n_g*p_s)) + (n_g-y_g) log((n_g-y_g)/(n_g*(1-p_s)))]; p_s is the pooled c-FOS/DAPI rate in stratum s",
            ),
        ):
            observed = float(sum(record[observed_key] for record in stratum_records))
            combined = np.asarray([0.0], dtype=float)
            for record in stratum_records:
                combined = np.add.outer(combined, record[statistics_key]).reshape(-1)
            if len(combined) != total:
                raise RuntimeError(
                    f"{region} {family_key}: evaluated {len(combined)} != {total}"
                )
            extreme = int(np.count_nonzero(combined >= observed - 1e-12))
            rows.append({
                "endpoint": "cfos_cells / dapi_nuclei",
                "region": region,
                "sensitivity_family": family_key,
                "test": test,
                "statistic": observed,
                "extreme_labelings": extreme,
                "p_value_raw": float(extreme / total),
                "p_value_holm_across_four_regions": np.nan,
                "permutation_unit": "whole animal (counts and denominator move together)",
                "stratification": "documented acquisition family: June/August/September",
                "stratum_condition_counts": count_receipt,
                "strata_informative": informative,
                "strata_no_information": no_information,
                "single_condition_stratum_policy": "recorded; zero statistic; one fixed labeling; no information",
                "labelings_by_stratum": labeling_receipt,
                "permutation_exhaustive": True,
                "permutations_evaluated": total,
                "permutations_total": total,
                "p_value_method": "exhaustive whole-animal label permutation within acquisition strata",
                "formula": formula,
                "cell_independence_claim": False,
            })
    result = pd.DataFrame(rows)
    for family_key in result["sensitivity_family"].unique():
        selected = result["sensitivity_family"].eq(family_key)
        result.loc[selected, "p_value_holm_across_four_regions"] = holm_adjust(
            result.loc[selected, "p_value_raw"].astype(float).tolist()
        )
    if not (
        result["permutation_exhaustive"].astype(bool).all()
        and (
            result["permutations_evaluated"].astype(int)
            == result["permutations_total"].astype(int)
        ).all()
    ):
        raise RuntimeError("A stratified sensitivity enumeration is incomplete")
    return result.sort_values(["sensitivity_family", "region"]).reset_index(drop=True)


def validate_paper_results(
    values: pd.DataFrame,
    stats: pd.DataFrame,
    outlier_audit: pd.DataFrame,
) -> None:
    """Fail closed if the fixed paper dataset no longer reproduces its audit."""
    counts = (
        values.groupby(["region", "condition"], observed=True)["animal"]
        .nunique().to_dict()
    )
    for region, expected_by_condition in EXPECTED_N.items():
        for condition, expected in expected_by_condition.items():
            observed = int(counts.get((region, condition), 0))
            if observed != expected:
                raise RuntimeError(
                    f"Unexpected {region}/{condition} n: {observed}; expected {expected}"
                )
    global_rows = stats.loc[stats["test"].eq("One_way_ANOVA")]
    pair_rows = stats.loc[stats["test"].eq("Mann_Whitney_U_two_sided")]
    descriptive_rows = stats.loc[stats["test"].eq("descriptive")]
    if len(global_rows) != 4 or len(pair_rows) != 12 or len(descriptive_rows) != 12:
        raise RuntimeError("Expected 4 global, 12 pairwise and 12 descriptive rows")
    for region, expected in EXPECTED_GLOBAL_P.items():
        row = global_rows.loc[global_rows["region"].eq(region)]
        if len(row) != 1:
            raise RuntimeError(f"Expected exactly one {region} one-way ANOVA result")
        result = row.iloc[0]
        if not np.isclose(float(result["p_value_raw"]), expected, rtol=1e-12, atol=1e-14):
            raise RuntimeError(f"{region} canonical one-way ANOVA p changed")
    expected_across_regions = holm_adjust(
        global_rows["p_value_raw"].astype(float).tolist()
    )
    if not np.allclose(
        global_rows["p_value_holm_across_four_regions"].astype(float),
        expected_across_regions, rtol=1e-12, atol=1e-14,
    ):
        raise RuntimeError("Figure-wide Holm adjustment of the four global ANOVA p values changed")
    if not stats.loc[
        stats["test"].eq("Mann_Whitney_U_two_sided"),
        "permutation_exhaustive",
    ].astype(bool).all():
        raise RuntimeError("Every inferential row must prove exhaustive enumeration")
    inferential = stats["test"].eq("Mann_Whitney_U_two_sided")
    if not (
        stats.loc[inferential, "permutations_evaluated"].astype(int)
        == stats.loc[inferential, "permutations_total"].astype(int)
    ).all():
        raise RuntimeError("Not every planned permutation was evaluated")
    expected_pairs = {
        "Water|Sucrose", "Water|Allulose", "Sucrose|Allulose",
    }
    if any(
        set(pair_rows.loc[pair_rows["region"].eq(region), "comparison"])
        != expected_pairs
        for region in REGIONS
    ):
        raise RuntimeError("Each region must contain the three required pairwise contrasts")
    for region in REGIONS:
        region_pairs = pair_rows.loc[pair_rows["region"].eq(region)]
        for comparison, (expected_raw, expected_delta) in (
            EXPECTED_PAIRWISE[region].items()
        ):
            row = region_pairs.loc[region_pairs["comparison"].eq(comparison)]
            if len(row) != 1:
                raise RuntimeError(f"Expected one {region} {comparison} pairwise row")
            observed = row.iloc[0]
            for column, expected in (
                ("p_value_raw", expected_raw),
                ("effect_size", expected_delta),
            ):
                if not np.isclose(
                    float(observed[column]), expected, rtol=1e-12, atol=1e-14,
                ):
                    raise RuntimeError(
                        f"{region} {comparison} canonical {column} changed"
                    )
    flagged = {
        (str(row.region), str(row.condition), str(row.animal))
        for row in outlier_audit.loc[outlier_audit["descriptive_iqr_flag"]].itertuples(index=False)
    }
    if flagged != EXPECTED_IQR_FLAGS:
        raise RuntimeError(
            f"Canonical descriptive 1.5 x IQR flags changed: {sorted(flagged)}; "
            f"expected {sorted(EXPECTED_IQR_FLAGS)}"
        )
    if not outlier_audit["inferential_inclusion"].all():
        raise RuntimeError("The canonical primary analysis must retain every IQR-flagged value")


def build_condition_provenance(values: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    unique = values[
        [
            "animal", "condition", "source_condition_labels", "condition_override_applied",
        ]
    ].drop_duplicates("animal").sort_values("animal")
    weights = CONDITION_EVIDENCE["behavior_weights"]
    bins = CONDITION_EVIDENCE["binning_metadata"]
    for row in unique.itertuples(index=False):
        if row.animal == "E3_FR2-1":
            status = "audited_override_resolved"
            basis = (
                "Sucrose in behavior weights rows 32-36 and assay-local binning rows 74-85 "
                "(sucrose2-rep1); replaces upstream REVIEW"
            )
        else:
            status = "upstream_label_retained"
            basis = "Upstream HIL condition label retained"
        rows.append({
            "animal": row.animal,
            "resolved_condition": str(row.condition),
            "upstream_condition_labels": row.source_condition_labels,
            "condition_override_applied": bool(row.condition_override_applied),
            "resolution_status": status,
            "resolution_basis": basis,
            "behavior_evidence_path": str(weights) if row.animal == "E3_FR2-1" else "",
            "behavior_evidence_sha256": sha256(weights) if row.animal == "E3_FR2-1" else "",
            "binning_evidence_path": str(bins) if row.animal == "E3_FR2-1" else "",
            "binning_evidence_sha256": sha256(bins) if row.animal == "E3_FR2-1" else "",
            "all_condition_evidence_paths": "",
            "all_condition_evidence_sha256": "",
        })
    return pd.DataFrame(rows)


def write_audit(
    outdir: Path,
    source: Path,
    values: pd.DataFrame,
    audit: pd.DataFrame,
    stats: pd.DataFrame,
    outlier_audit: pd.DataFrame,
    sensitivities: pd.DataFrame,
) -> None:
    n_by_group = (
        values.groupby(["region", "condition"], observed=True)["animal"]
        .nunique()
        .rename("n")
        .reset_index()
    )
    max_discrepancy = float(values["upstream_ratio_abs_difference"].max())
    pseudo = audit.loc[audit["is_exact_pseudoanimal_horizontal"], ["animal", "cond", "region", "source_dataset"]]
    true_horizontal = values.loc[
        values["source_dataset"].astype(str).str.contains("selected_horizontal_human_regions", regex=False),
        ["animal", "condition", "region", "source_dataset"],
    ]
    global_summary = stats.loc[
        stats["test"].eq("One_way_ANOVA"),
        [
            "region", "n", "p_value_raw", "p_value_holm_across_four_regions",
            "statistic", "df_between", "df_within",
        ],
    ]
    pairwise_summary = stats.loc[
        stats["test"].eq("Mann_Whitney_U_two_sided"),
        [
            "region", "comparison", "n_left", "n_right", "p_value_raw",
            "effect_size", "effect_size_definition", "permutations_evaluated",
        ],
    ]
    flagged = outlier_audit.loc[outlier_audit["descriptive_iqr_flag"]]
    lines = [
        "# Global c-FOS statistical and provenance audit",
        "",
        "## Decision",
        "",
        "The endpoint is the animal-level c-FOS+/DAPI ratio, recomputed after raw",
        "source rows are summed to one animal/condition/region value. All finite",
        "values from eligible, unambiguously labelled animals are retained; there",
        "is no IQR or other outcome-dependent exclusion.",
        "ARC and ME are the primary anatomical regions, VMN is exploratory, and",
        "OTHERS is unpainted residual tissue/QC rather than an anatomical ROI.",
        "",
        "The row whose exact animal ID is `horizontal` is a directory-derived",
        "pseudoanimal and is excluded. Real mice quantified from the dedicated",
        "horizontal-section run remain included (notably E7_FR7-5H and E9_FR8-1M).",
        "E3_FR2-1 is explicitly overridden from upstream REVIEW to Sucrose because",
        "both behavior records and assay-local bin metadata independently label it",
        "Sucrose. Exact evidence paths and hashes are in condition_provenance.csv.",
        "## Four-region sample sizes",
        "",
        "```",
        n_by_group.to_string(index=False),
        "```",
        "",
        "ME is absent for both E3_FR2-1 and E4_FR2-5 because their annotated ME",
        "denominators are 0/0. These two structural missing values are not imputed;",
        "therefore Sucrose ME n=4 rather than the six otherwise eligible animals.",
        "",
        "## Primary statistics",
        "",
        "Within each region, the global comparison uses ordinary equal-variance",
        "one-way ANOVA with animal as the experimental unit. These four",
        "global p values shown in the plots are raw, matching Figure 4. A figure-wide",
        "Holm adjustment across the four global regional p values is additionally",
        "reported in the CSV, audit and legend.",
        "Each of the three pairwise comparisons exhaustively permutes labels for the",
        "two-sided rank-sum/Mann-Whitney U statistic. Water-Sucrose, Water-Allulose",
        "and Sucrose-Allulose pairwise p values are reported unadjusted, as requested.",
        "Effect size is Cliff's delta for the second condition minus the first. Group",
        "summaries use arithmetic mean +/- sample SD.",
        "",
        "Global one-way ANOVA summary:",
        "",
        "```",
        global_summary.to_string(index=False),
        "```",
        "",
        "Pairwise exhaustive rank-sum/Mann-Whitney U summary:",
        "",
        "```",
        pairwise_summary.to_string(index=False),
        "```",
        "",
        "## Acquisition-stratified sensitivity analyses",
        "",
        "The count sensitivity exhaustively permutes whole-animal condition labels",
        "within the documented June, August and September acquisition families while",
        "preserving the observed condition counts in each family. September is Water-only,",
        "so it is recorded with one fixed labeling and contributes zero information.",
        "The count statistic is the sum of within-acquisition binomial deviances comparing",
        "condition-specific pooled c-FOS/DAPI rates with the acquisition-pooled rate.",
        "Counts and denominators move together with each animal; nuclei are not treated as",
        "independent experimental units. Raw p values and Holm adjustment across the four",
        "regions are reported for the count sensitivity family.",
        "",
        "```",
        sensitivities[[
            "region", "sensitivity_family", "statistic", "p_value_raw",
            "p_value_holm_across_four_regions", "strata_informative",
            "strata_no_information", "permutations_evaluated",
            "permutations_total",
        ]].to_string(index=False),
        "```",
        "",
        "## Descriptive outlier audit",
        "",
        "A 1.5 x IQR rule is evaluated independently within every treatment x region",
        "cell for transparent display/QC only. It flags two VMN observations:",
        "E3_FR2-1 (Sucrose; 332/4542) and E8_FR6-4H (Allulose; 189/1824).",
        "Both have substantial nonzero DAPI denominators, accepted HIL provenance and",
        "no identified measurement/QC error. They are plausible biological values and",
        "remain in the ANOVA/MWU inference; removing them because of their outcome would be",
        "an outcome-dependent exclusion. The figure should show them with x markers",
        "while retaining them in the mean, SD, one-way ANOVA and pairwise calculations.",
        "",
        "```",
        flagged[
            [
                "animal", "condition", "region", "value", "dapi_nuclei", "cfos_cells",
                "lower_fence_1_5_iqr", "upper_fence_1_5_iqr", "inferential_inclusion",
                "decision_reason",
            ]
        ].to_string(index=False),
        "```",
        "",
        "## Integrity checks",
        "",
        f"Maximum absolute difference between each raw row's upstream ratio and its counts-recomputed ratio: {max_discrepancy:.3g}.",
        f"Exact pseudoanimal rows identified upstream: {len(pseudo)}.",
        f"Primary rows retained from the separately quantified horizontal-section dataset: {len(true_horizontal)}.",
        "The Water E10_FR1-1 and E11_FR3-3 rows span two source analyses and are",
        "summed at the animal/condition/region level before ratio calculation.",
        "Sampling depth is unequal: Water E10_FR1-1 represents four images and",
        "E11_FR3-3 two images, whereas most other animal-region values represent",
        "one image. Summed counts weight every animal once statistically but yield",
        "different within-animal measurement depth.",
        "",
        "## Important limitation",
        "",
        "Animals were pooled across multiple staining/acquisition/annotation datasets.",
        "Dataset and condition are not fully crossed, image pixel sizes differ, and",
        "sex/genotype metadata are not available in this combined table. Consequently,",
        "treatment effects cannot be cleanly separated from cohort/batch effects. These",
        "tests are exploratory pooled comparisons, not a confirmatory randomized",
        "single-batch analysis. The one-way ANOVA does not address this",
        "batch/condition confounding. c-FOS density per annotated area is written only as a",
        "descriptive sensitivity/QC endpoint because sampled polygon area and acquisition",
        "batch vary; it is not used for inference.",
        "OTHERS includes all unpainted residual pixels and sometimes very small DAPI",
        "denominators (as low as 12 nuclei); it must not be interpreted as anatomy.",
        "",
        "## Source",
        "",
        f"- `{source}`",
        f"- SHA-256: `{sha256(source)}`",
    ]
    (outdir / "STATISTICAL_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_sensitivity_results(
    mapping: pd.DataFrame,
    counts: pd.DataFrame,
    sensitivities: pd.DataFrame,
) -> None:
    """Fail closed on the fixed acquisition structure and enumeration receipts."""

    if len(sensitivities) != 4 or set(sensitivities["region"]) != set(REGIONS):
        raise RuntimeError("Expected one acquisition-stratified count sensitivity for four regions")
    expected_families = {"acquisition_stratified_count"}
    if set(sensitivities["sensitivity_family"]) != expected_families:
        raise RuntimeError("Unexpected Figure 2 sensitivity family")
    if set(mapping["acquisition_family"].astype(str)) != set(ACQUISITION_FAMILIES):
        raise RuntimeError("Acquisition mapping does not cover June, August and September")
    expected_totals = {"ARC": 11_200, "ME": 1_200, "VMN": 11_200, "OTHERS": 11_200}
    for row in sensitivities.itertuples(index=False):
        expected = expected_totals[str(row.region)]
        if (
            not bool(row.permutation_exhaustive)
            or int(row.permutations_evaluated) != expected
            or int(row.permutations_total) != expected
            or int(row.extreme_labelings) < 1
            or not (0 < float(row.p_value_raw) <= 1)
            or str(row.strata_no_information) != "September"
            or bool(row.cell_independence_claim)
        ):
            raise RuntimeError(f"Incomplete/invalid sensitivity receipt for {row.region}")
    for family in expected_families:
        rows = sensitivities.loc[sensitivities["sensitivity_family"].eq(family)]
        recalculated = holm_adjust(rows["p_value_raw"].astype(float).tolist())
        if not np.allclose(
            rows["p_value_holm_across_four_regions"].astype(float),
            recalculated, rtol=1e-12, atol=1e-14,
        ):
            raise RuntimeError(f"Four-region Holm adjustment changed for {family}")
    september = counts.loc[counts["acquisition_family"].eq("September")]
    if (
        not september["condition"].eq("Water").groupby(september["region"]).any().all()
        or september["informative_for_condition_test"].astype(bool).any()
    ):
        raise RuntimeError("September must be explicitly recorded as Water-only/no-information")


def write_sensitivity_methods(
    outdir: Path,
    mapping: pd.DataFrame,
    counts: pd.DataFrame,
    sensitivities: pd.DataFrame,
) -> None:
    count = sensitivities.loc[
        sensitivities["sensitivity_family"].eq("acquisition_stratified_count")
    ]
    lines = [
        "# Figure 2 acquisition-stratified sensitivity methods",
        "",
        "The mapping is determined solely from source acquisition directories:",
        "June (`24_06_2025`/`2025.06.24`), August (`Aug_2025_NPY`/selected runs),",
        "and September (`29.09.2025`). Condition labels are permuted as whole-animal",
        "units only within acquisition family, with the observed condition counts fixed.",
        "September contains Water only; it is retained in the mapping but has one fixed",
        "labeling and a zero contribution to the count statistic.",
        "",
        "Count sensitivity: `T_count = sum_s D_s`, where",
        "`D_s = 2 sum_g [y_sg log(y_sg/(n_sg p_s)) + (n_sg-y_sg)",
        "log((n_sg-y_sg)/(n_sg(1-p_s)))]`; `y_sg` and `n_sg` are pooled c-FOS",
        "and DAPI counts for assigned condition `g`, and `p_s=sum_g y_sg/sum_g n_sg`.",
        "Each animal's c-FOS count and DAPI denominator move together. The exact p value",
        "is the fraction of exhaustive constrained labelings with `T >= T_observed`.",
        "This is an animal-label randomization test and makes no cell-independence claim.",
        "",
        "Holm adjustment is applied across ARC, ME, VMN and OTHERS for the",
        "count sensitivity family.",
        "",
        "## Exact results",
        "",
        "```",
        count[[
            "region", "sensitivity_family", "statistic", "extreme_labelings",
            "permutations_total", "p_value_raw", "p_value_holm_across_four_regions",
        ]].to_string(index=False),
        "```",
        "",
        "The complete animal mapping is in `acquisition_strata.csv`; condition counts",
        "and informative/no-information status are in `acquisition_strata_counts.csv`.",
    ]
    (outdir / "SENSITIVITY_METHODS.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=_first_existing(DEFAULT_INPUT, SHIPPED_INPUT))
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--strict-paper-validation", action="store_true",
        help="Apply fixed canonical n, exact-permutation and IQR-flag assertions even to a noncanonical --input.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    source = args.input.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(source)
    values, audit, excluded = audit_and_prepare(raw)
    stats = inferential_statistics(values)
    outlier_audit = build_outlier_audit(values)
    strata, acquisition_mapping, acquisition_counts = build_acquisition_strata(values)
    sensitivities = acquisition_stratified_sensitivities(strata)
    input_sha256 = sha256(source)
    canonical_input = input_sha256 == CANONICAL_INPUT_SHA256
    fixed_validation_applied = canonical_input or args.strict_paper_validation
    if fixed_validation_applied:
        validate_paper_results(values, stats, outlier_audit)
        validate_sensitivity_results(
            acquisition_mapping, acquisition_counts, sensitivities,
        )
    condition_provenance = build_condition_provenance(values)
    desc = descriptives(values, "cfos_over_dapi_recomputed")
    density_desc = descriptives(values, "cfos_density_per_100k_um2")

    values.to_csv(outdir / "animal_region_values.csv", index=False)
    audit.to_csv(outdir / "analysis_inclusion_audit.csv", index=False)
    excluded.to_csv(outdir / "excluded_rows.csv", index=False)
    stats.to_csv(outdir / "global_cfos_statistics.csv", index=False)
    outlier_audit.to_csv(outdir / "outlier_audit.csv", index=False)
    desc.to_csv(outdir / "group_descriptives.csv", index=False)
    density_desc.to_csv(outdir / "density_sensitivity_descriptives.csv", index=False)
    condition_provenance.to_csv(outdir / "condition_provenance.csv", index=False)
    acquisition_mapping.to_csv(outdir / "acquisition_strata.csv", index=False)
    acquisition_counts.to_csv(outdir / "acquisition_strata_counts.csv", index=False)
    sensitivities.to_csv(
        outdir / "global_cfos_sensitivity_statistics.csv", index=False,
    )
    write_sensitivity_methods(
        outdir, acquisition_mapping, acquisition_counts, sensitivities,
    )
    source_inventory = pd.DataFrame(
        [
            {
                "role": "upstream_raw_per_animal_region_rows_before_paper_collapse",
                "path": str(source),
                "size_bytes": source.stat().st_size,
                "sha256": input_sha256,
                "read_only_upstream": True,
            },
            *[
                {
                    "role": f"condition_resolution_evidence_{role}",
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "read_only_upstream": True,
                }
                for role, path in CONDITION_EVIDENCE.items()
            ],
        ]
    )
    source_inventory.to_csv(outdir / "source_inventory.csv", index=False)

    manifest = {
        "script": str(HERE),
        "input": str(source),
        "input_sha256": input_sha256,
        "output_directory": str(outdir),
        "primary_endpoint": "cfos_cells / dapi_nuclei at the animal x region level",
        "regions": list(REGIONS),
        "region_roles": REGION_ROLES,
        "condition_order": list(CONDITIONS),
        "audited_condition_overrides": CONDITION_OVERRIDES,
        "condition_override_evidence": {
            role: {"path": str(path), "sha256": sha256(path)}
            for role, path in CONDITION_EVIDENCE.items()
        },
        "condition_provenance_file": "condition_provenance.csv",
        "collapse_rule": "sum counts/area across raw rows to one animal x canonical condition x canonical region before calculating the ratio",
        "exclusion": "exact case-insensitive animal ID 'horizontal', out-of-scope condition/region, and nonpositive DAPI denominator after collapse",
        "values_policy": "all finite values; no outlier exclusion",
        "descriptive_basis": "animal-level arithmetic mean +/- sample SD",
        "inference": "ordinary one-way ANOVA global and unadjusted exact two-sided rank-sum/Mann-Whitney U pairwise tests separately within each region",
        "global_p_value_method": "ordinary equal-variance one-way ANOVA F; visible p is raw",
        "global_p_value_additional_adjustment": "Holm across the four regional global ANOVA p values, reported in CSV/audit/legend but not substituted for visible labels",
        "pairwise_p_value_method": "exhaustive two-sided independent-label permutation of rank-sum/U",
        "pairwise_adjustment": "none; W-S, W-A and S-A are reported as unadjusted exact two-sided MWU p values",
        "effect_size": "Cliff's delta, displayed second condition minus first",
        "permutation_proof_columns": [
            "permutation_exhaustive", "permutations_evaluated", "permutations_total"
        ],
        "outlier_audit": {
            "file": "outlier_audit.csv",
            "rule": "1.5 x IQR within each treatment x region cell; descriptive flag only",
            "flagged": [
                {"region": region, "condition": condition, "animal": animal}
                for region, condition, animal in sorted(EXPECTED_IQR_FLAGS)
            ],
            "primary_inference": "all flagged values retained after HIL/source validation; no outcome-dependent exclusion",
            "figure_marker": "x for flagged values; circles otherwise",
        },
        "global_multiplicity": {
            "visible_plot_contract": "raw global one-way ANOVA p, matching current Figure 4 grammar",
            "additional_figure_wide_report": "Holm across ARC, ME, VMN and OTHERS",
        },
        "acquisition_stratified_sensitivities": {
            "files": {
                "statistics": "global_cfos_sensitivity_statistics.csv",
                "animal_mapping": "acquisition_strata.csv",
                "stratum_counts": "acquisition_strata_counts.csv",
                "methods": "SENSITIVITY_METHODS.md",
            },
            "strata": list(ACQUISITION_FAMILIES),
            "mapping_basis": "source acquisition directories only; never animal name or outcome",
            "permutation_unit": "whole animal; c-FOS count and DAPI denominator move together",
            "label_constraint": "condition counts fixed separately within acquisition stratum",
            "single_condition_policy": "recorded with one fixed labeling and zero/no-information contribution",
            "count_statistic": "sum of condition-rate versus acquisition-only binomial deviances using actual c-FOS/DAPI counts",
            "cell_independence_claim": False,
            "multiplicity": "Holm across the four regions separately within each sensitivity family",
            "results": json.loads(sensitivities.to_json(orient="records")),
        },
        "fixed_dataset_validation": {
            "canonical_input_sha256": CANONICAL_INPUT_SHA256,
            "input_is_canonical": canonical_input,
            "strict_flag": bool(args.strict_paper_validation),
            "applied": fixed_validation_applied,
            "expected_n": EXPECTED_N,
            "expected_anova_global_p": EXPECTED_GLOBAL_P,
            "expected_exact_pairwise_raw_cliffs_delta": {
                region: {
                    comparison: {
                        "p_raw": result[0],
                        "cliffs_delta_second_minus_first": result[1],
                    }
                    for comparison, result in comparisons.items()
                }
                for region, comparisons in EXPECTED_PAIRWISE.items()
            },
            "expected_descriptive_iqr_flags": [
                {"region": region, "condition": condition, "animal": animal}
                for region, condition, animal in sorted(EXPECTED_IQR_FLAGS)
            ],
            "expected_result": "unadjusted exact pairwise p values reported directly; all IQR-flagged values retained",
            "status": "passed" if fixed_validation_applied else "not_applied_to_noncanonical_input",
        },
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    (outdir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_audit(
        outdir, source, values, audit, stats, outlier_audit, sensitivities,
    )
    # Remove a stale diagnostics file from the earlier, different ANOVA implementation
    # only after the current ANOVA/MWU results and manifest have succeeded.
    (outdir / "anova_assumption_diagnostics.csv").unlink(missing_ok=True)

    print(f"[OK] Wrote audited global c-FOS analysis to {outdir}")
    print(values.groupby(["region", "condition"], observed=True)["animal"].nunique())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
