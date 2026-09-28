#!/usr/bin/env python3
"""Auditable endpoint-specific IQR exclusions requested for Figure 3 main analyses.

The original nine-animal experiment is selected before calculating fences.
Keyence animals belong to a different experiment and remain outside all main
analyses. Raw observed values remain intact; inclusion flags define each endpoint.
"""
from __future__ import annotations
import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, rankdata
import fig3_keyence_permanova as spatial_stats
import shape_spatial_analysis as shape

ABUNDANCE_ENDPOINTS = ("cfos_over_dapi", "double_over_npy", "npy_over_dapi")
SPATIAL_ENDPOINTS = ("cfos", "marker_cfos")
SPATIAL_FEATURES = ("core_enrichment", "cluster_excess")
CONDITIONS = ("Water", "Sucrose", "Allulose")
CAVEAT = ("Exploratory after endpoint-specific within-condition IQR exclusion; "
          "filter held fixed during permutations, not selection-adjusted inference. "
          "Small groups give unstable quartiles; flags do not establish measurement error.")
PLAN = dict(schema="fig3-keyence-iqr-main-v1", rule="strictly outside Q1 - 1.5*IQR or Q3 + 1.5*IQR",
            quantile_method="numpy.quantile(method=linear)", stratification="within condition and endpoint/feature",
            replication_unit="biological animal", abundance_endpoints=list(ABUNDANCE_ENDPOINTS),
            spatial_endpoints=list(SPATIAL_ENDPOINTS), spatial_features=list(SPATIAL_FEATURES),
            spatial_exclusion="union of flags on core_enrichment and cluster_excess within the same joint endpoint",
            missingness="nonfinite values excluded from quartiles and never flagged as outliers",
            iteration="single pass; do not recompute fences after removing flags",
            selection_caveat=CAVEAT, cohort_selection="exclude Keyence before quartiles; original nine-animal experiment only",
            original_data="parent mixed-experiment results retained separately; not the main experiment")


def fence_flags(values):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    quartiles = np.quantile(values[finite], [.25, .75], method="linear") if finite.any() else [np.nan, np.nan]
    q1, q3 = map(float, quartiles)
    width = q3 - q1
    lower, upper = q1 - 1.5 * width, q3 + 1.5 * width
    flags = finite & ((values < lower) | (values > upper))
    return flags, dict(n_finite=int(finite.sum()), q1=q1, q3=q3, iqr=width, lower_fence=lower, upper_fence=upper)


def feature_ledger(frame, domain, endpoint, features):
    rows = []
    for condition, group in frame.groupby("condition", sort=True):
        for feature in features:
            flags, summary = fence_flags(group[feature].to_numpy(dtype=float))
            for (_, animal), flag in zip(group.iterrows(), flags):
                value = float(animal[feature])
                rows.append(dict(domain=domain, endpoint=endpoint, feature=feature, animal_id=animal.animal_id,
                                 condition=condition, value=value, **summary, is_outlier=bool(flag),
                                 is_missing=not np.isfinite(value), quantile_method="linear", multiplier=1.5))
    return pd.DataFrame(rows)


def endpoint_eligibility(ledger):
    rows = []
    for key, group in ledger.groupby(["domain", "endpoint", "animal_id", "condition"], sort=True):
        flagged = ";".join(sorted(group.loc[group.is_outlier, "feature"]))
        missing = ";".join(sorted(group.loc[group.is_missing, "feature"]))
        reason = "included"
        if flagged:
            reason = "IQR exclusion:" + flagged
        if missing:
            reason = (reason + "; " if flagged else "") + "nonestimable:" + missing
        rows.append(dict(zip(["domain", "endpoint", "animal_id", "condition"], key),
                         include_in_main=not bool(flagged or missing), iqr_excluded=bool(flagged),
                         nonestimable=bool(missing), flagged_features=flagged, missing_features=missing,
                         exclusion_reason=reason))
    return pd.DataFrame(rows)


def exact_mwu(left, right):
    ranks = rankdata(np.r_[left, right], method="average")
    n_left = len(left)
    expected = n_left * (len(ranks) + 1) / 2
    observed_sum = ranks[:n_left].sum()
    distance = abs(observed_sum - expected)
    total = math.comb(len(ranks), n_left)
    exceedances = sum(abs(ranks[list(indices)].sum() - expected) >= distance - 1e-12
                      for indices in combinations(range(len(ranks)), n_left))
    return dict(statistic=float(observed_sum - n_left * (n_left + 1) / 2),
                p_value=exceedances / total, enumerated_labelings=total, extreme_labelings=exceedances)


def abundance_statistics(wide):
    rows = []
    for endpoint in ABUNDANCE_ENDPOINTS:
        groups = {condition: wide.loc[wide.condition.eq(condition) & wide[endpoint + "_include_in_main"], endpoint].to_numpy(float)
                  for condition in CONDITIONS}
        counts = {"n_" + condition.lower(): len(groups[condition]) for condition in CONDITIONS}
        recorded = {"n_recorded_" + condition.lower(): int(wide.condition.eq(condition).sum()) for condition in CONDITIONS}
        base = dict(endpoint=endpoint, p_adjustment="none; nominal inherited abundance comparison",
                    analysis_set="endpoint-specific IQR included animals", selection_caveat=CAVEAT,
                    **counts, **recorded)
        anova = dict(statistic=np.nan, p_value=np.nan, enumerated_labelings=np.nan, extreme_labelings=np.nan,
                     status="fewer_than_two_included_animals_per_condition")
        if min(counts.values()) >= 2:
            result = f_oneway(*groups.values())
            anova.update(statistic=float(result.statistic), p_value=float(result.pvalue), status="estimable_after_IQR_exclusion")
        rows.append(dict(**base, test="one_way_ANOVA", group_a="", group_b="", **anova))
        for first, second in combinations(CONDITIONS, 2):
            result = dict(statistic=np.nan, p_value=np.nan, enumerated_labelings=0, extreme_labelings=0,
                          status="fewer_than_two_included_animals_in_a_group")
            if min(len(groups[first]), len(groups[second])) >= 2:
                result.update(exact_mwu(groups[first], groups[second]), status="estimable_after_IQR_exclusion")
            rows.append(dict(**base, test="exact_two_sided_Mann_Whitney", group_a=first, group_b=second, **result))
    return pd.DataFrame(rows)


def attach_spatial_flags(frame, eligibility, animal_column="animal_id"):
    flags = eligibility.loc[eligibility.domain.eq("spatial"),
                            ["endpoint", "animal_id", "include_in_main", "iqr_excluded", "flagged_features", "exclusion_reason"]].rename(
                                columns={"animal_id": animal_column, "flagged_features": "iqr_flagged_features", "exclusion_reason": "iqr_exclusion_reason"})
    result = frame.merge(flags, on=["endpoint", animal_column], how="left", validate="one_to_one")
    if result.include_in_main.isna().any():
        raise ValueError("Spatial table has animals/endpoints absent from the eligibility ledger")
    return result


def select_pairs(frame, pairs):
    return frame.merge(pairs, on=["endpoint", "animal_id"], how="inner", validate="many_to_one")


def run(parent: Path, output: Path | None = None):
    parent = parent.resolve()
    output = (output or parent / "iqr_main").resolve()
    if output == parent:
        raise ValueError("IQR results must not overwrite the unfiltered analysis directory")
    output.mkdir(parents=True, exist_ok=True)
    spatial_stats.json_write(output / "iqr_analysis_plan.json", PLAN)
    names = ("animal_plot_values.csv", "animal_metrics.csv", "unit_metrics.csv", "animal_roster.csv",
             "stratum_metrics.csv", "cluster_membership.csv.gz")
    source_hashes = {name: spatial_stats.sha256(parent / name) for name in names}
    wide = pd.read_csv(parent / "animal_plot_values.csv")
    animals = pd.read_csv(parent / "animal_metrics.csv")
    units = pd.read_csv(parent / "unit_metrics.csv")
    roster = pd.read_csv(parent / "animal_roster.csv")
    source_recorded_animals = len(roster)
    cohort_exclusions = roster.loc[roster.platform.eq("Keyence")].copy()
    cohort_exclusions["exclusion_stage"] = "before_IQR_fences"
    cohort_exclusions["exclusion_reason"] = "Different Keyence experiment; excluded from all main plots and inference by user instruction"
    cohort_exclusions.to_csv(output / "cohort_exclusions.csv", index=False)
    roster = roster.loc[~roster.platform.eq("Keyence")].copy()
    composition = roster.groupby("condition").size().to_dict()
    if set(composition) != set(CONDITIONS) or min(composition.values()) < 3:
        raise ValueError(f"Main analysis needs at least three animals per condition; found {composition}")
    original_ids = set(roster.animal_id)
    wide = wide.loc[wide.animal_id.isin(original_ids)].copy()
    animals = animals.loc[animals.animal_id.isin(original_ids)].copy()
    units = units.loc[units.unit_id.isin(original_ids)].copy()
    if wide.animal_id.duplicated().any() or animals.duplicated(["endpoint", "animal_id"]).any():
        raise ValueError("IQR filtering requires one value per biological animal and endpoint")
    if set(animals.cohort) != {"NPY"} or set(animals.region) != {"FIELD"}:
        raise ValueError("This IQR extension is scoped only to Figure 3 NPY FIELD endpoints")
    ledgers = [feature_ledger(wide, "abundance", endpoint, [endpoint]) for endpoint in ABUNDANCE_ENDPOINTS]
    ledgers += [feature_ledger(animals.loc[animals.endpoint.eq(endpoint)], "spatial", endpoint, SPATIAL_FEATURES)
                for endpoint in SPATIAL_ENDPOINTS]
    ledger = pd.concat(ledgers, ignore_index=True)
    eligibility = endpoint_eligibility(ledger)
    ledger.to_csv(output / "outlier_feature_ledger.csv", index=False)
    eligibility.to_csv(output / "eligibility.csv", index=False)
    bounds = ledger.drop_duplicates(["domain", "endpoint", "feature", "condition"])
    bounds[["domain", "endpoint", "feature", "condition", "n_finite", "q1", "q3", "iqr", "lower_fence", "upper_fence", "quantile_method", "multiplier"]].to_csv(output / "outlier_fences.csv", index=False)
    summary = eligibility.groupby(["domain", "endpoint", "condition"], sort=True).agg(
        n_recorded=("animal_id", "size"), n_included=("include_in_main", "sum"),
        n_iqr_excluded=("iqr_excluded", "sum"), n_nonestimable=("nonestimable", "sum")).reset_index()
    summary.to_csv(output / "inclusion_counts.csv", index=False)
    long_rows = []
    for endpoint in ABUNDANCE_ENDPOINTS:
        selection = eligibility.loc[eligibility.domain.eq("abundance") & eligibility.endpoint.eq(endpoint)].set_index("animal_id")
        wide[endpoint + "_iqr_excluded"] = wide.animal_id.map(selection.iqr_excluded)
        wide[endpoint + "_include_in_main"] = wide.animal_id.map(selection.include_in_main)
        for row in wide.itertuples():
            selected = selection.loc[row.animal_id]
            long_rows.append(dict(endpoint=endpoint, animal_id=row.animal_id, condition=row.condition,
                                  platform=row.platform, biological_cohort=row.biological_cohort,
                                  value=getattr(row, endpoint), value_percent=getattr(row, endpoint) * 100,
                                  include_in_main=bool(selected.include_in_main), iqr_excluded=bool(selected.iqr_excluded),
                                  exclusion_reason=selected.exclusion_reason))
    wide.to_csv(output / "animal_plot_values.csv", index=False)
    pd.DataFrame(long_rows).to_csv(output / "abundance_values.csv", index=False)
    abundance = abundance_statistics(wide)
    abundance.to_csv(output / "abundance_statistics.csv", index=False)
    animals = attach_spatial_flags(animals, eligibility)
    units = attach_spatial_flags(units, eligibility, "unit_id")
    animals.to_csv(output / "animal_metrics.csv", index=False)
    units.to_csv(output / "unit_metrics.csv", index=False)
    roster.to_csv(output / "animal_roster.csv", index=False)
    spatial_eligibility = eligibility.loc[eligibility.domain.eq("spatial")].drop(columns="domain")
    spatial_eligibility.to_csv(output / "spatial_eligibility.csv", index=False)
    included_animals = animals.loc[animals.include_in_main].copy()
    included_units = units.loc[units.include_in_main].copy()
    pairs = included_animals[["endpoint", "animal_id"]].drop_duplicates()
    condition_summary = shape.condition_summary(included_animals, included_units)
    cluster_counts = shape.cluster_counts(included_animals)
    recorded_counts = roster.groupby("condition").animal_id.nunique().to_dict()
    for table in (condition_summary, cluster_counts):
        table["recorded_animals"] = table.condition.map(recorded_counts)
        table["analysis_set"] = "endpoint-specific IQR included and estimable animals"
    condition_summary.to_csv(output / "condition_summary.csv", index=False)
    cluster_counts.to_csv(output / "cluster_counts_by_condition.csv", index=False)
    strata = select_pairs(pd.read_csv(parent / "stratum_metrics.csv"), pairs)
    members = select_pairs(pd.read_csv(parent / "cluster_membership.csv.gz"), pairs)
    strata.to_csv(output / "stratum_metrics.csv", index=False)
    members.to_csv(output / "cluster_membership.csv.gz", index=False)
    spatial_receipt = spatial_stats.run(output, output / "spatial_eligibility.csv")
    for name, expected in source_hashes.items():
        if spatial_stats.sha256(parent / name) != expected:
            raise RuntimeError("Unfiltered source changed during IQR analysis: " + name)
    for feature in SPATIAL_FEATURES:
        original = pd.read_csv(parent / "animal_metrics.csv").set_index(["endpoint", "animal_id"])[feature]
        current = animals.set_index(["endpoint", "animal_id"])[feature]
        if not np.allclose(original.loc[current.index], current, equal_nan=True, rtol=0, atol=0):
            raise AssertionError("Observed spatial values were altered")
    report = dict(schema=PLAN["schema"], status="PASS", unfiltered_inputs_unchanged=True,
                  condition_counts=composition,
                  unfiltered_source_directory=str(parent), input_sha256=source_hashes,
                  script_sha256=spatial_stats.sha256(Path(__file__)), recorded_animals=len(roster),
                  source_recorded_animals=source_recorded_animals, cohort_excluded_animals=cohort_exclusions.animal_id.tolist(),
                  cohort_selection_precedes_IQR=True,
                  raw_animal_rows_retained=len(animals), raw_unit_rows_retained=len(units),
                  missing_values_never_flagged=bool(not ledger.loc[ledger.is_missing, "is_outlier"].any()),
                  single_pass_fences=True, strict_fences=True, quantile_method="linear", selection_caveat=CAVEAT,
                  exclusions=eligibility.loc[eligibility.iqr_excluded, ["domain", "endpoint", "animal_id", "condition", "flagged_features"]].to_dict(orient="records"),
                  spatial_statistics_status=spatial_receipt["status"])
    spatial_stats.json_write(output / "iqr_validation.json", report)
    text = """# Figure 3 original-cohort main analysis with requested IQR rule

The main analysis first selects the original nine-animal experiment (three animals per condition). Keyence FR6-1, FR6-2 and FR722 belong to a different experiment and are excluded before any quartiles or inference. The parent directory retains the mixed-experiment results separately; they are not treated as a matched main-analysis comparator.

Within the original cohort, the main analysis excludes animal values strictly outside within-condition Q1 - 1.5 IQR and Q3 + 1.5 IQR. Quartiles use NumPy's linear interpolation, finite animal values, and one filtering pass. Values exactly on a fence stay included. No cells or sections are treated as independent replicates.

Abundance filtering is separate for c-FOS/DAPI, double-positive/NPY, and NPY/DAPI. Spatial filtering is separate for the all-DAPI c-FOS and NPY/c-FOS endpoints: a flag on either core enrichment or clustering excess excludes that animal from both features of the joint endpoint. Raw clustered percentage does not determine exclusion. Observed values and counts remain intact; explicit include_in_main flags identify the analysis rows. WATER_NPY3's unavailable NPY/c-FOS spatial metrics are structural missingness, not an outlier.

All tests are exploratory after endpoint-specific within-condition IQR exclusion. The observed inclusion set is held fixed during permutations. The p-values are conditional on this data-dependent selection and are not selection-adjusted randomization tests. Small groups give unstable quartiles, and an IQR flag does not establish measurement error. In particular, three finite values cannot be flagged by this linear-quartile 1.5-IQR rule; this does not establish that the Water distribution is free of unusual observations. The original observed values are retained in the main tables; all data from the separate Keyence experiment remain outside the main analysis.

The existing abundance test family is retained (ANOVA and exact two-sided Mann-Whitney comparisons, nominal p-values). Spatial tests retain the two declared features, pooled scaling, Allulose-Sucrose comparison within the original shared-platform cohort, with all 20 label allocations, exploratory three-condition PERMANOVA, dispersion diagnostic, and two-endpoint Holm adjustment within each family. Actual group sizes and available label allocations are recomputed after endpoint-specific exclusions. Water remains confounded with acquisition history; filtering does not resolve this or unknown legacy cage independence.

outlier_feature_ledger.csv records every animal/feature value, quartile, fence, and flag. outlier_fences.csv gives one row per condition/feature. eligibility.csv and inclusion_counts.csv distinguish included animals, IQR exclusions, and non-estimable observations. animal_plot_values.csv, animal_metrics.csv and unit_metrics.csv preserve all original rows with explicit flags; abundance_values.csv is the tidy abundance table. condition_summary.csv, cluster_counts_by_condition.csv, stratum_metrics.csv and cluster_membership.csv.gz contain only included endpoint/animal pairs. The nine-animal roster supplies original recorded counts. cohort_exclusions.csv records the Keyence experiment exclusions before IQR filtering. Parent contour geometries are unchanged; parent bandwidth sensitivity remains a mixed-experiment archival analysis and is not a main result.
"""
    (output / "README.md").write_text(text)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.analysis_dir, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
