#!/usr/bin/env python3
"""Unequal-variance abundance tests for Figure 3.

The Water group spans a far wider range than either sugar group, and its
spread coincides with acquisition batch, so the equal-variance assumption of
the classical one-way ANOVA does not hold. This module adds tests that remain
valid under unequal variances and reports effect sizes that carry more
information than a p-value at these group sizes:

* Brown-Forsythe test of equal variance (Levene's test on absolute deviations
  from the group median), which documents the heteroscedasticity instead of
  assuming it away.
* Welch's one-way ANOVA, the unequal-variance omnibus.
* An exact studentized permutation test per pair: the Welch t statistic
  recomputed over every distinct allocation of the pooled animal values.
  Unlike a rank test, it targets a difference in means and stays valid when the
  two groups have different variances.
* Hodges-Lehmann median shift with an exact distribution-free confidence
  interval, and Cliff's delta.

Every test uses the biological animal as the unit and the same endpoint-specific
included set as the rest of Figure 3. All results remain exploratory: the Water
spread is confounded with acquisition batch, which no test can remove.

References: Welch (1951) Biometrika 38:330; Brown & Forsythe (1974) JASA 69:364;
Janssen (1997) Statist. Probab. Lett. 36:9; Hodges & Lehmann (1963) Ann. Math.
Statist. 34:598.
"""
from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution
from scipy.stats import f_oneway

CONDITIONS = ("Water", "Sucrose", "Allulose")
ABUNDANCE_ENDPOINTS = ("cfos_over_dapi", "double_over_npy", "npy_over_dapi")
SPATIAL_FEATURES = ("core_enrichment", "cluster_excess")
ALPHA = 0.05
CAVEAT = ("Exploratory. The Water group's dispersion coincides with acquisition batch, "
          "so an unequal-variance test describes the heterogeneity but cannot separate "
          "a treatment effect from acquisition.")


def welch_anova(groups: list[np.ndarray]) -> dict:
    """Welch's heteroscedastic one-way F test."""
    counts = np.array([len(group) for group in groups], dtype=float)
    means = np.array([group.mean() for group in groups])
    variances = np.array([group.var(ddof=1) for group in groups])
    k = len(groups)
    if k < 2 or (counts < 2).any() or not np.isfinite(variances).all() or (variances <= 0).any():
        return dict(statistic=np.nan, df1=np.nan, df2=np.nan, p_value=np.nan,
                    status="not_estimable_zero_variance_or_singleton_group")
    weights = counts / variances
    total = weights.sum()
    grand = float((weights * means).sum() / total)
    lam = ((1 - weights / total) ** 2 / (counts - 1)).sum()
    numerator = float((weights * (means - grand) ** 2).sum() / (k - 1))
    denominator = 1 + 2 * (k - 2) / (k * k - 1) * lam
    statistic = numerator / denominator
    df2 = (k * k - 1) / (3 * lam)
    p_value = float(f_distribution.sf(statistic, k - 1, df2))
    return dict(statistic=float(statistic), df1=float(k - 1), df2=float(df2),
                p_value=p_value, status="estimable")


def brown_forsythe(groups: list[np.ndarray]) -> dict:
    """Levene's equal-variance test using absolute deviations from the median."""
    deviations = [np.abs(group - np.median(group)) for group in groups]
    counts = [len(group) for group in groups]
    if len(groups) < 2 or min(counts) < 2:
        return dict(statistic=np.nan, df1=np.nan, df2=np.nan, p_value=np.nan,
                    status="not_estimable_fewer_than_two_per_group")
    if all(np.allclose(deviation, deviation[0]) for deviation in deviations):
        return dict(statistic=np.nan, df1=np.nan, df2=np.nan, p_value=np.nan,
                    status="not_estimable_degenerate_deviations")
    result = f_oneway(*deviations)
    return dict(statistic=float(result.statistic), df1=float(len(groups) - 1),
                df2=float(sum(counts) - len(groups)), p_value=float(result.pvalue), status="estimable")


def welch_t(left: np.ndarray, right: np.ndarray) -> float:
    variance = left.var(ddof=1) / len(left) + right.var(ddof=1) / len(right)
    if not np.isfinite(variance) or variance <= 0:
        return math.inf if left.mean() != right.mean() else 0.0
    return float((left.mean() - right.mean()) / math.sqrt(variance))


def exact_studentized_permutation(left: np.ndarray, right: np.ndarray) -> dict:
    """Enumerate every allocation of the pooled values and studentize each one.

    Studentizing is what keeps a label permutation meaningful when the two
    groups have different variances: the permutation null of an unstudentized
    mean difference is not the null of interest under heteroscedasticity.
    """
    n_left, n_right = len(left), len(right)
    if min(n_left, n_right) < 2:
        return dict(statistic=np.nan, p_value=np.nan, enumerated_labelings=0, extreme_labelings=0,
                    minimum_attainable_p=np.nan, status="fewer_than_two_animals_in_a_group")
    pooled = np.concatenate([left, right])
    observed = abs(welch_t(left, right))
    candidates = []
    for chosen in combinations(range(len(pooled)), n_left):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(chosen)] = True
        candidates.append(abs(welch_t(pooled[mask], pooled[~mask])))
    candidates = np.asarray(candidates, dtype=float)
    total = len(candidates)
    extreme = int((candidates >= observed - 1e-12).sum())
    floor = int((candidates >= candidates.max() - 1e-12).sum()) / total
    return dict(statistic=float(welch_t(left, right)), p_value=extreme / total,
                enumerated_labelings=total, extreme_labelings=extreme,
                minimum_attainable_p=float(floor), status="estimable")


def exact_permutation_welch_anova(groups: list[np.ndarray]) -> dict:
    """Welch's F with its null obtained by enumerating every animal allocation.

    The F approximation to Welch's statistic relies on denominator degrees of
    freedom that are themselves estimated from three-animal variances, so at
    these group sizes it is anti-conservative. Enumerating all distinct
    allocations of the pooled animal values replaces that approximation with
    the finite-sample reference distribution of the same robust statistic.
    """
    counts = [len(group) for group in groups]
    if len(groups) < 2 or min(counts) < 2:
        return dict(statistic=np.nan, p_value=np.nan, enumerated_labelings=0, extreme_labelings=0,
                    minimum_attainable_p=np.nan, status="not_estimable_fewer_than_two_per_group")
    pooled = np.concatenate(groups)
    total_n = len(pooled)
    observed = welch_anova(groups)["statistic"]
    if not np.isfinite(observed):
        return dict(statistic=np.nan, p_value=np.nan, enumerated_labelings=0, extreme_labelings=0,
                    minimum_attainable_p=np.nan, status="not_estimable_zero_variance_group")
    candidates = []
    remaining_first = list(combinations(range(total_n), counts[0]))
    for first in remaining_first:
        rest = [index for index in range(total_n) if index not in set(first)]
        for second in combinations(rest, counts[1]):
            third = [index for index in rest if index not in set(second)]
            split = [pooled[list(first)], pooled[list(second)], pooled[third]]
            value = welch_anova(split)["statistic"]
            candidates.append(value if np.isfinite(value) else np.inf)
    candidates = np.asarray(candidates, dtype=float)
    finite = candidates[np.isfinite(candidates)]
    extreme = int((candidates >= observed - 1e-12).sum())
    floor = int((candidates >= (finite.max() if len(finite) else observed) - 1e-12).sum()) / len(candidates)
    return dict(statistic=float(observed), p_value=extreme / len(candidates),
                enumerated_labelings=int(len(candidates)), extreme_labelings=extreme,
                minimum_attainable_p=float(floor), status="estimable")


def mann_whitney_null_counts(n_left: int, n_right: int) -> np.ndarray:
    """Exact null distribution of the Mann-Whitney U statistic, as counts.

    Counts the arrangements of the two samples by the number of pairs in which
    a left observation exceeds a right one, which is U. Each of the n_left
    left observations contributes between 0 and n_right such pairs, subject to
    the ordering constraint, so the counts follow the Gaussian binomial
    coefficient recurrence C(u, m, n) = C(u - n, m - 1, n) + C(u, m, n - 1).
    """
    table = np.zeros((n_left + 1, n_right + 1, n_left * n_right + 1))
    table[0, :, 0] = 1
    table[:, 0, 0] = 1
    for left in range(1, n_left + 1):
        for right in range(1, n_right + 1):
            for value in range(left * right + 1):
                total = table[left, right - 1, value]
                if value >= right:
                    total += table[left - 1, right, value - right]
                table[left, right, value] = total
    return table[n_left, n_right]


def hodges_lehmann(left: np.ndarray, right: np.ndarray, alpha: float = ALPHA) -> dict:
    """Median pairwise difference with an exact rank-based confidence interval."""
    differences = np.sort((left[:, None] - right[None, :]).ravel())
    estimate = float(np.median(differences))
    n_left, n_right = len(left), len(right)
    total = n_left * n_right
    counts = mann_whitney_null_counts(n_left, n_right)
    probabilities = counts / counts.sum()
    cumulative = np.cumsum(probabilities)
    k = int(np.searchsorted(cumulative, alpha / 2, side="right"))
    if k <= 0 or 2 * k >= total:
        return dict(estimate=estimate, ci_low=np.nan, ci_high=np.nan,
                    ci_level=1 - alpha, ci_method="exact Mann-Whitney rank interval",
                    ci_status="not_resolvable_at_this_group_size")
    return dict(estimate=estimate, ci_low=float(differences[k - 1]), ci_high=float(differences[total - k]),
                ci_level=float(1 - 2 * cumulative[k - 1]),
                ci_method="exact Mann-Whitney rank interval", ci_status="estimable")


def cliffs_delta(left: np.ndarray, right: np.ndarray) -> float:
    comparison = np.sign(left[:, None] - right[None, :])
    return float(comparison.sum() / comparison.size)


def endpoint_groups(values: pd.DataFrame, endpoint: str) -> dict[str, np.ndarray]:
    include = endpoint + "_include_in_main"
    groups = {}
    for condition in CONDITIONS:
        selected = values.loc[values.condition.eq(condition)]
        if include in selected:
            selected = selected.loc[selected[include].astype(str).str.lower().isin(["true", "1"])]
        numbers = pd.to_numeric(selected[endpoint], errors="coerce").to_numpy(float)
        groups[condition] = numbers[np.isfinite(numbers)]
    return groups


def rows_for(domain: str, endpoint: str, groups: dict[str, np.ndarray]) -> list[dict]:
    ordered = [groups[condition] for condition in CONDITIONS]
    counts = {"n_" + condition.lower(): len(groups[condition]) for condition in CONDITIONS}
    base = dict(domain=domain, endpoint=endpoint, unit="biological animal",
                analysis_set="endpoint-specific IQR included and estimable animals",
                selection_caveat=CAVEAT, **counts)
    rows = []
    variance = brown_forsythe(ordered)
    rows.append(dict(**base, test="brown_forsythe_equal_variance", group_a="", group_b="",
                     role="assumption check", **variance,
                     note="Levene's test on absolute deviations from each group's median; a small p indicates unequal spread."))
    permuted = exact_permutation_welch_anova(ordered)
    rows.append(dict(**base, test="exact_permutation_welch_anova", group_a="", group_b="", role="primary omnibus",
                     **permuted,
                     note=("Welch's robust F with its null enumerated over every distinct animal allocation; "
                           "no degrees-of-freedom approximation. minimum_attainable_p is the smallest p this "
                           "number of animals can produce.")))
    welch = welch_anova(ordered)
    rows.append(dict(**base, test="welch_one_way_anova", group_a="", group_b="", role="parametric omnibus",
                     **welch,
                     note=("Unequal-variance omnibus using the Welch-Satterthwaite F approximation. Its denominator "
                           "degrees of freedom are estimated from three-animal variances, so it is anti-conservative "
                           "here; read the enumerated permutation p above instead.")))
    if min(len(group) for group in ordered) >= 2:
        classical = f_oneway(*ordered)
        rows.append(dict(**base, test="classical_one_way_anova", group_a="", group_b="", role="sensitivity",
                         statistic=float(classical.statistic), df1=float(len(ordered) - 1),
                         df2=float(sum(len(group) for group in ordered) - len(ordered)),
                         p_value=float(classical.pvalue), status="estimable",
                         note="Equal-variance omnibus retained for comparison with the previously reported analysis."))
    for first, second in combinations(CONDITIONS, 2):
        left, right = groups[first], groups[second]
        permutation = exact_studentized_permutation(left, right)
        shift = hodges_lehmann(left, right) if min(len(left), len(right)) >= 1 else {}
        rows.append(dict(**base, test="exact_studentized_permutation", group_a=first, group_b=second,
                         role="primary pairwise", **permutation,
                         effect_name="Hodges-Lehmann median difference", **{
                             "effect": shift.get("estimate", np.nan),
                             "effect_ci_low": shift.get("ci_low", np.nan),
                             "effect_ci_high": shift.get("ci_high", np.nan),
                             "effect_ci_level": shift.get("ci_level", np.nan),
                             "effect_ci_method": shift.get("ci_method", ""),
                             "effect_ci_status": shift.get("ci_status", ""),
                         },
                         cliffs_delta=cliffs_delta(left, right) if len(left) and len(right) else np.nan,
                         note="Welch t recomputed over every distinct allocation of the pooled animal values."))
    return rows


def run(analysis_dir: Path) -> pd.DataFrame:
    analysis_dir = Path(analysis_dir).resolve()
    values = pd.read_csv(analysis_dir / "animal_plot_values.csv")
    rows: list[dict] = []
    for endpoint in ABUNDANCE_ENDPOINTS:
        rows.extend(rows_for("abundance", endpoint, endpoint_groups(values, endpoint)))
    metrics = pd.read_csv(analysis_dir / "animal_metrics.csv")
    if "include_in_main" in metrics:
        metrics = metrics.loc[metrics.include_in_main.astype(str).str.lower().isin(["true", "1"])]
    for spatial_endpoint in sorted(set(metrics.endpoint)):
        block = metrics.loc[metrics.endpoint.eq(spatial_endpoint)]
        for feature in SPATIAL_FEATURES:
            groups = {}
            for condition in CONDITIONS:
                numbers = pd.to_numeric(block.loc[block.condition.eq(condition), feature], errors="coerce").to_numpy(float)
                groups[condition] = numbers[np.isfinite(numbers)]
            rows.extend(rows_for("spatial", f"{spatial_endpoint}:{feature}", groups))
    table = pd.DataFrame(rows)
    ordered = ["domain", "endpoint", "test", "role", "group_a", "group_b",
               "n_water", "n_sucrose", "n_allulose", "statistic", "df1", "df2", "p_value",
               "minimum_attainable_p", "enumerated_labelings", "extreme_labelings", "effect_name", "effect",
               "effect_ci_low", "effect_ci_high", "effect_ci_level", "effect_ci_method",
               "effect_ci_status", "cliffs_delta", "status", "unit", "analysis_set", "note", "selection_caveat"]
    table = table.reindex(columns=[name for name in ordered if name in table.columns])
    table.to_csv(analysis_dir / "abundance_robust_statistics.csv", index=False)
    plan = dict(schema="fig3-robust-abundance-v1",
                motivation=("The Water group's variance greatly exceeds the sugar groups' and coincides with "
                            "acquisition batch, so the classical equal-variance ANOVA is not appropriate."),
                equal_variance_check="Brown-Forsythe (Levene on absolute deviations from the group median)",
                omnibus="Welch's robust F with an exact enumerated permutation null",
                parametric_omnibus="Welch's unequal-variance one-way ANOVA, reported as a comparison only",
                resolution_limit=("minimum_attainable_p records the smallest p each exact test can return at these "
                                  "group sizes; a 3-versus-3 contrast cannot fall below 0.10"),
                pairwise="exact studentized (Welch t) permutation over all distinct animal-label allocations",
                effect_size="Hodges-Lehmann median difference with an exact Mann-Whitney rank interval; Cliff's delta",
                retained_for_comparison=["classical one-way ANOVA", "exact two-sided Mann-Whitney"],
                multiplicity="nominal, unadjusted, matching the established Figure 3 abundance family",
                alpha=ALPHA, unit="biological animal", caveat=CAVEAT,
                references=["Welch 1951 Biometrika 38:330", "Brown & Forsythe 1974 JASA 69:364",
                            "Janssen 1997 Statist. Probab. Lett. 36:9", "Hodges & Lehmann 1963 Ann. Math. Statist. 34:598"])
    (analysis_dir / "abundance_robust_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    args = parser.parse_args()
    table = run(args.analysis_dir)
    shown = table.loc[table.domain.eq("abundance"),
                      ["endpoint", "test", "group_a", "group_b", "statistic", "p_value", "effect", "cliffs_delta"]]
    print(shown.to_string(index=False))


if __name__ == "__main__":
    main()
