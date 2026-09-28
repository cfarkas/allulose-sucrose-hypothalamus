#!/usr/bin/env python3
"""Fixed, animal-level spatial PERMANOVA plan for the Figure 3 extension.

The primary identifiable contrast is Allulose versus Sucrose after platform
and biological cohort, with labels permuted only within their joint strata.
The unrestricted three-condition test is an exploratory condition/acquisition
association. Neither cells nor acquisitions are independent replicates.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ENDPOINTS = ("cfos", "marker_cfos")
CONDITIONS = ("Water", "Sucrose", "Allulose")
FEATURES = ("core_enrichment", "cluster_excess")
PLAN = {
    "schema": "fig3-keyence-spatial-permanova-v1",
    "features": list(FEATURES),
    "feature_handling": "each measure is tested separately; there is no joint two-measure test",
    "feature_question": "spatial organization conditional on positive-cell abundance",
    "scaling": "pooled mean and sample SD within the eligible analysis set; fixed under every permutation",
    "distance": "Euclidean on pooled standardized features",
    "missingness": "complete cases per endpoint; never impute conditional metrics for zero positive cells",
    "primary_contrast": "Allulose versus Sucrose per measure, condition after platform and biological cohort",
    "primary_statistic": "partial pseudo-F from reduced/full multivariate residual sums of squares",
    "single_shared_platform_cohort": "If sugar animals share one platform and biological cohort, the reduced design is intercept-only and this is an ordinary within-cohort comparison, not platform adjustment",
    "primary_effect_size": "partial R-squared = condition sum of squares / reduced-model residual sum of squares",
    "primary_permutations": "all unique label allocations preserving condition counts within platform by biological_cohort strata",
    "exploratory_omnibus": ("unrestricted exact three-condition animal-label permutation test per measure; with "
                            "one standardized feature the pseudo-F is the ordinary one-way F and the null is an "
                            "exhaustive enumeration of animal-label allocations"),
    "multiplicity": "Holm across the four single-measure tests, two endpoints by two measures, within each family",
    "per_measure": ("Core enrichment and clustering excess are tested separately. No joint two-measure test is "
                    "performed, so no multivariate dispersion diagnostic is needed: the single-measure analogue "
                    "is the Brown-Forsythe equal-variance test reported in abundance_robust_statistics.csv."),
    "experimental_unit": "biological animal; cage independence remains unresolved unless authoritative cage metadata establish it",
    "source_references": [
        "https://vegandevs.github.io/vegan/reference/adonis.html",
        "https://vegandevs.github.io/vegan/reference/betadisper.html",
        "https://vegandevs.github.io/vegan/reference/permutest.betadisper.html",
    ],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_write(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def allocations(labels: np.ndarray, blocks: np.ndarray | None = None):
    """All distinct label allocations, including the observed allocation once."""
    labels = np.asarray(labels, dtype=str)
    if blocks is not None:
        blocks = np.asarray(blocks, dtype=str)
        indices = [np.flatnonzero(blocks == value) for value in sorted(set(blocks))]
        choices = [list(allocations(labels[index])) for index in indices]
        for selected in product(*choices):
            result = labels.copy()
            for index, values in zip(indices, selected):
                result[index] = values
            yield result
        return
    levels, sizes = np.unique(labels, return_counts=True)
    result = np.empty(len(labels), dtype=object)

    def recurse(level: int, remaining: tuple[int, ...]):
        if level == len(levels) - 1:
            result[list(remaining)] = levels[level]
            yield result.astype(str).copy()
            return
        for selected in combinations(remaining, int(sizes[level])):
            result[list(selected)] = levels[level]
            chosen = set(selected)
            yield from recurse(level + 1, tuple(index for index in remaining if index not in chosen))

    if len(labels):
        yield from recurse(0, tuple(range(len(labels))))


def allocation_count(labels: np.ndarray, blocks: np.ndarray | None = None) -> int:
    labels = np.asarray(labels)
    if blocks is not None:
        return math.prod(allocation_count(labels[np.asarray(blocks) == block]) for block in set(blocks))
    return math.factorial(len(labels)) // math.prod(math.factorial(int(n)) for n in np.unique(labels, return_counts=True)[1])


def standardize(matrix: np.ndarray):
    means = np.mean(matrix, axis=0)
    scales = np.std(matrix, axis=0, ddof=1)
    retained = np.isfinite(scales) & (scales > 1e-12)
    return (matrix[:, retained] - means[retained]) / scales[retained], means, scales, retained


def pseudo_f(matrix: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    levels = np.unique(labels)
    overall = matrix.mean(axis=0)
    total = float(np.square(matrix - overall).sum())
    within = 0.0
    between = 0.0
    for level in levels:
        group = matrix[labels == level]
        center = group.mean(axis=0)
        within += float(np.square(group - center).sum())
        between += float(len(group) * np.square(center - overall).sum())
    if total <= 1e-24:
        return 0.0, 0.0
    if within <= 1e-24:
        return math.inf, 1.0
    return ((between / (len(levels) - 1)) / (within / (len(matrix) - len(levels))), between / total)


def nuisance_design(frame: pd.DataFrame) -> np.ndarray:
    columns = [np.ones(len(frame))]
    for name in ("platform", "biological_cohort"):
        levels = sorted(frame[name].astype(str).unique())
        columns.extend(frame[name].astype(str).eq(level).to_numpy(dtype=float) for level in levels[1:])
    return np.column_stack(columns)


def residual_ss(matrix: np.ndarray, design: np.ndarray) -> tuple[float, int]:
    coefficients, _, rank, _ = np.linalg.lstsq(design, matrix, rcond=None)
    return float(np.square(matrix - design @ coefficients).sum()), int(rank)


def partial_f(matrix: np.ndarray, reduced: np.ndarray, labels: np.ndarray) -> dict:
    labels = np.asarray(labels)
    levels = sorted(set(labels))
    treatment = [np.equal(labels, level).astype(float) for level in levels[1:]]
    full = np.column_stack([reduced, *treatment])
    reduced_ss, reduced_rank = residual_ss(matrix, reduced)
    full_ss, full_rank = residual_ss(matrix, full)
    df_effect, df_residual = full_rank - reduced_rank, len(matrix) - full_rank
    effect_ss = max(0.0, reduced_ss - full_ss)
    output = dict(reduced_rank=reduced_rank, full_rank=full_rank, df_condition=df_effect,
                  df_residual=df_residual, reduced_residual_ss=reduced_ss,
                  full_residual_ss=full_ss, condition_ss=effect_ss,
                  partial_pseudo_f=np.nan, partial_r_squared=np.nan)
    if df_effect <= 0 or df_residual <= 0 or reduced_ss <= 1e-24:
        return output
    output["partial_pseudo_f"] = (effect_ss / df_effect) / (full_ss / df_residual) if full_ss > 1e-24 else math.inf
    output["partial_r_squared"] = effect_ss / reduced_ss
    return output


def extreme(statistic: float, observed: float) -> bool:
    if math.isinf(observed):
        return math.isinf(statistic)
    return statistic >= observed - 1e-10 * max(1.0, abs(observed))


def holm(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    available = np.isfinite(values)
    working = np.where(available, values, 1.0)
    order = np.argsort(working, kind="stable")
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(1.0, np.maximum.accumulate(working[order] * np.arange(len(values), 0, -1)))
    adjusted[~available] = np.nan
    return adjusted


def geometric_median(points: np.ndarray) -> np.ndarray:
    """Modified Weiszfeld iteration, including medians equal to input points."""
    points = np.asarray(points, dtype=float)
    center = points.mean(axis=0)
    for _ in range(20000):
        distances = np.linalg.norm(points - center, axis=1)
        nonzero = distances > 1e-12
        if not nonzero.any():
            return center
        weights = 1.0 / distances[nonzero]
        candidate = np.sum(points[nonzero] * weights[:, None], axis=0) / weights.sum()
        coincident = int((~nonzero).sum())
        if coincident:
            score = np.linalg.norm(np.sum((points[nonzero] - center) / distances[nonzero, None], axis=0))
            if score <= coincident:
                return center
            fraction = min(1.0, coincident / score)
            candidate = fraction * center + (1.0 - fraction) * candidate
        if np.linalg.norm(candidate - center) < 1e-11 * (1.0 + np.linalg.norm(center)):
            return candidate
        center = candidate
    raise RuntimeError("Spatial median failed to converge")


def dispersion(matrix: np.ndarray, labels: np.ndarray):
    distances = np.empty(len(labels))
    for level in np.unique(labels):
        selected = labels == level
        count = int(selected.sum())
        distances[selected] = np.linalg.norm(matrix[selected] - geometric_median(matrix[selected]), axis=1) * math.sqrt(count / (count - 1))
    observed, _ = pseudo_f(distances, labels)
    residuals = distances.copy()
    for level in np.unique(labels):
        selected = labels == level
        residuals[selected] -= distances[selected].mean()
    statistics = [pseudo_f(residuals, assigned)[0] for assigned in allocations(labels)]
    exceedances = sum(extreme(value, observed) for value in statistics)
    return dict(dispersion_f=observed, p_value_permutation=(exceedances + 1) / (len(statistics) + 1),
                enumerated_residual_allocations=len(statistics), extreme_residual_allocations=exceedances,
                p_calibration="(extreme+1)/(residual allocations+1)",
                centering="spatial median", bias_correction="sqrt(n/(n-1))",
                permutation_scheme="full-model distance residuals; centers/distances fixed from observed groups"), distances


def read_inputs(directory: Path, eligibility_csv: Path | None = None):
    roster = pd.read_csv(directory / "animal_roster.csv", dtype=str).fillna("unknown")
    required = {"animal_id", "condition", "platform", "biological_cohort", "cage"}
    if not required <= set(roster):
        raise ValueError(f"Roster lacks {sorted(required - set(roster))}")
    if roster.animal_id.duplicated().any():
        raise ValueError("Roster contains repeated biological animal IDs")
    if not set(roster.condition) <= set(CONDITIONS):
        raise ValueError("Unrecognized condition in animal roster")
    if roster[["platform", "biological_cohort"]].isin(["unknown", "", "nan"]).any().any():
        raise ValueError("Platform and biological cohort must be known for all roster animals")
    metrics = pd.read_csv(directory / "animal_metrics.csv")
    metrics = metrics.loc[metrics.cohort.eq("NPY") & metrics.region.eq("FIELD") & metrics.endpoint.isin(ENDPOINTS)].copy()
    if metrics.duplicated(["endpoint", "animal_id"]).any():
        raise ValueError("Multiple rows for one endpoint/animal; aggregate acquisitions before inference")
    if not set(metrics.animal_id) <= set(roster.animal_id):
        raise ValueError("Metrics contain animals absent from the roster")
    mapping = roster.set_index("animal_id").condition
    if not metrics.condition.eq(metrics.animal_id.map(mapping)).all():
        raise ValueError("Metric condition conflicts with the authoritative animal roster")
    payload = metrics.drop(columns=[name for name in roster.columns if name != "animal_id" and name in metrics])
    frames = []
    for endpoint in ENDPOINTS:
        selected = payload.loc[payload.endpoint.eq(endpoint)].drop(columns="endpoint")
        merged = roster.merge(selected, on="animal_id", how="left", validate="one_to_one")
        merged["endpoint"] = endpoint
        finite = np.isfinite(merged[list(FEATURES)].to_numpy(dtype=float))
        merged["joint_eligible"] = finite.all(axis=1)
        merged["joint_exclusion_reason"] = ["eligible" if row.all() else "unavailable:" + ";".join(np.asarray(FEATURES)[~row]) for row in finite]
        frames.append(merged)
    frame = pd.concat(frames, ignore_index=True)
    if eligibility_csv is not None:
        eligibility = pd.read_csv(eligibility_csv)
        required = {"endpoint", "animal_id", "include_in_main", "exclusion_reason"}
        if not required <= set(eligibility):
            raise ValueError("Eligibility table lacks required columns")
        if eligibility.duplicated(["endpoint", "animal_id"]).any():
            raise ValueError("Repeated endpoint/animal in eligibility table")
        included = eligibility.include_in_main.astype(str).str.lower()
        if not included.isin(["true", "false", "1", "0"]).all():
            raise ValueError("Eligibility include_in_main must be boolean")
        eligibility["external_include_in_main"] = included.isin(["true", "1"])
        eligibility = eligibility.rename(columns={"exclusion_reason": "external_exclusion_reason"})
        frame = frame.merge(eligibility[["endpoint", "animal_id", "external_include_in_main", "external_exclusion_reason"]],
                            on=["endpoint", "animal_id"], how="left", validate="one_to_one")
        if frame.external_include_in_main.isna().any():
            raise ValueError("Eligibility table must cover every endpoint/roster animal")
        excluded = frame.joint_eligible & ~frame.external_include_in_main
        frame.loc[excluded, "joint_exclusion_reason"] = "excluded_by_eligibility:" + frame.loc[excluded, "external_exclusion_reason"].fillna("unspecified").astype(str)
        frame["joint_eligible"] &= frame.external_include_in_main
    return roster, frame


def group_counts(selected: pd.DataFrame, roster: pd.DataFrame) -> dict:
    counts = dict(n_animals=len(selected))
    for condition in CONDITIONS:
        counts["n_" + condition.lower()] = int(selected.condition.eq(condition).sum())
        counts["n_recorded_" + condition.lower()] = int(roster.condition.eq(condition).sum())
    return counts


def design_diagnostic(frame: pd.DataFrame) -> dict:
    reduced = nuisance_design(frame)
    indicators = [frame.condition.eq(condition).to_numpy(dtype=float) for condition in CONDITIONS[1:]]
    full = np.column_stack([reduced, *indicators])
    reduced_rank = int(np.linalg.matrix_rank(reduced))
    full_rank = int(np.linalg.matrix_rank(full))
    return dict(n_animals=len(frame), nuisance_rank=reduced_rank, nuisance_columns=int(reduced.shape[1]),
                full_rank=full_rank, full_columns=int(full.shape[1]),
                estimable_condition_degrees_of_freedom=full_rank - reduced_rank,
                requested_condition_degrees_of_freedom=2,
                all_three_condition_contrasts_identifiable=bool(full_rank - reduced_rank == 2),
                condition_platform_counts=pd.crosstab(frame.condition, frame.platform).astype(int).to_dict(),
                condition_cohort_counts=pd.crosstab(frame.condition, frame.biological_cohort).astype(int).to_dict())


def run(directory: Path, eligibility_csv: Path | None = None) -> dict:
    directory = directory.resolve()
    # The plan is written before any newly extended animal values are read.
    plan = dict(PLAN)
    selection_note = ""
    if eligibility_csv is not None:
        eligibility_csv = Path(eligibility_csv).resolve()
        plan["eligibility_source"] = str(eligibility_csv)
        plan["eligibility_sha256"] = sha256(eligibility_csv)
        selection_note = "; exploratory after endpoint-specific within-condition IQR exclusion; filter held fixed during permutations, not selection-adjusted inference"
        plan["selection_caveat"] = selection_note.lstrip("; ")
    json_write(directory / "spatial_permanova_plan.json", plan)
    roster, frame = read_inputs(directory, eligibility_csv)
    eligibility_columns = ["endpoint", "animal_id", "condition", "platform", "biological_cohort", "cage", "joint_eligible", "joint_exclusion_reason", *FEATURES]
    eligibility_columns += [name for name in ("core_reason", "cluster_reason", "m", "shape_m", "iqr_excluded", "iqr_flagged_features", "iqr_exclusion_reason", "external_include_in_main", "external_exclusion_reason") if name in frame]
    frame[eligibility_columns].to_csv(directory / "spatial_feature_eligibility.csv", index=False)
    omnibus_rows, sugar_rows = [], []
    features_rows, scales_rows, allocation_rows, sugar_allocation_rows = [], [], [], []
    diagnostics = {"roster_design": design_diagnostic(roster), "endpoint_designs": {},
                   "water_platform_adjustment": "A Water contrast is not identified when Water has no platform/cohort overlap with sugar animals; omitting Leica Water removes that contrast rather than estimating it.",
                   "cage_status": "Animal-level exploratory inference; recorded cage strings alone do not establish independent treatment assignment.",
                   "cage_counts": roster.groupby(["condition", "biological_cohort", "cage"], dropna=False).size().rename("animals").reset_index().to_dict(orient="records")}
    for endpoint in ENDPOINTS:
        selected = frame.loc[frame.endpoint.eq(endpoint) & frame.joint_eligible].sort_values("animal_id").copy()
        diagnostics["endpoint_designs"][endpoint] = design_diagnostic(selected)
        counts = group_counts(selected, roster)
        valid = all(counts["n_" + condition.lower()] >= 2 for condition in CONDITIONS)
        matrix = means = scales = retained = None
        if valid:
            matrix, means, scales, retained = standardize(selected[list(FEATURES)].to_numpy(dtype=float))
            for i, feature in enumerate(FEATURES):
                scales_rows.append(dict(endpoint=endpoint, analysis="three_condition", feature=feature,
                                        pooled_mean=means[i], pooled_sd=scales[i], retained=bool(retained[i]),
                                        n_animals=len(selected)))
            for _, animal in selected.iterrows():
                entry = dict(endpoint=endpoint, analysis="three_condition", animal_id=animal.animal_id,
                             condition=animal.condition, platform=animal.platform,
                             biological_cohort=animal.biological_cohort)
                for j, feature in enumerate(FEATURES):
                    entry[feature] = animal[feature]
                    entry["z_" + feature] = (animal[feature] - means[j]) / scales[j] if retained[j] else np.nan
                features_rows.append(entry)

        sugar = selected.loc[selected.condition.isin(["Sucrose", "Allulose"])].copy()
        sugar_valid = len(sugar) and min(sugar.condition.value_counts().reindex(["Sucrose", "Allulose"], fill_value=0)) >= 2
        sugar_matrix = sugar_retained = None
        if sugar_valid:
            sugar_matrix, sugar_means, sugar_scales, sugar_retained = standardize(sugar[list(FEATURES)].to_numpy(dtype=float))
            for i, feature in enumerate(FEATURES):
                scales_rows.append(dict(endpoint=endpoint, analysis="sugar_platform_adjusted", feature=feature,
                                        pooled_mean=sugar_means[i], pooled_sd=sugar_scales[i],
                                        retained=bool(sugar_retained[i]), n_animals=len(sugar)))
            for _, animal in sugar.iterrows():
                entry = dict(endpoint=endpoint, analysis="sugar_platform_adjusted", animal_id=animal.animal_id,
                             condition=animal.condition, platform=animal.platform,
                             biological_cohort=animal.biological_cohort)
                for j, feature in enumerate(FEATURES):
                    entry[feature] = animal[feature]
                    entry["z_" + feature] = (animal[feature] - sugar_means[j]) / sugar_scales[j] if sugar_retained[j] else np.nan
                features_rows.append(entry)

        # Each displayed measure is tested on its own; there is no joint test.
        for index, feature in enumerate(FEATURES):
            base = dict(endpoint=endpoint, feature=feature, distance=PLAN["distance"],
                        scaling=PLAN["scaling"], **counts)
            result = dict(**base, status="unavailable", pseudo_f=np.nan, r_squared=np.nan, p_value_exact=np.nan,
                          enumerated_labelings=0, extreme_labelings=0,
                          interpretation=("Exploratory single-measure condition/acquisition association; Water "
                                          "treatment and acquisition are not separable; animal/cage independence "
                                          "unresolved"))
            if not valid:
                result["status"] = "fewer_than_two_complete_animals_in_a_condition"
            elif not retained[index]:
                result["status"] = "no_variation_in_this_measure"
            else:
                column = matrix[:, [index]]
                labels = selected.condition.to_numpy(dtype=str)
                observed, rsquared = pseudo_f(column, labels)
                count, exceedances, identities = 0, 0, 0
                for count, assigned in enumerate(allocations(labels), 1):
                    statistic, _ = pseudo_f(column, assigned)
                    is_observed = bool(np.array_equal(assigned, labels))
                    is_extreme = extreme(statistic, observed)
                    identities += is_observed
                    exceedances += is_extreme
                    allocation_rows.append(dict(endpoint=endpoint, feature=feature, allocation_index=count,
                        pseudo_f=statistic,
                        water_animals=";".join(selected.animal_id[assigned == "Water"]),
                        sucrose_animals=";".join(selected.animal_id[assigned == "Sucrose"]),
                        allulose_animals=";".join(selected.animal_id[assigned == "Allulose"]),
                        is_observed_allocation=is_observed, is_as_or_more_extreme=is_extreme))
                assert identities == 1 and count == allocation_count(labels)
                result.update(status="estimable_exploratory", pseudo_f=observed, r_squared=rsquared,
                              p_value_exact=exceedances / count, enumerated_labelings=count,
                              extreme_labelings=exceedances,
                              minimum_attainable_p=1.0 / count)
            result["interpretation"] += selection_note
            omnibus_rows.append(result)

            sugar_result = dict(endpoint=endpoint, feature=feature, **group_counts(sugar, roster),
                                status="unavailable", partial_pseudo_f=np.nan, partial_r_squared=np.nan,
                                p_value_exact=np.nan, enumerated_labelings=0, extreme_labelings=0,
                                excluded_water_animals=";".join(roster.loc[roster.condition.eq("Water"), "animal_id"]),
                                nuisance_covariates="platform;biological_cohort",
                                permutation_blocks="platform x biological_cohort",
                                interpretation=("Allulose-Sucrose association for this measure, conditional on "
                                                "acquisition platform and biological cohort; no Water contrast; "
                                                "animal/cage independence unresolved"))
            if not sugar_valid:
                sugar_result["status"] = "fewer_than_two_complete_animals_in_a_sugar_condition"
            elif not sugar_retained[index]:
                sugar_result["status"] = "no_variation_in_this_measure"
            else:
                column = sugar_matrix[:, [index]]
                labels = sugar.condition.to_numpy(dtype=str)
                blocks = (sugar.platform.astype(str) + "::" + sugar.biological_cohort.astype(str)).to_numpy()
                reduced = nuisance_design(sugar)
                observed = partial_f(column, reduced, labels)
                sugar_result.update(observed)
                if np.isfinite(observed["partial_pseudo_f"]) or math.isinf(observed["partial_pseudo_f"]):
                    count, exceedances, identities = 0, 0, 0
                    for count, assigned in enumerate(allocations(labels, blocks), 1):
                        statistic = partial_f(column, reduced, assigned)["partial_pseudo_f"]
                        is_extreme = extreme(statistic, observed["partial_pseudo_f"])
                        is_observed = bool(np.array_equal(assigned, labels))
                        exceedances += is_extreme
                        identities += is_observed
                        sugar_allocation_rows.append(dict(endpoint=endpoint, feature=feature,
                            allocation_index=count, partial_pseudo_f=statistic,
                            sucrose_animals=";".join(sugar.animal_id[assigned == "Sucrose"]),
                            allulose_animals=";".join(sugar.animal_id[assigned == "Allulose"]),
                            is_observed_allocation=is_observed, is_as_or_more_extreme=is_extreme))
                    assert identities == 1 and count == allocation_count(labels, blocks)
                    sugar_result.update(status="estimable_platform_cohort_adjusted",
                                        p_value_exact=exceedances / count, enumerated_labelings=count,
                                        extreme_labelings=exceedances, minimum_attainable_p=1.0 / count)
                else:
                    sugar_result["status"] = "condition_not_estimable_after_platform_cohort_or_no_residual_variation"
            if len(sugar) and sugar.platform.nunique() == 1 and sugar.biological_cohort.nunique() == 1:
                sugar_result["nuisance_covariates"] = "intercept only; single shared platform and biological cohort"
                sugar_result["permutation_blocks"] = "one original-cohort block"
                sugar_result["interpretation"] = ("Allulose-Sucrose comparison for this measure within the original "
                                                  "shared-platform cohort; no platform-adjustment or Water contrast "
                                                  "is identified; animal/cage independence unresolved")
                if sugar_result["status"] == "estimable_platform_cohort_adjusted":
                    sugar_result["status"] = "estimable_within_original_cohort"
            sugar_result["interpretation"] += selection_note
            sugar_rows.append(sugar_result)

    omnibus = pd.DataFrame(omnibus_rows)
    sugar = pd.DataFrame(sugar_rows)
    # Four tests per family: two endpoints by two separately tested measures.
    omnibus["p_holm_four_measure_family"] = holm(omnibus.p_value_exact)
    omnibus["family_size"] = len(omnibus)
    sugar["p_holm_four_measure_family"] = holm(sugar.p_value_exact)
    sugar["family_size"] = len(sugar)
    for filename, table in (("spatial_permanova.csv", omnibus), ("spatial_sugar_platform_sensitivity.csv", sugar),
                            ("spatial_permanova_features.csv", pd.DataFrame(features_rows)),
                            ("spatial_permanova_scaling.csv", pd.DataFrame(scales_rows)),
                            ("spatial_permanova_permutations.csv.gz", pd.DataFrame(allocation_rows)),
                            ("spatial_sugar_platform_permutations.csv.gz", pd.DataFrame(sugar_allocation_rows))):
        table.to_csv(directory / filename, index=False)
    json_write(directory / "spatial_platform_diagnostics.json", diagnostics)
    receipt = dict(schema=PLAN["schema"], status="PASS", inputs={name: sha256(directory / name) for name in ("animal_roster.csv", "animal_metrics.csv")},
                   script_sha256=sha256(Path(__file__)), plan_sha256=sha256(directory / "spatial_permanova_plan.json"),
                   recorded_animal_count=len(roster), animal_ids_unique=True,
                   no_missing_features_imputed=True, permutation_unit="biological animal", endpoint_count=2,
                   primary_family_size=4, exploratory_family_size=4,
                   per_measure_tests=int(len(omnibus)),
                   per_measure_family="Holm across the four single-measure exact tests",
                   joint_two_measure_test="not performed; each displayed measure is tested on its own",
                   caveat="Cage independence unresolved; observational exchangeability assumptions remain.")
    if eligibility_csv is not None:
        key = str(eligibility_csv.relative_to(directory)) if eligibility_csv.is_relative_to(directory) else str(eligibility_csv)
        receipt["inputs"][key] = sha256(eligibility_csv)
        receipt["conditional_on_data_dependent_selection"] = True
        receipt["selection_caveat"] = selection_note.lstrip("; ")
    json_write(directory / "spatial_permanova_validation.json", receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--eligibility-csv", type=Path, help="Explicit per-endpoint inclusion flags; observed values are preserved")
    args = parser.parse_args()
    print(json.dumps(run(args.analysis_dir, args.eligibility_csv), indent=2))


if __name__ == "__main__":
    main()
