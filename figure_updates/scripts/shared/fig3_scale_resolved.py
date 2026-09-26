#!/usr/bin/env python3
"""Scale-resolved clustering curves and a global envelope test for Figure 3.

A single clustering number at one distance cap answers "are activated cells
closer together than chance?" only at that distance. This module answers it at
every distance, and tests whether conditions differ anywhere in that range with
one p-value that already controls the error across all distances.

For each animal and endpoint the observed count of positive-positive pairs
within r is compared with its exact mean and variance under conditional random
labelling: positions, stratum membership and the number of positives are held
fixed, and only which cells are positive is permuted. Writing T(r) for the
count, P(r) for the number of pairs within r and Q(r) for the number of pairs
of pairs sharing a vertex,

    E[T]  = P p2
    E[T2] = P p2 + 2 (Q p3 + R p4),   R = C(P,2) - Q

with p_k the probability that k named cells are all positive when m of N are.
These moments are exact, so no simulation is needed for the per-animal curve.
The null conditions on the observed geometry. It does not remove differences
in sampling window, tissue coverage, section depth or marker detection.

Strata are independent under the null, so counts and moments add within an
animal and z(r) = (sum T - sum E) / sqrt(sum Var) is that animal's clustering
in standard-deviation units, comparable across animals with different cell
counts, positive counts and field shapes.

Conditions are then compared with a studentized global envelope test
(Myllymaki et al. 2017, J. R. Stat. Soc. B 79:381): the one-way F across
conditions is computed at every radius, standardized by its own permutation
mean and spread, and the maximum over radii is referred to the exhaustive
enumeration of animal-label allocations. One p-value covers the whole curve,
and the accompanying critical curve shows at which distances the conditions
separate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import fig3_keyence_permanova as permutation

CONDITIONS = ("Water", "Sucrose", "Allulose")
ENDPOINTS = ("cfos", "marker_cfos")
# Each cohort's sampled regions and the code they carry in the geometry raster.
# ME carries too few marker-associated cells to support a curve, so Figure 4's
# scale-resolved analysis is restricted to ARC.
REGION_CODES = {"NPY": {"FIELD": 1}, "POMC": {"ARC": 1}}
# Figure 3 compares animals; Figure 4 compares equal-weight cage means.
COHORT_UNIT = {"NPY": "animal", "POMC": "cage"}
RADII_UM = tuple(float(r) for r in range(20, 151, 10))
ALPHA = 0.05
CAVEATS = {
    "NPY": ("Exploratory. The random-labelling null conditions on observed cell positions, not on "
            "unobserved tissue or detection differences. Two Water animals share the sugar microscope, "
            "but acquisition/staining batches remain confounded. Finite-sample label-permutation "
            "validity requires exchangeable animal curves; studentizing F over radii does not "
            "remove groupwise heteroscedasticity. Marker curves require two or more positive cells."),
    "POMC": ("Exploratory ARC analysis with biological cages as independent units. The random-labelling "
             "null conditions on observed cell positions; it does not remove acquisition or sampling "
             "differences. Finite-sample label-permutation validity requires exchangeable cage curves; "
             "studentizing F over radii does not remove groupwise heteroscedasticity.")}


def resolve_geometry_path(value: str, bundle: Path) -> Path:
    """Resolve archived Paper paths inside the current tree before old absolutes."""
    path = Path(value)
    paper = next((parent for parent in (bundle, *bundle.parents)
                  if (parent / "scripts/shared").is_dir()), None)
    candidates = [bundle / path] if not path.is_absolute() else []
    if paper is not None:
        if not path.is_absolute():
            candidates.append(paper / path)
        elif "Paper" in path.parts:
            relative = Path(*path.parts[path.parts.index("Paper") + 1:])
            candidates.insert(0, paper / relative)
    candidates.append(path)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Missing geometry for {value}; restore the geometry inputs in the current Paper tree")


def pair_statistics(xy: np.ndarray, labels: np.ndarray, radii: np.ndarray):
    """Observed count, exact null mean and variance of positive pairs within r."""
    n = len(xy)
    counts = np.zeros(len(radii))
    expected = np.zeros(len(radii))
    variance = np.zeros(len(radii))
    pairs_within = np.zeros(len(radii))
    if n < 2:
        return counts, expected, variance, pairs_within
    tree = cKDTree(xy)
    pairs = np.asarray(sorted(tree.query_pairs(float(radii[-1]))), dtype=np.int64).reshape(-1, 2)
    if not len(pairs):
        return counts, expected, variance, pairs_within
    distance = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1)
    order = np.argsort(distance, kind="stable")
    pairs, distance = pairs[order], distance[order]
    both = labels[pairs[:, 0]] & labels[pairs[:, 1]]
    cut = np.searchsorted(distance, radii, side="right")
    m = int(labels.sum())
    cumulative_positive = np.concatenate([[0], np.cumsum(both)])
    for index, stop in enumerate(cut):
        p = int(stop)
        pairs_within[index] = p
        counts[index] = cumulative_positive[p]
        if p == 0 or m < 2:
            continue
        degree = np.bincount(pairs[:p].ravel(), minlength=n)
        q = float(np.sum(degree * (degree - 1) / 2))
        r_disjoint = p * (p - 1) / 2 - q
        p2 = m * (m - 1) / (n * (n - 1))
        p3 = m * (m - 1) * (m - 2) / (n * (n - 1) * (n - 2)) if n > 2 and m > 2 else 0.0
        p4 = (m * (m - 1) * (m - 2) * (m - 3) / (n * (n - 1) * (n - 2) * (n - 3))
              if n > 3 and m > 3 else 0.0)
        mean = p * p2
        second = p * p2 + 2 * (q * p3 + r_disjoint * p4)
        expected[index] = mean
        variance[index] = max(second - mean * mean, 0.0)
    return counts, expected, variance, pairs_within


def animal_curves(bundle: Path, radii: np.ndarray, cohort: str = "NPY") -> pd.DataFrame:
    """One standardized clustering curve per animal, region and endpoint."""
    regions = REGION_CODES[cohort]
    cells = pd.read_csv(bundle / f"inputs/{cohort}_cells.csv.gz", low_memory=False)
    manifest = pd.read_csv(bundle / f"inputs/{cohort}_acquisition_manifest.csv")
    if "include_acquisition" in manifest:
        manifest = manifest.loc[manifest.include_acquisition.astype(str).str.lower().isin(("true", "1"))]
    for column in ("is_cfos", "accepted_legacy", "marker_channel_present", "in_tissue"):
        cells[column] = cells[column].astype(str).str.lower().isin(("true", "1"))
    rows = []
    for record in manifest.to_dict("records"):
        with np.load(resolve_geometry_path(record["geometry_path"], bundle)) as archive:
            side_map, roi = archive["side"], archive["roi"]
        for region, code in regions.items():
            frame = cells.loc[cells.acquisition_id.eq(record["acquisition_id"]) & cells.region.eq(region)].copy()
            if not len(frame):
                continue
            iy = np.clip(frame.y_px.round().astype(int), 0, side_map.shape[0] - 1)
            ix = np.clip(frame.x_px.round().astype(int), 0, side_map.shape[1] - 1)
            frame["geometry_side"] = side_map[iy, ix]
            for side in np.unique(side_map[roi == code]):
                if side <= 0:
                    continue
                base = frame.loc[frame.geometry_side.eq(side)].sort_values("cell_id")
                base = base.loc[base.in_tissue]
                marker = base.loc[base.accepted_legacy & base.marker_channel_present]
                for endpoint in ENDPOINTS:
                    selected = base if endpoint == "cfos" else marker
                    xy = selected[["x_um", "y_um"]].to_numpy(float)
                    labels = selected.is_cfos.to_numpy(bool)
                    observed, expected, variance, pairs = pair_statistics(xy, labels, radii)
                    for index, radius in enumerate(radii):
                        rows.append(dict(animal_id=record["animal_id"], condition=record["condition"],
                                         cage=str(record.get("cage", "unknown")),
                                         acquisition_id=record["acquisition_id"], region=region,
                                         hemifield=f"SIDE_{int(side)}",
                                         endpoint=endpoint, radius_um=float(radius), cells=len(selected),
                                         positives=int(labels.sum()), pairs_within=float(pairs[index]),
                                         observed_pairs=float(observed[index]),
                                         expected_pairs=float(expected[index]),
                                         variance_pairs=float(variance[index])))
    strata = pd.DataFrame(rows)
    grouped = strata.groupby(["animal_id", "condition", "cage", "region", "endpoint", "radius_um"], as_index=False).agg(
        cells=("cells", "sum"), positives=("positives", "sum"), pairs_within=("pairs_within", "sum"),
        observed_pairs=("observed_pairs", "sum"), expected_pairs=("expected_pairs", "sum"),
        variance_pairs=("variance_pairs", "sum"), strata=("hemifield", "size"))
    spread = np.sqrt(grouped.variance_pairs.to_numpy(float))
    grouped["z"] = np.where(spread > 0, (grouped.observed_pairs - grouped.expected_pairs) / np.where(spread > 0, spread, 1), np.nan)
    grouped["enrichment"] = np.where(grouped.expected_pairs > 0,
                                     grouped.observed_pairs / np.where(grouped.expected_pairs > 0, grouped.expected_pairs, 1),
                                     np.nan)
    return strata, grouped


def one_way_f(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """One-way F at every radius; values is animals by radii."""
    groups = [values[labels == condition] for condition in CONDITIONS if (labels == condition).any()]
    n = len(values)
    k = len(groups)
    if k < 2 or n <= k:
        return np.full(values.shape[1], np.nan)
    grand = values.mean(axis=0)
    between = sum(len(g) * (g.mean(axis=0) - grand) ** 2 for g in groups) / (k - 1)
    within = sum(((g - g.mean(axis=0)) ** 2).sum(axis=0) for g in groups) / (n - k)
    return np.where(within > 0, between / np.where(within > 0, within, 1), np.where(between > 0, np.inf, 0.0))


def envelope_test(values: np.ndarray, labels: np.ndarray, radii: np.ndarray, alpha: float = ALPHA) -> dict:
    """Studentized maximum-deviation global envelope over an exact enumeration."""
    values, labels, radii = np.asarray(values, float), np.asarray(labels), np.asarray(radii, float)
    if (values.ndim != 2 or values.shape != (len(labels), len(radii))
            or not len(radii) or not np.isfinite(values).all() or not 0 < alpha < 1):
        raise ValueError("The envelope needs finite unit-by-radius values, matching labels/radii and 0 < alpha < 1")
    if not set(labels).issubset(CONDITIONS) or len(set(labels)) < 2:
        raise ValueError("At least two known treatment groups are required")
    statistics = []
    observed_index = None
    for index, assigned in enumerate(permutation.allocations(labels)):
        statistics.append(one_way_f(values, assigned))
        if observed_index is None and np.array_equal(assigned, labels):
            observed_index = index
    statistics = np.asarray(statistics, dtype=float)
    if observed_index is None:
        raise ValueError("The observed allocation was absent from the enumeration")
    if not np.isfinite(statistics).all():
        raise ValueError("Undefined or infinite F in the permutation distribution; the studentized envelope is not estimable")
    centre = statistics.mean(axis=0)
    spread = statistics.std(axis=0, ddof=1)
    informative = spread > 0
    deviation = np.zeros_like(statistics)
    deviation[:, informative] = (statistics[:, informative] - centre[informative]) / spread[informative]
    # Invariant radii provide no evidence. An entirely invariant dataset must
    # return p=1, not p=0 from comparisons against an all-NaN maximum.
    extremes = deviation[:, informative].max(axis=1) if informative.any() else np.zeros(len(statistics))
    observed_extreme = extremes[observed_index]
    p_value = float(np.mean(extremes >= observed_extreme - 1e-12))
    ordered = np.sort(extremes)[::-1]
    position = min(int(np.floor(alpha * len(ordered))), len(ordered) - 1)
    critical_extreme = float(ordered[position])
    return dict(p_value=p_value, allocations=int(len(extremes)),
                extreme_allocations=int(np.sum(extremes >= observed_extreme - 1e-12)),
                minimum_attainable_p=float(np.mean(extremes >= ordered[0] - 1e-12)),
                observed_statistic=float(observed_extreme), critical_statistic=critical_extreme,
                observed_f=statistics[observed_index], centre=centre, spread=spread,
                critical_curve=centre + critical_extreme * spread, alpha=alpha)


def run(bundle: Path, radii=RADII_UM, alpha: float = ALPHA, cohort: str = "NPY") -> dict:
    bundle = Path(bundle).resolve()
    radii = np.asarray(radii, dtype=float)
    unit = COHORT_UNIT[cohort]
    strata, curves = animal_curves(bundle, radii, cohort)
    strata.to_csv(bundle / "scale_resolved_strata.csv.gz", index=False, compression="gzip")
    eligibility = pd.read_csv(bundle / "spatial_eligibility.csv") if (bundle / "spatial_eligibility.csv").is_file() else None
    results, envelope_rows = [], []
    for region in REGION_CODES[cohort]:
        for endpoint in ENDPOINTS:
            block = curves.loc[curves.endpoint.eq(endpoint) & curves.region.eq(region)]
            wide = block.pivot_table(index=["animal_id", "condition", "cage"], columns="radius_um", values="z")
            # pivot_table can drop a radius with no finite values. Reindex so
            # that the claimed full distance range is actually required.
            usable = wide.reindex(columns=radii).dropna()
            if eligibility is not None and cohort == "NPY":
                allowed = eligibility.loc[eligibility.endpoint.eq(endpoint) & eligibility.include_in_main.astype(str).str.lower().isin(["true", "1"]), "animal_id"]
                usable = usable.loc[usable.index.get_level_values("animal_id").isin(set(allowed))]
            if unit == "cage":
                # Cages, not animals, are the independent units for this cohort.
                usable = usable.groupby(level=["condition", "cage"]).mean()
                labels = np.asarray([index[0] for index in usable.index], dtype=str)
            else:
                labels = np.asarray([index[1] for index in usable.index], dtype=str)
            counts = {f"n_{condition.lower()}": int((labels == condition).sum()) for condition in CONDITIONS}
            base = dict(cohort=cohort, region=region, endpoint=endpoint,
                        radii_um=";".join(f"{r:g}" for r in radii), **counts,
                        test="studentized_maximum_deviation_global_envelope",
                        unit=f"biological {unit}", selection_caveat=CAVEATS[cohort])
            if min(counts.values()) < 2 or len(usable) < 4:
                results.append(dict(**base, status=f"fewer_than_two_estimable_{unit}s_in_a_condition",
                                    p_value=np.nan, allocations=0, extreme_allocations=0,
                                    minimum_attainable_p=np.nan, observed_statistic=np.nan,
                                    critical_statistic=np.nan))
                continue
            outcome = envelope_test(usable.to_numpy(float), labels, radii, alpha)
            results.append(dict(**base, status="estimable_exploratory",
                                **{key: outcome[key] for key in ("p_value", "allocations", "extreme_allocations",
                                                                 "minimum_attainable_p", "observed_statistic",
                                                                 "critical_statistic")}))
            for index, radius in enumerate(radii):
                envelope_rows.append(dict(cohort=cohort, region=region, endpoint=endpoint,
                                          radius_um=float(radius),
                                          observed_f=float(outcome["observed_f"][index]),
                                          permutation_mean_f=float(outcome["centre"][index]),
                                          permutation_sd_f=float(outcome["spread"][index]),
                                          critical_f=float(outcome["critical_curve"][index]),
                                          exits_envelope=bool(outcome["observed_f"][index] > outcome["critical_curve"][index])))
    table = pd.DataFrame(results)
    if len(table):
        # Holm covers the two endpoints within each region.
        table["p_holm_two_endpoint_family"] = np.nan
        for _, index in table.groupby("region").groups.items():
            table.loc[index, "p_holm_two_endpoint_family"] = permutation.holm(table.loc[index, "p_value"])
    table.to_csv(bundle / "scale_resolved_test.csv", index=False)
    pd.DataFrame(envelope_rows).to_csv(bundle / "scale_resolved_envelope.csv", index=False)
    curves.to_csv(bundle / "scale_resolved_curves.csv", index=False)
    plan = dict(schema="fig3-scale-resolved-envelope-v1",
                question="Do conditions differ in how activated cells are spatially arranged, at any distance in the examined range?",
                curve="count of positive-positive pairs within r, standardized by its exact conditional random-labelling mean and standard deviation",
                moments="exact; E[T]=P*p2 and E[T^2]=P*p2+2(Q*p3+R*p4) with p_k the probability that k named cells are all positive",
                geometry="conditional on observed positions and field boundaries; no inference to unobserved tissue and no removal of cross-acquisition differences",
                strata="acquisition by tissue side; independent under the null, so counts and moments add within an animal",
                radii_um=[float(r) for r in radii],
                test="one-way F at each radius, studentized by its permutation mean and standard deviation, maximum taken over radii",
                null=f"exhaustive enumeration of {unit}-label allocations, the observed allocation included; exact test validity requires exchangeability",
                validity="Studentization is across radii in the permutation distribution, not a heteroscedastic groupwise Welch correction.",
                multiplicity="the maximum over radii controls the error across the whole curve; Holm then covers the two endpoints",
                alpha=alpha, cohort=cohort, unit=f"biological {unit}",
                regions=list(REGION_CODES[cohort]), caveat=CAVEATS[cohort],
                references=["Myllymaki et al. 2017 J. R. Stat. Soc. B 79:381",
                            "Baddeley, Rubak & Turner 2015, Spatial Point Patterns, chapters 7 and 14"])
    (bundle / "scale_resolved_plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    return dict(plan=plan, results=table)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--cohort", default="NPY", choices=sorted(REGION_CODES))
    args = parser.parse_args()
    outcome = run(args.analysis_dir, cohort=args.cohort)
    columns = ["region", "endpoint", "n_water", "n_sucrose", "n_allulose", "observed_statistic",
               "critical_statistic", "p_value", "p_holm_two_endpoint_family", "minimum_attainable_p",
               "allocations", "status"]
    print(outcome["results"][[c for c in columns if c in outcome["results"]]].to_string(index=False))


if __name__ == "__main__":
    main()
