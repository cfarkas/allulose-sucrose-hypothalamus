#!/usr/bin/env python3
"""Prepare Figure 3 from the extended cohort with the 2026 Water controls.

The cohort is the original nine Figure 3 animals plus the Water animals
acquired on 24 September 2026 on the same Zeiss LSM 780 and the same
LD LCI Plan-Apochromat 25x/0.8 objective as the Sucrose and Allulose animals.
The original nine animals' abundance counts and spatial measurements are copied
verbatim from the immutable nine-animal bundle and are verified by hash on
every rebuild; only the added animals are new.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import fig3_keyence_iqr as iqr
import fig3_keyence_permanova as statistics
import fig3_robust_abundance as robust
import fig3_scale_resolved as scale
import shape_spatial_analysis as shape

PAPER = Path(__file__).resolve().parents[2]
DEFAULT = PAPER / "analyses/fig3_water_extended_20260924"
ORIGINAL_NINE = PAPER / "analyses/fig3_original_cohort_20260920"
CONDITIONS = ("Water", "Sucrose", "Allulose")
ADDED_ROSTER = {
    "WATER_NPY4": dict(animal_id="WATER_NPY4", condition="Water", platform="Zeiss_LSM780",
                       biological_cohort="Water_2026_LSM780", cage="unknown", sections=1),
    "WATER_NPY5": dict(animal_id="WATER_NPY5", condition="Water", platform="Zeiss_LSM780",
                       biological_cohort="Water_2026_LSM780", cage="unknown", sections=1),
}
ADDED_ANIMALS = tuple(ADDED_ROSTER)
# Every ordered pair shown as a bracket in panels H and I.
PAIRWISE_COMPARISONS = (("Allulose", "Water"), ("Allulose", "Sucrose"), ("Water", "Sucrose"))
CONTRAST_PLAN = dict(schema="fig3-water-extended-pairwise-spatial-v1", experimental_unit="biological animal",
    comparisons=["Allulose - Water", "Allulose - Sucrose", "Water - Sucrose"], endpoints=["cfos", "marker_cfos"],
    metrics=["core_enrichment", "cluster_excess"], statistic="two-sided absolute difference of animal means",
    permutations="all distinct pairwise animal-label allocations, conditional on fixed eligible animals",
    multiplicity="Holm separately across six Figure 3 comparisons per metric (two endpoints by three comparisons)",
    other_figures_contribute=False,
    graph_panel_annotation="Observed clustering fraction is shown, but its comparison uses observed-minus-expected clustering excess",
    selection_caveat=iqr.CAVEAT)
SEED_SCHEMA = "fig3-water-extended-immutable-inputs-v1"


def copy(source: Path, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def added_abundance(cells: pd.DataFrame, animal: str) -> dict:
    rows = cells.loc[cells.animal_id.eq(animal)]
    if rows.empty:
        raise ValueError(f"The added animal is absent from the cell ledger: {animal}")
    roster = ADDED_ROSTER[animal]
    dapi = int(len(rows))
    npy = int(rows.accepted_legacy.sum())
    cfos = int(rows.is_cfos.sum())
    double = int((rows.accepted_legacy & rows.is_cfos).sum())
    return dict(animal=animal, animal_id=animal, condition="Water", sample=animal,
                platform=roster["platform"], biological_cohort=roster["biological_cohort"], sections=1,
                dapi_nuclei=dapi, npy_positive_nuclei=npy, cfos_positive_nuclei=cfos, double_positive_nuclei=double,
                cfos_over_dapi=cfos / dapi, double_over_npy=(double / npy) if npy else np.nan,
                npy_over_dapi=npy / dapi)


def seed(output: Path, spatial: Path, methodology_change: str | None = None):
    immutable = output / "original_inputs"
    manifest_path = immutable / "seed_manifest.json"
    if manifest_path.exists():
        raise ValueError("Immutable inputs already exist; seed only once")
    nine = ORIGINAL_NINE / "original_inputs"
    if not (nine / "seed_manifest.json").is_file():
        raise ValueError(f"The immutable nine-animal bundle is required: {nine}")
    nine_manifest = json.loads((nine / "seed_manifest.json").read_text())
    for name, expected in nine_manifest["files"].items():
        if statistics.sha256(nine / name) != expected:
            raise ValueError("The immutable nine-animal input changed: " + name)
    immutable.mkdir(parents=True, exist_ok=True)

    roster = pd.read_csv(nine / "animal_roster.csv")
    roster = pd.concat([roster, pd.DataFrame(list(ADDED_ROSTER.values()))[roster.columns]], ignore_index=True)
    composition = roster.groupby("condition").size().to_dict()
    if set(composition) != set(CONDITIONS) or min(composition.values()) < 3:
        raise ValueError(f"Unexpected cohort composition: {composition}")
    ids = set(roster.animal_id)
    roster.to_csv(immutable / "animal_roster.csv", index=False)

    source_cells = spatial.parent / "source/data/NPY_cells.csv.gz"
    cells = pd.read_csv(source_cells, low_memory=False)
    abundance = pd.read_csv(nine / "animal_plot_values.csv")
    original_ids = sorted(set(abundance.animal_id))
    added_rows = [added_abundance(cells, animal) for animal in ADDED_ANIMALS]
    abundance = pd.concat([abundance, pd.DataFrame(added_rows)[abundance.columns]], ignore_index=True)
    abundance.to_csv(immutable / "animal_plot_values.csv", index=False)

    history = [dict(source=str(nine / "animal_plot_values.csv"), source_sha256=statistics.sha256(nine / "animal_plot_values.csv"),
                    selected_animals=original_ids, output="animal_plot_values.csv",
                    note="Original nine abundance rows copied verbatim; " + ", ".join(ADDED_ANIMALS) + " appended from their own cell ledgers")]
    for name in ("animal_metrics.csv", "unit_metrics.csv", "stratum_metrics.csv", "cluster_membership.csv.gz"):
        frame = pd.read_csv(spatial / name, low_memory=False)
        identity = "animal_id" if "animal_id" in frame else "unit_id"
        frame = frame.loc[frame[identity].isin(ids)].copy()
        remove = [column for column in frame if column.startswith("iqr_") or column.endswith("_iqr_excluded")
                  or column.endswith("_include_in_main") or column == "include_in_main"]
        frame = frame.drop(columns=remove)
        frame.to_csv(immutable / name, index=False)
        history.append(dict(source=str(spatial / name), source_sha256=statistics.sha256(spatial / name),
                            selected_animals=sorted(ids), output=name))

    # The nine original animals keep the values the immutable nine-animal bundle
    # recorded, unless a measurement definition was deliberately changed. A
    # declared change is recorded value by value rather than waved through.
    reference = pd.read_csv(nine / "animal_metrics.csv").set_index(["endpoint", "animal_id"])
    checked = pd.read_csv(immutable / "animal_metrics.csv").set_index(["endpoint", "animal_id"])
    shared = checked.index.intersection(reference.index)
    if len(shared) != len(reference.index):
        raise ValueError("The recomputed cohort does not cover every original endpoint/animal")
    moved = []
    for feature in shape.METRICS:
        left = reference.loc[shared, feature].to_numpy(float)
        right = checked.loc[shared, feature].to_numpy(float)
        if np.allclose(left, right, equal_nan=True, rtol=0, atol=0):
            continue
        if methodology_change is None:
            raise ValueError("An original animal's spatial feature changed: " + feature
                             + "; declare the change with --methodology-change if it is intended")
        for (endpoint, animal), before, after in zip(shared, left, right):
            if not (np.isnan(before) and np.isnan(after)) and before != after:
                moved.append(dict(endpoint=endpoint, animal_id=animal, feature=feature,
                                  before=None if np.isnan(before) else float(before),
                                  after=None if np.isnan(after) else float(after),
                                  difference=None if np.isnan(after - before) else float(after - before)))
    # Counts are observations, not a modelling choice; they may never move.
    original_counts = pd.read_csv(nine / "animal_plot_values.csv").set_index("animal_id")
    current_counts = pd.read_csv(immutable / "animal_plot_values.csv").set_index("animal_id")
    for column in ("dapi_nuclei", "npy_positive_nuclei", "cfos_positive_nuclei", "double_positive_nuclei"):
        if not (current_counts.loc[original_counts.index, column] == original_counts[column]).all():
            raise ValueError("An original animal's observed counts changed: " + column)

    runtime = immutable / "runtime"
    (runtime / "inputs").mkdir(parents=True, exist_ok=True)
    cells.loc[cells.animal_id.isin(ids)].to_csv(runtime / "inputs/NPY_cells.csv.gz", index=False)
    acquisition = pd.read_csv(spatial.parent / "source/audit/NPY_acquisition_manifest.csv")
    acquisition = acquisition.loc[acquisition.animal_id.isin(ids)].copy()
    if len(acquisition) != len(ids) or acquisition.animal_id.nunique() != len(ids):
        raise ValueError("Every cohort animal needs exactly one accepted acquisition")
    acquisition = acquisition.merge(roster[["animal_id", "platform", "biological_cohort"]], on="animal_id", validate="one_to_one")
    geometry_records = []
    for index, row in acquisition.iterrows():
        source = Path(row.geometry_path)
        relative = Path("geometry") / (row.acquisition_id + "_geometry.npz")
        copy(source, runtime / relative)
        geometry_records.append(dict(acquisition_id=row.acquisition_id, source=str(source),
                                     source_sha256=statistics.sha256(source), local_path=str(relative)))
        acquisition.loc[index, "geometry_path"] = str(output / relative)
    acquisition.to_csv(runtime / "inputs/NPY_acquisition_manifest.csv", index=False)

    contours = pd.read_csv(spatial / "contour_manifest.csv")
    contours = contours.loc[contours.cohort.eq("NPY") & contours.animal_id.isin(ids)].copy()
    contours.to_csv(runtime / "contour_manifest.csv", index=False)
    for row in contours.itertuples():
        if isinstance(row.contour_path, str) and row.contour_path:
            copy(spatial / row.contour_path, runtime / row.contour_path)
    sensitivity = pd.read_csv(spatial / "bandwidth_sensitivity_condition_summary.csv")
    sensitivity.loc[sensitivity.cohort.eq("NPY")].to_csv(runtime / "bandwidth_sensitivity_condition_summary.csv", index=False)
    copy(ORIGINAL_NINE / "original_inputs/runtime/example_npy.npz", runtime / "example_npy.npz")

    added_receipt = json.loads((spatial.parent / "source/new_water_cells_receipt.json").read_text())
    provenance = dict(schema="fig3-water-extended-input-provenance-v1", cohort_animal_ids=sorted(ids),
        declared_methodology_change=methodology_change,
        original_spatial_values_changed_by_that_declaration=moved,
        conditions=composition, acquisition_count=len(ids),
        original_nine=dict(bundle=str(ORIGINAL_NINE), seed_manifest_sha256=statistics.sha256(nine / "seed_manifest.json"),
                           animal_ids=original_ids, values_copied_verbatim=True),
        added_animals={animal: added_receipt["per_animal"][animal] for animal in ADDED_ANIMALS},
        added_cells_receipt_sha256=statistics.sha256(spatial.parent / "source/new_water_cells_receipt.json"),
        cell_ledger_source=dict(path=str(source_cells), sha256=statistics.sha256(source_cells)),
        acquisition_manifest_source=dict(path=str(spatial.parent / "source/audit/NPY_acquisition_manifest.csv"),
                                         sha256=statistics.sha256(spatial.parent / "source/audit/NPY_acquisition_manifest.csv")),
        spatial_analysis_source=dict(path=str(spatial), metadata_sha256=statistics.sha256(spatial / "analysis_metadata.json")),
        original_geometry=geometry_records, table_seed_history=history,
        original_spatial_values_independently_matched=not moved,
        original_observed_counts_unchanged=True,
        accepted_associations="accepted_legacy nuclear marker identities; original Figure 3 masks unchanged",
        runtime_policy="Only this cohort's fixed copies are read during preparation; no other experiment is extracted or analyzed")
    statistics.json_write(immutable / "source_provenance.json", provenance)
    files = {str(path.relative_to(immutable)): statistics.sha256(path) for path in sorted(immutable.rglob("*")) if path.is_file()}
    statistics.json_write(manifest_path, dict(schema=SEED_SCHEMA, animal_ids=sorted(ids), animal_count=len(ids),
                                              condition_counts=composition, added_animals=list(ADDED_ANIMALS), files=files))
    for name in ("before_change_hashes.json", "figure4_before_hashes.json"):
        source = ORIGINAL_NINE / name
        if source.is_file():
            copy(source, output / name)
    # The publisher compares the incoming animal table with the currently published one.
    archive = output / "archive/Fig3/source_data"
    for name in ("Figure_3_cFos_NPY_spatial_animal_counts.csv", "Figure_3_cFos_NPY_animal_plot_values.csv",
                 "Figure_3_cFos_NPY_spatial_segmentation_model_setting_audit.csv",
                 "Figure_3_cFos_NPY_legacy_nine_segmentation_model_setting_audit.csv"):
        for candidate in (PAPER / "Fig3/source_data" / name, ORIGINAL_NINE / "archive/Fig3/source_data" / name):
            if candidate.is_file():
                copy(candidate, archive / name)
                break


def verify_rederived(output: Path, spatial: Path):
    """Confirm a from-raw rebuild reproduces the hash-bound cohort inputs.

    Seeding happens once; every later rebuild reads the frozen tables. When the
    two added animals are re-derived from their CZI files, that frozen copy is
    only meaningful if the fresh computation still lands on the same numbers,
    so this compares them column by column instead of assuming it.
    """
    output, spatial = output.resolve(), spatial.resolve()
    immutable = output / "original_inputs"
    manifest = json.loads((immutable / "seed_manifest.json").read_text())
    ids = set(manifest["animal_ids"])
    cells = pd.read_csv(spatial.parent / "source/NPY_cells.csv.gz", low_memory=False) if (
        spatial.parent / "source/NPY_cells.csv.gz").is_file() else pd.read_csv(
        spatial.parent / "source/data/NPY_cells.csv.gz", low_memory=False)
    differences = []
    frozen_abundance = pd.read_csv(immutable / "animal_plot_values.csv").set_index("animal_id")
    for animal in ADDED_ANIMALS:
        fresh = added_abundance(cells, animal)
        recorded = frozen_abundance.loc[animal]
        for column, value in fresh.items():
            if column in ("animal", "animal_id", "sample"):
                continue
            previous = recorded[column]
            same = (np.isclose(float(value), float(previous), rtol=0, atol=1e-12)
                    if isinstance(value, (int, float)) and not isinstance(value, bool)
                    else str(value) == str(previous))
            if not same:
                differences.append(f"{animal}.{column}: {previous!r} -> {value!r}")
    for name in ("animal_metrics.csv", "unit_metrics.csv", "stratum_metrics.csv"):
        fresh = pd.read_csv(spatial / name, low_memory=False)
        identity = "animal_id" if "animal_id" in fresh else "unit_id"
        fresh = fresh.loc[fresh[identity].isin(ids)].reset_index(drop=True)
        frozen = pd.read_csv(immutable / name, low_memory=False).reset_index(drop=True)
        fresh = fresh[frozen.columns]
        if len(fresh) != len(frozen):
            differences.append(f"{name}: {len(frozen)} rows recorded, {len(fresh)} recomputed")
            continue
        for column in frozen.columns:
            left, right = frozen[column], fresh[column]
            if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
                if not np.allclose(left.to_numpy(float), right.to_numpy(float), equal_nan=True, rtol=0, atol=1e-12):
                    differences.append(f"{name}:{column}")
            elif not left.astype(str).equals(right.astype(str)):
                differences.append(f"{name}:{column}")
    report = dict(schema="fig3-water-extended-rederivation-check-v1",
                  status="PASS" if not differences else "FAIL",
                  spatial_directory=str(spatial), differences=differences,
                  added_animals=list(ADDED_ANIMALS), cohort_animals=sorted(ids))
    statistics.json_write(output / "rederivation_check.json", report)
    if differences:
        raise ValueError("A from-raw rebuild did not reproduce the recorded cohort inputs: "
                         + "; ".join(differences[:8]))
    print(json.dumps(report, indent=2))
    return report


def verify_inputs(output: Path):
    immutable = output / "original_inputs"
    manifest = json.loads((immutable / "seed_manifest.json").read_text())
    if manifest.get("schema") != SEED_SCHEMA:
        raise ValueError("Invalid immutable input manifest for this cohort")
    for name, expected in manifest["files"].items():
        if statistics.sha256(immutable / name) != expected:
            raise ValueError("Immutable input changed: " + name)
    roster = pd.read_csv(immutable / "animal_roster.csv")
    if roster.platform.eq("Keyence").any() or roster.animal_id.str.contains("Keyence", case=False).any():
        raise ValueError("Different-experiment animals entered the cohort inputs")
    if set(roster.animal_id) != set(manifest["animal_ids"]):
        raise ValueError("Immutable roster does not match the seeded animal identities")
    if not set(manifest["added_animals"]) <= set(roster.animal_id):
        raise ValueError("An added Water animal is missing from the roster")
    return manifest


def pairwise_contrasts(output: Path):
    statistics.json_write(output / "treatment_contrast_plan.json", CONTRAST_PLAN)
    units = pd.read_csv(output / "unit_metrics.csv")
    units = units.loc[units.include_in_main].copy()
    rows = []
    for metric in ("core_enrichment", "cluster_excess"):
        for endpoint in ("cfos", "marker_cfos"):
            selected = units.loc[units.endpoint.eq(endpoint)]
            for first, second in PAIRWISE_COMPARISONS:
                result = shape.exact_comparison(selected.loc[selected.condition.eq(first), metric],
                                                selected.loc[selected.condition.eq(second), metric])
                result.update(cohort="NPY", endpoint=endpoint, region="FIELD", metric=metric,
                    group_a=first, group_b=second, comparison=f"{first} - {second}",
                    family=metric, planned_family_size=2 * len(PAIRWISE_COMPARISONS), unit_type="animal",
                    interpretation="Exploratory Figure 3 cohort only; " + iqr.CAVEAT)
                rows.append(result)
    table = pd.DataFrame(rows)
    expected = 2 * len(PAIRWISE_COMPARISONS)
    for _, indices in table.groupby("family").groups.items():
        if len(indices) != expected:
            raise ValueError(f"Each Figure 3 metric must have {expected} planned comparisons")
        table.loc[indices, "p_holm"] = statistics.holm(table.loc[indices, "p_raw"])
    table.to_csv(output / "treatment_contrasts.csv", index=False)


def prepare(output: Path):
    output = output.resolve()
    manifest = verify_inputs(output)
    immutable = output / "original_inputs"
    saved_readme = (output / "README.md").read_text() if (output / "README.md").is_file() else None
    for source in sorted((immutable / "runtime").rglob("*")):
        if source.is_file():
            copy(source, output / source.relative_to(immutable / "runtime"))
    report = iqr.run(immutable, output)
    pairwise_contrasts(output)
    robust_table = robust.run(output)
    scale_outcome = scale.run(output)
    if saved_readme is not None:
        (output / "README.md").write_text(saved_readme)
    counts = manifest["condition_counts"]
    provenance = json.loads((immutable / "source_provenance.json").read_text())
    declared = provenance.get("declared_methodology_change")
    moved = provenance.get("original_spatial_values_changed_by_that_declaration") or []
    metadata = dict(schema="fig3-water-extended-cohort-analysis-v1", cohort_animal_ids=manifest["animal_ids"],
        condition_counts=counts, added_animals=manifest["added_animals"], independent_of_other_experiments=True,
        immutable_input_manifest=dict(path=str(immutable / "seed_manifest.json"), sha256=statistics.sha256(immutable / "seed_manifest.json")),
        input_files=[dict(path=str(immutable / name), sha256=value) for name, value in manifest["files"].items()],
        analysis_script_sha256=statistics.sha256(Path(__file__)), iqr_script_sha256=statistics.sha256(Path(iqr.__file__)),
        statistics_script_sha256=statistics.sha256(Path(statistics.__file__)),
        robust_abundance_script_sha256=statistics.sha256(Path(robust.__file__)),
        scale_resolved_script_sha256=statistics.sha256(Path(scale.__file__)),
        scale_resolved=dict(radii_um=scale_outcome["plan"]["radii_um"],
                            test=scale_outcome["plan"]["test"],
                            endpoints=int(len(scale_outcome["results"]))),
        joint_two_measure_test="not performed; core enrichment and clustering excess are tested separately",
        unequal_variance_tests=dict(
            equal_variance_check="Brown-Forsythe", omnibus="Welch one-way ANOVA",
            pairwise="exact studentized (Welch t) permutation",
            effect_size="Hodges-Lehmann shift with exact rank interval; Cliff's delta",
            rows=int(len(robust_table)),
            reason="The Water group's variance greatly exceeds the sugar groups' and coincides with acquisition batch"),
        parameters=dict(k=6, edge_cap_um=75., bandwidth_um=40., grid_um=5., contour_mass=.8, minimum_marker_anchors=5),
        interpretation=("Original nine animals plus further Water controls on the sugar cohort's own microscope; "
                        "abundance and spatial IQR rules applied per endpoint; acquisition session and staining batch "
                        "still differ between Water and the sugar groups"),
        spatial_recomputation=("Spatial features recomputed for every animal; the original nine reproduce their "
                               "recorded values exactly" if not moved else
                               "Spatial features recomputed for every animal under a declared measurement change; "
                               "see declared_methodology_change"),
        declared_methodology_change=declared,
        original_spatial_values_changed=len(moved),
        native_masks_modified=False)
    statistics.json_write(output / "analysis_metadata.json", metadata)
    validation = dict(status="PASS", animals=manifest["animal_count"], conditions=counts,
        original_nine_observed_counts_unchanged=True,
        original_nine_spatial_values_unchanged=not moved,
        declared_methodology_change=declared,
        original_spatial_values_changed_by_declaration=len(moved),
        cohort_only=True,
        no_other_experiment_rows=True, immutable_input_hashes_verified=True, acquisitions=manifest["animal_count"],
        missing_conditional_features_not_imputed=True, original_images_and_masks_modified=False,
        iqr_outlier_count=len(report["exclusions"]),
        pairwise_figure3_only_family_size=2 * len(PAIRWISE_COMPARISONS),
        unequal_variance_tests_written=bool(len(robust_table)),
        scale_resolved_envelope_tests=int(len(scale_outcome["results"])),
        joint_two_measure_permanova_performed=False,
        no_figure4_data_in_statistics=True)
    statistics.json_write(output / "analysis_validation.json", validation)
    verify_inputs(output)
    print(json.dumps(validation, indent=2))
    return validation


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT)
    parser.add_argument("--seed-from-spatial", type=Path, help="One-time seeding source: this cohort's spatial analysis directory")
    parser.add_argument("--verify-rederived", type=Path, help="Check a from-raw rebuild against the hash-bound inputs and return")
    parser.add_argument("--methodology-change", help="Declare, and record per animal, a deliberate change to a spatial measurement definition")
    parser.add_argument("--snapshot-figure4", type=Path, help="Write the current Figure 4 publication guard and return")
    args = parser.parse_args()
    if args.verify_rederived is not None:
        verify_rederived(args.analysis_dir.resolve(), args.verify_rederived.resolve())
        return
    if args.snapshot_figure4 is not None:
        hashes = {name: statistics.sha256(PAPER / name) for name in ("Fig4/Figure_4.pdf", "Fig4/Figure_4.png")}
        args.snapshot_figure4.parent.mkdir(parents=True, exist_ok=True)
        statistics.json_write(args.snapshot_figure4, hashes)
        print(json.dumps(dict(status="PASS", figure4_snapshot=str(args.snapshot_figure4)), indent=2))
        return
    output = args.analysis_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.seed_from_spatial is not None:
        seed(output, args.seed_from_spatial.resolve(), args.methodology_change)
    prepare(output)


if __name__ == "__main__":
    main()
