#!/usr/bin/env python3
"""Prepare Figure 3 solely from the immutable original nine-animal cohort."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import fig3_keyence_iqr as iqr
import fig3_keyence_permanova as statistics
import shape_spatial_analysis as shape

PAPER = Path(__file__).resolve().parents[2]
DEFAULT = PAPER / "analyses/fig3_original_cohort_20260920"
RECOVERED = PAPER / "thesis_and_manuscript/Scientific_Reports_Segura_et_al_2026/spatial_reanalysis_20260919"
ORIGINAL_SPATIAL = PAPER / "analyses/spatial_shape_20260920"
CONDITIONS = ("Water", "Sucrose", "Allulose")
CONTRAST_PLAN = dict(schema="fig3-original-pairwise-spatial-v1", experimental_unit="biological animal",
    comparisons=["Allulose - Water", "Allulose - Sucrose"], endpoints=["cfos", "marker_cfos"],
    metrics=["core_enrichment", "cluster_excess"], statistic="two-sided absolute difference of animal means",
    permutations="all distinct pairwise animal-label allocations, conditional on fixed eligible animals",
    multiplicity="Holm separately across four Figure 3 comparisons per metric (two endpoints by two comparisons)",
    other_figures_contribute=False,
    graph_panel_annotation="Observed clustering fraction is shown, but its comparison uses observed-minus-expected clustering excess",
    selection_caveat=iqr.CAVEAT)


def copy(source: Path, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def seed(output: Path, legacy: Path):
    immutable = output / "original_inputs"
    manifest_path = immutable / "seed_manifest.json"
    if manifest_path.exists():
        raise ValueError("Immutable original inputs already exist; seed only once")
    immutable.mkdir(parents=True, exist_ok=True)
    roster = pd.read_csv(legacy / "animal_roster.csv")
    roster = roster.loc[roster.platform.ne("Keyence")].copy()
    if roster.groupby("condition").size().to_dict() != {"Water": 3, "Sucrose": 3, "Allulose": 3}:
        raise ValueError("Seeding requires the nine original animals")
    ids = set(roster.animal_id)
    roster.to_csv(immutable / "animal_roster.csv", index=False)
    history = []
    for name in ("animal_plot_values.csv", "animal_metrics.csv", "unit_metrics.csv", "stratum_metrics.csv", "cluster_membership.csv.gz"):
        source = legacy / name
        frame = pd.read_csv(source, low_memory=False)
        identity = "animal_id" if "animal_id" in frame else "unit_id"
        frame = frame.loc[frame[identity].isin(ids)].copy()
        remove = [column for column in frame if column.startswith("iqr_") or column.endswith("_iqr_excluded") or column.endswith("_include_in_main") or column == "include_in_main"]
        frame = frame.drop(columns=remove)
        frame.to_csv(immutable / name, index=False)
        history.append(dict(source=str(source), source_sha256=statistics.sha256(source), selected_animals=sorted(ids), output=name))
    # Independently confirm the spatial features against the earlier original-cohort analysis.
    recovered_metrics = pd.read_csv(ORIGINAL_SPATIAL / "animal_metrics.csv")
    recovered_metrics = recovered_metrics.loc[recovered_metrics.cohort.eq("NPY") & recovered_metrics.animal_id.isin(ids)]
    checked = pd.read_csv(immutable / "animal_metrics.csv").set_index(["endpoint", "animal_id"])
    reference = recovered_metrics.set_index(["endpoint", "animal_id"])
    for feature in shape.METRICS:
        if not np.allclose(checked[feature], reference.loc[checked.index, feature], equal_nan=True, rtol=0, atol=1e-12):
            raise ValueError("Original spatial feature differs from recovered cohort: " + feature)
    runtime = immutable / "runtime"
    cells = pd.read_csv(RECOVERED / "data/NPY_cells.csv.gz", low_memory=False)
    cells = cells.loc[cells.animal_id.isin(ids)].copy()
    (runtime / "inputs").mkdir(parents=True, exist_ok=True)
    cells.to_csv(runtime / "inputs/NPY_cells.csv.gz", index=False)
    acquisition = pd.read_csv(RECOVERED / "audit/NPY_acquisition_manifest.csv")
    acquisition = acquisition.loc[acquisition.animal_id.isin(ids)].copy()
    if len(acquisition) != 9 or acquisition.animal_id.nunique() != 9:
        raise ValueError("Original cohort must have nine accepted acquisitions")
    acquisition = acquisition.merge(roster[["animal_id", "platform", "biological_cohort"]], on="animal_id", validate="one_to_one")
    original_geometry = []
    for index, row in acquisition.iterrows():
        source = Path(row.geometry_path)
        relative = Path("geometry") / (row.acquisition_id + "_geometry.npz")
        copy(source, runtime / relative)
        original_geometry.append(dict(acquisition_id=row.acquisition_id, source=str(source), source_sha256=statistics.sha256(source), local_path=str(relative)))
        acquisition.loc[index, "geometry_path"] = str(output / relative)
    acquisition.to_csv(runtime / "inputs/NPY_acquisition_manifest.csv", index=False)
    contours = pd.read_csv(ORIGINAL_SPATIAL / "contour_manifest.csv")
    contours = contours.loc[contours.cohort.eq("NPY") & contours.animal_id.isin(ids)].copy()
    contours.to_csv(runtime / "contour_manifest.csv", index=False)
    for row in contours.itertuples():
        if isinstance(row.contour_path, str) and row.contour_path:
            copy(ORIGINAL_SPATIAL / row.contour_path, runtime / row.contour_path)
    sensitivity = pd.read_csv(ORIGINAL_SPATIAL / "bandwidth_sensitivity_condition_summary.csv")
    sensitivity.loc[sensitivity.cohort.eq("NPY")].to_csv(runtime / "bandwidth_sensitivity_condition_summary.csv", index=False)
    copy(ORIGINAL_SPATIAL / "example_npy.npz", runtime / "example_npy.npz")
    source_provenance = dict(schema="fig3-original-input-provenance-v1", original_animal_ids=sorted(ids),
        conditions={condition: 3 for condition in CONDITIONS}, acquisition_count=9,
        cell_ledger_source=dict(path=str(RECOVERED / "data/NPY_cells.csv.gz"), sha256=statistics.sha256(RECOVERED / "data/NPY_cells.csv.gz")),
        acquisition_manifest_source=dict(path=str(RECOVERED / "audit/NPY_acquisition_manifest.csv"), sha256=statistics.sha256(RECOVERED / "audit/NPY_acquisition_manifest.csv")),
        original_geometry=original_geometry, table_seed_history=history,
        original_spatial_values_independently_matched=True,
        accepted_associations="accepted_legacy nuclear marker identities; original Figure 3 masks unchanged",
        runtime_policy="Only immutable original-cohort copies are read during preparation; no other experiment is extracted or analyzed")
    statistics.json_write(immutable / "source_provenance.json", source_provenance)
    files = {str(path.relative_to(immutable)): statistics.sha256(path) for path in sorted(immutable.rglob("*")) if path.is_file()}
    statistics.json_write(manifest_path, dict(schema="fig3-original-nine-immutable-inputs-v1", animal_ids=sorted(ids),
        animal_count=9, condition_counts={condition: 3 for condition in CONDITIONS}, files=files))
    # Publication guards are independent of the fixed scientific inputs.
    for name in ("before_change_hashes.json", "figure4_before_hashes.json"):
        copy(legacy / name, output / name)
    archive = output / "archive/Fig3/source_data"
    for name in ("Figure_3_cFos_NPY_spatial_animal_counts.csv", "Figure_3_cFos_NPY_animal_plot_values.csv", "Figure_3_cFos_NPY_spatial_segmentation_model_setting_audit.csv"):
        source = legacy / "archive/Fig3/source_data" / name
        if source.is_file():
            copy(source, archive / name)


def verify_inputs(output: Path):
    immutable = output / "original_inputs"
    path = immutable / "seed_manifest.json"
    manifest = json.loads(path.read_text())
    if manifest.get("schema") != "fig3-original-nine-immutable-inputs-v1" or manifest.get("animal_count") != 9:
        raise ValueError("Invalid original-cohort immutable input manifest")
    for name, expected in manifest["files"].items():
        if statistics.sha256(immutable / name) != expected:
            raise ValueError("Immutable original input changed: " + name)
    roster = pd.read_csv(immutable / "animal_roster.csv")
    if roster.platform.eq("Keyence").any() or roster.animal_id.str.contains("Keyence", case=False).any():
        raise ValueError("Different-experiment animals entered the original-cohort inputs")
    if set(roster.animal_id) != set(manifest["animal_ids"]):
        raise ValueError("Immutable roster does not match seed animal identities")
    return manifest


def pairwise_contrasts(output: Path):
    statistics.json_write(output / "treatment_contrast_plan.json", CONTRAST_PLAN)
    units = pd.read_csv(output / "unit_metrics.csv")
    units = units.loc[units.include_in_main].copy()
    rows = []
    for metric in ("core_enrichment", "cluster_excess"):
        for endpoint in ("cfos", "marker_cfos"):
            selected = units.loc[units.endpoint.eq(endpoint)]
            for comparison in ("Water", "Sucrose"):
                result = shape.exact_comparison(selected.loc[selected.condition.eq("Allulose"), metric],
                                                selected.loc[selected.condition.eq(comparison), metric])
                result.update(cohort="NPY", endpoint=endpoint, region="FIELD", metric=metric,
                    group_a="Allulose", group_b=comparison, comparison="Allulose - " + comparison,
                    family=metric, planned_family_size=4, unit_type="animal",
                    interpretation="Exploratory original Figure 3 cohort only; " + iqr.CAVEAT)
                rows.append(result)
    table = pd.DataFrame(rows)
    for _, indices in table.groupby("family").groups.items():
        if len(indices) != 4:
            raise ValueError("Each Figure 3 metric must have four planned comparisons")
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
    if saved_readme is not None:
        (output / "README.md").write_text(saved_readme)
    else:
        (output / "README.md").write_text("""# Figure 3: original nine-animal cohort

This bundle uses only the original Figure 3 experiment: three Water, three Sucrose and three Allulose animals. No animals from another experiment enter the data, plots, inference, or multiplicity families. Fixed original inputs and their provenance are in original_inputs/; every rebuild verifies their SHA256 hashes.

The specified within-condition 1.5-IQR rule is applied separately by endpoint after selecting the original cohort. No original animals are flagged. WATER_NPY3 has zero detected NPY/c-FOS double-positive nuclei and remains in abundance counts; its conditional NPY spatial features are unavailable, yielding two of three estimable Water animals.

The density contour and spatial cluster definitions are unchanged. Exact pairwise differences of animal means provide panel annotations, with Holm adjustment across four Figure 3 comparisons per metric. Clustering comparisons use observed-minus-expected excess. PERMANOVA and its dispersion diagnostic remain in source tables and legends. All inference is exploratory; the within-cohort Allulose-Sucrose comparison does not resolve the Water acquisition confounding. Any IQR selection is held fixed during permutation calibration.

Rebuild and publish from Paper root: ./scripts/run_fig3_original_cohort.sh
""")
    metadata = dict(schema="fig3-original-cohort-analysis-v1", original_animal_ids=manifest["animal_ids"],
        condition_counts={condition: 3 for condition in CONDITIONS}, independent_of_other_experiments=True,
        immutable_input_manifest=dict(path=str(immutable / "seed_manifest.json"), sha256=statistics.sha256(immutable / "seed_manifest.json")),
        input_files=[dict(path=str(immutable / name), sha256=value) for name, value in manifest["files"].items()],
        analysis_script_sha256=statistics.sha256(Path(__file__)), iqr_script_sha256=statistics.sha256(Path(iqr.__file__)),
        statistics_script_sha256=statistics.sha256(Path(statistics.__file__)),
        parameters=dict(k=6, edge_cap_um=75., bandwidth_um=40., grid_um=5., contour_mass=.8, minimum_marker_anchors=5),
        interpretation="Original experiment only; abundance and spatial IQR rules applied per endpoint; Water acquisition confounding remains",
        spatial_recomputation="Fixed, previously validated original-cohort animal/stratum metrics; current tests and summaries rebuilt from immutable tables",
        native_masks_modified=False)
    statistics.json_write(output / "analysis_metadata.json", metadata)
    validation = dict(status="PASS", animals=9, conditions={condition: 3 for condition in CONDITIONS},
        original_nine_abundance_and_spatial_values_unchanged=True, original_cohort_only=True,
        no_other_experiment_rows=True, immutable_input_hashes_verified=True, acquisitions=9,
        missing_conditional_features_not_imputed=True, original_images_and_masks_modified=False,
        iqr_outlier_count=len(report["exclusions"]), pairwise_figure3_only_family_size=4,
        no_figure4_data_in_statistics=True)
    statistics.json_write(output / "analysis_validation.json", validation)
    verify_inputs(output)
    print(json.dumps(validation, indent=2))
    return validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT)
    parser.add_argument("--seed-from", type=Path, help="One-time read-only source for the validated original-animal subset")
    parser.add_argument("--snapshot-figure4", type=Path, help="Write the current Figure 4 publication guard and return without recomputation")
    args = parser.parse_args()
    if args.snapshot_figure4 is not None:
        hashes = {name: statistics.sha256(PAPER / name) for name in ("Fig4/Figure_4.pdf", "Fig4/Figure_4.png")}
        args.snapshot_figure4.parent.mkdir(parents=True, exist_ok=True)
        statistics.json_write(args.snapshot_figure4, hashes)
        print(json.dumps(dict(status="PASS", figure4_snapshot=str(args.snapshot_figure4)), indent=2))
        return
    output = args.analysis_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.seed_from is not None:
        seed(output, args.seed_from.resolve())
    prepare(output)


if __name__ == "__main__":
    main()
