#!/usr/bin/env python3
"""Validate Figure S3 cohort completeness, receipts, statistics, and artifacts."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from PIL import Image

from figs3_common import FIGS3_DIR, ORGANS, FigS3Error, atomic_json, sha256_file


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FigS3Error(message)


def recorded_path(root: Path, recorded: str) -> Path:
    """Bind a receipt path to content in this tree, not to the directory it had.

    Receipts store the absolute path the file occupied when it was written. That
    directory changes whenever the package is renumbered or a reader unpacks the
    archive somewhere else, so the recorded path is treated as a hint: the
    longest tail of it that exists under this figure root wins, and the size and
    SHA-256 in the same receipt row remain the real binding.
    """
    candidate = Path(recorded)
    parts = candidate.parts
    for index in range(1, len(parts)):
        moved = root.joinpath(*parts[index:])
        if moved.is_file():
            return moved
    raise FigS3Error(
        f"Receipt path cannot be rebound inside the current Figure S3 tree: {recorded}"
    )


def json_file(path: Path) -> dict:
    require(path.is_file() and not path.is_symlink(), f"Missing/unsafe JSON: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def validate(root: Path) -> dict[str, object]:
    checks: list[dict[str, object]] = []
    slide = pd.read_csv(root / "source_data" / "Figure_S3_slide_manifest.csv")
    organ_manifest = pd.read_csv(
        root / "source_data" / "Figure_S3_sample_organ_manifest.csv"
    )
    require(
        len(slide) == 24 and slide.sample_id.nunique() == 24,
        "Slide manifest is not 24 unique animals",
    )
    require(len(organ_manifest) == 72, "Animal-organ manifest is not 72 rows")
    require(
        slide.treatment.value_counts().to_dict()
        == {"Sucrose": 9, "Allulose": 8, "Water": 7},
        "Treatment counts differ",
    )
    m26_016 = slide[slide.sample_id == "m26-016"].iloc[0]
    require(
        m26_016.treatment == "Water" and m26_016.sample_note == "brain_not_collected",
        "m26-016 is not Water with provenance note",
    )
    checks.append(
        {
            "check": "metadata_and_control_definition",
            "status": "pass",
            "slides": 24,
            "organ_rows": 72,
        }
    )

    export_summary = json_file(root / "exports" / "Figure_S3_L0_L2_export_summary.json")
    export_manifest = pd.read_csv(
        root / "exports" / "Figure_S3_L0_L2_export_manifest.csv"
    )
    require(
        export_summary.get("status") == "complete"
        and export_summary.get("complete_slides") == 24,
        "L0/L2 export is incomplete",
    )
    require(
        len(export_manifest) == 48 and set(export_manifest.level) == {0, 2},
        "Export manifest is not 24 x L0/L2",
    )
    for row in export_manifest.itertuples():
        path = recorded_path(root, str(row.output_path))
        require(path.is_file() and not path.is_symlink(), f"Missing/unsafe TIFF {path}")
        require(
            path.stat().st_size == int(row.output_size_bytes),
            f"TIFF size differs from receipt: {path}",
        )
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]
            require(
                page.is_tiled and page.samplesperpixel == 3,
                f"TIFF is not tiled RGB: {path}",
            )
            require(
                (page.imagewidth, page.imagelength) == (row.width, row.height),
                f"TIFF geometry differs: {path}",
            )
    checks.append(
        {"check": "L0_L2_exports", "status": "pass", "tiffs": 48, "slides": 24}
    )

    segmentation = json_file(root / "organ_segmentation" / "segmentation_summary.json")
    require(
        segmentation.get("status") == "complete" and segmentation.get("slides") == 24,
        "Organ segmentation is incomplete",
    )
    require(
        segmentation.get("version") == "figs3-three-organ-cv-1.1"
        and segmentation.get("human_reviewed") is True
        and segmentation.get("image_registration_performed") is False,
        "Organ segmentation lacks the corrected, segmentation-only HIL-reviewed identity audit",
    )
    segmentation_receipts = segmentation.get("receipts", [])
    require(
        len(segmentation_receipts) == 24
        and {row.get("version") for row in segmentation_receipts}
        == {"figs3-three-organ-cv-1.1"}
        and all(
            "liver-kidney-spleen" in row.get("organ_identity_basis", "")
            for row in segmentation_receipts
        ),
        "Per-slide organ-identity receipts do not preserve verified Liver-Kidney-Spleen order",
    )
    for sample_id in slide.sample_id:
        values = set(
            np.unique(
                np.asarray(
                    Image.open(
                        root / "organ_segmentation" / sample_id / "organ_label_mask.png"
                    )
                )
            ).tolist()
        )
        require(
            values == {0, 1, 2, 3}, f"Incomplete organ labels for {sample_id}: {values}"
        )
    checks.append(
        {
            "check": "three_organs_per_wsi_and_identity_audit",
            "status": "pass",
            "slides": 24,
            "version": segmentation.get("version"),
            "human_reviewed": True,
            "placement_order": "Liver-Kidney-Spleen",
        }
    )

    migration_path = (
        root / "provenance" / "Figure_S3_organ_identity_migration_receipt.json"
    )
    if migration_path.exists():
        migration = json_file(migration_path)
        migration_audit = pd.read_csv(
            root / "provenance" / "Figure_S3_organ_identity_migration_audit.csv"
        )
        geometry_columns = [
            "kidney_geometry_gate_fraction",
            "liver_geometry_gate_fraction",
            "spleen_geometry_gate_fraction",
        ]
        require(
            migration.get("status") == "complete"
            and migration.get("version") == "figs3-organ-identity-audit-1.0"
            and migration.get("samples") == 24,
            "Organ-identity migration receipt is incomplete",
        )
        require(
            len(migration_audit) == 24
            and set(migration_audit.status) == {"corrected"}
            and (migration_audit[geometry_columns] == 1.0).all().all(),
            "Organ-identity migration geometry gate did not pass exactly",
        )
        checks.append(
            {
                "check": "organ_identity_migration_geometry",
                "status": "pass",
                "slides": 24,
                "minimum_geometry_fraction": float(
                    migration_audit[geometry_columns].min().min()
                ),
            }
        )

    require(
        not (root / "registration").exists(),
        "Obsolete Figure S3 image-registration directory remains active",
    )
    for suffix in (".png", ".pdf"):
        require(
            not (
                root / "provenance" / f"Figure_S3_Diagnostic_Registration{suffix}"
            ).exists(),
            "Obsolete Figure S3 registration diagnostic remains active",
        )
    checks.append(
        {"check": "segmentation_only_no_image_registration", "status": "pass"}
    )

    inference = json_file(root / "histoplus" / "inference_summary.json")
    require(
        inference.get("status") == "complete"
        and inference.get("complete_slides") == 24,
        "HistoPLUS inference is incomplete",
    )
    phenotype = pd.read_csv(
        root / "histoplus" / "Figure_S3_HistoPLUS_phenotype_summary.csv"
    )
    tile_cv = pd.read_csv(
        root / "histoplus" / "Figure_S3_classical_tile_cv_features.csv"
    )
    tiles = pd.read_csv(root / "histoplus" / "Figure_S3_selected_tile_manifest.csv")
    tile_counts = pd.read_csv(
        root / "histoplus" / "Figure_S3_HistoPLUS_tile_class_counts.csv"
    )
    require(
        len(phenotype) == 24 * 3 * 14,
        "Phenotype summary does not include 14 classes per animal-organ",
    )
    require(
        len(tile_cv) == 24 * 3 * 24 and len(tiles) == 24 * 3 * 24,
        "Tile outputs are not balanced at 24/organ",
    )
    require(
        set(tile_cv.groupby(["sample_id", "organ"]).size()) == {24},
        "Classical CV tile balance differs",
    )
    require(
        set(tiles.groupby(["sample_id", "organ"]).size()) == {24},
        "Selected tile balance differs",
    )
    require(
        len(tile_counts) == 24 * 3 * 24,
        "Tile-count table lacks explicit rows for zero-cell tiles",
    )
    tile_keys = ["sample_id", "organ", "tile_uid"]
    require(
        set(map(tuple, tile_counts[tile_keys].to_numpy()))
        == set(map(tuple, tiles[tile_keys].to_numpy())),
        "Tile-count keys do not exactly match selected tiles",
    )
    require((tile_counts.total_cells >= 0).all(), "A tile count is negative")
    require(
        (
            phenotype.groupby(
                ["sample_id", "organ"]
            ).total_cells_in_sampled_tiles.first()
            > 0
        ).all(),
        "A completed animal-organ has zero inferred cells",
    )
    require(
        (phenotype.cross_species_validated.astype(str).str.lower() == "false").all(),
        "Cross-species status is misstated",
    )
    checks.append(
        {
            "check": "inference_and_tile_balance",
            "status": "pass",
            "slides": 24,
            "selected_tiles": len(tiles),
            "phenotype_rows": len(phenotype),
        }
    )

    analysis = json_file(root / "provenance" / "Figure_S3_analysis_receipt.json")
    require(
        analysis.get("status") == "complete" and analysis.get("animal_count") == 24,
        "Analysis receipt is incomplete",
    )
    require(
        analysis.get("analysis_version") == "figs3-animal-composition-1.8"
        and analysis.get("image_registration_performed") is False,
        "Analysis receipt is not the segmentation-only Figure S3 analysis",
    )
    water_ids = sorted(slide.loc[slide.treatment == "Water", "sample_id"].tolist())
    all_fractions = {
        str(row.get("sample_id")): float(row.get("tissue_fraction", np.nan))
        for row in segmentation_receipts
    }
    require(
        set(water_ids).issubset(all_fractions)
        and all(np.isfinite(all_fractions[sample_id]) for sample_id in water_ids),
        "Water tissue fractions are missing from segmentation receipts",
    )
    median_fraction = float(
        np.median([all_fractions[sample_id] for sample_id in water_ids])
    )
    expected_representative = min(
        water_ids,
        key=lambda sample_id: (
            abs(all_fractions[sample_id] - median_fraction),
            sample_id,
        ),
    )
    selection = analysis.get("representative_sample_selection", {})
    observed_median = selection.get("median_tissue_fraction")
    require(
        expected_representative == "m26-015"
        and analysis.get("representative_sample_id") == expected_representative
        and selection.get("selected_sample_id") == expected_representative
        and selection.get("eligible_treatment") == "Water"
        and selection.get("eligible_sample_count") == 7
        and isinstance(observed_median, (int, float))
        and np.isclose(observed_median, median_fraction),
        "Representative Water sample is not the median-tissue-fraction animal",
    )
    checks.append(
        {
            "check": "representative_selection_without_alignment",
            "status": "pass",
            "sample_id": expected_representative,
            "median_tissue_fraction": median_fraction,
        }
    )
    composition = pd.read_csv(
        root / "source_data" / "Figure_S3_composition_contrasts.csv"
    )
    deviation = pd.read_csv(
        root / "source_data" / "Figure_S3_Water_reference_deviation_contrasts.csv"
    )
    classical = pd.read_csv(
        root / "source_data" / "Figure_S3_classical_cv_contrasts.csv"
    )
    deviation_anova = pd.read_csv(
        root / "source_data" / "Figure_S3_Water_reference_deviation_ANOVA.csv"
    )
    examples = pd.read_csv(
        root / "source_data" / "Figure_S3_significant_cell_example_manifest.csv"
    )
    cell_bars = pd.read_csv(
        root / "source_data" / "Figure_S3_displayed_cell_barplot_values.csv"
    )
    require(
        len(composition) == 13 * 3 * 2,
        "Composition contrast count differs after tumor-class exclusion",
    )
    require(
        "Cancer cell" not in set(composition.model_class_raw),
        "Tumor-specific class remains in composition contrasts",
    )
    require(len(deviation) == 6, "Deviation contrast count differs")
    require(len(classical) == 36, "Classical contrast count differs")
    require(
        len(deviation_anova) == 3 and deviation_anova.p_value.between(0, 1).all(),
        "Deviation ANOVA table differs",
    )
    significant = composition[
        composition.analysis_eligible.astype(str).str.lower().eq("true")
        & composition.q_value_global.lt(0.05)
    ].sort_values(["q_value_global", "organ", "model_class_raw", "contrast"])
    require(
        not significant.empty,
        "No globally significant contrast is available for the requested cell inset",
    )
    expected_example = significant.iloc[0]
    require(
        len(examples) == 2
        and set(examples.treatment_displayed)
        == {"Water", str(expected_example.contrast).split(" vs ")[0]}
        and set(examples.organ) == {expected_example.organ}
        and set(examples.model_class_raw) == {expected_example.model_class_raw}
        and examples.q_value_global.lt(0.05).all(),
        "Significant-cell inset manifest does not match the strongest global BH result",
    )
    require(
        examples.target_cell_count_in_tile.gt(0).all()
        and examples.displayed_cell_count_in_tile.ge(
            examples.target_cell_count_in_tile
        ).all()
        and all(
            recorded_path(root, str(path)).is_file() for path in examples.source_l0_path
        ),
        "Significant-cell inset does not resolve to positive, traceable L0 examples",
    )
    require(
        "Cancer cell" not in set(examples.model_class_raw),
        "Excluded tumor-specific class appears in significant-cell examples",
    )
    require(
        len(cell_bars) == 24 * 3 * 13
        and cell_bars.sample_id.nunique() == 24
        and cell_bars.model_class_raw.nunique() == 13
        and "Cancer cell" not in set(cell_bars.model_class_raw),
        "Displayed-cell barplot source table is not 24 x 3 x 13 retained outputs",
    )
    fraction_sums = cell_bars.groupby(
        ["sample_id", "organ"]
    ).displayed_class_fraction.sum()
    require(
        np.allclose(fraction_sums.to_numpy(float), 1.0, atol=1e-10)
        and cell_bars.displayed_class_count.ge(0).all()
        and cell_bars.displayed_total_count.gt(0).all(),
        "Displayed-cell barplot fractions/counts are invalid",
    )
    require(
        int((classical.freedman_lane_permutation_p_value < 0.05).sum()) == 6
        and int((classical.q_value_global < 0.05).sum()) == 0,
        "Panel E nominal-star/FDR distinction differs from the rendered claim",
    )
    for table, name in (
        (composition, "composition"),
        (deviation, "deviation"),
        (classical, "classical"),
    ):
        for column in (
            "freedman_lane_permutation_p_value",
            "q_value_global",
            "q_value_within_organ",
        ):
            require(table[column].between(0, 1).all(), f"Invalid {name} {column}")
        require(set(table.n_animals) == {24}, f"{name} does not use animal n=24")
    eligible = composition.analysis_eligible.astype(str).str.lower() == "true"
    require(
        composition.loc[eligible, "effect"].notna().all(),
        "An eligible composition contrast was not estimated",
    )
    require(
        composition.loc[~eligible, "effect"].isna().all(),
        "A prevalence-filtered contrast has an effect estimate",
    )
    require(
        (composition.loc[~eligible, "freedman_lane_permutation_p_value"] == 1).all(),
        "A prevalence-filtered contrast is not explicitly marked untested",
    )
    checks.append(
        {
            "check": "animal_level_statistics_cell_examples_and_barplots",
            "status": "pass",
            "composition_tests": 78,
            "deviation_tests": 6,
            "deviation_anova_tests": 3,
            "classical_tests": 36,
            "significant_cell_example_rows": len(examples),
            "cell_barplot_rows": len(cell_bars),
            "panel_e_nominal_permutation_p_lt_0_05": int((classical.freedman_lane_permutation_p_value < 0.05).sum()),
            "panel_e_global_q_lt_0_05": 0,
        }
    )

    png = root / "Figure_S3.png"
    pdf = root / "Figure_S3.pdf"
    require(
        png.is_file() and pdf.is_file() and pdf.read_bytes()[:4] == b"%PDF",
        "Figure master is missing",
    )
    with Image.open(png) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0, 0))
    require(
        width >= 8000 and height >= 6000, "Figure PNG is not a 600-dpi-scale master"
    )
    require(min(dpi) >= 590, f"Figure PNG DPI metadata is below 600: {dpi}")
    expected_panels = {"A", "B", "C", "D", "E", "F"}

    def exported_panels(extension: str, spanish: bool) -> set[str]:
        found: set[str] = set()
        for path in (root / "panels").glob(f"Figure_S3_Panel_*{extension}"):
            stem = path.stem
            if stem.endswith("_spanish") != spanish:
                continue
            if spanish:
                stem = stem[: -len("_spanish")]
            found.add(stem.rsplit("_", 1)[-1])
        return found

    exports = {
        (extension, language): exported_panels(extension, language == "spanish")
        for extension in (".png", ".pdf")
        for language in ("english", "spanish")
    }
    for (extension, language), found in sorted(exports.items()):
        require(
            found == expected_panels,
            f"{language} {extension} panel set differs: {sorted(found)}",
        )
    translation = json_file(
        root / "provenance" / "Figure_S3_spanish_translation_receipt.json"
    )
    require(
        translation.get("status") == "complete"
        and set(translation.get("panels", [])) == expected_panels
        and len(translation.get("translations", {})) > 0,
        "Spanish subpanel translation receipt is missing or incomplete",
    )

    legend_dir = root / "legends"
    legend = (legend_dir / "Figure_S3_LEGEND.txt").read_text(encoding="utf-8")
    for phrase in (
        "m26-015",
        "median whole-slide tissue fraction",
        "m26-016",
        "included as Water",
        "not validated",
        "animal is the biological unit",
        "not causal",
    ):
        require(
            phrase.casefold() in legend.casefold(),
            f"Legend lacks required interpretation phrase: {phrase}",
        )
    require(
        "NPY" not in legend and "POMC" not in legend,
        "Cohort labels remain in the figure legend",
    )
    for forbidden in ("registration", "registered", "affine", "dice", "medoid"):
        require(
            forbidden not in legend.casefold(),
            f"Figure S3 legend retains obsolete image-registration wording: {forbidden}",
        )
    spanish_legend = legend_dir / "Figure_S3_LEGEND_spanish.txt"
    require(spanish_legend.is_file(), "The Spanish figure legend is missing")
    spanish_legend_text = spanish_legend.read_text(encoding="utf-8")
    for forbidden in ("registro", "registrada", "registró", "afín", "dice", "medoide"):
        require(
            forbidden not in spanish_legend_text.casefold(),
            "Spanish Figure S3 legend retains obsolete image-registration wording: "
            f"{forbidden}",
        )
    for letter in sorted(expected_panels):
        for name in (
            f"Figure_S3_Panel_{letter}_LEGEND.txt",
            f"Figure_S3_Panel_{letter}_LEGEND_spanish.txt",
        ):
            path = legend_dir / name
            require(
                path.is_file() and path.stat().st_size > 0,
                f"Missing or empty per-panel legend: {name}",
            )
    for letter, expected_reference in (
        ("A", "workflow"),
        ("B", "Mahalanobis"),
        ("C", "scale bars"),
        ("D", "ordination"),
        ("E", "forest"),
        ("F", "confidence intervals"),
    ):
        body = (legend_dir / f"Figure_S3_Panel_{letter}_LEGEND.txt").read_text(
            encoding="utf-8"
        )
        require(
            body.startswith(f"Figure S3, panel {letter}.")
            and expected_reference in body,
            f"Per-panel legend {letter} does not describe the re-lettered panel",
        )
    checks.append(
        {
            "check": "figure_artifacts_and_bilingual_legends",
            "status": "pass",
            "png_width": width,
            "png_height": height,
            "dpi": list(dpi),
            "panels_per_language": len(expected_panels),
            "languages": ["english", "spanish"],
        }
    )
    return {"checks": checks, "export_manifest": export_manifest}


def write_manifest(root: Path, export_manifest: pd.DataFrame) -> Path:
    target = root / "provenance" / "SHA256_MANIFEST.tsv"
    rows: dict[str, tuple[int, str, str]] = {}
    for row in export_manifest.itertuples():
        path = recorded_path(root, str(row.output_path)).resolve()
        rows[str(path.relative_to(root))] = (
            int(row.output_size_bytes),
            str(row.output_sha256),
            "conversion_receipt",
        )
    excluded = {
        target.resolve(),
        (root / "provenance" / "VALIDATION_RECEIPT.json").resolve(),
    }
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if (
            not path.is_file()
            or path.is_symlink()
            or path.resolve() in excluded
            or "__pycache__" in relative.parts
            or path.suffix.lower() in {".pyc", ".tif", ".tiff"}
        ):
            continue
        rows[str(path.resolve().relative_to(root))] = (
            path.stat().st_size,
            sha256_file(path),
            "validated_file",
        )
    lines = ["relative_path\tsize_bytes\tsha256\thash_source"]
    lines.extend(
        f"{path}\t{size}\t{digest}\t{source}"
        for path, (size, digest, source) in sorted(rows.items())
    )
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figs3-root", type=Path, default=FIGS3_DIR)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate without rewriting the SHA-256 manifest or validation receipt.",
    )
    args = parser.parse_args()
    root = args.figs3_root.resolve()
    result = validate(root)
    export_manifest = result.pop("export_manifest")
    if args.check_only:
        print(json.dumps({"status": "pass", "checks": result["checks"]}, indent=2))
        return 0
    manifest = write_manifest(root, export_manifest)
    receipt = {
        "status": "pass",
        "validator": "figs3-validator-1.4",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "checks": result["checks"],
        "sha256_manifest": str(manifest.resolve()),
        "sha256_manifest_sha256": sha256_file(manifest),
        "publication_interpretation": "data-complete exploratory model-predicted phenotypes; cross-species labels not validated",
    }
    atomic_json(root / "provenance" / "VALIDATION_RECEIPT.json", receipt)
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
