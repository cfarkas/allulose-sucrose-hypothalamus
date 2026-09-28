#!/usr/bin/env python3
"""Fail-closed readiness audit for the independent Figure 2 raw-data bundle.

This validator never edits raw payload files.  It reconciles the five completed
HIL runs, verifies the copied-file hash receipt, maps every top-level Cellpose
mask to its microscopy plane/acquisition, and checks native-image nucleus
coordinates against the saved HIL-defined region masks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


HERE = Path(__file__).resolve()
PAPER = next(path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
FIG2 = PAPER / "Fig2"
RAW_DATA = FIG2 / "raw_data"
APOTOME = RAW_DATA / "Apotome"
# Acquisition-time prefix recorded inside discovered_images.csv. This is a
# fixed property of the archived CSVs, not of the machine: it is only used to
# rebase those recorded paths onto raw_data/Apotome, so it stays correct on
# any host and in any relocated copy of the tree.
SOURCE_APOTOME = Path("/media/server/STORAGE/Apotome")
RAW_MANIFEST = RAW_DATA / "RAW_DATA_MANIFEST.csv"
# This active maintenance helper verifies the accepted HIL bundle without
# making region decisions; the analysis engine remains at the front of Fig2.
CANONICAL_HIL = FIG2 / "09_rebuild_hil_inputs.py"
CANONICAL_ANALYSIS = FIG2 / "01_analyze_global_cfos.py"
CANONICAL_COMBINED_SHA256 = (
    "d80bfe4dd5b259f1d309dcf59f54cb0effd62d7911f72de769cfc2b9b8170ebf"
)


RUNS = (
    {
        "run_id": "june_2025_tiled",
        "dataset": APOTOME / "2025.06.24",
        "hil": APOTOME / "2025.06.24_human_regions_v9_2",
    },
    {
        "run_id": "june_2025_primary",
        "dataset": APOTOME / "24_06_2025" / "analysis",
        "hil": APOTOME / "24_06_2025" / "analysis" / "human_regions_v9",
    },
    {
        "run_id": "september_2025_water",
        "dataset": APOTOME / "29.09.2025",
        "hil": APOTOME / "29.09.2025_human_regions_v9_2",
    },
    {
        "run_id": "august_2025_horizontal",
        "dataset": APOTOME / "Aug_2025_NPY" / "selected" / "horizontal",
        "hil": APOTOME / "Aug_2025_NPY" / "selected_horizontal_human_regions_v9_2",
    },
    {
        "run_id": "august_2025_selected",
        "dataset": APOTOME / "Aug_2025_NPY" / "selected",
        "hil": APOTOME / "Aug_2025_NPY" / "selected_human_regions_v9_2",
    },
)


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load canonical module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    temporary.replace(path)


def atomic_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty receipt: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent,
        prefix=f".{path.name}.", delete=False,
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def require_manifest() -> tuple[pd.DataFrame, dict[str, pd.Series], dict[str, pd.Series]]:
    manifest = pd.read_csv(RAW_MANIFEST, dtype=str, keep_default_na=False)
    if len(manifest) != 560:
        raise RuntimeError(f"Raw manifest row count changed: {len(manifest)} != 560")
    sizes = pd.to_numeric(manifest["size_bytes"], errors="raise")
    if int(sizes.sum()) != 108_000_611_886:
        raise RuntimeError("Raw manifest payload byte total changed")
    if not manifest["validation_status"].eq("PASS").all():
        raise RuntimeError("Raw manifest contains a non-PASS row")
    if not manifest["distinct_inode"].str.casefold().eq("true").all():
        raise RuntimeError("Raw manifest contains a non-distinct source/destination inode")
    if not manifest["metadata_match"].str.casefold().eq("true").all():
        raise RuntimeError("Raw manifest contains a metadata mismatch")
    if not manifest["source_sha256"].eq(manifest["destination_sha256"]).all():
        raise RuntimeError("Raw manifest contains a source/destination hash mismatch")
    if "destination_relative_to_figure" not in manifest:
        raise RuntimeError("Raw manifest lacks portable destination-relative paths")
    if (
        manifest["source_path"].duplicated().any()
        or manifest["destination_path"].duplicated().any()
        or manifest["destination_relative_to_figure"].duplicated().any()
    ):
        raise RuntimeError("Raw manifest contains a duplicate source or destination path")
    figure_root = FIG2.resolve()
    resolved_destinations: list[str] = []
    for row in manifest.itertuples(index=False):
        relative = Path(str(row.destination_relative_to_figure))
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError(f"Unsafe manifest-relative destination: {relative}")
        rebased = FIG2 / relative
        resolved = rebased.resolve()
        if not resolved.is_relative_to(figure_root):
            raise RuntimeError(f"Manifest destination escapes Figure 2: {relative}")
        if not rebased.is_file() or rebased.is_symlink():
            raise RuntimeError(
                f"Manifest destination is missing, nonregular, or a link: {rebased}"
            )
        if rebased.stat().st_size != int(row.size_bytes):
            raise RuntimeError(
                f"Manifest destination size differs: {rebased}: "
                f"{rebased.stat().st_size} != {row.size_bytes}"
            )
        resolved_destinations.append(str(resolved))
    by_source = {str(row.source_path): row for _, row in manifest.iterrows()}
    by_destination = {
        destination: manifest.iloc[index]
        for index, destination in enumerate(resolved_destinations)
    }
    return manifest, by_source, by_destination


def manifest_record(
    path: Path, by_destination: dict[str, pd.Series], expected_type: str | None = None
) -> pd.Series:
    row = by_destination.get(str(path.resolve()))
    if row is None:
        raise RuntimeError(f"Copied path is absent from RAW_DATA_MANIFEST.csv: {path}")
    if expected_type is not None and row.file_type != expected_type:
        raise RuntimeError(f"Unexpected manifest type for {path}: {row.file_type}")
    if row.source_sha256 != row.destination_sha256 or row.validation_status != "PASS":
        raise RuntimeError(f"Hash receipt is not PASS for {path}")
    return row


def companion_for_ancillary(mask: Path) -> Path:
    base = mask.name[: -len("_seg.npy")]
    candidates: list[Path] = []
    for extension in (".tif", ".tiff"):
        candidate = mask.parent / f"{base}{extension}"
        if candidate.is_file():
            candidates.append(candidate)
    lower = base.casefold()
    if lower.endswith("_masks"):
        prefix = base[: -len("_masks")]
        for extension in (".tif", ".tiff"):
            candidate = mask.parent / f"{prefix}_dapi{extension}"
            if candidate.is_file():
                candidates.append(candidate)
    if lower.endswith("_overlay"):
        prefix = base[: -len("_overlay")]
        candidates.extend(sorted(mask.parent.glob(f"{prefix}*dapi.tif*")))
    if not candidates:
        normalized = re.sub(r"[-_ ]+", "-", base.casefold())
        for candidate in (*mask.parent.glob("*.tif"), *mask.parent.glob("*.tiff")):
            candidate_normalized = re.sub(
                r"[-_ ]+", "-", candidate.stem.casefold()
            )
            if candidate_normalized == normalized:
                candidates.append(candidate)
    candidates = sorted(set(path.resolve() for path in candidates))
    if len(candidates) != 1:
        raise RuntimeError(
            f"Ancillary Cellpose mask has {len(candidates)} candidate microscopy planes: "
            f"{mask}; {candidates}"
        )
    return candidates[0]


def channel_from_mask(mask: Path) -> str:
    match = re.search(
        r"(?:^|[-_ ])(dapi|fos|gfp|neun|masks|overlay)_seg\.npy$",
        mask.name,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"Cannot classify Cellpose mask channel: {mask}")
    return match.group(1).upper()


def validate_coordinates(run: dict[str, Path | str]) -> dict[str, int]:
    hil = Path(run["hil"])
    status = pd.read_csv(hil / "annotation_status.csv")
    nuclei = pd.read_csv(hil / "per_nucleus_human_region_assignments.csv")
    required = {"tile_key", "nucleus_label", "centroid_y", "centroid_x", "region"}
    missing = required - set(nuclei.columns)
    if missing:
        raise RuntimeError(f"Missing spatial columns in {hil}: {sorted(missing)}")
    if nuclei.duplicated(["tile_key", "nucleus_label"]).any():
        raise RuntimeError(f"Duplicate tile/nucleus labels in {hil}")
    x = pd.to_numeric(nuclei["centroid_x"], errors="coerce")
    y = pd.to_numeric(nuclei["centroid_y"], errors="coerce")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or (x < 0).any() or (y < 0).any():
        raise RuntimeError(f"Invalid native-image nucleus coordinates in {hil}")
    coordinate_rows = 0
    polygon_points = 0
    polygon_points_outside_native_bounds = 0
    annotated = status.loc[status["status"].eq("annotated")]
    for row in annotated.itertuples(index=False):
        tile = str(row.tile_key)
        mask = hil / "region_masks" / Path(str(row.region_mask)).name
        annotation = hil / "annotations" / Path(str(row.annotation_json)).name
        with tifffile.TiffFile(mask) as tif:
            shape = tuple(int(value) for value in tif.pages[0].shape)
        if len(shape) < 2:
            raise RuntimeError(f"Region mask is not two-dimensional: {mask}; shape={shape}")
        height, width = shape[-2], shape[-1]
        subset = nuclei.loc[nuclei["tile_key"].astype(str).eq(tile)]
        if subset.empty:
            raise RuntimeError(f"Annotated tile has no per-nucleus coordinates: {tile}")
        sx = pd.to_numeric(subset["centroid_x"], errors="raise")
        sy = pd.to_numeric(subset["centroid_y"], errors="raise")
        if (sx >= width).any() or (sy >= height).any():
            raise RuntimeError(
                f"Nucleus coordinate escapes region-mask bounds for {run['run_id']}/{tile}"
            )
        coordinate_rows += len(subset)
        payload = json.loads(annotation.read_text(encoding="utf-8"))
        for polygon in payload.get("polygons", []):
            for point in polygon.get("points", []):
                if not isinstance(point, list) or len(point) != 2:
                    raise RuntimeError(f"Malformed polygon point in {annotation}")
                px, py = float(point[0]), float(point[1])
                if not math.isfinite(px) or not math.isfinite(py):
                    raise RuntimeError(f"Nonfinite polygon point in {annotation}")
                if px < 0 or py < 0 or px > width or py > height:
                    polygon_points_outside_native_bounds += 1
                polygon_points += 1
    if coordinate_rows != len(nuclei):
        raise RuntimeError(
            f"Not every per-nucleus row belongs to an annotated tile in {hil}: "
            f"{coordinate_rows} != {len(nuclei)}"
        )
    return {
        "annotated_tiles": len(annotated),
        "nucleus_coordinate_rows": coordinate_rows,
        "polygon_points": polygon_points,
        "polygon_points_outside_native_bounds": polygon_points_outside_native_bounds,
    }


def targeted_registration_inventory() -> list[str]:
    pattern = re.compile(
        r"(^|[/_.-])(atlas|registration|registered|transform|affine|landmark|brainreg|quicknii)([/_.-]|$)",
        flags=re.IGNORECASE,
    )
    return [
        str(path.relative_to(RAW_DATA))
        for path in RAW_DATA.rglob("*")
        if path.is_file() and pattern.search(str(path.relative_to(RAW_DATA)))
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-receipts",
        action="store_true",
        help="Write CELLPOSE_INPUT_MAP.csv and readiness receipts under Fig2/raw_data.",
    )
    args = parser.parse_args()

    outputs = (
        RAW_DATA / "CELLPOSE_INPUT_MAP.csv",
        RAW_DATA / "SPATIAL_READINESS.json",
        RAW_DATA / "RECONSTRUCTION_READINESS.md",
    )
    if args.write_receipts:
        existing = [path for path in outputs if path.exists()]
        if existing:
            raise FileExistsError(f"Refusing to overwrite readiness receipt(s): {existing}")

    manifest, _, by_destination = require_manifest()
    hil_module = load_module("fig2_hil_rebuilder_for_readiness", CANONICAL_HIL)
    analysis_module = load_module("fig2_analysis_for_readiness", CANONICAL_ANALYSIS)

    staged_frames: list[pd.DataFrame] = []
    required_masks: dict[str, dict[str, object]] = {}
    hil_summary: list[dict[str, object]] = []
    coordinate_summary: list[dict[str, object]] = []
    original_specs = {run.run_id: run for run in hil_module.HIL_RUNS}
    for index, run in enumerate(RUNS, start=1):
        original = original_specs[str(run["run_id"])]
        relocated = hil_module.HilRun(
            original.run_id,
            str(Path(run["hil"]).resolve()),
            original.source_dataset,
            original.annotated,
            original.excluded,
            original.pending,
            original.per_animal_rows,
        )
        staged, _, record = hil_module.validate_run(index, relocated, Path("/"))
        staged_frames.append(staged)
        hil_summary.append(
            {
                "run_id": run["run_id"],
                **record["status_counts"],
                **record["quantitative_chain"],
            }
        )
        coordinate_summary.append(
            {"run_id": run["run_id"], **validate_coordinates(run)}
        )

        discovered = pd.read_csv(Path(run["hil"]) / "discovered_images.csv")
        expected_rows = original.annotated + original.excluded + original.pending
        if len(discovered) != expected_rows:
            raise RuntimeError(
                f"Discovery row count changed for {run['run_id']}: {len(discovered)} != {expected_rows}"
            )
        for pair in discovered.itertuples(index=False):
            for channel, raw_column, mask_column in (
                ("DAPI", "dapi_raw", "dapi_seg"),
                ("FOS", "fos_raw", "fos_seg"),
            ):
                source_raw = Path(str(getattr(pair, raw_column)))
                source_mask = Path(str(getattr(pair, mask_column)))
                raw_destination = APOTOME / source_raw.relative_to(SOURCE_APOTOME)
                mask_destination = APOTOME / source_mask.relative_to(SOURCE_APOTOME)
                key = str(mask_destination.resolve())
                required = required_masks.setdefault(key, {
                    "run_ids": [],
                    "acquisition_groups": [],
                    "channel": channel,
                    "input_plane": raw_destination.resolve(),
                    "hil_required": True,
                    "annotation_statuses": [],
                    "logical_uses": 0,
                })
                if required["channel"] != channel:
                    raise RuntimeError(f"Cellpose mask changes channel role: {mask_destination}")
                if Path(required["input_plane"]) != raw_destination.resolve():
                    raise RuntimeError(f"Cellpose mask changes input plane: {mask_destination}")
                required["run_ids"].append(str(run["run_id"]))
                required["acquisition_groups"].append(str(pair.tile_key))
                required["annotation_statuses"].append(str(pair.annotation_status))
                required["logical_uses"] += 1

    combined = pd.concat(staged_frames, ignore_index=True, sort=False)
    combined_sha = hashlib.sha256(combined.to_csv(index=False).encode("utf-8")).hexdigest()
    direct_combined = APOTOME / "combined_human_region_plots" / "combined_raw_per_animal_human_region_summary.csv"
    direct_frame = pd.read_csv(direct_combined)
    direct_combined_sha = sha256(direct_combined)
    if direct_combined_sha != CANONICAL_COMBINED_SHA256:
        raise RuntimeError("Copied canonical combined CSV hash changed")
    if list(combined.columns) != list(direct_frame.columns) or len(combined) != len(direct_frame):
        raise RuntimeError("Relocated HIL rebuild schema/row count differs from audited direct CSV")
    semantic_columns = [
        column for column in combined.columns
        if column not in {"source_csv", "source_dataset"}
    ]
    for column in semantic_columns:
        left, right = combined[column], direct_frame[column]
        numeric = pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right)
        if numeric:
            if not np.allclose(left.to_numpy(float), right.to_numpy(float), rtol=1e-12, atol=1e-12, equal_nan=True):
                raise RuntimeError(f"Relocated HIL numeric column changed: {column}")
        elif not left.astype(str).equals(right.astype(str)):
            raise RuntimeError(f"Relocated HIL semantic column changed: {column}")


    values, _, _ = analysis_module.audit_and_prepare(pd.read_csv(direct_combined))
    n_table = (
        values.groupby(["region", "condition"], observed=True)["animal"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(index=["ARC", "ME", "VMN", "OTHERS"], columns=["Water", "Sucrose", "Allulose"])
    )
    anatomical = n_table.loc[["ARC", "ME", "VMN"]]
    if anatomical.isna().any().any() or (anatomical < 3).any().any():
        raise RuntimeError(f"Independent-animal inference gate failed:\n{anatomical}")

    top_level_masks: list[Path] = []
    for run in RUNS:
        top_level_masks.extend(sorted(Path(run["dataset"]).glob("*_seg.npy")))
    top_level_masks = sorted(set(path.resolve() for path in top_level_masks))
    all_npy = sorted(set(path.resolve() for path in APOTOME.rglob("*.npy")))
    if len(all_npy) != 113:
        raise RuntimeError(f"NPY inventory changed: {len(all_npy)} != 113")
    if len(top_level_masks) != 66:
        raise RuntimeError(f"Top-level Cellpose mask count changed: {len(top_level_masks)} != 66")
    if len(required_masks) != 50:
        raise RuntimeError(f"HIL-required unique mask count changed: {len(required_masks)} != 50")

    dataset_to_run = {str(Path(run["dataset"]).resolve()): str(run["run_id"]) for run in RUNS}
    mask_rows: list[dict[str, object]] = []
    for mask in top_level_masks:
        required = required_masks.get(str(mask))
        if required is None:
            input_plane = companion_for_ancillary(mask)
            channel = channel_from_mask(mask)
            run_id = dataset_to_run[str(mask.parent.resolve())]
            acquisition = re.sub(
                r"(?i)(?:[-_ ](?:000[1-4][-_])?)?(dapi|fos|gfp|neun|masks|overlay)$",
                "",
                mask.name[: -len("_seg.npy")],
            ).strip("-_ ")
            required = {
                "run_ids": [run_id],
                "acquisition_groups": [acquisition],
                "channel": channel,
                "input_plane": input_plane,
                "hil_required": False,
                "annotation_statuses": ["not_in_hil_pair_ancillary_preserved_mask"],
                "logical_uses": 0,
            }
        input_plane = Path(required["input_plane"]).resolve()
        raw_row = manifest_record(input_plane, by_destination, "microscopy_image")
        mask_row = manifest_record(mask, by_destination, "cellpose_mask")
        mask_rows.append(
            {
                "run_ids": ";".join(sorted(set(required["run_ids"]))),
                "acquisition_groups": ";".join(sorted(set(required["acquisition_groups"]))),
                "channel": required["channel"],
                "input_plane_source": raw_row.source_path,
                "input_plane_destination": raw_row.destination_path,
                "input_plane_sha256": raw_row.destination_sha256,
                "cellpose_mask_source": mask_row.source_path,
                "cellpose_mask_destination": mask_row.destination_path,
                "cellpose_mask_sha256": mask_row.destination_sha256,
                "hil_required": str(bool(required["hil_required"])).lower(),
                "annotation_statuses": ";".join(sorted(set(required["annotation_statuses"]))),
                "hil_logical_uses": int(required["logical_uses"]),
                "bundle_contract_required": "true",
                "mapping_status": "PASS",
            }
        )
    hil_unique = sum(row["hil_required"] == "true" for row in mask_rows)
    hil_logical_uses = sum(int(row["hil_logical_uses"]) for row in mask_rows)
    if len(mask_rows) != 66 or hil_unique != 50 or hil_logical_uses != 54:
        raise RuntimeError("Cellpose mapping cardinality failed")

    registration_assets = targeted_registration_inventory()
    if registration_assets:
        raise RuntimeError(
            "Unexpected atlas/registration-named assets require manual classification: "
            + repr(registration_assets)
        )

    coordinates_total = sum(int(row["nucleus_coordinate_rows"]) for row in coordinate_summary)
    polygons_total = sum(int(row["polygon_points"]) for row in coordinate_summary)
    polygons_outside = sum(int(row["polygon_points_outside_native_bounds"]) for row in coordinate_summary)
    completed_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    readiness = {
        "status": "PASS_WITH_DOCUMENTED_NONRAW_BLOCKERS",
        "completed_utc": completed_utc,
        "raw_payload": {
            "files": len(manifest),
            "bytes": int(pd.to_numeric(manifest["size_bytes"]).sum()),
            "all_source_destination_sha256_equal": True,
            "all_distinct_inodes": True,
            "all_metadata_match": True,
        },
        "hil": {
            "runs": hil_summary,
            "relocated_hil_rebuild_sha256": combined_sha,
            "audited_direct_combined_sha256": direct_combined_sha,
            "audited_direct_combined_matches_canonical": True,
            "relocated_rebuild_semantics_match_except_provenance_fields": True,
        },
        "cellpose": {
            "all_npy_files": len(all_npy),
            "top_level_cellpose_masks": len(mask_rows),
            "hil_required_unique_masks": 50,
            "hil_required_logical_mask_uses": 54,
            "ancillary_preserved_masks": 16,
            "hil_dapi_fos_logical_pairs": 27,
            "hil_dapi_fos_unique_acquisitions": 25,
            "mapping_receipt": str(outputs[0]),
            "all_mapped_and_hash_verified": True,
        },
        "spatial": {
            "coordinate_system": "native_image_pixel_yx",
            "nucleus_coordinate_rows": coordinates_total,
            "annotation_polygon_points": polygons_total,
            "annotation_polygon_points_outside_native_bounds_before_raster_clipping": polygons_outside,
            "all_nucleus_coordinates_within_region_mask_bounds": True,
            "saved_region_masks_are_authoritative_clipped_annotation_geometry": True,
            "region_mask_tiffs_present": True,
            "atlas_registration_assets_found": 0,
            "section_local_spatial_shell_ready": True,
            "atlas_normalized_spatial_shell_ready": False,
        },
        "inference_gate": {
            "unit": "independent_animal",
            "minimum_per_W_S_A_condition": 3,
            "n_by_region_condition": {
                region: {condition: int(n_table.loc[region, condition]) for condition in n_table.columns}
                for region in n_table.index
            },
            "anatomical_endpoints_pass": True,
        },
        "nonraw_blockers": [
            {
                "path": str(PAPER / "reference_assets" / "Fig2" / "Fig2.pdf"),
                "status": "RENDERER_ASSET_BLOCKER_SOURCE_UNBOUND",
            },
            {
                "path": str(
                    PAPER / "reference_assets" / "Fig2"
                    / "ChatGPT Image Jun 5, 2026, 01_51_01 PM.png"
                ),
                "status": "RENDERER_ASSET_BLOCKER_SOURCE_UNBOUND",
            },
            {
                "path": str(
                    APOTOME / "24_06_2025" / "analysis"
                    / "apotome_cfos_human_in_loop_v9.py"
                ),
                "status": "MISSING_HISTORICAL_CODE_NONRAW",
            },
            {
                "path": "atlas_registration_or_transform_assets",
                "status": "ABSENT_NOT_REQUIRED_FOR_SECTION_LOCAL_SPATIAL_PANELS",
            },
        ],
    }

    table_lines = [
        "| Region | Water n | Sucrose n | Allulose n | >=3 each |",
        "|---|---:|---:|---:|:---:|",
    ]
    for region in n_table.index:
        counts = [int(n_table.loc[region, condition]) for condition in n_table.columns]
        table_lines.append(
            f"| {region} | {counts[0]} | {counts[1]} | {counts[2]} | "
            f"{'PASS' if min(counts) >= 3 else 'DESCRIPTIVE ONLY'} |"
        )
    report = "\n".join(
        [
            "# Figure 2 reconstruction readiness",
            "",
            f"Audit status: **{readiness['status']}**  ",
            f"Completed: {completed_utc}",
            "",
            "## Verified raw and HIL chain",
            "",
            "- 560 physically independent copied files; 108,000,611,886 bytes.",
            "- Every recorded source/destination SHA-256 is equal; metadata and distinct-inode gates pass.",
            "- All five completed HIL runs reconcile nucleus -> image -> animal without pending annotations.",
            f"- Audited direct combined HIL SHA-256: `{direct_combined_sha}` (canonical match).",
            f"- Relocated HIL rebuild SHA-256: `{combined_sha}`; quantitative fields match, while absolute source paths and eight historical source-dataset labels are intentionally provenance-different.",
            "- 113 NPY files are preserved. The explicit map has 66 unique top-level masks: 50 HIL-required unique masks used 54 times by 27 logical DAPI/FOS pairs (25 unique acquisitions), plus 16 ancillary channel/alternate masks.",
            "",
            "## Spatial readiness",
            "",
            f"- All {coordinates_total:,} per-nucleus native-image coordinates are finite and inside their saved region-mask TIFF bounds.",
            f"- The {polygons_total:,} saved annotation vertices are finite; {polygons_outside:,} pre-rasterization vertices extend outside native bounds and are clipped in the authoritative saved region masks.",
            "- Section-local animal-level spatial shells can be constructed from these coordinates, Cellpose labels, region assignments, c-FOS positivity, and region masks.",
            "- No atlas-registration, affine, landmark, BrainGlobe/brainreg, or QuickNII asset is present. Atlas-normalized cross-section spatial shells therefore remain blocked unless a proven registration source is supplied; this does not block native-section panels.",
            "",
            "## Independent-animal inference gate",
            "",
            *table_lines,
            "",
            "Cells and sections are not n. Anatomical inference may run only where every Water/Sucrose/Allulose cell has at least three independent animals; otherwise the endpoint must remain descriptive with no p-value.",
            "",
            "## Remaining non-raw/final-render blockers",
            "",
            "- The legacy Fig2 PDF and exact timeline PNG have no proven surviving source; both remain renderer-asset blockers, not raw-data omissions.",
            "- The broken historical HIL engine link was correctly excluded as `MISSING_HISTORICAL_CODE_NONRAW`; completed HIL reuse/validation is unaffected, but interactive re-annotation would require a separately audited engine.",
            "- The current canonical renderer does not yet implement the requested full-width Fig3-style spatial row, enhanced-red-only display policy, bilingual isolated panels, paired-language count audit, >=96%/94% occupancy gate, or label/data overlap/clearance receipt. Final reconstruction is not publication-ready until those code/output gates are implemented.",
            "- Final statistical graphics must reuse the Figure 4 mean +/- SD plus individual-animal grammar (not quartile boxes): width .52, edge 1.7, error/cap 1.6, s=32, jitter +/- .065, W #b9e3f2, S #e31a1c, A #2ecc71, DejaVu Sans, axes 1.6.",
            "",
            "## Runnable staging command",
            "",
            "Run `Fig2/RUN_RECONSTRUCTION_STAGED.sh /absolute/fresh/stage_directory`. It revalidates this bundle and writes canonical animal-level analysis only to the new stage. It deliberately stops before native-CZI rendering/final layout while the documented renderer assets and new layout/export contract remain unresolved.",
            "",
        ]
    )

    if args.write_receipts:
        atomic_csv(outputs[0], mask_rows)
        atomic_text(outputs[1], json.dumps(readiness, indent=2, sort_keys=True) + "\n")
        atomic_text(outputs[2], report)

    print(
        json.dumps(
            {
                "status": readiness["status"],
                "payload_files": len(manifest),
                "payload_bytes": int(pd.to_numeric(manifest["size_bytes"]).sum()),
                "hil_runs": len(hil_summary),
                "audited_direct_combined_sha256": direct_combined_sha,
                "relocated_hil_rebuild_sha256": combined_sha,
                "npy_files": len(all_npy),
                "cellpose_masks": len(mask_rows),
                "hil_required_unique_masks": 50,
                "hil_required_logical_uses": 54,
                "nucleus_coordinate_rows": coordinates_total,
                "atlas_registration_assets": 0,
                "receipts_written": bool(args.write_receipts),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
