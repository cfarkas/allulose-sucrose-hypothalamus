#!/usr/bin/env python3
"""Rebuild Figure 2's canonical input from completed manual HIL runs.

This is a fail-closed staging step between the five native Apotome
HIL region-annotation runs and the Figure 2 statistical analysis. It does
not infer or edit any region. It verifies the saved completion receipts,
annotation decisions and masks, reconciles nucleus/image/animal tables, and
then concatenates the five per-animal HIL summaries in their audited order.

Outputs are the Figure 2 raw combined CSV consumed by
01_analyze_global_cfos.py and file-level CSV/JSON provenance inventories for
every completion marker, status/quantitative table, annotation JSON and
corresponding region-mask TIFF.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
_paper_override = os.environ.get("FIG2_PAPER_ROOT")
PAPER = (
    Path(_paper_override).expanduser().resolve()
    if _paper_override
    else next(
        (
            path
            for path in HERE.parents
            if path.name == "Paper"
            or (
                (path / "scripts" / "setup").is_dir()
                and (path / "README.txt").is_file()
            )
        ),
        None,
    )
)
if PAPER is None:  # pragma: no cover - protects relocated standalone copies
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
FIG2 = PAPER / "analyses" / "Fig2"
DEFAULT_HIL_ROOT = FIG2 / "hil" / "runs"

DEFAULT_OUTPUT = FIG2 / "raw" / "combined_raw_per_animal_human_region_summary.csv"
DEFAULT_PROVENANCE = FIG2 / "provenance"
AUDITED_COMBINED_SHA256 = "d80bfe4dd5b259f1d309dcf59f54cb0effd62d7911f72de769cfc2b9b8170ebf"
EXPECTED_PIPELINE_VERSION = "9.2.0"
EXPECTED_REGIONS = ("ME", "ARC", "VMN", "OTHERS")


@dataclass(frozen=True)
class HilRun:
    run_id: str
    relative_root: str
    source_dataset: str
    annotated: int
    excluded: int
    pending: int
    per_animal_rows: int


# This order reproduces the audited upstream concatenation byte for byte.
HIL_RUNS = (
    HilRun(
        "june_2025_tiled",
        "june_2025_tiled",
        "2025.06.24_human_regions_v9_2",
        6,
        4,
        0,
        20,
    ),
    HilRun(
        "june_2025_primary",
        "june_2025_primary",
        "human_regions_v9",
        7,
        0,
        0,
        20,
    ),
    HilRun(
        "september_2025_water",
        "september_2025_water",
        "september_2025_water",
        2,
        0,
        0,
        8,
    ),
    HilRun(
        "august_2025_horizontal",
        "august_2025_horizontal",
        "selected_horizontal_human_regions_v9_2",
        2,
        0,
        0,
        8,
    ),
    HilRun(
        "august_2025_selected",
        "august_2025_selected",
        "selected_human_regions_v9_2",
        6,
        0,
        0,
        20,
    ),
)


REQUIRED_ANIMAL_COLUMNS = {
    "animal",
    "cond",
    "region",
    "sample",
    "dapi_nuclei",
    "cfos_cells",
    "area_px",
    "area_um2",
    "n_images",
    "cfos_over_dapi",
}
REQUIRED_IMAGE_COLUMNS = {
    "tile_key",
    "animal",
    "cond",
    "region",
    "dapi_nuclei",
    "cfos_cells",
    "area_px",
    "area_um2",
}
REQUIRED_NUCLEUS_COLUMNS = {
    "tile_key",
    "region",
    "nucleus_label",
    "cfos_positive",
}


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_completion_marker(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "=" not in line:
            raise RuntimeError(f"Malformed completion-marker line in {path}: {raw_line!r}")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    required = {
        "completed_at_utc",
        "pipeline_version",
        "annotated",
        "excluded",
        "pending",
        "input",
        "output",
    }
    missing = sorted(required - values.keys())
    if missing:
        raise RuntimeError(f"Completion marker {path} is missing: {missing}")
    return values


def require_columns(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise RuntimeError(f"{path} is missing required columns: {missing}")


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def artifact_row(
    run_index: int,
    run: HilRun,
    run_root: Path,
    role: str,
    path: Path,
    *,
    tile_key: str = "",
    status: str = "",
    rows: int | None = None,
    schema_version: object = "",
    pipeline_version: object = "",
) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required HIL artifact is missing: {path}")
    return {
        "run_index": run_index,
        "run_id": run.run_id,
        "source_dataset": run.source_dataset,
        "run_root": str(run_root),
        "artifact_role": role,
        "tile_key": tile_key,
        "decision_status": status,
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "rows": "" if rows is None else rows,
        "schema_version": schema_version,
        "pipeline_version": pipeline_version,
        "expected_annotated": run.annotated,
        "expected_excluded": run.excluded,
        "expected_pending": run.pending,
    }


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    mapped = series.astype(str).str.strip().str.casefold().map(
        {"true": True, "1": True, "yes": True, "false": False, "0": False, "no": False}
    )
    if mapped.isna().any():
        bad = sorted(series.loc[mapped.isna()].astype(str).unique())
        raise RuntimeError(f"Unrecognised cfos_positive values: {bad}")
    return mapped.astype(bool)


def relocated_run_artifact(run_root: Path, stored: object, directory: str) -> Path:
    """Resolve copied HIL artifacts from their internal run before old absolute paths."""
    original = Path(str(stored)).expanduser()
    local = (run_root / directory / original.name).resolve()
    if local.is_file():
        return local
    raise RuntimeError(
        f"Manifest-bound HIL artifact is missing inside the current tree: {local}"
    )


def validate_quantitative_chain(
    run_root: Path,
    status: pd.DataFrame,
    per_nucleus: pd.DataFrame,
    per_image: pd.DataFrame,
    per_animal: pd.DataFrame,
) -> dict[str, int]:
    """Reconcile nucleus -> image -> animal counts without changing source data."""

    require_columns(
        per_nucleus,
        REQUIRED_NUCLEUS_COLUMNS,
        run_root / "per_nucleus_human_region_assignments.csv",
    )
    require_columns(
        per_image,
        REQUIRED_IMAGE_COLUMNS,
        run_root / "per_image_human_region_summary.csv",
    )
    require_columns(
        per_animal,
        REQUIRED_ANIMAL_COLUMNS,
        run_root / "per_animal_human_region_summary.csv",
    )

    annotated_tiles = set(status.loc[status["status"].eq("annotated"), "tile_key"].astype(str))
    nonannotated_tiles = set(status.loc[~status["status"].eq("annotated"), "tile_key"].astype(str))
    image_tiles = set(per_image["tile_key"].astype(str))
    nucleus_tiles = set(per_nucleus["tile_key"].astype(str))
    if image_tiles != annotated_tiles or nucleus_tiles != annotated_tiles:
        raise RuntimeError(
            f"HIL quantitative tile coverage mismatch under {run_root}: "
            f"status={len(annotated_tiles)}, image={len(image_tiles)}, nucleus={len(nucleus_tiles)}"
        )
    if image_tiles & nonannotated_tiles:
        raise RuntimeError(f"Excluded/pending tiles appear in quantitative tables under {run_root}")

    per_image_region_counts = per_image.groupby("tile_key", observed=True)["region"].agg(
        lambda values: tuple(sorted(set(values.astype(str))))
    )
    expected_region_set = tuple(sorted(EXPECTED_REGIONS))
    bad_tiles = per_image_region_counts.loc[
        per_image_region_counts.map(lambda observed: observed != expected_region_set)
    ]
    if not bad_tiles.empty:
        raise RuntimeError(
            f"Per-image region coverage is incomplete under {run_root}: {bad_tiles.to_dict()}"
        )

    nuclei = per_nucleus.copy()
    nuclei["cfos_positive_bool"] = bool_series(nuclei["cfos_positive"])
    nucleus_counts = (
        nuclei.groupby(["tile_key", "region"], observed=True)
        .agg(
            dapi_from_nuclei=("nucleus_label", "size"),
            cfos_from_nuclei=("cfos_positive_bool", "sum"),
        )
        .reset_index()
    )
    image_check = per_image.merge(nucleus_counts, on=["tile_key", "region"], how="left")
    image_check[["dapi_from_nuclei", "cfos_from_nuclei"]] = image_check[
        ["dapi_from_nuclei", "cfos_from_nuclei"]
    ].fillna(0)
    dapi_match = image_check["dapi_nuclei"].astype(int).eq(
        image_check["dapi_from_nuclei"].astype(int)
    )
    cfos_match = image_check["cfos_cells"].astype(int).eq(
        image_check["cfos_from_nuclei"].astype(int)
    )
    if not bool((dapi_match & cfos_match).all()):
        raise RuntimeError(
            f"Nucleus-to-image DAPI/c-FOS reconciliation failed for "
            f"{int((~(dapi_match & cfos_match)).sum())} rows under {run_root}"
        )

    image_aggregate = (
        per_image.groupby(["animal", "cond", "region"], observed=True, as_index=False)
        .agg(
            dapi_nuclei=("dapi_nuclei", "sum"),
            cfos_cells=("cfos_cells", "sum"),
            area_px=("area_px", "sum"),
            area_um2=("area_um2", "sum"),
            n_images=("tile_key", "nunique"),
        )
    )
    keys = ["animal", "cond", "region"]
    animal_check = per_animal.merge(
        image_aggregate,
        on=keys,
        how="outer",
        suffixes=("_saved", "_from_images"),
        indicator=True,
    )
    if not animal_check["_merge"].eq("both").all():
        raise RuntimeError(f"Image-to-animal key reconciliation failed under {run_root}")
    for column in ("dapi_nuclei", "cfos_cells", "area_px", "area_um2", "n_images"):
        saved = pd.to_numeric(animal_check[f"{column}_saved"], errors="coerce")
        rebuilt = pd.to_numeric(animal_check[f"{column}_from_images"], errors="coerce")
        if not np.allclose(saved, rebuilt, rtol=1e-12, atol=1e-9, equal_nan=True):
            raise RuntimeError(f"Image-to-animal {column} reconciliation failed under {run_root}")

    return {
        "annotated_tiles_reconciled": len(annotated_tiles),
        "nucleus_rows_reconciled": len(per_nucleus),
        "image_region_rows_reconciled": len(per_image),
        "animal_region_rows_reconciled": len(per_animal),
    }


def validate_run(
    run_index: int, run: HilRun, hil_root: Path
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    run_root = (hil_root / run.relative_root).resolve()
    if not run_root.is_dir():
        raise FileNotFoundError(f"Required HIL run directory is missing: {run_root}")

    marker_path = run_root / "ANNOTATION_COMPLETE.txt"
    status_path = run_root / "annotation_status.csv"
    per_nucleus_path = run_root / "per_nucleus_human_region_assignments.csv"
    per_image_path = run_root / "per_image_human_region_summary.csv"
    per_animal_path = run_root / "per_animal_human_region_summary.csv"
    marker = parse_completion_marker(marker_path)
    if marker["pipeline_version"] != EXPECTED_PIPELINE_VERSION:
        raise RuntimeError(
            f"Unexpected HIL pipeline version for {run.run_id}: {marker['pipeline_version']}"
        )
    marker_counts = {
        name: int(marker[name]) for name in ("annotated", "excluded", "pending")
    }
    expected_counts = {
        "annotated": run.annotated,
        "excluded": run.excluded,
        "pending": run.pending,
    }
    if marker_counts["pending"] != 0:
        raise RuntimeError(
            f"HIL run {run.run_id} is incomplete: pending={marker_counts['pending']}"
        )
    marker["relocated_input"] = str((FIG2 / "raw" / "datasets" / run.run_id).resolve())
    marker["relocated_output"] = str(run_root)

    status = pd.read_csv(status_path)
    required_status = {
        "tile_key",
        "sample",
        "animal",
        "cond",
        "status",
        "annotation_json",
        "region_mask",
    }
    require_columns(status, required_status, status_path)
    if status["tile_key"].astype(str).duplicated().any():
        raise RuntimeError(f"Duplicate tile keys in {status_path}")
    status_counts = status["status"].value_counts().to_dict()
    observed_counts = {
        "annotated": int(status_counts.get("annotated", 0)),
        "excluded": int(status_counts.get("excluded", 0)),
        "pending": int(status_counts.get("pending", 0)),
    }
    unknown_statuses = sorted(set(status["status"].astype(str)) - set(observed_counts))
    if unknown_statuses or observed_counts != marker_counts:
        raise RuntimeError(
            f"Status table does not match completion marker for {run.run_id}: "
            f"counts={observed_counts}, unknown={unknown_statuses}"
        )

    per_nucleus = pd.read_csv(per_nucleus_path)
    per_image = pd.read_csv(per_image_path)
    per_animal = pd.read_csv(per_animal_path)
    chain = validate_quantitative_chain(
        run_root, status, per_nucleus, per_image, per_animal
    )

    inventory: list[dict[str, Any]] = [
        artifact_row(
            run_index,
            run,
            run_root,
            "completion_marker",
            marker_path,
            pipeline_version=marker["pipeline_version"],
        ),
        artifact_row(
            run_index,
            run,
            run_root,
            "annotation_status",
            status_path,
            rows=len(status),
        ),
        artifact_row(
            run_index,
            run,
            run_root,
            "per_nucleus_assignments",
            per_nucleus_path,
            rows=len(per_nucleus),
        ),
        artifact_row(
            run_index,
            run,
            run_root,
            "per_image_summary",
            per_image_path,
            rows=len(per_image),
        ),
        artifact_row(
            run_index,
            run,
            run_root,
            "per_animal_summary",
            per_animal_path,
            rows=len(per_animal),
        ),
    ]

    annotation_records: list[dict[str, Any]] = []
    for status_row in status.sort_values("tile_key").itertuples(index=False):
        annotation_path = relocated_run_artifact(
            run_root, status_row.annotation_json, "annotations"
        )
        mask_path = relocated_run_artifact(
            run_root, status_row.region_mask, "region_masks"
        )
        if annotation_path.parent != run_root / "annotations":
            raise RuntimeError(f"Annotation path escapes its HIL run: {annotation_path}")
        if mask_path.parent != run_root / "region_masks":
            raise RuntimeError(f"Mask path escapes its HIL run: {mask_path}")
        if not annotation_path.is_file() or not mask_path.is_file():
            raise FileNotFoundError(
                f"Missing annotation/mask for {run.run_id}/{status_row.tile_key}: "
                f"{annotation_path}, {mask_path}"
            )
        payload = json.loads(annotation_path.read_text(encoding="utf-8"))
        expected_status = {"include": "annotated", "exclude": "excluded"}.get(
            str(payload.get("decision", ""))
        )
        if expected_status != status_row.status:
            raise RuntimeError(
                f"Annotation decision/status mismatch for "
                f"{run.run_id}/{status_row.tile_key}: "
                f"decision={payload.get('decision')!r}, status={status_row.status!r}"
            )
        for key in ("tile_key", "sample", "animal", "cond"):
            if str(payload.get(key, "")) != str(getattr(status_row, key)):
                raise RuntimeError(
                    f"Annotation metadata mismatch for "
                    f"{run.run_id}/{status_row.tile_key}: {key}"
                )
        if str(payload.get("pipeline_version")) != EXPECTED_PIPELINE_VERSION:
            raise RuntimeError(
                f"Annotation pipeline version mismatch for "
                f"{run.run_id}/{status_row.tile_key}"
            )
        if int(payload.get("schema_version", -1)) != 1:
            raise RuntimeError(
                f"Annotation schema version mismatch for "
                f"{run.run_id}/{status_row.tile_key}"
            )
        ann_row = artifact_row(
            run_index,
            run,
            run_root,
            "annotation_json",
            annotation_path,
            tile_key=str(status_row.tile_key),
            status=str(status_row.status),
            schema_version=payload["schema_version"],
            pipeline_version=payload["pipeline_version"],
        )
        mask_row = artifact_row(
            run_index,
            run,
            run_root,
            "region_mask_tiff",
            mask_path,
            tile_key=str(status_row.tile_key),
            status=str(status_row.status),
            schema_version=payload["schema_version"],
            pipeline_version=payload["pipeline_version"],
        )
        inventory.extend((ann_row, mask_row))
        annotation_records.append(
            {
                "tile_key": str(status_row.tile_key),
                "sample": str(status_row.sample),
                "animal": str(status_row.animal),
                "condition": str(status_row.cond),
                "status": str(status_row.status),
                "decision": str(payload["decision"]),
                "saved_at_utc": str(payload.get("saved_at_utc", "")),
                "annotation_json": str(annotation_path),
                "annotation_json_sha256": ann_row["sha256"],
                "region_mask": str(mask_path),
                "region_mask_sha256": mask_row["sha256"],
            }
        )

    staged = per_animal.copy()
    staged["source_csv"] = str(per_animal_path)
    staged["source_dataset"] = run.source_dataset
    run_record = {
        "run_index": run_index,
        "run_id": run.run_id,
        "run_root": str(run_root),
        "source_dataset": run.source_dataset,
        "completion": marker,
        "status_counts": observed_counts,
        "accepted_2026_08_16_snapshot_counts": expected_counts,
        "counts_match_accepted_2026_08_16_snapshot": observed_counts == expected_counts,
        "accepted_2026_08_16_per_animal_rows": run.per_animal_rows,
        "per_animal_rows_match_accepted_2026_08_16_snapshot": (
            len(per_animal) == run.per_animal_rows
        ),
        "quantitative_chain": chain,
        "per_animal_source": str(per_animal_path),
        "per_animal_source_sha256": sha256(per_animal_path),
        "annotations": annotation_records,
        "validation_status": "passed",
    }
    return staged, inventory, run_record


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hil-root",
        type=Path,
        default=DEFAULT_HIL_ROOT,
        help="Directory containing the five named HIL run directories.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--provenance-dir", type=Path, default=DEFAULT_PROVENANCE)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate and reconcile all inputs without writing any output.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace only the three named output files if they already exist.",
    )
    parser.add_argument(
        "--allow-live-output",
        action="store_true",
        help=(
            "Explicitly permit writes inside this Paper tree. This does not "
            "imply --force and is never needed for --check-only."
        ),
    )
    return parser


def is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def validate_write_targets(
    output: Path, provenance: Path, *, force: bool, allow_live_output: bool
) -> tuple[Path, Path]:
    """Fail closed before validation if a requested publish target is unsafe."""
    inventory_path = provenance / "hil_input_inventory.csv"
    manifest_path = provenance / "hil_input_inventory.json"
    targets = (output, inventory_path, manifest_path)
    if len(set(targets)) != len(targets):
        raise RuntimeError(f"Output targets must be distinct: {targets}")
    if any(is_within(path, PAPER) for path in targets) and not allow_live_output:
        raise RuntimeError(
            "Refusing to write inside the live Paper tree. Use a fresh staging "
            "path outside Paper, or pass --allow-live-output after review."
        )
    existing = [path for path in targets if path.exists()]
    if existing and not force:
        rendered = ", ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Refusing to replace existing output(s) without --force: {rendered}"
        )
    for path in targets:
        if path.exists() and not path.is_file():
            raise RuntimeError(f"Output target exists but is not a regular file: {path}")
    return inventory_path, manifest_path


def main() -> int:
    args = build_parser().parse_args()
    hil_root = args.hil_root.expanduser().resolve()
    output = args.output.expanduser().resolve()
    provenance = args.provenance_dir.expanduser().resolve()
    inventory_path = provenance / "hil_input_inventory.csv"
    manifest_path = provenance / "hil_input_inventory.json"
    if not args.check_only:
        inventory_path, manifest_path = validate_write_targets(
            output,
            provenance,
            force=args.force,
            allow_live_output=args.allow_live_output,
        )

    staged_frames: list[pd.DataFrame] = []
    inventory_rows: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []
    for run_index, run in enumerate(HIL_RUNS, start=1):
        staged, inventory, record = validate_run(run_index, run, hil_root)
        staged_frames.append(staged)
        inventory_rows.extend(inventory)
        run_records.append(record)

    combined = pd.concat(staged_frames, ignore_index=True, sort=False)
    combined_bytes = combined.to_csv(index=False).encode("utf-8")
    combined_sha = hashlib.sha256(combined_bytes).hexdigest()

    inventory = pd.DataFrame(inventory_rows)
    inventory = inventory.sort_values(
        ["run_index", "artifact_role", "tile_key", "path"], kind="stable"
    ).reset_index(drop=True)
    inventory_bytes = inventory.to_csv(index=False).encode("utf-8")
    inventory_sha = hashlib.sha256(inventory_bytes).hexdigest()

    manifest = {
        "schema_version": 1,
        "purpose": (
            "Figure 2 input rebuilt exclusively from five completed "
            "HIL Apotome region-annotation runs"
        ),
        "script": str(HERE),
        "script_sha256": sha256(HERE),
        "root": str(hil_root),
        "hil_pipeline_version_required": EXPECTED_PIPELINE_VERSION,
        "validation_policy": {
            "pending_annotations_allowed": False,
            "current_complete_human_decisions_may_differ_from_accepted_snapshot": True,
            "accepted_snapshot_count_and_row_matches_are_recorded_not_enforced": True,
            "annotation_decision_must_match_status": True,
            "annotation_json_and_region_mask_required": True,
            "all_input_artifacts_sha256_hashed": True,
            "nucleus_to_image_counts_reconciled": True,
            "image_to_animal_counts_and_areas_reconciled": True,
            "statistical_exclusions_or_condition_overrides_applied_here": False,
        },
        "runs": run_records,
        "totals": {
            "runs": len(run_records),
            "annotated_receipts": sum(
                int(record["status_counts"]["annotated"]) for record in run_records
            ),
            "excluded_receipts": sum(
                int(record["status_counts"]["excluded"]) for record in run_records
            ),
            "pending_receipts": sum(
                int(record["status_counts"]["pending"]) for record in run_records
            ),
            "nucleus_rows_reconciled": sum(
                int(record["quantitative_chain"]["nucleus_rows_reconciled"])
                for record in run_records
            ),
            "combined_per_animal_region_source_rows": len(combined),
        },
        "output": {
            "path": str(output),
            "rows": len(combined),
            "columns": list(combined.columns),
            "size_bytes": len(combined_bytes),
            "sha256": combined_sha,
            "matches_audited_2026_08_16_snapshot": (
                combined_sha == AUDITED_COMBINED_SHA256
            ),
            "audited_2026_08_16_sha256": AUDITED_COMBINED_SHA256,
        },
        "inventory_csv": str(inventory_path),
        "inventory_csv_sha256": inventory_sha,
        "validation_status": "passed",
    }
    if not args.check_only:
        atomic_csv(combined, output)
        atomic_csv(inventory, inventory_path)
        atomic_json(manifest, manifest_path)

    print(f"[OK] Validated {len(HIL_RUNS)} completed HIL runs")
    print(
        "[OK] Receipts: "
        f"annotated={manifest['totals']['annotated_receipts']}, "
        f"excluded={manifest['totals']['excluded_receipts']}, "
        f"pending={manifest['totals']['pending_receipts']}"
    )
    print(
        f"[OK] Reconciled nucleus rows: "
        f"{manifest['totals']['nucleus_rows_reconciled']}"
    )
    if args.check_only:
        print(f"[OK] Check-only: reconciled {len(combined)} rows; wrote no files")
    else:
        print(f"[OK] Wrote {len(combined)} rows to {output}")
    print(f"[OK] SHA-256: {combined_sha}")
    if combined_sha != AUDITED_COMBINED_SHA256:
        print(
            "[NOTE] HIL inputs differ from the audited 2026-08-16 snapshot; "
            "rerun and audit statistics."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
