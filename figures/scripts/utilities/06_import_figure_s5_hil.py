#!/usr/bin/env python3
"""Validate and safely import the frozen, accepted Supplementary Figure S5 HIL annotations.

The importer is intentionally fail-closed: it validates the accepted receipt,
manifest, source DAPI stack and segmentation hashes, annotation/mask hashes,
region codes, dimensions, polygon counts, and each freshly prepared DAPI MIP
against either its frozen TIFF hash or its canonical reviewed-pixel hash before writing. It only stages into absent/empty destinations,
copies atomically, verifies every copy, and never deletes or overwrites.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import tifffile


HERE = Path(__file__).resolve()
PAPER = next(
    (
        parent
        for parent in HERE.parents
        if parent.name == "Paper"
        or (
            (parent / "scripts" / "setup").is_dir()
            and (parent / "README.txt").is_file()
        )
    ),
    None,
)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

DEFAULT_HIL = PAPER / "FigS5" / "hil_review" / "human_final_20260828_v1"
EXPECTED_CODES = {0, 1, 2, 3}
EXPECTED_REGIONS = {"ARC": 1, "ME": 2, "VMN": 3, "outside": 0}
EXPECTED_REVIEW_MIP_PIXEL_SHA256 = {
    "2026-07-09_FR4-1_Allulose": "9b57da6d433cac9d20ae884a41a65f8d562dd418d0f6f5f6a7993dcbd004aec9",
    "2026-07-09_FR5-1_Allulose": "a62b70496992dd96280b343b5921b043493bae2e722a0f6a4473db097d97c328",
    "2026-07-09_FR5-2_Allulose": "947123b99140e98bb83845e269c390cc8f36fea5e029f4d76330f57e9c604b8d",
    "2026-07-09_FR5-3_Sucrose": "2883e699f5a96e1fc501d803035536e123648bd3411c8dcd97f58d5c34a0e536",
    "2026-07-09_FR5-4_Sucrose": "6e2e64fafb224c7ae504dd9b629032c25a009a4ebe380ac37f76f658fcb41cfc",
    "2026-07-10_AGUA-M_Water": "f1bd18a6835c05c6ba7787fb4b614252864eab0228711b5d50526347b4691df7",
    "2026-07-10_FR4-3_Sucrose": "39a76fe442d192f37773c07770d538d3e7eb903f5f9d87aff147bed7773f6758",
    "2026-07-10_FR5-5_Water": "89a70d2ead3b5974c9804202259d020290e604244fe350a5b67959247a01f1d5",
    "2026-07-10_NPY-M_Water": "c5629347d9333cd5e62f77b8458e59d239fbf0a1558e724f6ee8375c36ffc66f",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_record_path(value: str, source: Path) -> Path:
    """Resolve frozen paths, including the documented Fig5-to-FigS5 move."""
    path = Path(str(value)).expanduser()
    candidates = [path] if path.is_absolute() else [PAPER / path, source / path]
    if not path.is_absolute() and path.parts and path.parts[0] == "Fig5":
        candidates.append(PAPER / "FigS5" / Path(*path.parts[1:]))
    if path.is_absolute():
        try:
            relative = path.relative_to(PAPER / "Fig5")
        except ValueError:
            pass
        else:
            candidates.append(PAPER / "FigS5" / relative)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Manifest-bound file is missing after the Fig5-to-FigS5 relocation: {value}"
    )


def require_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256(path)
    if actual != str(expected):
        raise RuntimeError(
            f"{label} SHA-256 mismatch for {path}: expected {expected}, got {actual}"
        )


def pixel_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def validate_review_mip(path: Path, row: dict[str, str]) -> str:
    """Accept the frozen TIFF hash or the exact reviewed decoded pixels.

    The Fig5-to-FigS5 move changed only TIFF description metadata (source path
    and pipeline label). The decoded pixels remain bound to both the historical
    reviewed-pixel hash and the hash-validated source-stack projection.
    """
    accepted_file_hash = str(row["review_dapi_mip_sha256"])
    if sha256(path) == accepted_file_hash:
        return "accepted_file_sha256"
    acquisition = str(row["acquisition_id"])
    prepared = np.asarray(tifffile.imread(path))
    expected_shape = (int(row["image_height"]), int(row["image_width"]))
    if prepared.shape != expected_shape or prepared.dtype != np.uint16:
        raise RuntimeError(
            f"{acquisition}: prepared DAPI MIP shape/dtype changed: "
            f"{prepared.shape}/{prepared.dtype} vs {expected_shape}/uint16"
        )
    expected_pixel_hash = EXPECTED_REVIEW_MIP_PIXEL_SHA256.get(acquisition)
    actual_pixel_hash = pixel_sha256(prepared)
    if actual_pixel_hash != expected_pixel_hash:
        raise RuntimeError(
            f"{acquisition}: DAPI MIP file metadata and decoded pixels differ "
            f"from the accepted review (pixels expected {expected_pixel_hash}, "
            f"got {actual_pixel_hash})"
        )
    source_stack = resolve_record_path(row["source_dapi_stack"], DEFAULT_HIL)
    stack = np.squeeze(np.asarray(tifffile.imread(source_stack)))
    while stack.ndim > 2:
        stack = np.max(stack, axis=0)
    if not np.array_equal(prepared, stack):
        raise RuntimeError(
            f"{acquisition}: prepared DAPI MIP does not equal the validated "
            "source-stack maximum projection"
        )
    return "canonical_reviewed_pixel_sha256"


def load_and_validate_source(
    source: Path,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    receipt_path = source / "HIL_RECEIPT.json"
    manifest_path = source / "HIL_MANIFEST.csv"
    if not receipt_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(f"Accepted HIL receipt/manifest missing under {source}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema") not in {"fig5_human_hil_receipt_v1", "figs5_human_hil_receipt_v2"}:
        raise RuntimeError("Unrecognized frozen Figure 5/S5 HIL receipt schema")
    if receipt.get("status") != "human_accepted":
        raise RuntimeError("Supplementary Figure S5 HIL receipt is not human_accepted")
    if receipt.get("raw_images_modified") or receipt.get("raw_cellpose_masks_modified"):
        raise RuntimeError("Accepted HIL receipt reports modified raw inputs")
    if receipt.get("region_codes") != EXPECTED_REGIONS:
        raise RuntimeError(
            f"Unexpected HIL region codes: {receipt.get('region_codes')}"
        )
    require_hash(manifest_path, receipt.get("manifest_sha256", ""), "HIL manifest")

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "acquisition_id",
        "sample_id",
        "animal_id",
        "condition",
        "decision",
        "schema_version",
        "annotation_path",
        "annotation_sha256",
        "region_mask_path",
        "region_mask_sha256",
        "image_height",
        "image_width",
        "mask_codes",
        "arc_polygons",
        "me_polygons",
        "vmn_polygons",
        "arc_pixels",
        "me_pixels",
        "vmn_pixels",
        "source_dapi_stack",
        "source_dapi_stack_sha256",
        "source_dapi_segmentation",
        "source_dapi_segmentation_sha256",
        "review_dapi_mip_sha256",
        "anatomy_channels_used",
    }
    if not rows or required - set(rows[0]):
        raise RuntimeError(
            f"HIL manifest missing columns: {sorted(required - set(rows[0] if rows else []))}"
        )
    acquisitions = [row["acquisition_id"] for row in rows]
    if len(acquisitions) != len(set(acquisitions)):
        raise RuntimeError("HIL manifest contains duplicate acquisition IDs")
    expected_ids = set(EXPECTED_REVIEW_MIP_PIXEL_SHA256)
    if len(acquisitions) == 8:
        expected_ids.remove("2026-07-10_NPY-M_Water")
    if set(acquisitions) != expected_ids:
        raise RuntimeError(
            "Canonical reviewed DAPI MIP pixel-hash inventory differs from HIL manifest"
        )
    expected_count = int(receipt.get("accepted_pairs", -1))
    if (
        len(rows) != expected_count
        or int(receipt.get("excluded_pairs", -1)) != 0
        or int(receipt.get("pending_pairs", -1)) != 0
    ):
        raise RuntimeError("HIL receipt counts do not describe a complete accepted set")

    for row in rows:
        acq = row["acquisition_id"]
        if row["decision"] != "include" or int(row["schema_version"]) != 2:
            raise RuntimeError(f"{acq}: expected included schema-v2 annotation")
        if row["anatomy_channels_used"] != "DAPI only":
            raise RuntimeError(f"{acq}: anatomy was not recorded as DAPI only")
        if row["mask_codes"] != "0=outside;1=ARC;2=ME;3=VMN":
            raise RuntimeError(f"{acq}: unexpected mask-code declaration")

        annotation_path = resolve_record_path(row["annotation_path"], source)
        mask_path = resolve_record_path(row["region_mask_path"], source)
        dapi_stack = resolve_record_path(row["source_dapi_stack"], source)
        dapi_seg = resolve_record_path(row["source_dapi_segmentation"], source)
        require_hash(annotation_path, row["annotation_sha256"], f"{acq} annotation")
        require_hash(mask_path, row["region_mask_sha256"], f"{acq} region mask")
        require_hash(dapi_stack, row["source_dapi_stack_sha256"], f"{acq} DAPI stack")
        require_hash(
            dapi_seg,
            row["source_dapi_segmentation_sha256"],
            f"{acq} DAPI segmentation",
        )

        annotation = json.loads(annotation_path.read_text(encoding="utf-8"))
        shape = (int(row["image_height"]), int(row["image_width"]))
        if (
            annotation.get("acquisition_id") != acq
            or annotation.get("decision") != "include"
            or int(annotation.get("schema_version", -1)) != 2
            or annotation.get("anatomy_channels_used") != ["DAPI"]
            or (
                int(annotation.get("image_height", -1)),
                int(annotation.get("image_width", -1)),
            )
            != shape
        ):
            raise RuntimeError(
                f"{acq}: annotation identity/channel/shape contract failed"
            )
        polygon_counts = {
            region: sum(
                str(poly.get("region", "")).upper() == region
                for poly in annotation.get("polygons", [])
            )
            for region in ("ARC", "ME", "VMN")
        }
        expected_polygons = {
            "ARC": int(row["arc_polygons"]),
            "ME": int(row["me_polygons"]),
            "VMN": int(row["vmn_polygons"]),
        }
        if polygon_counts != expected_polygons:
            raise RuntimeError(
                f"{acq}: polygon counts differ: {polygon_counts} vs {expected_polygons}"
            )

        mask = np.asarray(tifffile.imread(mask_path))
        if mask.shape != shape:
            raise RuntimeError(f"{acq}: region-mask shape {mask.shape} != {shape}")
        codes = {int(value) for value in np.unique(mask)}
        if not codes.issubset(EXPECTED_CODES):
            raise RuntimeError(f"{acq}: unexpected region-mask codes {sorted(codes)}")
        expected_pixels = {
            1: int(row["arc_pixels"]),
            2: int(row["me_pixels"]),
            3: int(row["vmn_pixels"]),
        }
        actual_pixels = {
            code: int(np.count_nonzero(mask == code)) for code in (1, 2, 3)
        }
        if actual_pixels != expected_pixels:
            raise RuntimeError(
                f"{acq}: region pixel counts differ: {actual_pixels} vs {expected_pixels}"
            )
    return receipt, rows


def validate_prepared_work(
    rows: list[dict[str, str]], workdir: Path
) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    sheet_path = workdir / "analysis_samplesheet.csv"
    if not sheet_path.is_file():
        raise FileNotFoundError(
            f"Run Supplementary Figure S5 prepare first; missing {sheet_path}"
        )
    with sheet_path.open(newline="", encoding="utf-8") as handle:
        sheet_rows = list(csv.DictReader(handle))
    by_id = {row["acquisition_id"]: row for row in sheet_rows}
    expected = {row["acquisition_id"] for row in rows}
    if set(by_id) != expected or len(by_id) != len(sheet_rows):
        raise RuntimeError(
            f"Prepared/HIL acquisition mismatch: prepared={sorted(by_id)}, HIL={sorted(expected)}"
        )
    mip_validation_modes: dict[str, str] = {}
    for row in rows:
        acq = row["acquisition_id"]
        prepared = by_id[acq]
        for column in ("sample_id", "animal_id", "condition"):
            if str(prepared.get(column, "")) != str(row[column]):
                raise RuntimeError(
                    f"{acq}: prepared {column}={prepared.get(column)!r} "
                    f"does not match accepted HIL {row[column]!r}"
                )
        mip = Path(prepared["dapi_mip"]).expanduser()
        if not mip.is_absolute():
            mip = (workdir / mip).resolve()
        if not mip.is_file():
            raise FileNotFoundError(f"{acq}: prepared DAPI MIP missing: {mip}")
        mip_validation_modes[acq] = validate_review_mip(mip, row)
    return by_id, mip_validation_modes


def copy_verified(source: Path, target: Path, expected_hash: str) -> None:
    if target.exists() or target.is_symlink():
        raise FileExistsError(f"Refusing to overwrite HIL target: {target}")
    incoming = target.with_name(target.name + ".incoming")
    if incoming.exists() or incoming.is_symlink():
        raise FileExistsError(f"Stale incoming path blocks safe HIL import: {incoming}")
    shutil.copy2(source, incoming)
    require_hash(incoming, expected_hash, "copied HIL file")
    os.replace(incoming, target)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from",
        dest="source",
        type=Path,
        default=DEFAULT_HIL,
        help="Frozen accepted Supplementary Figure S5 HIL directory.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        help="Fresh prepared Supplementary Figure S5 work directory.",
    )
    parser.add_argument(
        "--source-only",
        action="store_true",
        help="Validate the frozen source without requiring or modifying a work directory.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Validate source and prepared work but write nothing.",
    )
    args = parser.parse_args()

    source = args.source.expanduser().resolve()
    receipt, rows = load_and_validate_source(source)
    if args.source_only:
        print(
            f"[PASS] accepted Supplementary Figure S5 HIL source: {len(rows)} included DAPI-only fields"
        )
        return 0
    if args.workdir is None:
        parser.error("--workdir is required unless --source-only is used")
    workdir = args.workdir.expanduser().resolve()
    _, mip_validation_modes = validate_prepared_work(rows, workdir)
    if args.verify_only:
        print(
            f"[PASS] accepted Supplementary Figure S5 HIL matches prepared work: {len(rows)} fields"
        )
        return 0

    annotations_target = workdir / "annotations"
    masks_target = workdir / "region_masks"
    for directory in (annotations_target, masks_target):
        if directory.exists() and any(directory.iterdir()):
            raise FileExistsError(f"Refusing non-empty HIL destination: {directory}")
    annotations_target.mkdir(parents=True, exist_ok=True)
    masks_target.mkdir(parents=True, exist_ok=True)

    copied: list[dict[str, Any]] = []
    for row in rows:
        annotation_source = resolve_record_path(row["annotation_path"], source)
        mask_source = resolve_record_path(row["region_mask_path"], source)
        annotation_target = annotations_target / annotation_source.name
        mask_target = masks_target / mask_source.name
        copy_verified(annotation_source, annotation_target, row["annotation_sha256"])
        copy_verified(mask_source, mask_target, row["region_mask_sha256"])
        for target in (annotation_target, mask_target):
            copied.append(
                {
                    "path": target.relative_to(workdir).as_posix(),
                    "bytes": target.stat().st_size,
                    "sha256": sha256(target),
                }
            )

    complete_source = source / "ANNOTATION_COMPLETE.txt"
    complete_target = workdir / "ANNOTATION_COMPLETE.txt"
    if complete_source.is_file():
        copy_verified(complete_source, complete_target, sha256(complete_source))

    import_receipt = {
        "schema": "figs5_hil_import_receipt_v1",
        "imported_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_hil": str(source),
        "source_hil_receipt_sha256": sha256(source / "HIL_RECEIPT.json"),
        "source_hil_manifest_sha256": sha256(source / "HIL_MANIFEST.csv"),
        "source_status": receipt["status"],
        "workdir": ".",
        "analysis_samplesheet_sha256": sha256(workdir / "analysis_samplesheet.csv"),
        "accepted_fields": len(rows),
        "prepared_dapi_mip_validation_modes": mip_validation_modes,
        "copied_files": copied,
        "raw_images_modified": False,
        "raw_cellpose_masks_modified": False,
    }
    receipt_target = workdir / "HIL_IMPORT_RECEIPT.json"
    if receipt_target.exists() or receipt_target.is_symlink():
        raise FileExistsError(f"Refusing to overwrite import receipt: {receipt_target}")
    incoming = receipt_target.with_name(receipt_target.name + ".incoming")
    if incoming.exists() or incoming.is_symlink():
        raise FileExistsError(f"Stale incoming path blocks import receipt: {incoming}")
    incoming.write_text(
        json.dumps(import_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(incoming, receipt_target)
    print(
        f"[PASS] imported {len(rows)} accepted Supplementary Figure S5 HIL pairs into {workdir}; "
        f"receipt: {receipt_target}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
