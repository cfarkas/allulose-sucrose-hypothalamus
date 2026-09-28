#!/usr/bin/env python3
"""Apply the audited Liver/Kidney identity correction without rerunning inference.

The neural model was run on fixed pixel tiles and does not consume organ names.
This migration swaps only Kidney/Liver annotations and tile identifiers, verifies
their geometry against the corrected masks, preserves cell polygons/probabilities,
and refreshes every affected receipt hash. It is idempotent.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from PIL import Image, ImageDraw

from figs3_common import FIGS3_DIR, FigS3Error, atomic_csv, atomic_json, sha256_file


MIGRATION_VERSION = "figs3-organ-identity-audit-1.0"
CORRECTED_SEGMENTATION_VERSION = "figs3-three-organ-cv-1.1"
CSV_NAMES = (
    "candidate_tile_manifest.csv",
    "selected_tile_manifest.csv",
    "provisional_tile_pixel_qc.csv",
    "tile_cv_features.csv",
    "phenotype_summary.csv",
    "tile_class_counts.csv",
)
ORGAN_COLORS = {"Kidney": "#0097ff", "Liver": "#22b14c", "Spleen": "#b446c8"}
EXPECTED_NEW_MASK_VALUE = {"Kidney": 2, "Liver": 1, "Spleen": 3}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figs3-root", type=Path, default=FIGS3_DIR)
    return parser.parse_args()


def swap_organ(value: str) -> str:
    return {"Kidney": "Liver", "Liver": "Kidney"}.get(value, value)


def swap_tile_uid(value: str) -> str:
    if not value:
        return value
    return value.replace("_kidney_", "_temporary_").replace("_liver_", "_kidney_").replace("_temporary_", "_liver_")


def atomic_frame(frame: pd.DataFrame, path: Path) -> None:
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    temp = Path(temp_name)
    try:
        frame.to_csv(temp, index=False, lineterminator="\n")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def swap_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "organ" in frame:
        frame["organ"] = frame["organ"].map(swap_organ)
    if "tile_uid" in frame:
        frame["tile_uid"] = frame["tile_uid"].map(swap_tile_uid)
    atomic_frame(frame, path)
    return frame


def swap_cells_parquet(path: Path) -> int:
    table = pq.read_table(path)
    for column_name, transform in (("organ", swap_organ), ("tile_uid", swap_tile_uid)):
        index = table.schema.get_field_index(column_name)
        if index < 0:
            raise FigS3Error(f"Missing {column_name} in {path}")
        field = table.schema.field(index)
        values = [transform(value) if value is not None else None for value in table[column_name].to_pylist()]
        table = table.set_column(index, field, pa.array(values, type=field.type))
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        pq.write_table(table, temp, compression="zstd", write_statistics=True)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()
    return table.num_rows


def geometry_gate(selected: pd.DataFrame, corrected_mask: np.ndarray, true_width: int, true_height: int) -> dict[str, float]:
    fractions: dict[str, float] = {}
    for old_organ in ("Kidney", "Liver", "Spleen"):
        rows = selected[selected.organ == old_organ]
        matches = 0
        for row in rows.itertuples():
            x0 = max(0, int(np.floor(float(row.x_l0) * corrected_mask.shape[1] / true_width)))
            y0 = max(0, int(np.floor(float(row.y_l0) * corrected_mask.shape[0] / true_height)))
            x1 = min(corrected_mask.shape[1], int(np.ceil((float(row.x_l0) + float(row.width_l0)) * corrected_mask.shape[1] / true_width)))
            y1 = min(corrected_mask.shape[0], int(np.ceil((float(row.y_l0) + float(row.height_l0)) * corrected_mask.shape[0] / true_height)))
            values = corrected_mask[y0:y1, x0:x1]
            nonzero = values[values > 0]
            if len(nonzero) and int(np.bincount(nonzero).argmax()) == EXPECTED_NEW_MASK_VALUE[old_organ]:
                matches += 1
        fraction = matches / len(rows) if len(rows) else 0.0
        fractions[old_organ] = fraction
        if fraction < 0.90:
            raise FigS3Error(
                f"Only {fraction:.1%} of old {old_organ} tiles map to the expected corrected mask"
            )
    return fractions


def redraw_tile_overview(
    overview_path: Path,
    selected: pd.DataFrame,
    true_width: int,
    true_height: int,
    output: Path,
) -> None:
    with Image.open(overview_path) as source:
        canvas = source.convert("RGB")
    draw = ImageDraw.Draw(canvas)
    for row in selected.itertuples():
        x0 = float(row.x_l0) * canvas.width / true_width
        y0 = float(row.y_l0) * canvas.height / true_height
        x1 = (float(row.x_l0) + float(row.width_l0)) * canvas.width / true_width
        y1 = (float(row.y_l0) + float(row.height_l0)) * canvas.height / true_height
        draw.rectangle((x0, y0, x1, y1), outline=ORGAN_COLORS[row.organ], width=3)
    temp = output.with_name(f".{output.name}.{os.getpid()}.tmp.png")
    try:
        canvas.save(temp)
        os.replace(temp, output)
    finally:
        if temp.exists():
            temp.unlink()


def output_identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def main() -> int:
    args = parse_args()
    root = args.figs3_root.resolve()
    segmentation_summary_path = root / "organ_segmentation" / "segmentation_summary.json"
    segmentation_summary = json.loads(segmentation_summary_path.read_text(encoding="utf-8"))
    if segmentation_summary.get("version") != CORRECTED_SEGMENTATION_VERSION:
        raise FigS3Error("Corrected segmentation must complete before inference relabeling")
    slide_manifest = pd.read_csv(root / "source_data" / "Figure_S3_slide_manifest.csv").set_index("sample_id")
    inference_summary_path = root / "histoplus" / "inference_summary.json"
    backup = root / "provenance" / "inference_summary_before_organ_identity_audit.json"
    if inference_summary_path.is_file() and not backup.exists():
        shutil.copy2(inference_summary_path, backup)

    audit_rows: list[dict[str, object]] = []
    for sample_id in slide_manifest.index:
        sample_dir = root / "histoplus" / sample_id
        receipt_path = sample_dir / "inference_receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("organ_identity_migration", {}).get("version") == MIGRATION_VERSION:
            audit_rows.append({"sample_id": sample_id, "status": "already_corrected"})
            continue
        selected_old = pd.read_csv(sample_dir / "selected_tile_manifest.csv")
        conversion = json.loads((root / "exports" / sample_id / "conversion_receipt.json").read_text(encoding="utf-8"))
        corrected_mask_path = root / "organ_segmentation" / sample_id / "organ_label_mask.png"
        corrected_mask = np.asarray(Image.open(corrected_mask_path), dtype=np.uint8)
        gate = geometry_gate(
            selected_old,
            corrected_mask,
            int(conversion["true_width"]),
            int(conversion["true_height"]),
        )
        phenotype_old = pd.read_csv(sample_dir / "phenotype_summary.csv")
        totals_old = phenotype_old.groupby("organ").cell_count_in_sampled_tiles.sum().to_dict()

        corrected_frames: dict[str, pd.DataFrame] = {}
        for name in CSV_NAMES:
            corrected_frames[name] = swap_csv(sample_dir / name)
        cells_count = swap_cells_parquet(sample_dir / "cells.parquet")
        selected_new = corrected_frames["selected_tile_manifest.csv"]
        redraw_tile_overview(
            Path(slide_manifest.loc[sample_id, "source_overview_path"]),
            selected_new,
            int(conversion["true_width"]),
            int(conversion["true_height"]),
            sample_dir / "selected_tiles_overview.png",
        )
        phenotype_new = corrected_frames["phenotype_summary.csv"].copy()
        phenotype_new["cell_count_in_sampled_tiles"] = pd.to_numeric(
            phenotype_new["cell_count_in_sampled_tiles"]
        )
        totals_new = phenotype_new.groupby("organ").cell_count_in_sampled_tiles.sum().to_dict()
        if not (
            int(totals_new["Kidney"]) == int(totals_old["Liver"])
            and int(totals_new["Liver"]) == int(totals_old["Kidney"])
            and int(totals_new["Spleen"]) == int(totals_old["Spleen"])
        ):
            raise FigS3Error(f"Count-preservation check failed for {sample_id}")

        receipt["organ_label_mask"] = str(corrected_mask_path.resolve())
        receipt["organ_label_mask_sha256"] = sha256_file(corrected_mask_path)
        receipt["representative_tile_uid"] = swap_tile_uid(str(receipt["representative_tile_uid"]))
        receipt["run_status"] = "organ_identity_relabelled_without_model_recompute"
        receipt["organ_identity_migration"] = {
            "version": MIGRATION_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "correction": "old Kidney -> Liver; old Liver -> Kidney; Spleen unchanged",
            "basis": "full-resolution histology audit; renal glomeruli/tubules and hepatic parenchyma",
            "pixel_tiles_changed": False,
            "model_predictions_changed": False,
            "cell_geometry_probabilities_changed": False,
            "geometry_gate_fraction_by_old_label": gate,
        }
        for key, identity in receipt["outputs"].items():
            path = Path(identity["path"])
            receipt["outputs"][key] = output_identity(path)
        atomic_json(receipt_path, receipt)
        audit_rows.append(
            {
                "sample_id": sample_id,
                "status": "corrected",
                "cells_preserved": cells_count,
                "old_kidney_new_liver_count": int(totals_old["Kidney"]),
                "old_liver_new_kidney_count": int(totals_old["Liver"]),
                "spleen_count": int(totals_old["Spleen"]),
                "kidney_geometry_gate_fraction": gate["Liver"],
                "liver_geometry_gate_fraction": gate["Kidney"],
                "spleen_geometry_gate_fraction": gate["Spleen"],
            }
        )
        print(json.dumps({"sample_id": sample_id, "status": "corrected", "cells": cells_count}), flush=True)

    atomic_csv(root / "provenance" / "Figure_S3_organ_identity_migration_audit.csv", audit_rows)
    atomic_json(
        root / "provenance" / "Figure_S3_organ_identity_migration_receipt.json",
        {
            "status": "complete",
            "version": MIGRATION_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "samples": len(audit_rows),
            "correction": "Liver/Kidney labels swapped; pixels and neural predictions preserved",
            "audit_table": str((root / "provenance" / "Figure_S3_organ_identity_migration_audit.csv").resolve()),
        },
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
