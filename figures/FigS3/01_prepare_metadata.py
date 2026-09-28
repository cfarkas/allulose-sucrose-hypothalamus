#!/usr/bin/env python3
"""Convert the source Excel nomenclature into strict, analysis-ready CSV files."""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from figs3_common import (
    DEFAULT_SOURCE_ROOT,
    FIGS3_DIR,
    ORGANS,
    FigS3Error,
    atomic_csv,
    atomic_json,
    file_identity,
    require_directory,
    require_regular_file,
    sha256_file,
    validate_sample_id,
)


EXPECTED_COLUMNS = ("GRUPO", "CONDICIÓN", "CÓDIGO", "ÓRGANO", "ID")
ORGAN_MAP = {"RINON": "Kidney", "HIGADO": "Liver", "BAZO": "Spleen"}
TREATMENT_MAP = {"SACAROSA": "Sucrose", "ALULOSA": "Allulose", "AGUA": "Water"}


def ascii_upper(value: object) -> str:
    text = unicodedata.normalize("NFKD", str(value).strip())
    return "".join(char for char in text if not unicodedata.combining(char)).upper()


def safe_code(value: object) -> str:
    text = ascii_upper(value).lower().replace("_", "-")
    text = re.sub(r"[^a-z0-9-]+", "-", text).strip("-")
    if not text:
        raise FigS3Error(f"Cannot sanitize experimental code: {value!r}")
    return text


def normalize_treatment(raw: object) -> tuple[str, str]:
    value = ascii_upper(raw)
    base = value.split("_", 1)[0]
    if base not in TREATMENT_MAP:
        raise FigS3Error(f"Unknown treatment in Excel: {raw!r}")
    note = "brain_not_collected" if "SIN_CEREBRO" in value else ""
    return TREATMENT_MAP[base], note


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=FIGS3_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_root = require_directory(args.source_root, "source WSI directory")
    output_root = args.output_root.expanduser().resolve()
    excel = require_regular_file(source_root / "EXCEL_ORGANOS.xlsx", "metadata Excel")
    table = pd.read_excel(excel, sheet_name="Hoja1", dtype=str)
    columns = tuple(str(item).strip() for item in table.columns)
    if columns != EXPECTED_COLUMNS:
        raise FigS3Error(f"Unexpected Excel columns: {columns!r}")
    table.columns = EXPECTED_COLUMNS
    for column in ("GRUPO", "CONDICIÓN", "CÓDIGO", "ID"):
        table[column] = table[column].ffill()
    if table.isna().any().any():
        raise FigS3Error("Excel still contains missing values after grouped-row forward fill")

    organ_rows: list[dict[str, object]] = []
    for source_row, row in table.iterrows():
        excel_id = str(row["ID"]).strip()
        sample_id = validate_sample_id(excel_id)
        folder_name = sample_id.upper()
        sample_dir = require_directory(source_root / folder_name, "sample directory")
        treatment, sample_note = normalize_treatment(row["GRUPO"])
        organ_key = ascii_upper(row["ÓRGANO"])
        if organ_key not in ORGAN_MAP:
            raise FigS3Error(f"Unknown organ in Excel row {source_row + 2}: {row['ÓRGANO']!r}")
        cohort = ascii_upper(row["CONDICIÓN"])
        if cohort not in {"NPY", "POMC"}:
            raise FigS3Error(f"Unexpected cohort/condition: {cohort!r}")
        organ_rows.append(
            {
                "sample_id": sample_id,
                "source_folder": folder_name,
                "excel_id_raw": excel_id,
                "experimental_code": safe_code(row["CÓDIGO"]),
                "experimental_code_raw": str(row["CÓDIGO"]).strip(),
                "treatment": treatment,
                "treatment_raw": str(row["GRUPO"]).strip(),
                "cohort": cohort,
                "cohort_raw": str(row["CONDICIÓN"]).strip(),
                "organ": ORGAN_MAP[organ_key],
                "organ_raw": str(row["ÓRGANO"]).strip(),
                "sample_note": sample_note,
                "biological_unit": "animal",
                "source_excel_row": int(source_row + 2),
                "source_mds_path": str((sample_dir / "1.mds").resolve()),
                "source_overview_path": str((sample_dir / "3.jpg").resolve()),
                "source_info_ini_path": str((sample_dir / "info.ini").resolve()),
            }
        )
        for required_name in ("1.mds", "3.jpg", "info.ini"):
            require_regular_file(sample_dir / required_name, required_name)

    by_sample: dict[str, list[dict[str, object]]] = {}
    for row in organ_rows:
        by_sample.setdefault(str(row["sample_id"]), []).append(row)
    if len(by_sample) != 24:
        raise FigS3Error(f"Expected 24 animals, found {len(by_sample)}")
    expected_ids = {f"m26-{number:03d}" for number in range(1, 25)}
    if set(by_sample) != expected_ids:
        raise FigS3Error("Sample IDs are not the complete m26-001 through m26-024 series")

    slide_rows: list[dict[str, object]] = []
    for sample_id in sorted(by_sample):
        rows = by_sample[sample_id]
        if len(rows) != 3 or {str(row["organ"]) for row in rows} != set(ORGANS):
            raise FigS3Error(f"{sample_id} does not have exactly Kidney, Liver, and Spleen")
        invariant_fields = (
            "source_folder",
            "excel_id_raw",
            "experimental_code",
            "experimental_code_raw",
            "treatment",
            "treatment_raw",
            "cohort",
            "cohort_raw",
            "sample_note",
            "source_mds_path",
            "source_overview_path",
            "source_info_ini_path",
        )
        for field in invariant_fields:
            if len({str(row[field]) for row in rows}) != 1:
                raise FigS3Error(f"Slide-level field {field} varies within {sample_id}")
        first = rows[0]
        mds_identity = file_identity(Path(str(first["source_mds_path"])))
        slide_rows.append(
            {
                "sample_id": sample_id,
                "source_folder": first["source_folder"],
                "excel_id_raw": first["excel_id_raw"],
                "experimental_code": first["experimental_code"],
                "experimental_code_raw": first["experimental_code_raw"],
                "treatment": first["treatment"],
                "treatment_raw": first["treatment_raw"],
                "cohort": first["cohort"],
                "cohort_raw": first["cohort_raw"],
                "sample_note": first["sample_note"],
                "biological_unit": "animal",
                "organs_in_same_wsi": "Kidney|Liver|Spleen",
                "source_mds_path": first["source_mds_path"],
                "source_mds_size_bytes": mds_identity["size_bytes"],
                "source_overview_path": first["source_overview_path"],
                "source_info_ini_path": first["source_info_ini_path"],
            }
        )

    treatment_counts = Counter(str(row["treatment"]) for row in slide_rows)
    if treatment_counts != Counter({"Sucrose": 9, "Allulose": 8, "Water": 7}):
        raise FigS3Error(f"Unexpected treatment counts: {dict(treatment_counts)}")

    source_data = output_root / "source_data"
    raw_data = output_root / "raw_data"
    provenance = output_root / "provenance"
    atomic_csv(source_data / "Figure_S3_slide_manifest.csv", slide_rows)
    atomic_csv(source_data / "Figure_S3_sample_organ_manifest.csv", organ_rows)
    audit_rows = [
        {
            "field": "treatment",
            "raw_value": raw,
            "standardized_value": standardized,
            "rule": "trim, remove accents, map Spanish label; AGUA_Sin_cerebro retains note",
        }
        for raw, standardized in sorted(
            {(str(row["treatment_raw"]), str(row["treatment"])) for row in slide_rows}
        )
    ]
    audit_rows.extend(
        {
            "field": "organ",
            "raw_value": raw,
            "standardized_value": standardized,
            "rule": "trim, remove accents, map Spanish organ label",
        }
        for raw, standardized in sorted(
            {(str(row["organ_raw"]), str(row["organ"])) for row in organ_rows}
        )
    )
    atomic_csv(source_data / "Figure_S3_nomenclature_audit.csv", audit_rows)

    raw_manifest = []
    for row in slide_rows:
        for kind, key in (
            ("motic_mds_wsi", "source_mds_path"),
            ("scanner_overview", "source_overview_path"),
            ("scanner_info", "source_info_ini_path"),
        ):
            identity = file_identity(Path(str(row[key])))
            raw_manifest.append(
                {
                    "sample_id": row["sample_id"],
                    "asset_type": kind,
                    "source_path": identity["path"],
                    "size_bytes": identity["size_bytes"],
                    "sha256": "computed_by_export_step" if kind == "motic_mds_wsi" else sha256_file(Path(identity["path"])),
                    "copied_into_figs3": "false",
                }
            )
    raw_manifest.append(
        {
            "sample_id": "cohort",
            "asset_type": "metadata_excel",
            "source_path": str(excel),
            "size_bytes": excel.stat().st_size,
            "sha256": sha256_file(excel),
            "copied_into_figs3": "false",
        }
    )
    atomic_csv(raw_data / "RAW_DATA_MANIFEST.csv", raw_manifest)
    atomic_json(
        provenance / "metadata_receipt.json",
        {
            "status": "complete",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_root": str(source_root),
            "source_excel": str(excel),
            "source_excel_sha256": sha256_file(excel),
            "slide_count": len(slide_rows),
            "sample_organ_row_count": len(organ_rows),
            "treatment_counts": dict(sorted(treatment_counts.items())),
            "cohort_counts": dict(sorted(Counter(str(row["cohort"]) for row in slide_rows).items())),
            "organs_per_wsi": list(ORGANS),
            "normalization_note": "Raw labels are retained beside standardized English labels.",
        },
    )
    print(json.dumps({"slides": 24, "organ_rows": 72, "treatments": dict(treatment_counts)}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
