#!/usr/bin/env python3
"""Copy the 2026 Water controls and their companions into Fig3/raw_data.

``Fig3/raw_data`` is the figure's raw-data authority: every payload there is an
independent copy with its own inode, a destination link count of one, and a
recorded SHA-256.  This command reproduces that guarantee for the 24 September
2026 Zeiss LSM 780 acquisitions (specimens N132-5 and N132-3), their analytical
channel TIFFs, their Cellpose masks and their receipts, and it refreshes the
matching manifest rows.

Only this dataset's manifest rows are rewritten; every other row is carried
over untouched, and the previous manifest text is preserved beside the copies
before anything is replaced.  Nothing in ``Fig3/raw`` is moved or modified.
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

PAPER_ROOT = next(
    path
    for path in Path(__file__).resolve().parents
    if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir() and (path / "README.txt").is_file()))
)
DATASET = "new_water_lsm780_2026_09_24"
FOLDER = "new_water_control_LSM780_24_SEPT_2026"
RAW = PAPER_ROOT / "Fig3/raw" / FOLDER
RAW_DATA = PAPER_ROOT / "Fig3/raw_data"
ANIMALS = {
    "WATER_NPY4": {"label": "N132-5", "czi": "N132-5 CTRL NPY GFP 555 647 25X - 001.czi"},
    "WATER_NPY5": {"label": "N132-3", "czi": "N132-3 CTRL NPY GFP 555 647 25X - 001.czi"},
}
ANALYTICAL = ("dapi", "GFP", "fos")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def independent_copy(source: Path, destination: Path) -> dict:
    """Copy with a fresh inode and verify the payload byte-for-byte."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected = sha256_file(source)
    staged = destination.with_name(destination.name + ".copying")
    shutil.copy2(source, staged)
    os.replace(staged, destination)
    source_stat, destination_stat = source.stat(), destination.stat()
    observed = sha256_file(destination)
    if observed != expected:
        raise RuntimeError(f"Copy hash mismatch: {destination}")
    if source_stat.st_size != destination_stat.st_size:
        raise RuntimeError(f"Copy size mismatch: {destination}")
    if (source_stat.st_dev, source_stat.st_ino) == (destination_stat.st_dev, destination_stat.st_ino):
        raise RuntimeError(f"Copy shares the source inode: {destination}")
    if destination_stat.st_nlink != 1:
        raise RuntimeError(f"Copy is not a single-link regular file: {destination}")
    return dict(
        dataset=DATASET,
        source_absolute_path=str(source),
        figure_relative_path=str(destination.relative_to(RAW_DATA)),
        bytes=source_stat.st_size,
        source_device=source_stat.st_dev,
        source_inode=source_stat.st_ino,
        destination_device=destination_stat.st_dev,
        destination_inode=destination_stat.st_ino,
        destination_link_count=destination_stat.st_nlink,
        independent_inode_verified="yes",
        sha256=observed,
    )


def refresh_rows(path: Path, rows: list[dict], key: str, owned_prefix: str, preserve_into: Path) -> dict:
    """Replace this dataset's rows and carry every other manifest row over."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        existing = list(reader)
    if not fieldnames:
        raise RuntimeError(f"Manifest has no header: {path}")
    for row in rows:
        unknown = set(row) - set(fieldnames)
        if unknown:
            raise RuntimeError(f"Unknown manifest columns for {path.name}: {sorted(unknown)}")
    kept = [row for row in existing if not str(row.get(key, "")).startswith(owned_prefix)]
    preserve_into.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    shutil.copy2(path, preserve_into / f"{path.name}.before_{stamp}")
    staged = path.with_name(path.name + ".refreshing")
    with staged.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in kept + rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
    os.replace(staged, path)
    return dict(rows_written=len(rows), rows_replaced=len(existing) - len(kept), rows_carried_over=len(kept))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-dir", type=Path, default=RAW)
    parser.add_argument("--raw-data-dir", type=Path, default=RAW_DATA)
    args = parser.parse_args()
    raw, raw_data = args.raw_dir, args.raw_data_dir
    derived = raw / "derived"
    for animal, spec in ANIMALS.items():
        if not (raw / spec["czi"]).is_file():
            raise RuntimeError(f"Source acquisition is missing for {animal}: {raw / spec['czi']}")
    if not derived.is_dir():
        raise RuntimeError(f"Run 00c/00d first; derived outputs are missing: {derived}")

    payloads: list[tuple[Path, Path]] = [
        (raw / spec["czi"], raw_data / FOLDER / spec["czi"]) for spec in ANIMALS.values()
    ]
    for source in sorted(derived.rglob("*")):
        if source.is_file():
            payloads.append((source, raw_data / FOLDER / "derived" / source.relative_to(derived)))
    records = [independent_copy(source, destination) for source, destination in payloads]
    by_path = {record["figure_relative_path"]: record for record in records}
    preserve_into = raw_data / FOLDER / "superseded_manifests"
    raw_data_summary = refresh_rows(raw_data / "RAW_DATA_MANIFEST.csv", records,
                                    "figure_relative_path", FOLDER + "/", preserve_into)

    relative = f"{FOLDER}/derived"
    companions: list[dict] = []
    for animal, spec in ANIMALS.items():
        sample_id = animal + "_S01"
        czi_relative = f"{FOLDER}/{spec['czi']}"
        mask_paths = [f"{relative}/{sample_id}_{suffix}_seg.npy" for suffix in ANALYTICAL]
        companions.append(dict(
            microscopy_relative_path=czi_relative,
            microscopy_sha256=by_path[czi_relative]["sha256"],
            microscopy_role="native_2026_water_acquisition",
            cellpose_required_for_analysis="yes",
            companion_mask_relative_paths=";".join(mask_paths),
            companion_mask_sha256=";".join(by_path[name]["sha256"] for name in mask_paths),
            companion_mask_count=len(mask_paths),
            status="PASS_2026_WATER_DERIVED_MASKS_PRESENT",
            note=(f"{animal} (specimen {spec['label']}): DAPI, NPY/GFP and c-FOS masks and their analytical "
                  "TIFFs are present as independent figure-local files; extraction and segmentation receipts "
                  "are hash-bound in RAW_DATA_MANIFEST."),
        ))
        for suffix in ANALYTICAL:
            image = f"{relative}/{sample_id}_{suffix}.tif"
            mask = f"{relative}/{sample_id}_{suffix}_seg.npy"
            companions.append(dict(
                microscopy_relative_path=image,
                microscopy_sha256=by_path[image]["sha256"],
                microscopy_role="analytical_channel",
                cellpose_required_for_analysis="yes",
                companion_mask_relative_paths=mask,
                companion_mask_sha256=by_path[mask]["sha256"],
                companion_mask_count=1,
                status="PASS_EXACT_SAME_STEM_MASK",
                note="Cellpose mask produced with the models and settings recorded for the original nine animals.",
            ))
        display = f"{relative}/{sample_id}_NeuN_display_only.tif"
        companions.append(dict(
            microscopy_relative_path=display,
            microscopy_sha256=by_path[display]["sha256"],
            microscopy_role="display_only_channel",
            cellpose_required_for_analysis="no",
            companion_mask_relative_paths="",
            companion_mask_sha256="",
            companion_mask_count=0,
            status="PASS_DISPLAY_ONLY_NO_MASK_REQUIRED",
            note="NeuN is display-only and excluded from Figure 3 quantification.",
        ))
    companion_summary = refresh_rows(raw_data / "CELLPOSE_COMPANION_MANIFEST.csv", companions,
                                     "microscopy_relative_path", FOLDER + "/", preserve_into)

    receipt = dict(
        schema="figure3-new-water-raw-data-mirror-v1",
        dataset=DATASET,
        animals={animal: spec["label"] for animal, spec in ANIMALS.items()},
        mirrored_utc=datetime.now(timezone.utc).isoformat(),
        payload_count=len(records),
        payload_bytes=sum(record["bytes"] for record in records),
        raw_data_manifest=raw_data_summary,
        cellpose_companion_manifest=companion_summary,
        every_copy_independent_inode=True,
        every_copy_link_count_one=True,
        sources_unmodified=True,
        payloads=records,
        raw_data_manifest_sha256=sha256_file(raw_data / "RAW_DATA_MANIFEST.csv"),
        cellpose_companion_manifest_sha256=sha256_file(raw_data / "CELLPOSE_COMPANION_MANIFEST.csv"),
    )
    (raw_data / FOLDER / "NEW_WATER_RAW_DATA_MIRROR_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in receipt.items() if key != "payloads"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
