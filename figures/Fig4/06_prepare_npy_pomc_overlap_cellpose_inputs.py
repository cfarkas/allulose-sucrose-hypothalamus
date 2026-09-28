#!/usr/bin/env python3
"""Prepare the two July 2026 Figure 4 NPY/POMC overlap animals for Cellpose.

The source CZI files and exported OME-TIFF Z-stacks are immutable inputs. This
script validates the user-confirmed channel map, makes one native-resolution
maximum-intensity projection per channel, and writes a SHA-256-bound manifest.
It never overwrites an existing output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile


SCRIPT_VERSION = "1.0.0"
MARKER_ORDER = ("DAPI", "NPY", "POMC", "cFOS")
EXPECTED_ANIMALS = ("FR7-5", "FR8-1")


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def safe_token(value: Any) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return re.sub(r"_+", "_", token).strip("_") or "unnamed"


def maximum_projection(path: Path) -> np.ndarray:
    array = np.squeeze(np.asarray(tifffile.imread(str(path))))
    if array.ndim == 2:
        return array
    while array.ndim > 2:
        array = np.max(array, axis=0)
    if array.ndim != 2:
        raise ValueError(f"Could not make a 2-D projection from {path}: {array.shape}")
    return np.asarray(array)


def write_csv_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to replace manifest: {path}")
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    temporary = path.with_name(path.name + ".part")
    if temporary.exists():
        raise FileExistsError(f"Stale partial manifest blocks safe write: {temporary}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main() -> int:
    figure_dir = Path(__file__).resolve().parent
    extension = figure_dir / "raw_data" / "npy_pomc_overlap_extension_20260706"
    parser = argparse.ArgumentParser(
        description="Create native-resolution Cellpose MIPs for the Figure 4 NPY/POMC extension."
    )
    parser.add_argument("--channel-dir", type=Path, default=extension / "channel_tiffs")
    parser.add_argument("--output-dir", type=Path, default=extension / "cellpose_inputs")
    args = parser.parse_args()

    channel_dir = args.channel_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    channel_map_path = channel_dir / "channel_map.csv"
    if not channel_map_path.is_file():
        raise FileNotFoundError(channel_map_path)
    channel_map = pd.read_csv(channel_map_path)
    required = {
        "acquisition_id", "animal_id", "genotype", "condition", "marker",
        "channel_map_status", "output_tiff", "include_channel",
    }
    missing = required - set(channel_map.columns)
    if missing:
        raise ValueError(f"channel_map.csv missing columns: {sorted(missing)}")

    selected = channel_map[
        channel_map["include_channel"].astype(str).str.upper().isin({"TRUE", "1"})
    ].copy()
    if set(selected["animal_id"].astype(str)) != set(EXPECTED_ANIMALS):
        raise ValueError(
            f"Expected exactly {EXPECTED_ANIMALS}; found {sorted(set(selected['animal_id'].astype(str)))}"
        )
    if set(selected["genotype"].astype(str)) != {"NPY-Tg"}:
        raise ValueError("Both extension animals must be confirmed NPY-Tg")
    if set(selected["marker"].astype(str)) != set(MARKER_ORDER):
        raise ValueError(
            f"Expected DAPI/NPY/POMC/cFOS for each animal; found {sorted(set(selected['marker']))}"
        )
    if not selected["channel_map_status"].astype(str).str.startswith("CONFIRMED").all():
        raise ValueError("Every selected channel-map row must be confirmed")
    counts = selected.groupby("animal_id")["marker"].nunique()
    if not (counts == len(MARKER_ORDER)).all():
        raise ValueError(f"Incomplete marker set per animal: {counts.to_dict()}")

    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    order = {marker: index for index, marker in enumerate(MARKER_ORDER)}
    selected["_marker_order"] = selected["marker"].map(order)
    selected = selected.sort_values(["animal_id", "_marker_order"])
    for record in selected.to_dict("records"):
        source = channel_dir / str(record["output_tiff"])
        if not source.is_file():
            raise FileNotFoundError(source)
        marker = str(record["marker"])
        acquisition = str(record["acquisition_id"])
        output = output_dir / f"{safe_token(acquisition)}_{safe_token(marker)}_mip.tif"
        if output.exists():
            raise FileExistsError(f"Refusing to replace Cellpose input: {output}")
        temporary = output.with_name(output.name + ".part.tif")
        if temporary.exists():
            raise FileExistsError(f"Stale partial image blocks safe write: {temporary}")
        mip = maximum_projection(source)
        tifffile.imwrite(
            str(temporary),
            mip,
            photometric="minisblack",
            metadata={
                "axes": "YX",
                "source_zstack": str(source),
                "marker": marker,
                "pipeline": SCRIPT_VERSION,
            },
        )
        check = np.squeeze(np.asarray(tifffile.imread(str(temporary))))
        if check.shape != mip.shape or check.dtype != mip.dtype or not np.array_equal(check, mip):
            raise RuntimeError(f"Projection validation failed: {temporary}")
        os.replace(temporary, output)
        rows.append({
            "pipeline_version": SCRIPT_VERSION,
            "acquisition_id": acquisition,
            "sample_id": record["sample_id"],
            "animal_id": record["animal_id"],
            "genotype": record["genotype"],
            "condition": record["condition"],
            "analysis_scope": "solution_blind_NPY_POMC_overlap_extension",
            "marker": marker,
            "source_zstack": str(source),
            "source_zstack_sha256": sha256_file(source),
            "source_shape_zyx": "x".join(map(str, tifffile.imread(str(source)).shape)),
            "cellpose_input": str(output),
            "cellpose_input_sha256": sha256_file(output),
            "cellpose_input_shape_yx": "x".join(map(str, mip.shape)),
            "dtype": str(mip.dtype),
            "expected_cellpose_sidecar": str(output.with_name(output.name[:-4] + "_seg.npy")),
        })
        print(f"[MIP] {record['animal_id']} {marker}: {output}", flush=True)

    manifest = output_dir / "cellpose_input_manifest.csv"
    write_csv_atomic(manifest, rows)
    print(f"[DONE] {len(rows)}/8 Cellpose inputs; manifest={manifest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
