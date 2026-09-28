#!/usr/bin/env python3
"""Build an organised, condition-keyed index of the Figure 2 HIL QC plates.

The accepted v9.2 QC plates live in five separate acquisition runs under
Fig2/raw_data/Apotome. This collects them into one tree keyed by condition and
animal, so the panel B inputs are easy to find and audit.

Originals are copied, never moved or modified, and every copy carries the
SHA-256 of its source in the index. Re-running is safe: the output directory is
rebuilt only when --refresh is given, and the run otherwise refuses to write
over an existing collection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
# See the note in 02_stage_analysis.py: the root is found by marker so this runs
# from Fig2/ and from its scripts/Fig2 mirror alike.
PAPER = next((path for path in HERE.parents
              if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                  and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
FIGURE = PAPER / "Fig2"
RAW = FIGURE / "raw_data" / "Apotome"
DEFAULT_OUT = FIGURE / "raw_data" / "human_in_loop_qc"

# Every accepted HIL region-annotation run that contributes to Figure 2.
RUNS = (
    RAW / "24_06_2025/analysis/human_regions_v9",
    RAW / "2025.06.24_human_regions_v9_2",
    RAW / "29.09.2025_human_regions_v9_2",
    RAW / "Aug_2025_NPY/selected_human_regions_v9_2",
    RAW / "Aug_2025_NPY/selected_horizontal_human_regions_v9_2",
)
PANEL_B_TILES = {
    "e10-fr1-1-agua-dapi-pomc488-cfos647-001",
    "e7fr7-1h-sacarosa",
    "e8fr6-4h-alulosa",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def collect() -> list[dict[str, object]]:
    # A few runs nest a copy of another run's plates, so resolve each tile
    # against the union of every run's per-image summary.
    catalogue: dict[str, pd.Series] = {}
    statuses: dict[str, str] = {}
    for run in RUNS:
        if not run.is_dir():
            raise FileNotFoundError(run)
        summary = pd.read_csv(run / "per_image_human_region_summary.csv").drop_duplicates("tile_key")
        for key, row in summary.set_index("tile_key").iterrows():
            catalogue.setdefault(str(key), row)
        status = pd.read_csv(run / "annotation_status.csv")
        for key, value in zip(status.tile_key, status.status):
            statuses.setdefault(str(key), str(value))

    rows: list[dict[str, object]] = []
    for run in RUNS:
        for plate in sorted((run / "qc_human_regions").rglob("*_human_regions_qc.*")):
            tile = plate.name.split("_human_regions_qc")[0]
            tile = tile.split("__")[-1] if "__" in tile else tile
            if tile not in catalogue:
                raise RuntimeError(f"QC plate has no per-image summary row: {plate}")
            row = catalogue[tile]
            rows.append(
                {
                    "condition": str(row.cond),
                    "animal": str(row.animal),
                    "tile_key": tile,
                    "sample": str(row["sample"]),
                    "run": run.relative_to(FIGURE).as_posix(),
                    "annotation_status": statuses.get(tile, "unknown"),
                    "used_in_panel_B": tile in PANEL_B_TILES,
                    "suffix": plate.suffix.lstrip("."),
                    "source_path": str(plate),
                    "source_bytes": plate.stat().st_size,
                    "source_sha256": sha256(plate),
                }
            )
    if not rows:
        raise RuntimeError("No HIL QC plates found")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--refresh", action="store_true", help="Rebuild an existing collection")
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    if outdir.exists():
        if not args.refresh:
            raise FileExistsError(f"Collection already exists; pass --refresh to rebuild: {outdir}")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    rows = collect()
    for row in rows:
        destination = outdir / str(row["condition"]) / f"{row['animal']}__{row['tile_key']}.{row['suffix']}"
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if sha256(destination) == row["source_sha256"]:
                # The same plate is nested in a second run; keep one copy.
                row["organised_path"] = str(destination.relative_to(FIGURE))
                row["duplicate_source"] = True
                continue
            # Two runs rendered the same tile differently; keep both, run-tagged.
            token = Path(str(row["run"])).name
            destination = destination.with_name(
                f"{row['animal']}__{row['tile_key']}__{token}.{row['suffix']}"
            )
            if destination.exists():
                raise RuntimeError(f"Unresolved collision in the organised collection: {destination}")
        shutil.copy2(Path(str(row["source_path"])), destination)
        row["organised_path"] = str(destination.relative_to(FIGURE))
        row["duplicate_source"] = False
        if sha256(destination) != row["source_sha256"]:
            raise RuntimeError(f"Copy verification failed: {destination}")

    index = outdir / "HUMAN_IN_LOOP_QC_INDEX.csv"
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with index.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)
    (outdir / "README.txt").write_text(
        "Figure 2 HIL QC plates, organised by condition and animal.\n"
        "\n"
        "Every file here is a byte-exact copy of an accepted v9.2 QC plate that lives in\n"
        "one of the five acquisition runs under raw_data/Apotome. Nothing was moved or\n"
        "modified; HUMAN_IN_LOOP_QC_INDEX.csv records each source path and its SHA-256,\n"
        "the recorded condition and animal, the annotation status, and whether the tile\n"
        "is one of the three drawn in Figure 2 panel B.\n"
        "\n"
        "Rebuild with:  Fig2/08_organize_raw_qc.py --refresh\n",
        encoding="utf-8",
    )
    conditions = sorted({str(row["condition"]) for row in rows})
    print(f"[OK] {len(rows)} QC plates organised into {outdir} across {conditions}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
