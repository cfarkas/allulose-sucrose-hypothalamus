#!/usr/bin/env python3
"""Export Motic MDS L0/L2 pixels as tiled BigTIFF with resumable receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from typing import Iterator

import numpy as np
import olefile
import tifffile
from PIL import Image

from figs3_common import (
    FIGS3_DIR,
    FigS3Error,
    atomic_csv,
    atomic_json,
    read_csv,
    require_regular_file,
    sha256_file,
    validate_sample_id,
)


CONVERTER_VERSION = "figs3-mds-export-1.0"


@dataclass(frozen=True)
class Level:
    index: int
    name: str
    rows: int
    columns: int
    tile_width: int
    tile_height: int

    @property
    def width(self) -> int:
        return self.columns * self.tile_width

    @property
    def height(self) -> int:
        return self.rows * self.tile_height


def tile_coordinate(name: str) -> tuple[int, int] | None:
    parts = name.split("_")
    if len(parts) != 2:
        return None
    try:
        row, column = map(int, parts)
    except ValueError:
        return None
    return (row, column) if row >= 0 and column >= 0 else None


def level_key(name: str) -> tuple[int, float | str]:
    try:
        return 0, -float(name)
    except ValueError:
        return 1, name.casefold()


class MdsPixels:
    def __init__(self, path: Path) -> None:
        self.path = require_regular_file(path, "MDS input")
        try:
            self.ole = olefile.OleFileIO(str(self.path))
        except Exception as exc:
            raise FigS3Error(f"Cannot open MDS OLE file: {self.path}") from exc
        self._streams: dict[str, dict[tuple[int, int], tuple[str, ...]]] = {}
        self.levels = self._discover()

    def close(self) -> None:
        self.ole.close()

    def __enter__(self) -> "MdsPixels":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _discover(self) -> tuple[Level, ...]:
        grouped: dict[str, dict[tuple[int, int], tuple[str, ...]]] = {}
        for parts in self.ole.listdir(streams=True, storages=False):
            if len(parts) != 3 or parts[0] != "DSI0":
                continue
            coordinate = tile_coordinate(parts[2])
            if coordinate is not None:
                grouped.setdefault(parts[1], {})[coordinate] = tuple(parts)
        if not grouped:
            raise FigS3Error(f"No DSI0 pixel levels in {self.path}")
        discovered: list[Level] = []
        for index, name in enumerate(sorted(grouped, key=level_key)):
            streams = grouped[name]
            encoded = self.ole.openstream(list(streams[min(streams)])).read()
            with Image.open(BytesIO(encoded)) as image:
                width, height = image.size
            rows = max(item[0] for item in streams) + 1
            columns = max(item[1] for item in streams) + 1
            self._streams[name] = streams
            discovered.append(Level(index, name, rows, columns, width, height))
        return tuple(discovered)

    def tiles(self, level: Level) -> Iterator[np.ndarray]:
        streams = self._streams[level.name]
        shape = (level.tile_height, level.tile_width, 3)
        for row in range(level.rows):
            for column in range(level.columns):
                parts = streams.get((row, column))
                if parts is None:
                    yield np.full(shape, 255, dtype=np.uint8)
                    continue
                try:
                    encoded = self.ole.openstream(list(parts)).read()
                    with Image.open(BytesIO(encoded)) as image:
                        array = np.asarray(image.convert("RGB"), dtype=np.uint8)
                except Exception as exc:
                    raise FigS3Error(
                        f"Cannot decode {level.name}/{row}_{column} in {self.path}"
                    ) from exc
                if array.shape == shape:
                    yield array
                else:
                    padded = np.full(shape, 255, dtype=np.uint8)
                    height = min(shape[0], array.shape[0])
                    width = min(shape[1], array.shape[1])
                    padded[:height, :width] = array[:height, :width, :3]
                    yield padded


def scanner_values(path: Path) -> dict[str, float | int | str]:
    result: dict[str, str] = {}
    text = require_regular_file(path, "scanner info.ini").read_text(
        encoding="utf-8", errors="replace"
    )
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip().casefold()] = value.strip()
    try:
        return {
            "source_mpp": float(result["scale"]),
            "true_width": int(result["mifwidth"]),
            "true_height": int(result["mifheight"]),
            "objective": float(result["scanlens"]),
            "scanner": result.get("scanmachineinfo", ""),
            "acquired": result.get("createtimetext", ""),
        }
    except (KeyError, ValueError) as exc:
        raise FigS3Error(f"Missing scanner geometry in {path}") from exc


def digest_with_stability(path: Path) -> tuple[int, str]:
    before = path.stat()
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    after = path.stat()
    identity = lambda stat: (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )
    if identity(before) != identity(after):
        raise FigS3Error(f"Source changed while hashing: {path}")
    return size, digest.hexdigest()


def rational_float(value: object) -> float:
    if isinstance(value, tuple):
        return float(value[0]) / float(value[1])
    return float(value)


def validate_tiff(path: Path, width: int, height: int, mpp: float) -> None:
    if path.is_symlink() or not path.is_file():
        raise FigS3Error(f"Missing/unsafe TIFF: {path}")
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        if (int(page.imagewidth), int(page.imagelength)) != (width, height):
            raise FigS3Error(f"TIFF dimension mismatch: {path}")
        if not page.is_tiled or int(page.samplesperpixel) != 3:
            raise FigS3Error(f"TIFF is not tiled RGB: {path}")
        unit = int(page.tags["ResolutionUnit"].value)
        xres = rational_float(page.tags["XResolution"].value)
        measured = 10_000.0 / xres if unit == 3 else math.nan
        if not math.isclose(measured, mpp, rel_tol=3e-6, abs_tol=3e-6):
            raise FigS3Error(f"TIFF MPP mismatch: {path}: {measured} != {mpp}")


def export_level(
    slide: MdsPixels,
    level: Level,
    path: Path,
    source_mpp: float,
    compression_level: int,
) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        raise FigS3Error(f"Refusing symlink output: {path}")
    scale = float(slide.levels[0].name) / float(level.name)
    output_mpp = source_mpp * scale
    resolution = 10_000.0 / output_mpp
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        tifffile.imwrite(
            tmp,
            data=slide.tiles(level),
            shape=(level.height, level.width, 3),
            dtype=np.uint8,
            tile=(level.tile_height, level.tile_width),
            photometric="rgb",
            planarconfig="contig",
            compression="deflate",
            compressionargs={"level": compression_level},
            bigtiff=True,
            metadata=None,
            resolution=(resolution, resolution),
            resolutionunit="CENTIMETER",
            maxworkers=4,
        )
        validate_tiff(tmp, level.width, level.height, output_mpp)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    size, sha = digest_with_stability(path)
    return {
        "level": level.index,
        "internal_level": level.name,
        "output_path": str(path.resolve()),
        "width": level.width,
        "height": level.height,
        "tile_width": level.tile_width,
        "tile_height": level.tile_height,
        "output_mpp": output_mpp,
        "output_size_bytes": size,
        "output_sha256": sha,
    }


def receipt_valid(
    receipt: dict[str, object], source: Path, output_dir: Path, verify_sha: bool
) -> bool:
    if receipt.get("status") != "complete" or receipt.get("converter_version") != CONVERTER_VERSION:
        return False
    if receipt.get("source_mds_path") != str(source.resolve()):
        return False
    if int(receipt.get("source_mds_size_bytes", -1)) != source.stat().st_size:
        return False
    levels = receipt.get("levels")
    if not isinstance(levels, list) or {int(row["level"]) for row in levels} != {0, 2}:
        return False
    for row in levels:
        path = Path(str(row["output_path"]))
        try:
            path.resolve().relative_to(output_dir.resolve())
            validate_tiff(path, int(row["width"]), int(row["height"]), float(row["output_mpp"]))
        except (FigS3Error, ValueError, OSError, KeyError):
            return False
        if path.stat().st_size != int(row["output_size_bytes"]):
            return False
        if verify_sha and sha256_file(path) != str(row["output_sha256"]):
            return False
    if verify_sha and sha256_file(source) != str(receipt.get("source_mds_sha256")):
        return False
    return True


def worker(job: dict[str, object]) -> dict[str, object]:
    sample_id = validate_sample_id(str(job["sample_id"]))
    source = Path(str(job["source_mds_path"]))
    info_path = Path(str(job["source_info_ini_path"]))
    output_dir = Path(str(job["output_dir"]))
    sample_dir = output_dir / sample_id
    receipt_path = sample_dir / "conversion_receipt.json"
    resume = bool(job["resume"])
    verify_sha = bool(job["verify_sha"])
    try:
        if receipt_path.is_file() and resume:
            existing = json.loads(receipt_path.read_text(encoding="utf-8"))
            if receipt_valid(existing, source, output_dir, verify_sha):
                existing["run_status"] = "verified_existing"
                return existing
            raise FigS3Error(f"Untrusted/incompatible prior export for {sample_id}")
        if sample_dir.exists() and any(sample_dir.iterdir()):
            raise FigS3Error(f"Output exists without a reusable receipt: {sample_dir}")
        sample_dir.mkdir(parents=True, exist_ok=True)
        scanner = scanner_values(info_path)
        source_size, source_sha = digest_with_stability(source)
        with MdsPixels(source) as slide:
            if len(slide.levels) <= 2:
                raise FigS3Error(f"{sample_id} has fewer than three MDS levels")
            level_rows = [
                export_level(
                    slide,
                    slide.levels[index],
                    sample_dir / f"1_L{index}_rgb.tif",
                    float(scanner["source_mpp"]),
                    int(job["compression_level"]),
                )
                for index in (0, 2)
            ]
            level_count = len(slide.levels)
            all_geometry = [
                {"level": level.index, "name": level.name, "width": level.width, "height": level.height}
                for level in slide.levels
            ]
        receipt = {
            "status": "complete",
            "run_status": "exported",
            "converter_version": CONVERTER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_id": sample_id,
            "source_mds_path": str(source.resolve()),
            "source_mds_size_bytes": source_size,
            "source_mds_sha256": source_sha,
            "scanner_info_path": str(info_path.resolve()),
            "scanner_info_sha256": sha256_file(info_path),
            "source_mpp": scanner["source_mpp"],
            "objective": scanner["objective"],
            "true_width": scanner["true_width"],
            "true_height": scanner["true_height"],
            "padded_width": level_rows[0]["width"],
            "padded_height": level_rows[0]["height"],
            "scanner": scanner["scanner"],
            "acquired": scanner["acquired"],
            "source_level_count": level_count,
            "all_level_geometry": all_geometry,
            "levels": level_rows,
            "compression": "deflate",
            "compression_level": int(job["compression_level"]),
        }
        atomic_json(receipt_path, receipt)
        return receipt
    except Exception as exc:
        failure = {
            "status": "failed",
            "converter_version": CONVERTER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_id": sample_id,
            "source_mds_path": str(source),
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        sample_dir.mkdir(parents=True, exist_ok=True)
        atomic_json(sample_dir / "conversion_failure.json", failure)
        return failure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-manifest",
        type=Path,
        default=FIGS3_DIR / "source_data" / "Figure_S3_slide_manifest.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=FIGS3_DIR / "exports")
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--compression-level", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fast-resume", action="store_true", help="Skip full SHA verification only on resume")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers < 1 or args.workers > 4:
        raise FigS3Error("--workers must be between 1 and 4")
    if not 0 <= args.compression_level <= 9:
        raise FigS3Error("--compression-level must be between 0 and 9")
    rows = read_csv(require_regular_file(args.slide_manifest, "slide manifest"))
    requested = {validate_sample_id(value) for value in args.sample_id}
    if requested:
        rows = [row for row in rows if row["sample_id"] in requested]
        missing = requested - {row["sample_id"] for row in rows}
        if missing:
            raise FigS3Error(f"Requested sample IDs absent from manifest: {sorted(missing)}")
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        {
            **row,
            "output_dir": str(output_dir),
            "resume": args.resume,
            "verify_sha": not args.fast_resume,
            "compression_level": args.compression_level,
        }
        for row in rows
    ]
    if args.dry_run:
        print(json.dumps({"status": "planned", "slides": [row["sample_id"] for row in rows]}, indent=2))
        return 0
    results: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(worker, job): str(job["sample_id"]) for job in jobs}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps({"sample_id": result["sample_id"], "status": result["status"], "run_status": result.get("run_status")}), flush=True)
    results.sort(key=lambda item: str(item["sample_id"]))
    manifest_rows: list[dict[str, object]] = []
    for result in results:
        if result["status"] != "complete":
            continue
        for level in result["levels"]:
            manifest_rows.append(
                {
                    "sample_id": result["sample_id"],
                    "level": level["level"],
                    "internal_level": level["internal_level"],
                    "output_path": level["output_path"],
                    "width": level["width"],
                    "height": level["height"],
                    "output_mpp": level["output_mpp"],
                    "output_size_bytes": level["output_size_bytes"],
                    "output_sha256": level["output_sha256"],
                    "source_mds_sha256": result["source_mds_sha256"],
                    "true_width_l0": result["true_width"],
                    "true_height_l0": result["true_height"],
                    "run_status": result["run_status"],
                }
            )
    atomic_csv(output_dir / "Figure_S3_L0_L2_export_manifest.csv", manifest_rows)
    atomic_json(
        output_dir / "Figure_S3_L0_L2_export_summary.json",
        {
            "status": "complete" if len(results) == len(rows) and all(row["status"] == "complete" for row in results) else "failed",
            "converter_version": CONVERTER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "requested_slides": len(rows),
            "complete_slides": sum(row["status"] == "complete" for row in results),
            "failed_slides": [row["sample_id"] for row in results if row["status"] != "complete"],
            "levels": [0, 2],
            "results": results,
        },
    )
    failed = [row for row in results if row["status"] != "complete"]
    if failed:
        for row in failed:
            print(f"ERROR {row['sample_id']}: {row['error']}", file=sys.stderr)
        return 2
    print(json.dumps({"status": "complete", "slides": len(results), "levels": [0, 2]}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
