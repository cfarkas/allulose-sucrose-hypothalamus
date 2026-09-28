#!/usr/bin/env python3
"""Stage the Figure 2 HIL-defined anatomy panel from Fig2/raw_data.

For one representative section per condition this rebuilds the accepted v9.2
HIL-defined ARC / ME / VMN annotation over the Cellpose DAPI
segmentation, with strictly DAPI-associated c-FOS-positive nuclei in red.

The HIL result is taken from the archived pipeline outputs rather than
re-derived: the anatomy comes from the accepted annotation polygons and region
mask, and the nucleus-to-region assignment and c-FOS positivity come from
per_nucleus_human_region_assignments.csv. Those three artefacts are
cross-checked against per_image_human_region_summary.csv and against the
Cellpose label image itself, and the run fails closed on any disagreement.

Output is a fresh directory with one RGB base image per condition, all sharing
one display aspect so the rendered panels are equal in width, plus a JSON
carrying the annotated polygons in display coordinates so the figure renderer
can draw the anatomy outlines and labels as vectors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from PIL import Image

HERE = Path(__file__).resolve()
# The paper root is found by walking up to the directory that holds README.txt
# and scripts/setup, so this works from Fig2/ and from its scripts/Fig2 mirror
# alike. Deriving it as HERE.parent.parent silently resolved to Paper/scripts
# when run from the mirror, and every raw path under it then pointed at
# scripts/Fig2/raw_data, which does not exist.
PAPER = next((path for path in HERE.parents
              if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                  and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
FIGURE = PAPER / "Fig2"
RAW = FIGURE / "raw_data" / "Apotome"

CONDITIONS = ("Water", "Sucrose", "Allulose")
REGION_CODES = {0: "OTHERS", 1: "ARC", 2: "ME", 3: "VMN"}
REGION_ORDER = ("ME", "ARC", "VMN")
DISPLAY_MAX_DIMENSION_PX = 2400
CROP_MARGIN_FRACTION = 0.06
GRAY_PERCENTILES = (0.5, 99.5)
GRAY_SAMPLE_STRIDE = 8
FOS_RED = np.array([236, 26, 26], dtype=np.float32)
FOS_ALPHA = 0.95

# One representative accepted section per condition, matching the panel C animals.
TILES = {
    "Water": {
        "tile_key": "e10-fr1-1-agua-dapi-pomc488-cfos647-001",
        "run": RAW / "24_06_2025/analysis/human_regions_v9",
        "dapi_tif": RAW / "24_06_2025/analysis/E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001-0001_dapi.tif",
        "dapi_seg": RAW / "24_06_2025/analysis/E10_FR1-1_AGUA_DAPI_POMC488_CFOS647_001-0001_dapi_seg.npy",
    },
    "Sucrose": {
        "tile_key": "e7fr7-1h-sacarosa",
        "run": RAW / "Aug_2025_NPY/selected_human_regions_v9_2",
        "dapi_tif": RAW / "Aug_2025_NPY/selected/E7FR7-1H-sacarosa-DAPI.tif",
        "dapi_seg": RAW / "Aug_2025_NPY/selected/E7FR7-1H-sacarosa-DAPI_seg.npy",
    },
    "Allulose": {
        "tile_key": "e8fr6-4h-alulosa",
        "run": RAW / "Aug_2025_NPY/selected_human_regions_v9_2",
        "dapi_tif": RAW / "Aug_2025_NPY/selected/E8FR6-4H-alulosa-DAPI.tif",
        "dapi_seg": RAW / "Aug_2025_NPY/selected/E8FR6-4H-alulosa-DAPI_seg.npy",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_seg_labels(path: Path) -> np.ndarray:
    obj = np.load(path, allow_pickle=True)
    if isinstance(obj, np.ndarray) and obj.dtype == object and obj.shape == ():
        obj = obj.item()
    if isinstance(obj, dict) and "masks" in obj:
        return np.asarray(obj["masks"])
    raise RuntimeError(f"Unrecognised Cellpose segmentation payload: {path}")


def read_gray(path: Path) -> np.ndarray:
    """Grey DAPI plate on one uniform percentile stretch for every condition."""
    array = np.asarray(tifffile.imread(path)).astype(np.float32)
    if array.ndim == 3:
        array = array[..., :3].mean(axis=2) if array.shape[-1] in (3, 4) else array[..., 0]
    sample = array[::GRAY_SAMPLE_STRIDE, ::GRAY_SAMPLE_STRIDE]
    low, high = (float(v) for v in np.percentile(sample, GRAY_PERCENTILES))
    if not high > low:
        raise RuntimeError(f"Degenerate DAPI display range for {path}")
    return np.clip((array - low) / (high - low), 0.0, 1.0)


def tiff_pixel_size_um(path: Path) -> float:
    with tifffile.TiffFile(path) as handle:
        page = handle.pages[0]
        meta = handle.imagej_metadata or {}
        if str(meta.get("unit", "")).lower() not in {"micron", "um", "µm"}:
            raise RuntimeError(f"{path} is not calibrated in microns")
        numerator, denominator = page.tags["XResolution"].value
        pixels_per_um = float(numerator) / float(denominator)
    if pixels_per_um <= 0:
        raise RuntimeError(f"Non-positive calibration in {path}")
    return 1.0 / pixels_per_um


def build(condition: str, spec: dict) -> dict:
    run = Path(spec["run"])
    tile = spec["tile_key"]
    annotation_path = run / "annotations" / f"{tile}.json"
    mask_path = run / "region_masks" / f"{tile}_regions.tif"
    summary_path = run / "per_image_human_region_summary.csv"
    nucleus_path = run / "per_nucleus_human_region_assignments.csv"
    status_path = run / "annotation_status.csv"
    complete_path = run / "ANNOTATION_COMPLETE.txt"
    for path in (annotation_path, mask_path, summary_path, nucleus_path, status_path, complete_path, spec["dapi_tif"], spec["dapi_seg"]):
        if not Path(path).is_file():
            raise FileNotFoundError(path)
    qc_png = next((run / "qc_human_regions").rglob(f"*{tile}_human_regions_qc.png"), None)
    if qc_png is None:
        raise FileNotFoundError(f"No accepted QC plate for {tile}")

    annotation = json.loads(annotation_path.read_text(encoding="utf-8"))
    if annotation.get("decision") != "include" or annotation.get("pipeline_version") != "9.2.0":
        raise RuntimeError(f"{condition}: annotation is not an accepted v9.2 include decision")
    if annotation.get("cond") != condition:
        raise RuntimeError(f"{condition}: annotation records condition {annotation.get('cond')!r}")
    polygons = annotation["polygons"]
    if {str(item["region"]) for item in polygons} != set(REGION_ORDER):
        raise RuntimeError(f"{condition}: annotation lacks a complete ARC/ME/VMN set")

    status = pd.read_csv(status_path)
    status_row = status.loc[status.tile_key.eq(tile)]
    if len(status_row) != 1 or str(status_row.iloc[0]["status"]) != "annotated":
        raise RuntimeError(f"{condition}: HIL status is not uniquely 'annotated'")

    regions = np.asarray(tifffile.imread(mask_path))
    height, width = regions.shape
    if (height, width) != (int(annotation["image_height"]), int(annotation["image_width"])):
        raise RuntimeError(f"{condition}: region mask shape disagrees with the annotation")
    if set(np.unique(regions).tolist()) - set(REGION_CODES) != set():
        raise RuntimeError(f"{condition}: unexpected region mask codes")

    summary = pd.read_csv(summary_path)
    summary = summary.loc[summary.tile_key.eq(tile)].set_index("region")
    nuclei = pd.read_csv(nucleus_path)
    nuclei = nuclei.loc[nuclei.tile_key.eq(tile)]
    if nuclei.empty:
        raise RuntimeError(f"{condition}: no per-nucleus assignments for {tile}")

    # Gate 1: region mask areas agree with the recorded per-image summary.
    for code, name in REGION_CODES.items():
        recorded = int(summary.loc[name, "area_px"])
        computed = int((regions == code).sum())
        if recorded != computed:
            raise RuntimeError(f"{condition}/{name}: mask area {computed} != recorded {recorded}")

    # Gate 2: per-nucleus assignments aggregate to the recorded per-region counts.
    grouped = nuclei.groupby("region")
    for name in REGION_CODES.values():
        block = grouped.get_group(name) if name in grouped.groups else nuclei.iloc[:0]
        if len(block) != int(summary.loc[name, "dapi_nuclei"]):
            raise RuntimeError(f"{condition}/{name}: {len(block)} nuclei != recorded {int(summary.loc[name, 'dapi_nuclei'])}")
        positives = int(block.cfos_positive.astype(bool).sum())
        if positives != int(summary.loc[name, "cfos_cells"]):
            raise RuntimeError(f"{condition}/{name}: {positives} c-FOS+ != recorded {int(summary.loc[name, 'cfos_cells'])}")

    labels = load_seg_labels(Path(spec["dapi_seg"]))
    if labels.shape != (height, width):
        raise RuntimeError(f"{condition}: Cellpose label image shape disagrees with the region mask")
    # Gate 3: the assignment table covers exactly the segmented nuclei.
    recorded_labels = set(int(value) for value in nuclei.nucleus_label)
    if recorded_labels != set(range(1, int(labels.max()) + 1)):
        raise RuntimeError(f"{condition}: per-nucleus labels do not cover the Cellpose label image")

    dapi_gray = read_gray(Path(spec["dapi_tif"]))
    if dapi_gray.shape != (height, width):
        raise RuntimeError(f"{condition}: DAPI image shape disagrees with the region mask")

    annotated = regions > 0
    rows = np.flatnonzero(annotated.any(axis=1))
    cols = np.flatnonzero(annotated.any(axis=0))
    pad_y = int(round((rows[-1] - rows[0] + 1) * CROP_MARGIN_FRACTION))
    pad_x = int(round((cols[-1] - cols[0] + 1) * CROP_MARGIN_FRACTION))
    box = {
        "y0": max(0, int(rows[0]) - pad_y),
        "y1": min(height, int(rows[-1]) + 1 + pad_y),
        "x0": max(0, int(cols[0]) - pad_x),
        "x1": min(width, int(cols[-1]) + 1 + pad_x),
    }
    positive_labels = np.array(sorted(int(v) for v in nuclei.loc[nuclei.cfos_positive.astype(bool), "nucleus_label"]), dtype=np.int64)
    return {
        "condition": condition,
        "animal": annotation["animal"],
        "sample": annotation["sample"],
        "tile_key": tile,
        "native_height_px": height,
        "native_width_px": width,
        "pixel_size_um": tiff_pixel_size_um(Path(spec["dapi_tif"])),
        "box": box,
        "polygons": [
            {"region": str(item["region"]), "points": [[float(x), float(y)] for x, y in item["points"]]}
            for item in polygons
        ],
        "region_counts": {
            name: {
                "dapi_nuclei": int(summary.loc[name, "dapi_nuclei"]),
                "cfos_cells": int(summary.loc[name, "cfos_cells"]),
                "area_um2": float(summary.loc[name, "area_um2"]),
            }
            for name in REGION_CODES.values()
        },
        "total_dapi_nuclei": int(labels.max()),
        "total_cfos_positive_nuclei": int(positive_labels.size),
        "accepted_fos_rois": int(summary.iloc[0]["accepted_fos_rois"]),
        "total_fos_rois": int(summary.iloc[0]["total_fos_rois"]),
        "annotation_saved_at_utc": annotation["saved_at_utc"],
        "sources": {
            role: {"path": str(path), "sha256": sha256(Path(path))}
            for role, path in (
                ("annotation_json", annotation_path),
                ("region_mask", mask_path),
                ("accepted_qc_plate", qc_png),
                ("per_image_summary", summary_path),
                ("per_nucleus_assignments", nucleus_path),
                ("annotation_status", status_path),
                ("annotation_complete_marker", complete_path),
                ("dapi_tif", spec["dapi_tif"]),
                ("dapi_seg", spec["dapi_seg"]),
            )
        },
        "recorded_counts_match": True,
        "_arrays": (dapi_gray, labels, positive_labels),
    }


def fit_box(record: dict, target_aspect: float) -> tuple[dict, tuple[int, int], tuple[int, int]]:
    """Grow the crop to the shared aspect and report its letterboxed canvas."""
    height, width = record["native_height_px"], record["native_width_px"]
    box = dict(record["box"])
    box_h, box_w = box["y1"] - box["y0"], box["x1"] - box["x0"]
    if box_w / box_h < target_aspect:
        grow = (box_h * target_aspect - box_w) / 2.0
        box["x0"] = max(0, box["x0"] - int(np.floor(grow)))
        box["x1"] = min(width, box["x1"] + int(np.ceil(grow)))
    else:
        grow = (box_w / target_aspect - box_h) / 2.0
        box["y0"] = max(0, box["y0"] - int(np.floor(grow)))
        box["y1"] = min(height, box["y1"] + int(np.ceil(grow)))
    box_h, box_w = box["y1"] - box["y0"], box["x1"] - box["x0"]
    canvas_w = max(box_w, int(round(box_h * target_aspect)))
    canvas_h = max(box_h, int(round(canvas_w / target_aspect)))
    offset = ((canvas_w - box_w) // 2, (canvas_h - box_h) // 2)
    return box, (canvas_w, canvas_h), offset


def rasterise(record: dict, box: dict, canvas: tuple[int, int], offset: tuple[int, int], target_size: tuple[int, int], outdir: Path) -> None:
    dapi_gray, labels, positive_labels = record.pop("_arrays")
    canvas_w, canvas_h = canvas
    box_h, box_w = box["y1"] - box["y0"], box["x1"] - box["x0"]

    window = labels[box["y0"] : box["y1"], box["x0"] : box["x1"]]
    base = np.repeat(dapi_gray[box["y0"] : box["y1"], box["x0"] : box["x1"], None] * 255.0, 3, axis=2)
    positive = np.isin(window, positive_labels)
    base[positive] = base[positive] * (1.0 - FOS_ALPHA) + FOS_RED * FOS_ALPHA
    image = Image.fromarray(np.rint(np.clip(base, 0, 255)).astype(np.uint8))

    if (canvas_w, canvas_h) != (box_w, box_h):
        plate = Image.new("RGB", (canvas_w, canvas_h), (0, 0, 0))
        plate.paste(image, offset)
        image = plate
    image = image.resize(target_size, Image.Resampling.LANCZOS)
    png = outdir / f"{record['condition']}_human_regions_base.png"
    image.save(png, optimize=True)

    scale_x = image.width / canvas_w
    scale_y = image.height / canvas_h
    for polygon in record["polygons"]:
        polygon["display_points"] = [
            [(x - box["x0"] + offset[0]) * scale_x, (y - box["y0"] + offset[1]) * scale_y]
            for x, y in polygon["points"]
        ]
        polygon.pop("points")
    record.update(
        {
            "box": box,
            "letterbox_offset_px": list(offset),
            "canvas_width_px": canvas_w,
            "canvas_height_px": canvas_h,
            "display_width_px": image.width,
            "display_height_px": image.height,
            "display_scale_x": scale_x,
            "display_scale_y": scale_y,
            "display_pixel_size_um": record["pixel_size_um"] / scale_x,
            "png": png.name,
            "png_sha256": sha256(png),
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Fresh absolute output directory")
    args = parser.parse_args()
    outdir = args.outdir.expanduser()
    if not outdir.is_absolute():
        raise ValueError("Output directory must be absolute")
    if outdir.exists() or outdir.is_symlink():
        raise FileExistsError(f"Refusing to replace an existing HIL anatomy stage: {outdir}")
    outdir.mkdir(parents=True)

    records = [build(condition, TILES[condition]) for condition in CONDITIONS]
    aspects = [(r["box"]["x1"] - r["box"]["x0"]) / (r["box"]["y1"] - r["box"]["y0"]) for r in records]
    target_aspect = float(np.mean(aspects))
    fitted = [fit_box(record, target_aspect) for record in records]
    # One display size for all three, so the rendered panels are equal in width.
    shortest = min(canvas[1] for _, canvas, _ in fitted)
    target_h = min(shortest, int(round(DISPLAY_MAX_DIMENSION_PX)))
    target_w = int(round(target_h * target_aspect))
    if target_w > DISPLAY_MAX_DIMENSION_PX:
        target_w = DISPLAY_MAX_DIMENSION_PX
        target_h = int(round(target_w / target_aspect))
    for record, (box, canvas, offset) in zip(records, fitted):
        rasterise(record, box, canvas, offset, (target_w, target_h), outdir)
    if len({(r["display_width_px"], r["display_height_px"]) for r in records}) != 1:
        raise RuntimeError("HIL anatomy panels did not converge on one display size")

    manifest = {
        "status": "passed",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "panel": "B accepted HIL-defined ARC/ME/VMN anatomy over Cellpose DAPI segmentation",
        "pipeline_version": "9.2.0",
        "nucleus_assignment_source": "per_nucleus_human_region_assignments.csv (accepted HIL run)",
        "crop_rule": "bounding box of the annotated ARC/ME/VMN union, expanded 6%, then grown symmetrically to the shared panel aspect",
        "shared_display_aspect": target_aspect,
        "equal_panel_width": True,
        "display_only": True,
        "manual_scope": "HIL anatomical annotation only; DAPI and c-FOS object segmentation and their association remain automated",
        "conditions": records,
    }
    (outdir / "human_regions_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for record in records:
        counts = record["region_counts"]
        print(
            f"[OK] {record['condition']} {record['animal']}: "
            + ", ".join(f"{name} {counts[name]['cfos_cells']}/{counts[name]['dapi_nuclei']}" for name in REGION_ORDER)
            + f", {record['display_width_px']}x{record['display_height_px']} px"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
