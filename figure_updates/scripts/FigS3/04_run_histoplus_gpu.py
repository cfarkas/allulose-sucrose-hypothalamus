#!/usr/bin/env python3
"""Run organ-stratified LazySlide/HistoPLUS inference and classical H&E CV."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import os
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from types import MethodType

import cv2
import geopandas as gpd
import lazyslide as zs
import numpy as np
import pandas as pd
import torch
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from shapely.affinity import scale as scale_geometry, translate as translate_geometry
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from skimage import color, feature, filters, measure, morphology, segmentation
from spatialdata.models import ShapesModel
from torch.utils.data import DataLoader

from figs3_common import (
    FIGS3_DIR,
    ORGANS,
    FigS3Error,
    atomic_csv,
    atomic_json,
    read_csv,
    require_regular_file,
    sha256_file,
    validate_sample_id,
)


RUNNER_VERSION = "figs3-lazyslide-histoplus-1.1"
MODEL_NAME = "HistoPLUS-20x"
MODEL_LICENSE = "CC-BY-NC-ND-4.0"
TARGET_MPP = 0.5
TILE_SIZE = 840
ORGAN_VALUE = {"Kidney": 1, "Liver": 2, "Spleen": 3}
CLASS_NAMES = {
    1: "Cancer cell",
    2: "Lymphocytes",
    3: "Fibroblasts",
    4: "Plasmocytes",
    5: "Eosinophils",
    6: "Neutrophils",
    7: "Macrophages",
    8: "Muscle Cell",
    9: "Endothelial Cell",
    10: "Red blood cell",
    11: "Epithelial",
    12: "Apoptotic Body",
    13: "Mitotic Figures",
    14: "Minor Stromal Cell",
}
CLASS_COLORS = {
    name: color_value
    for name, color_value in zip(
        CLASS_NAMES.values(),
        [
            "#d73027", "#4575b4", "#fdae61", "#8073ac", "#fee090", "#91bfdb",
            "#1a9850", "#a6611a", "#018571", "#b2182b", "#66bd63", "#762a83",
            "#e08214", "#999999",
        ],
    )
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-manifest",
        type=Path,
        default=FIGS3_DIR / "source_data" / "Figure_S3_slide_manifest.csv",
    )
    parser.add_argument("--export-dir", type=Path, default=FIGS3_DIR / "exports")
    parser.add_argument("--segmentation-dir", type=Path, default=FIGS3_DIR / "organ_segmentation")
    parser.add_argument("--output-dir", type=Path, default=FIGS3_DIR / "histoplus")
    parser.add_argument(
        "--weight-file",
        type=Path,
        default=Path("/home/server/.cache/histoplus/histoplus_cellvit_segmentor_20x.pt"),
    )
    parser.add_argument("--tiles-per-organ", type=int, default=24)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=1707)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--prepare-only", action="store_true", help="Create/validate organ-stratified tile manifests without model inference")
    return parser.parse_args()


def scanner_values(path: Path) -> tuple[float, int, int]:
    values: dict[str, str] = {}
    for line in require_regular_file(path, "scanner info.ini").read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip().casefold()] = value.strip()
    try:
        return float(values["scale"]), int(values["mifwidth"]), int(values["mifheight"])
    except (KeyError, ValueError) as exc:
        raise FigS3Error(f"Incomplete scanner geometry in {path}") from exc


def mask_geometry(mask: np.ndarray, true_width: int, true_height: int):
    contours, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    polygons = []
    for contour in contours:
        points = contour[:, 0, :]
        if len(points) < 3 or cv2.contourArea(contour) < 25:
            continue
        polygon = Polygon(points).buffer(0)
        if polygon.is_empty:
            continue
        polygons.append(
            scale_geometry(
                polygon,
                xfact=true_width / mask.shape[1],
                yfact=true_height / mask.shape[0],
                origin=(0, 0),
            )
        )
    if not polygons:
        raise FigS3Error("Organ mask has no valid polygon")
    geometry = unary_union(polygons).buffer(0)
    if geometry.is_empty:
        raise FigS3Error("Organ polygon union is empty")
    return geometry


def attach_organ_shapes(wsi, label_mask: np.ndarray, true_width: int, true_height: int) -> None:
    rows = []
    for index, organ in enumerate(ORGANS):
        rows.append(
            {
                "tissue_id": index,
                "organ": organ,
                "geometry": mask_geometry(label_mask == ORGAN_VALUE[organ], true_width, true_height),
            }
        )
    wsi.shapes["organs"] = ShapesModel.parse(gpd.GeoDataFrame(rows, geometry="geometry"))


def stable_seed(sample_id: str, organ: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}|{sample_id}|{organ}".encode()).digest()
    return int.from_bytes(digest[:4], "little")


def farthest_point_sample(frame: gpd.GeoDataFrame, count: int, seed: int) -> gpd.GeoDataFrame:
    if len(frame) < count:
        raise FigS3Error(f"Only {len(frame)} eligible tiles are available; {count} required")
    centers = np.column_stack([frame.geometry.centroid.x, frame.geometry.centroid.y]).astype(float)
    minimum = centers.min(axis=0)
    spread = np.maximum(centers.max(axis=0) - minimum, 1.0)
    points = (centers - minimum) / spread
    rng = np.random.default_rng(seed)
    centroid = points.mean(axis=0)
    distances_to_centroid = np.sum((points - centroid) ** 2, axis=1)
    tied = np.flatnonzero(np.isclose(distances_to_centroid, distances_to_centroid.min()))
    first = int(rng.choice(tied))
    selected = [first]
    minimum_distance = np.sum((points - points[first]) ** 2, axis=1)
    while len(selected) < count:
        minimum_distance[selected] = -1.0
        maximum = minimum_distance.max()
        tied = np.flatnonzero(np.isclose(minimum_distance, maximum))
        choice = int(rng.choice(tied))
        selected.append(choice)
        minimum_distance = np.minimum(
            minimum_distance, np.sum((points - points[choice]) ** 2, axis=1)
        )
    output = frame.iloc[selected].copy().reset_index(drop=True)
    output["sampling_rank"] = np.arange(1, len(output) + 1)
    return output


def choose_tiles(
    candidates: gpd.GeoDataFrame,
    sample_id: str,
    count: int,
    seed: int,
    reserve: int = 0,
) -> gpd.GeoDataFrame:
    tissue_to_organ = {0: "Kidney", 1: "Liver", 2: "Spleen"}
    candidates = candidates.copy()
    candidates["organ"] = candidates["tissue_id"].map(tissue_to_organ)
    if candidates["organ"].isna().any():
        raise FigS3Error("LazySlide returned an unknown tissue_id")
    selected_frames = []
    for organ in ORGANS:
        pool = candidates[candidates["organ"] == organ].reset_index(drop=True)
        target = min(len(pool), count + reserve)
        if target < count:
            raise FigS3Error(f"Only {target} eligible {organ} tiles are available; {count} required")
        chosen = farthest_point_sample(pool, target, stable_seed(sample_id, organ, seed))
        chosen["candidate_tiles_in_organ"] = len(pool)
        selected_frames.append(chosen)
    selected = gpd.GeoDataFrame(pd.concat(selected_frames, ignore_index=True), geometry="geometry")
    selected["tile_uid"] = [f"{sample_id}_{organ.lower()}_{rank:03d}" for organ, rank in zip(selected["organ"], selected["sampling_rank"])]
    selected["selected_tile_id"] = np.arange(len(selected))
    return selected


def install_selected_tiles(wsi, selected: gpd.GeoDataFrame) -> None:
    wsi.shapes["selected_tiles"] = ShapesModel.parse(selected)
    specs = wsi.attrs[wsi.TILE_SPEC_KEY]
    specs["selected_tiles"] = dict(specs["candidate_tiles"])


def tile_manifest_rows(sample_id: str, candidates: gpd.GeoDataFrame, selected: gpd.GeoDataFrame, source_mpp: float) -> list[dict[str, object]]:
    selected_ids = set(int(value) for value in selected["tile_id"])
    selected_by_id = {int(row.tile_id): row for row in selected.itertuples()}
    tissue_to_organ = {0: "Kidney", 1: "Liver", 2: "Spleen"}
    rows = []
    for row in candidates.itertuples():
        bounds = row.geometry.bounds
        chosen = int(row.tile_id) in selected_ids
        chosen_row = selected_by_id.get(int(row.tile_id))
        rows.append(
            {
                "sample_id": sample_id,
                "candidate_tile_id": int(row.tile_id),
                "tile_uid": chosen_row.tile_uid if chosen else "",
                "organ": tissue_to_organ[int(row.tissue_id)],
                "selected": str(chosen).lower(),
                "sampling_rank": int(chosen_row.sampling_rank) if chosen else "",
                "x_l0": int(round(bounds[0])),
                "y_l0": int(round(bounds[1])),
                "width_l0": int(round(bounds[2] - bounds[0])),
                "height_l0": int(round(bounds[3] - bounds[1])),
                "source_mpp": source_mpp,
                "target_mpp": TARGET_MPP,
                "tile_px_at_target_mpp": TILE_SIZE,
                "physical_tile_area_mm2": (TILE_SIZE * TARGET_MPP / 1000.0) ** 2,
                "sampling_design": "deterministic spatial farthest-point within organ",
            }
        )
    return rows


def classical_tile_features(rgb: np.ndarray) -> dict[str, float | int]:
    rgb = np.asarray(rgb, dtype=np.uint8)[..., :3]
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    tissue = ((hsv[..., 1] >= 18) & (hsv[..., 2] <= 250)) | (hsv[..., 2] <= 205)
    tissue = morphology.remove_small_objects(tissue, max_size=63)
    hed = color.rgb2hed(rgb)
    h_channel = hed[..., 0]
    e_channel = hed[..., 1]
    gray = color.rgb2gray(rgb)
    if tissue.sum() < 256:
        raise FigS3Error("Selected tile contains insufficient tissue")
    threshold = filters.threshold_otsu(h_channel[tissue])
    nuclei = (h_channel > threshold) & tissue
    nuclei = morphology.opening(nuclei, morphology.disk(1))
    nuclei = morphology.remove_small_objects(nuclei, max_size=9)
    nuclei = morphology.remove_small_holes(nuclei, max_size=15)
    distance = ndi.distance_transform_edt(nuclei)
    coordinates = feature.peak_local_max(
        distance, min_distance=4, threshold_abs=1.5, labels=nuclei, exclude_border=False
    )
    markers = np.zeros(nuclei.shape, dtype=np.int32)
    if len(coordinates):
        markers[tuple(coordinates.T)] = np.arange(1, len(coordinates) + 1)
        labels = segmentation.watershed(-distance, markers, mask=nuclei)
    else:
        labels = measure.label(nuclei)
    regions = [region for region in measure.regionprops(labels) if 12 <= region.area <= 1200]
    areas = np.array([region.area for region in regions], dtype=float)
    eccentricity = np.array([region.eccentricity for region in regions], dtype=float)
    solidity = np.array([region.solidity for region in regions], dtype=float)
    tissue_area_mm2 = float(tissue.sum() * TARGET_MPP**2 / 1_000_000.0)
    return {
        "tissue_fraction": float(tissue.mean()),
        "analyzed_tissue_area_mm2": tissue_area_mm2,
        "hematoxylin_mean": float(np.mean(h_channel[tissue])),
        "hematoxylin_p95": float(np.quantile(h_channel[tissue], 0.95)),
        "eosin_mean": float(np.mean(e_channel[tissue])),
        "gray_mean": float(np.mean(gray[tissue])),
        "gray_std": float(np.std(gray[tissue])),
        "laplacian_variance": float(cv2.Laplacian((gray * 255).astype(np.uint8), cv2.CV_64F).var()),
        "classical_nuclei_count": int(len(regions)),
        "classical_nuclei_density_mm2": float(len(regions) / tissue_area_mm2) if tissue_area_mm2 > 0 else math.nan,
        "classical_nucleus_area_median_um2": float(np.median(areas) * TARGET_MPP**2) if len(areas) else math.nan,
        "classical_nucleus_eccentricity_median": float(np.median(eccentricity)) if len(eccentricity) else math.nan,
        "classical_nucleus_solidity_median": float(np.median(solidity)) if len(solidity) else math.nan,
    }


def as_rgb_array(image) -> np.ndarray:
    """Normalize the wsidata reader's PIL-or-ndarray return type to RGB uint8."""
    if isinstance(image, Image.Image):
        return np.asarray(image.convert("RGB"), dtype=np.uint8)
    array = np.asarray(image)
    if array.ndim != 3 or array.shape[2] < 3:
        raise FigS3Error(f"Unexpected tile image shape: {array.shape}")
    return np.asarray(array[..., :3], dtype=np.uint8)


def direct_tissue_fraction(rgb: np.ndarray) -> float:
    hsv = cv2.cvtColor(as_rgb_array(rgb), cv2.COLOR_RGB2HSV)
    tissue = ((hsv[..., 1] >= 18) & (hsv[..., 2] <= 250)) | (hsv[..., 2] <= 205)
    tissue = morphology.remove_small_objects(tissue, max_size=63)
    return float(tissue.mean())


@contextlib.contextmanager
def local_weight_override(weight_file: Path):
    import huggingface_hub

    original = huggingface_hub.hf_hub_download

    def patched(*args, **kwargs):
        filename = kwargs.get("filename")
        if filename is None and len(args) >= 2:
            filename = args[1]
        if filename != weight_file.name:
            raise FigS3Error(f"Unexpected gated-model filename request: {filename}")
        return str(weight_file)

    huggingface_hub.hf_hub_download = patched
    try:
        yield
    finally:
        huggingface_hub.hf_hub_download = original


def build_model(weight_file: Path):
    with local_weight_override(weight_file):
        return zs.models.segmentation.HistoPLUS(tile_size=TILE_SIZE, magnification="20x")


def enable_parallel_postprocessing(model, workers: int) -> None:
    """Parallelize the unchanged HoVer postprocessor across a model batch."""
    from lazyslide.models.segmentation.cellvit_family.postprocess import np_hv_postprocess

    workers = max(1, int(workers))

    def segment_parallel(self, images):
        with torch.inference_mode():
            output = self.model(images)
        flattened = [dict(zip(output.keys(), values)) for values in zip(*output.values())]
        cpu_batches = [
            (
                batch["np"].softmax(0).detach().float().cpu().numpy()[1],
                batch["hv"].detach().float().cpu().numpy(),
                batch["tp"].softmax(0).detach().float().cpu().numpy(),
            )
            for batch in flattened
        ]

        def postprocess(item):
            nuclei, hv_map, probability = item
            return np_hv_postprocess(nuclei, hv_map, variant=self.variant), probability

        with ThreadPoolExecutor(max_workers=min(workers, len(cpu_batches))) as pool:
            processed = list(pool.map(postprocess, cpu_batches))
        return {
            "instance_map": np.asarray([item[0] for item in processed]),
            "class_map": np.asarray([item[1] for item in processed]),
        }

    model.segment = MethodType(segment_parallel, model)


def vectorize_instance_tile(
    instance_map: np.ndarray,
    probability_map: np.ndarray,
    pos_x: float,
    pos_y: float,
    downsample: float,
) -> gpd.GeoDataFrame:
    """Vectorize tight instance crops with LazySlide-equivalent CV operations."""
    rows: list[dict[str, object]] = []
    edge_boundary = box(0, 0, TILE_SIZE, TILE_SIZE).buffer(-2).boundary
    objects = ndi.find_objects(instance_map)
    for offset, slices in enumerate(objects, start=1):
        if slices is None:
            continue
        y_slice, x_slice = slices
        local_mask = np.asarray(instance_map[y_slice, x_slice] == offset, dtype=np.uint8)
        # One-pixel zero padding preserves the same exterior contour that the
        # full-tile LazySlide implementation obtains, but avoids allocating an
        # 840x840 mask for every nucleus.
        local_mask = np.pad(local_mask, 1, mode="constant")
        local_probability = np.pad(
            probability_map[:, y_slice, x_slice],
            ((0, 0), (1, 1), (1, 1)),
            mode="constant",
        )
        # These are the operations used by LazySlide's
        # binary_mask_to_polygons_with_prob. Performing them directly avoids
        # constructing one GeoDataFrame per nucleus while preserving contours,
        # filled-polygon class probabilities, and class selection.
        contours, hierarchy = cv2.findContours(
            local_mask, mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE
        )
        if hierarchy is None:
            continue
        polygon = None
        for contour in contours:
            if cv2.contourArea(contour) <= 0:
                continue
            points = np.squeeze(contour, axis=1)
            if len(points) >= 4:
                polygon = Polygon(shell=points, holes=[])
                break
        if polygon is None:
            continue
        polygon_mask = np.zeros_like(local_mask, dtype=np.uint8)
        points = np.asarray(polygon.exterior.coords, dtype=np.int32)
        cv2.drawContours(polygon_mask, [points], -1, 1, thickness=cv2.FILLED)
        pixel_count = int(np.sum(polygon_mask))
        if pixel_count == 0:
            continue
        class_probabilities = np.sum(
            local_probability * polygon_mask, axis=(1, 2)
        ) / pixel_count
        class_index = int(np.argmax(class_probabilities))
        probability = float(class_probabilities[class_index])
        geometry = translate_geometry(
            polygon,
            xoff=float(x_slice.start - 1),
            yoff=float(y_slice.start - 1),
        )
        if edge_boundary.intersects(geometry):
            continue
        geometry = scale_geometry(
            geometry,
            xfact=downsample,
            yfact=downsample,
            origin=(0, 0),
        )
        geometry = translate_geometry(geometry, xoff=float(pos_x), yoff=float(pos_y)).buffer(0)
        if geometry.is_empty or not geometry.is_valid:
            continue
        rows.append(
            {
                "geometry": geometry,
                "prob": probability,
                "class": "Background" if class_index == 0 else CLASS_NAMES[class_index],
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry") if rows else gpd.GeoDataFrame(columns=["geometry", "prob", "class"], geometry="geometry")


def selected_tiles_overlap(selected: gpd.GeoDataFrame) -> bool:
    index = selected.sindex
    for row_index, geometry in enumerate(selected.geometry):
        for other in index.query(geometry, predicate="intersects"):
            if int(other) != row_index and geometry.intersection(selected.geometry.iloc[int(other)]).area > 0:
                return True
    return False


def run_histoplus_cells(
    wsi,
    selected: gpd.GeoDataFrame,
    model,
    device: str,
    batch_size: int,
    num_workers: int,
) -> gpd.GeoDataFrame:
    """Run HistoPLUS with parallel postprocessing and crop-local vectorization."""
    from lazyslide.cv import nms

    model.to(device)
    dataset = wsi.ds.tile_images(
        tile_key="selected_tiles", transform=model.get_transform()
    )
    loader = DataLoader(dataset, batch_size=batch_size, num_workers=num_workers)
    spec = wsi.tile_spec("selected_tiles")
    frames: list[gpd.GeoDataFrame] = []
    processed_tiles = 0
    with torch.inference_mode():
        for chunk in loader:
            images = chunk["image"].to(device)
            output = model.segment(images)
            tasks = [
                (
                    output["instance_map"][index],
                    output["class_map"][index],
                    float(chunk["x"][index]),
                    float(chunk["y"][index]),
                    float(spec.base_downsample),
                )
                for index in range(len(images))
            ]
            with ThreadPoolExecutor(max_workers=min(batch_size, len(tasks))) as pool:
                frames.extend(pool.map(lambda values: vectorize_instance_tile(*values), tasks))
            processed_tiles += len(tasks)
            print(
                json.dumps(
                    {
                        "stage": "histoplus_gpu",
                        "processed_tiles": processed_tiles,
                        "total_tiles": len(dataset),
                    }
                ),
                flush=True,
            )
    frames = [frame for frame in frames if len(frame)]
    if not frames:
        raise FigS3Error("HistoPLUS returned no vectorized nuclei")
    cells = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry")
    cells = cells[cells["class"] != "Background"].reset_index(drop=True)
    if selected_tiles_overlap(selected):
        cells = nms(cells, "prob").reset_index(drop=True)
    organ_union = wsi.shapes["organs"].geometry.union_all()
    cells = cells[cells.intersects(organ_union)].reset_index(drop=True)
    return cells


def cell_to_tile(cells: gpd.GeoDataFrame, selected: gpd.GeoDataFrame) -> pd.DataFrame:
    points = gpd.GeoDataFrame(
        {"cell_index": np.arange(len(cells))},
        geometry=cells.geometry.centroid,
        index=cells.index,
    )
    joined = gpd.sjoin(
        points,
        selected[["tile_uid", "organ", "geometry"]],
        how="left",
        predicate="within",
    )
    joined = joined[~joined.index.duplicated(keep="first")].reindex(cells.index)
    if joined["tile_uid"].isna().any():
        raise FigS3Error(f"{joined['tile_uid'].isna().sum()} cells cannot be assigned to selected tiles")
    return joined[["tile_uid", "organ"]]


def cell_rows(cells: gpd.GeoDataFrame, assigned: pd.DataFrame, sample_id: str, source_mpp: float) -> gpd.GeoDataFrame:
    output = cells.copy().reset_index(drop=True)
    assigned = assigned.reset_index(drop=True)
    output["sample_id"] = sample_id
    output["tile_uid"] = assigned["tile_uid"].to_numpy()
    output["organ"] = assigned["organ"].to_numpy()
    output["model_class_raw"] = output["class"].astype(str)
    output["phenotype_label"] = "model_predicted_" + output["model_class_raw"].str.lower().str.replace(" ", "_", regex=False)
    centroids = output.geometry.centroid
    output["centroid_x_l0"] = centroids.x
    output["centroid_y_l0"] = centroids.y
    output["nuclear_area_um2"] = output.geometry.area * source_mpp**2
    output["nuclear_perimeter_um"] = output.geometry.length * source_mpp
    output["nuclear_circularity"] = np.clip(
        4.0 * np.pi * output.geometry.area / np.maximum(output.geometry.length**2, 1e-12), 0, 1
    )
    return output


def summary_rows(
    cells: gpd.GeoDataFrame,
    tile_features: list[dict[str, object]],
    sample_row: dict[str, str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    tile_frame = pd.DataFrame(tile_features)
    phenotype_rows: list[dict[str, object]] = []
    tile_count_rows: list[dict[str, object]] = []
    for organ in ORGANS:
        organ_cells = cells[cells["organ"] == organ]
        organ_tiles = tile_frame[tile_frame["organ"] == organ]
        area = float(organ_tiles["analyzed_tissue_area_mm2"].sum())
        total = len(organ_cells)
        for class_name in CLASS_NAMES.values():
            count = int((organ_cells["model_class_raw"] == class_name).sum())
            phenotype_rows.append(
                {
                    "sample_id": sample_row["sample_id"],
                    "treatment": sample_row["treatment"],
                    "cohort": sample_row["cohort"],
                    "organ": organ,
                    "model": MODEL_NAME,
                    "model_class_raw": class_name,
                    "phenotype_label": "model_predicted_" + class_name.lower().replace(" ", "_"),
                    "cell_count_in_sampled_tiles": count,
                    "cell_fraction_in_sampled_tiles": count / total if total else math.nan,
                    "cell_density_per_analyzed_mm2": count / area if area else math.nan,
                    "total_cells_in_sampled_tiles": total,
                    "selected_tile_count": len(organ_tiles),
                    "analyzed_tissue_area_mm2": area,
                    "biological_unit": "animal",
                    "extrapolated_to_whole_organ": "false",
                    "cross_species_validated": "false",
                }
            )
        for tile_uid in organ_tiles["tile_uid"]:
            tile_cells = organ_cells[organ_cells["tile_uid"] == tile_uid]
            base = {
                "sample_id": sample_row["sample_id"],
                "treatment": sample_row["treatment"],
                "cohort": sample_row["cohort"],
                "organ": organ,
                "tile_uid": tile_uid,
                "total_cells": len(tile_cells),
            }
            counts = Counter(tile_cells["model_class_raw"])
            for class_name in CLASS_NAMES.values():
                base["count_" + class_name.lower().replace(" ", "_")] = counts[class_name]
            tile_count_rows.append(base)
    return phenotype_rows, tile_count_rows


def selected_tile_overview(
    overview_path: Path,
    selected: gpd.GeoDataFrame,
    true_width: int,
    true_height: int,
    output: Path,
) -> None:
    with Image.open(overview_path) as image:
        canvas = image.convert("RGB")
    draw = ImageDraw.Draw(canvas)
    colors = {"Kidney": "#0097ff", "Liver": "#22b14c", "Spleen": "#b446c8"}
    for row in selected.itertuples():
        x0, y0, x1, y1 = row.geometry.bounds
        box = (
            x0 * canvas.width / true_width,
            y0 * canvas.height / true_height,
            x1 * canvas.width / true_width,
            y1 * canvas.height / true_height,
        )
        draw.rectangle(box, outline=colors[row.organ], width=3)
    canvas.save(output)


def representative_cell_overlay(
    wsi,
    selected: gpd.GeoDataFrame,
    cells: gpd.GeoDataFrame,
    output: Path,
) -> str:
    counts = cells.groupby("tile_uid").size()
    ordered = selected.copy()
    ordered["cell_count"] = ordered["tile_uid"].map(counts).fillna(0)
    target = ordered.iloc[(ordered["cell_count"] - ordered["cell_count"].median()).abs().argsort().iloc[0]]
    tile_index = int(ordered.index[ordered["tile_uid"] == target["tile_uid"]][0])
    dataset = wsi.ds.tile_images(tile_key="selected_tiles")
    rgb = as_rgb_array(dataset[tile_index]["image"])
    image = Image.fromarray(rgb)
    draw = ImageDraw.Draw(image)
    bounds = target.geometry.bounds
    spec = wsi.tile_spec("selected_tiles")
    tile_cells = cells[cells["tile_uid"] == target["tile_uid"]]
    for cell in tile_cells.itertuples():
        center = cell.geometry.centroid
        x = (center.x - bounds[0]) / spec.base_downsample
        y = (center.y - bounds[1]) / spec.base_downsample
        color_value = CLASS_COLORS.get(cell.model_class_raw, "#000000")
        draw.ellipse((x - 3, y - 3, x + 3, y + 3), outline=color_value, width=2)
    image.save(output)
    return str(target["tile_uid"])


def existing_result(
    sample_dir: Path,
    args: argparse.Namespace,
    sample_id: str,
    weight_sha: str,
) -> dict[str, object] | None:
    receipt_path = sample_dir / "inference_receipt.json"
    if not receipt_path.is_file():
        return None
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "complete" or receipt.get("runner_version") != RUNNER_VERSION:
        return None
    if (
        receipt.get("sample_id") != sample_id
        or int(receipt.get("tiles_per_organ", -1)) != args.tiles_per_organ
        or float(receipt.get("target_mpp", -1)) != TARGET_MPP
        or int(receipt.get("tile_size_px", -1)) != TILE_SIZE
        or receipt.get("model_weight_sha256") != weight_sha
        or receipt.get("model") != MODEL_NAME
        or receipt.get("amp") is not False
    ):
        return None
    conversion = args.export_dir / sample_id / "conversion_receipt.json"
    labels = args.segmentation_dir / sample_id / "organ_label_mask.png"
    if (
        not conversion.is_file()
        or not labels.is_file()
        or sha256_file(conversion) != receipt.get("conversion_receipt_sha256")
        or sha256_file(labels) != receipt.get("organ_label_mask_sha256")
    ):
        return None
    for key in ("phenotype_summary", "tile_cv_features", "selected_tile_manifest", "tile_class_counts", "cells_parquet"):
        path = Path(str(receipt["outputs"][key]["path"]))
        if not path.is_file() or path.stat().st_size != int(receipt["outputs"][key]["size_bytes"]):
            return None
        if sha256_file(path) != str(receipt["outputs"][key]["sha256"]):
            return None
    return receipt


def output_identity(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def process_slide(
    row: dict[str, str],
    args: argparse.Namespace,
    model,
    weight_sha: str,
) -> dict[str, object]:
    sample_id = validate_sample_id(row["sample_id"])
    sample_dir = args.output_dir.resolve() / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    if args.resume:
        receipt = existing_result(sample_dir, args, sample_id, weight_sha)
        if receipt is not None:
            receipt["run_status"] = "verified_existing"
            return receipt
    started = time.time()
    try:
        source_mpp, true_width, true_height = scanner_values(Path(row["source_info_ini_path"]))
        slide_path = require_regular_file(args.export_dir / sample_id / "1_L0_rgb.tif", "L0 TIFF")
        conversion_receipt = require_regular_file(args.export_dir / sample_id / "conversion_receipt.json", "conversion receipt")
        label_path = require_regular_file(args.segmentation_dir / sample_id / "organ_label_mask.png", "organ label mask")
        label_mask = np.asarray(Image.open(label_path), dtype=np.uint8)
        wsi = zs.open_wsi(
            slide_path,
            store=None,
            reader="tiffslide",
            attach_thumbnail=False,
            save_thumbnail=False,
        )
        wsi.set_mpp(source_mpp)
        attach_organ_shapes(wsi, label_mask, true_width, true_height)
        result = zs.pp.tile_tissues(
            wsi,
            TILE_SIZE,
            mpp=TARGET_MPP,
            slide_mpp=source_mpp,
            background_filter=True,
            background_fraction=0.20,
            background_filter_mode="exact",
            tissue_key="organs",
            key_added="candidate_tiles",
            return_tiles=True,
        )
        if result is None:
            raise FigS3Error("LazySlide produced no eligible tiles")
        candidates, _ = result
        reserve = max(12, args.tiles_per_organ // 2)
        provisional = choose_tiles(
            candidates,
            sample_id,
            args.tiles_per_organ,
            args.seed,
            reserve=reserve,
        )
        install_selected_tiles(wsi, provisional)
        provisional_dataset = wsi.ds.tile_images(
            tile_key="selected_tiles", target_key="organ"
        )
        provisional_qc: list[dict[str, object]] = []
        rgb_by_candidate: dict[int, np.ndarray] = {}
        for index in range(len(provisional_dataset)):
            item = provisional_dataset[index]
            tile = provisional.iloc[index]
            rgb_tile = as_rgb_array(item["image"])
            fraction = direct_tissue_fraction(rgb_tile)
            accepted = fraction >= 0.35
            provisional_qc.append(
                {
                    "sample_id": sample_id,
                    "candidate_tile_id": int(tile["tile_id"]),
                    "organ": tile["organ"],
                    "provisional_sampling_rank": int(tile["sampling_rank"]),
                    "direct_tissue_fraction": fraction,
                    "direct_pixel_qc": "pass" if accepted else "reject",
                    "threshold": 0.35,
                }
            )
            if accepted:
                rgb_by_candidate[int(tile["tile_id"])] = rgb_tile
        qc_frame = pd.DataFrame(provisional_qc)
        accepted_ids = set(
            qc_frame.loc[qc_frame["direct_pixel_qc"] == "pass", "candidate_tile_id"]
            .astype(int)
            .tolist()
        )
        final_frames = []
        for organ in ORGANS:
            organ_pool = provisional[
                (provisional["organ"] == organ)
                & provisional["tile_id"].astype(int).isin(accepted_ids)
            ].sort_values("sampling_rank")
            if len(organ_pool) < args.tiles_per_organ:
                raise FigS3Error(
                    f"Only {len(organ_pool)} direct-pixel-QC {organ} tiles remain; "
                    f"{args.tiles_per_organ} required"
                )
            final_frames.append(organ_pool.head(args.tiles_per_organ).copy())
        selected = gpd.GeoDataFrame(
            pd.concat(final_frames, ignore_index=True), geometry="geometry"
        )
        selected["provisional_sampling_rank"] = selected["sampling_rank"].astype(int)
        selected["sampling_rank"] = selected.groupby("organ").cumcount() + 1
        selected["tile_uid"] = [
            f"{sample_id}_{organ.lower()}_{rank:03d}"
            for organ, rank in zip(selected["organ"], selected["sampling_rank"])
        ]
        selected["selected_tile_id"] = np.arange(len(selected))
        install_selected_tiles(wsi, selected)
        all_tile_rows = tile_manifest_rows(sample_id, candidates, selected, source_mpp)
        selected_rows = [tile for tile in all_tile_rows if tile["selected"] == "true"]
        selected_manifest = sample_dir / "selected_tile_manifest.csv"
        candidate_manifest = sample_dir / "candidate_tile_manifest.csv"
        atomic_csv(selected_manifest, selected_rows)
        atomic_csv(candidate_manifest, all_tile_rows)
        pixel_qc_path = sample_dir / "provisional_tile_pixel_qc.csv"
        atomic_csv(pixel_qc_path, provisional_qc)
        selected_tile_overview(
            Path(row["source_overview_path"]), selected, true_width, true_height, sample_dir / "selected_tiles_overview.png"
        )
        selected_rgbs = [rgb_by_candidate[int(tile.tile_id)] for tile in selected.itertuples()]
        with ThreadPoolExecutor(
            max_workers=min(max(1, args.num_workers), len(selected_rgbs))
        ) as pool:
            selected_features = list(pool.map(classical_tile_features, selected_rgbs))
        del rgb_by_candidate, selected_rgbs
        tile_features: list[dict[str, object]] = []
        for tile, features in zip(selected.itertuples(), selected_features):
            tile_features.append(
                {
                    "sample_id": sample_id,
                    "treatment": row["treatment"],
                    "cohort": row["cohort"],
                    "organ": tile.organ,
                    "tile_uid": tile.tile_uid,
                    "sampling_rank": int(tile.sampling_rank),
                    "provisional_sampling_rank": int(tile.provisional_sampling_rank),
                    **features,
                }
            )
        tile_cv_path = sample_dir / "tile_cv_features.csv"
        atomic_csv(tile_cv_path, tile_features)
        if args.prepare_only:
            return {
                "status": "prepared",
                "sample_id": sample_id,
                "candidate_tiles": len(candidates),
                "selected_tiles": len(selected),
            }
        torch.cuda.reset_peak_memory_stats() if args.device.startswith("cuda") else None
        # Full precision avoids the installed OpenCV float16 minMaxIdx failure.
        # The crop-local vectorizer is algorithmically equivalent to LazySlide's
        # InstanceMap path while avoiding a full-tile mask per nucleus.
        cells = run_histoplus_cells(
            wsi,
            selected,
            model,
            args.device,
            args.batch_size,
            args.num_workers,
        )
        if cells is None or len(cells) == 0:
            raise FigS3Error("HistoPLUS returned no cells; result is not encoded as zero")
        assigned = cell_to_tile(cells, selected)
        cells = cell_rows(cells, assigned, sample_id, source_mpp)
        phenotype_rows, tile_count_rows = summary_rows(cells, tile_features, row)
        cells_path = sample_dir / "cells.parquet"
        phenotype_path = sample_dir / "phenotype_summary.csv"
        tile_counts_path = sample_dir / "tile_class_counts.csv"
        cells.to_parquet(cells_path, index=False)
        atomic_csv(phenotype_path, phenotype_rows)
        atomic_csv(tile_counts_path, tile_count_rows)
        representative_tile = representative_cell_overlay(
            wsi, selected, cells, sample_dir / "representative_histoplus_overlay.png"
        )
        outputs = {
            "phenotype_summary": output_identity(phenotype_path),
            "tile_cv_features": output_identity(tile_cv_path),
            "tile_pixel_qc": output_identity(pixel_qc_path),
            "selected_tile_manifest": output_identity(selected_manifest),
            "candidate_tile_manifest": output_identity(candidate_manifest),
            "tile_class_counts": output_identity(tile_counts_path),
            "cells_parquet": output_identity(cells_path),
            "selected_tiles_overview": output_identity(sample_dir / "selected_tiles_overview.png"),
            "representative_overlay": output_identity(sample_dir / "representative_histoplus_overlay.png"),
        }
        receipt = {
            "status": "complete",
            "run_status": "inferred",
            "runner_version": RUNNER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_id": sample_id,
            "treatment": row["treatment"],
            "cohort": row["cohort"],
            "source_l0_tiff": str(slide_path),
            "conversion_receipt": str(conversion_receipt),
            "conversion_receipt_sha256": sha256_file(conversion_receipt),
            "organ_label_mask": str(label_path),
            "organ_label_mask_sha256": sha256_file(label_path),
            "source_mpp": source_mpp,
            "target_mpp": TARGET_MPP,
            "tile_size_px": TILE_SIZE,
            "tiles_per_organ": args.tiles_per_organ,
            "candidate_tiles": len(candidates),
            "selected_tiles": len(selected),
            "selected_tiles_by_organ": dict(Counter(selected["organ"])),
            "sampling_seed": args.seed,
            "sampling_design": "deterministic spatial farthest-point within each CV organ mask",
            "model": MODEL_NAME,
            "model_weight_file": str(args.weight_file.resolve()),
            "model_weight_sha256": weight_sha,
            "model_license": MODEL_LICENSE,
            "device": args.device,
            "gpu_name": torch.cuda.get_device_name(0) if args.device.startswith("cuda") else "",
            "gpu_peak_memory_bytes": int(torch.cuda.max_memory_allocated()) if args.device.startswith("cuda") else 0,
            "amp": False,
            "precision": "float32",
            "vectorization": "LazySlide-equivalent contour and probability operations on tight padded instance crops",
            "batch_size": args.batch_size,
            "postprocess_workers": args.batch_size,
            "cell_count": len(cells),
            "representative_tile_uid": representative_tile,
            "elapsed_seconds": time.time() - started,
            "biological_unit": "animal",
            "extrapolated_to_whole_organ": False,
            "cross_species_validated": False,
            "interpretation_label": "exploratory model-predicted nuclear phenotypes",
            "outputs": outputs,
        }
        atomic_json(sample_dir / "inference_receipt.json", receipt)
        return receipt
    except Exception as exc:
        failure = {
            "status": "failed",
            "runner_version": RUNNER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_id": sample_id,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
            "failed_result_encoded_as_zero": False,
        }
        atomic_json(sample_dir / "inference_failure.json", failure)
        return failure


def aggregate_results(receipts: list[dict[str, object]], output_dir: Path) -> None:
    phenotype_frames = []
    cv_frames = []
    tile_frames = []
    selected_frames = []
    for receipt in receipts:
        if receipt["status"] != "complete":
            continue
        outputs = receipt["outputs"]
        phenotype_frames.append(pd.read_csv(outputs["phenotype_summary"]["path"]))
        cv_frames.append(pd.read_csv(outputs["tile_cv_features"]["path"]))
        tile_frames.append(pd.read_csv(outputs["tile_class_counts"]["path"]))
        selected_frames.append(pd.read_csv(outputs["selected_tile_manifest"]["path"]))
    if phenotype_frames:
        phenotype = pd.concat(phenotype_frames, ignore_index=True)
        tile_counts = pd.concat(tile_frames, ignore_index=True)
        selected = pd.concat(selected_frames, ignore_index=True)
        key_columns = ["sample_id", "organ", "tile_uid"]
        expected = selected[key_columns].drop_duplicates()
        observed = tile_counts[key_columns].drop_duplicates()
        missing = expected.merge(observed, on=key_columns, how="left", indicator=True)
        missing = missing[missing["_merge"] == "left_only"].drop(columns="_merge")
        if len(missing):
            metadata = phenotype[["sample_id", "treatment", "cohort"]].drop_duplicates()
            missing = missing.merge(metadata, on="sample_id", how="left")
            missing["total_cells"] = 0
            for class_name in CLASS_NAMES.values():
                missing["count_" + class_name.lower().replace(" ", "_")] = 0
            missing = missing[tile_counts.columns]
            tile_counts = pd.concat([tile_counts, missing], ignore_index=True)
        tile_counts = tile_counts.sort_values(key_columns).reset_index(drop=True)
        if len(tile_counts) != len(selected):
            raise FigS3Error(
                f"Aggregate tile-count table has {len(tile_counts)} rows; {len(selected)} expected"
            )
        phenotype.to_csv(output_dir / "Figure_S3_HistoPLUS_phenotype_summary.csv", index=False)
        pd.concat(cv_frames, ignore_index=True).to_csv(output_dir / "Figure_S3_classical_tile_cv_features.csv", index=False)
        tile_counts.to_csv(output_dir / "Figure_S3_HistoPLUS_tile_class_counts.csv", index=False)
        selected.to_csv(output_dir / "Figure_S3_selected_tile_manifest.csv", index=False)


def main() -> int:
    args = parse_args()
    args.output_dir = args.output_dir.expanduser().resolve()
    args.export_dir = args.export_dir.expanduser().resolve()
    args.segmentation_dir = args.segmentation_dir.expanduser().resolve()
    args.weight_file = require_regular_file(args.weight_file, "HistoPLUS weight")
    if args.tiles_per_organ < 8:
        raise FigS3Error("--tiles-per-organ must be at least 8")
    if args.batch_size < 1:
        raise FigS3Error("--batch-size must be positive")
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        raise FigS3Error("CUDA was requested but torch.cuda.is_available() is false")
    rows = read_csv(require_regular_file(args.slide_manifest, "slide manifest"))
    requested = {validate_sample_id(value) for value in args.sample_id}
    if requested:
        rows = [row for row in rows if row["sample_id"] in requested]
        if requested != {row["sample_id"] for row in rows}:
            raise FigS3Error("One or more requested samples are absent")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    weight_sha = sha256_file(args.weight_file)
    model = None if args.prepare_only else build_model(args.weight_file)
    if model is not None:
        enable_parallel_postprocessing(model, args.batch_size)
    receipts: list[dict[str, object]] = []
    for index, row in enumerate(rows, start=1):
        receipt = process_slide(row, args, model, weight_sha)
        receipts.append(receipt)
        print(
            json.dumps(
                {
                    "progress": f"{index}/{len(rows)}",
                    "sample_id": row["sample_id"],
                    "status": receipt["status"],
                    "cells": receipt.get("cell_count"),
                    "elapsed_seconds": receipt.get("elapsed_seconds"),
                }
            ),
            flush=True,
        )
        if receipt["status"] == "failed" and args.device.startswith("cuda"):
            torch.cuda.empty_cache()
    if args.prepare_only:
        status = "prepared" if all(row["status"] == "prepared" for row in receipts) else "failed"
    else:
        aggregate_results(receipts, args.output_dir)
        status = "complete" if all(row["status"] == "complete" for row in receipts) else "failed"
    summary = {
        "status": status,
        "runner_version": RUNNER_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "requested_slides": len(rows),
        "complete_slides": sum(row["status"] == "complete" for row in receipts),
        "prepared_slides": sum(row["status"] == "prepared" for row in receipts),
        "failed_slides": [row["sample_id"] for row in receipts if row["status"] == "failed"],
        "model": MODEL_NAME,
        "weight_sha256": weight_sha,
        "model_license": MODEL_LICENSE,
        "cross_species_validated": False,
        "publication_interpretation": "exploratory model-predicted nuclear phenotypes",
        "receipts": receipts,
    }
    atomic_json(args.output_dir / "inference_summary.json", summary)
    print(json.dumps({key: summary[key] for key in ("status", "requested_slides", "complete_slides", "failed_slides")}, indent=2))
    return 0 if status in {"complete", "prepared"} else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
