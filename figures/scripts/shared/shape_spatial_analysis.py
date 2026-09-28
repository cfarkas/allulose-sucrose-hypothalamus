#!/usr/bin/env python3
"""Closed marker footprints and c-FOS-positive cell clusters for Figures 3/4.

Uses the recovered, accepted-legacy nuclear associations. No source masks are
changed. Coordinates never cross acquisition, region, or tissue-side boundaries.
The graph construction follows spatial_reanalysis_20260919/src/spatial.py, with
k_eff=min(k,N-1) so true pairs can be evaluated in sparse strata. KDE clipping is
within tissue/region/side; unlike graph edges, individual kernel contributions
are not line-of-sight filtered. There are no individual-animal significance tests.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
import time

# Each acquisition is a process-level task; prevent nested BLAS thread pools.
for _thread_key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_key] = "1"

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.special import gammaln

SCHEMA = "marker-shape-clusters-v1"
CONDITIONS = ("Water", "Sucrose", "Allulose")
ENDPOINTS = ("cfos", "marker_cfos")
REGIONS = {"NPY": ("FIELD",), "POMC": ("ARC", "ME")}
METRICS = ("core_enrichment", "cluster_fraction", "cluster_expected_fraction", "cluster_excess")
MEMBERSHIP_COLUMNS = ["cohort", "endpoint", "condition", "animal_id", "cage", "acquisition_id", "region", "hemifield", "cell_id", "x_um", "y_um", "cluster_id", "cluster_size"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: dict) -> None:
    def convert(item):
        if isinstance(item, np.generic):
            return item.item()
        if isinstance(item, Path):
            return str(item)
        raise TypeError(type(item).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False) + "\n")


def bool_column(frame: pd.DataFrame, name: str) -> pd.Series:
    values = frame[name]
    if values.dtype == bool:
        return values
    if pd.api.types.is_numeric_dtype(values):
        return values.fillna(0).ne(0)
    return values.fillna(False).astype(str).str.lower().isin(("true", "1"))


def graph_edges(xy: np.ndarray, k: int = 6, cap: float = 75.0) -> np.ndarray:
    """Symmetric, unique kNN edges with an adaptive neighbor count for N<7."""
    n = len(xy)
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    distances, indices = cKDTree(xy).query(xy, k=min(k, n - 1) + 1)
    edges = set()
    for i, (neighbors, lengths) in enumerate(zip(indices, distances)):
        for j, length in zip(neighbors, lengths):
            if j != i and length <= cap:
                edges.add(tuple(sorted((i, int(j)))))
    return np.asarray(sorted(edges), dtype=np.int64).reshape(-1, 2)


def segment_valid(xy: np.ndarray, edges: np.ndarray, support: np.ndarray,
                  scale: tuple[float, float], step: float = 1.0) -> np.ndarray:
    """Check every edge including endpoints at <=1 µm steps in native tissue."""
    sx, sy = scale
    height, width = support.shape
    accepted = np.zeros(len(edges), dtype=bool)
    for i, (u, v) in enumerate(edges):
        a, b = xy[u], xy[v]
        t = np.linspace(0, 1, max(2, int(np.ceil(np.linalg.norm(a - b) / step)) + 1))
        points = a[None, :] + t[:, None] * (b - a)
        xx = np.rint(points[:, 0] / sx).astype(int)
        yy = np.rint(points[:, 1] / sy).astype(int)
        inside = (xx >= 0) & (xx < width) & (yy >= 0) & (yy < height)
        accepted[i] = bool(inside.all() and support[yy[inside], xx[inside]].all())
    return accepted


def expected_clustered_cells(n: int, m: int, edges: np.ndarray) -> float:
    """Analytic E[number of positive vertices with >=1 positive neighbor].

    Given a positive vertex i, absence of another positive in its degree d_i
    neighborhood has probability C(N-1-d_i,m-1)/C(N-1,m-1).
    """
    if n < 2 or m < 2 or len(edges) == 0:
        return 0.0
    if not 0 <= m <= n:
        raise ValueError("Invalid positive count")
    degree = np.bincount(edges.ravel(), minlength=n)
    available = n - 1 - degree
    choose = m - 1
    log_ratio = np.full(n, -np.inf)
    valid = available >= choose
    # Difference form avoids the common factorial term and improves precision.
    log_ratio[valid] = (gammaln(available[valid] + 1) - gammaln(available[valid] - choose + 1)
                        - gammaln(n) + gammaln(n - choose))
    probability_clustered = -np.expm1(np.minimum(log_ratio, 0.0))
    probability_clustered[~valid] = 1.0
    probability_clustered[degree == 0] = 0.0
    result = (m / n) * probability_clustered.sum()
    return float(np.clip(result, 0, m))


def observed_clusters(n: int, labels: np.ndarray, edges: np.ndarray):
    """Return component labels/sizes for activated components with >=2 cells."""
    positive_edges = edges[labels[edges[:, 0]] & labels[edges[:, 1]]]
    ids = np.full(n, -1, dtype=int)
    sizes = np.zeros(n, dtype=int)
    if len(positive_edges) == 0:
        return ids, sizes, positive_edges, []
    graph = coo_matrix((np.ones(len(positive_edges), dtype=np.uint8),
                       (positive_edges[:, 0], positive_edges[:, 1])), shape=(n, n))
    _, components = connected_components(graph, directed=False)
    counts = np.bincount(components[labels], minlength=n)
    found = []
    for component in np.flatnonzero(counts >= 2):
        take = labels & (components == component)
        cluster_id = len(found)
        ids[take] = cluster_id
        sizes[take] = int(take.sum())
        found.append(int(take.sum()))
    assert int(np.sum(sizes > 0)) == len(np.unique(positive_edges))
    return ids, sizes, positive_edges, found


def grid_geometry(shape: tuple[int, int], scale: tuple[float, float], grid: float):
    """Physical raster centers include an explicit exterior false padding row."""
    height, width = shape
    sx, sy = scale
    xs = np.arange(-grid, (width - 1) * sx + 2 * grid, grid, dtype=float)
    ys = np.arange(-grid, (height - 1) * sy + 2 * grid, grid, dtype=float)
    xx, yy = np.meshgrid(xs, ys)
    ix = np.rint(xx / sx).astype(int)
    iy = np.rint(yy / sy).astype(int)
    inside = (xx >= 0) & (yy >= 0) & (xx <= (width - 1) * sx) & (yy <= (height - 1) * sy)
    return xs, ys, xx, yy, ix, iy, inside


def edge_weights(raster_support: np.ndarray, bandwidth: float, grid: float, truncate: float = 3.0) -> np.ndarray:
    """Fraction of a kernel centred at each raster centre that lies in the support.

    Without this the estimate is the raw kernel sum, which is biased downwards
    wherever a kernel overlaps the tissue boundary. The arcuate marker band sits
    against that boundary, and the severity differs several-fold between
    acquisitions, so the bias is not a constant the contour can absorb.
    """
    weights = ndimage.gaussian_filter(raster_support.astype(float), sigma=bandwidth / grid,
                                      mode="constant", cval=0.0, truncate=truncate)
    return np.clip(weights, 0.0, 1.0)


def marker_footprint(marker_xy: np.ndarray, support: np.ndarray, scale: tuple[float, float],
                     bandwidth: float = 40, grid: float = 5, mass: float = .8,
                     min_anchors: int = 5, edge_correction: bool = True,
                     min_edge_weight: float = .15) -> dict:
    xs, ys, xx, yy, ix, iy, inside = grid_geometry(support.shape, scale, grid)
    raster_support = np.zeros(xx.shape, dtype=bool)
    raster_support[inside] = support[iy[inside], ix[inside]]
    density = np.zeros(xx.shape, dtype=float)
    core = np.zeros(xx.shape, dtype=bool)
    result = dict(x_um=xs, y_um=ys, density=density, support_mask=raster_support,
                  core_mask=core, available=False, reason="fewer_than_five_marker_anchors",
                  anchors=len(marker_xy), threshold=np.nan, achieved_mass=np.nan,
                  core_area_um2=0.0, component_count=0, edge_correction=bool(edge_correction),
                  min_edge_weight=float(min_edge_weight), edge_corrected_centers=0,
                  edge_weight_floor_hits=0, minimum_edge_weight_observed=np.nan)
    if len(marker_xy) < min_anchors:
        return result
    queries = np.c_[xx[raster_support], yy[raster_support]]
    if not len(queries):
        result["reason"] = "no_supported_raster_centers"
        return result
    # Exact fixed-width Gaussian sums at supported raster centers; no point/bin
    # snapping of marker positions. Kernels truncated at 3 sigma, as prior maps.
    tree = cKDTree(marker_xy)
    values = np.zeros(len(queries), dtype=float)
    neighbors = tree.query_ball_point(queries, 3 * bandwidth)
    for i, nearby in enumerate(neighbors):
        if nearby:
            distance2 = np.sum((marker_xy[nearby] - queries[i]) ** 2, axis=1)
            values[i] = np.exp(-distance2 / (2 * bandwidth * bandwidth)).sum()
    if edge_correction:
        weights = edge_weights(raster_support, bandwidth, grid)
        supported = weights[raster_support]
        floored = supported < min_edge_weight
        result.update(edge_corrected_centers=int(raster_support.sum()),
                      edge_weight_floor_hits=int(floored.sum()),
                      minimum_edge_weight_observed=float(supported.min()) if len(supported) else np.nan)
        values = values / np.maximum(supported, min_edge_weight)
    density[raster_support] = values
    positive = values[values > 0]
    if not len(positive):
        result["reason"] = "zero_supported_kernel_mass"
        return result
    ordered = np.sort(positive)[::-1]
    threshold = ordered[min(int(np.searchsorted(np.cumsum(ordered), mass * ordered.sum(), side="left")), len(ordered) - 1)]
    core[:] = raster_support & (density >= threshold)
    achieved = float(density[core].sum() / density.sum())
    _, component_count = ndimage.label(core, structure=np.ones((3, 3), dtype=int))
    assert not (core[0].any() or core[-1].any() or core[:, 0].any() or core[:, -1].any())
    assert achieved + 1e-12 >= mass
    result.update(available=True, reason="estimable", threshold=float(threshold),
                  achieved_mass=achieved, core_area_um2=float(core.sum() * grid * grid),
                  component_count=int(component_count))
    return result


def points_in_core(xy: np.ndarray, footprint: dict, grid: float) -> np.ndarray:
    """Nearest raster center, without interpolation across support boundaries."""
    if not len(xy):
        return np.zeros(0, dtype=bool)
    ix = np.rint((xy[:, 0] - footprint["x_um"][0]) / grid).astype(int)
    iy = np.rint((xy[:, 1] - footprint["y_um"][0]) / grid).astype(int)
    mask = footprint["core_mask"]
    valid = (ix >= 0) & (ix < mask.shape[1]) & (iy >= 0) & (iy < mask.shape[0])
    result = np.zeros(len(xy), dtype=bool)
    result[valid] = mask[iy[valid], ix[valid]]
    return result


def process_acquisition(task: dict) -> dict:
    frame = task["cells"]
    row = task["manifest"]
    cfg = task["config"]
    cohort = row["cohort"]
    acquisition = row["acquisition_id"]
    output = Path(task["output"])
    with np.load(row["geometry_path"]) as archive:
        geometry = {key: archive[key] for key in ("tissue", "side", "roi")}
    sx, sy = float(row["pixel_x_um"]), float(row["pixel_y_um"])
    scale = (sx, sy)
    common = {key: row[key] for key in ("cohort", "animal_id", "condition", "cage", "acquisition_id")}
    channel_present = bool(frame["marker_channel_present"].any()) if len(frame) else bool(row.get("pomc_channel_complete", True))
    strata = []
    memberships = []
    contours = []
    sensitivity = []
    example = {"xy": [], "labels": [], "edges": [], "positive_edges": [], "inside": []}
    example_offset = 0
    example_mask = example_density = example_support = None
    frame = frame.copy()
    if len(frame):
        ix = np.rint(frame.x_px).astype(int).clip(0, geometry["roi"].shape[1] - 1)
        iy = np.rint(frame.y_px).astype(int).clip(0, geometry["roi"].shape[0] - 1)
        frame["geometry_side"] = geometry["side"][iy, ix]
    else:
        frame["geometry_side"] = pd.Series(dtype=int)
    for region in REGIONS[cohort]:
        region_frame = frame.loc[frame.region.eq(region)]
        region_code = 1 if region in ("FIELD", "ARC") else 2
        sides = np.unique(geometry["side"][geometry["roi"] == region_code])
        sides = sides[sides > 0]
        if not len(sides):
            sides = [1]
        for side in sides:
            support = geometry["tissue"] & (geometry["side"] == side) & (geometry["roi"] == region_code)
            base = region_frame.loc[region_frame.geometry_side.eq(side)].sort_values("cell_id")
            base = base.loc[base.in_tissue]
            marker = base.loc[base.accepted_legacy & base.marker_channel_present]
            marker_xy = marker[["x_um", "y_um"]].to_numpy(float)
            footprint = marker_footprint(marker_xy, support, scale, cfg["bandwidth_um"], cfg["grid_um"], cfg["contour_mass"], cfg["minimum_marker_anchors"])
            sensitivity_footprints = {}
            if cfg.get("bandwidth_sensitivity", False):
                for bandwidth in sorted(set((20.0, float(cfg["bandwidth_um"]), 80.0))):
                    sensitivity_footprints[bandwidth] = (footprint if bandwidth == cfg["bandwidth_um"] else
                        marker_footprint(marker_xy, support, scale, bandwidth, cfg["grid_um"], cfg["contour_mass"], cfg["minimum_marker_anchors"]))
            hemi = str(base.hemifield.iloc[0]) if len(base) else f"SIDE_{int(side)}"
            stratum_name = f"{cohort}|{acquisition}|{region}|{hemi}"
            cache_name = f"{cohort}_{acquisition}_{region}_side{int(side)}.npz"
            np.savez_compressed(output / "contours" / cache_name,
                                x_um=footprint["x_um"], y_um=footprint["y_um"],
                                density=footprint["density"], core_mask=footprint["core_mask"],
                                support_mask=footprint["support_mask"],
                                threshold=np.array(footprint["threshold"]),
                                bandwidth_um=np.array(cfg["bandwidth_um"]), grid_um=np.array(cfg["grid_um"]),
                                contour_mass=np.array(cfg["contour_mass"]))
            contours.append(dict(**common, region=region, hemifield=hemi, side_code=int(side),
                                 contour_path=f"contours/{cache_name}", available=footprint["available"],
                                 reason=footprint["reason"], anchors=footprint["anchors"],
                                 threshold=footprint["threshold"], achieved_mass=footprint["achieved_mass"],
                                 core_area_um2=footprint["core_area_um2"], component_count=footprint["component_count"]))
            for endpoint in ENDPOINTS:
                cells = base if endpoint == "cfos" else marker
                xy = cells[["x_um", "y_um"]].to_numpy(float)
                labels = cells.is_cfos.to_numpy(bool)
                n, m = len(cells), int(labels.sum())
                edges = graph_edges(xy, cfg["k"], cfg["edge_cap_um"])
                edges = edges[segment_valid(xy, edges, support, scale)]
                expected_clustered = expected_clustered_cells(n, m, edges)
                cluster_ids, cluster_sizes, positive_edges, sizes = observed_clusters(n, labels, edges)
                missing_marker = endpoint == "marker_cfos" and not channel_present
                inside = points_in_core(xy, footprint, cfg["grid_um"])
                shape_available = bool(footprint["available"] and n)
                shape_n = n if shape_available else 0
                shape_m = m if shape_available else 0
                ncore = int(inside.sum()) if shape_available else 0
                observed_core = int(np.sum(inside & labels)) if shape_available else 0
                expected_core = m * ncore / n if shape_available else 0.0
                clustered = int(np.sum(cluster_sizes > 0))
                result = dict(**common, endpoint=endpoint, region=region, hemifield=hemi,
                              stratum_id=stratum_name + "|" + endpoint, marker_channel_present=channel_present,
                              N=np.nan if missing_marker else n, m=np.nan if missing_marker else m,
                              shape_N=shape_n, shape_m=shape_m, shape_core_N=ncore,
                              shape_available=shape_available, shape_reason=footprint["reason"],
                              core_observed=observed_core, core_expected=expected_core,
                              core_enrichment=observed_core / expected_core if expected_core > 0 else np.nan,
                              graph_N=n, graph_m=m, graph_edges=len(edges), k_effective=min(cfg["k"], max(0, n - 1)),
                              cluster_count=np.nan if missing_marker else len(sizes),
                              clustered_cells=np.nan if missing_marker else clustered,
                              cluster_expected_cells=np.nan if missing_marker else expected_clustered,
                              cluster_fraction=clustered / m if m else np.nan,
                              cluster_expected_fraction=expected_clustered / m if m else np.nan,
                              cluster_excess=(clustered - expected_clustered) / m if expected_clustered > 0 and m else np.nan,
                              largest_cluster=max(sizes, default=0),
                              graph_boundary_pass=True, graph_unique_pass=(len(edges) == len(np.unique(edges, axis=0))))
                strata.append(result)
                for bandwidth, varied in sensitivity_footprints.items():
                    varied_inside = points_in_core(xy, varied, cfg["grid_um"])
                    available = bool(varied["available"] and n)
                    varied_ncore = int(varied_inside.sum()) if available else 0
                    varied_observed = int(np.sum(varied_inside & labels)) if available else 0
                    varied_expected = m * varied_ncore / n if available else 0.0
                    # Graph and cluster statistics are copied exactly, never recomputed.
                    variant = dict(result)
                    variant.update(bandwidth_um=bandwidth, shape_N=n if available else 0,
                                   shape_m=m if available else 0, shape_core_N=varied_ncore,
                                   shape_available=available, shape_reason=varied["reason"],
                                   core_observed=varied_observed, core_expected=varied_expected,
                                   core_enrichment=varied_observed / varied_expected if varied_expected > 0 else np.nan)
                    sensitivity.append(variant)
                for index in np.flatnonzero(cluster_sizes > 0):
                    memberships.append(dict(**common, endpoint=endpoint, region=region, hemifield=hemi,
                                            cell_id=cells.iloc[index].cell_id, x_um=xy[index, 0], y_um=xy[index, 1],
                                            cluster_id=f"{stratum_name}|{endpoint}|C{cluster_ids[index] + 1:03d}",
                                            cluster_size=int(cluster_sizes[index])))
                if cohort == "NPY" and row["animal_id"] == "WATER_NPY1" and endpoint == "marker_cfos":
                    example["xy"].append(xy)
                    example["labels"].append(labels)
                    example["inside"].append(inside)
                    example["edges"].append(edges + example_offset)
                    example["positive_edges"].append(positive_edges + example_offset)
                    example_offset += n
                    if example_mask is None:
                        example_mask = footprint["core_mask"].copy()
                        example_density = footprint["density"].copy()
                        example_support = footprint["support_mask"].copy()
                    else:
                        example_mask |= footprint["core_mask"]
                        example_density += footprint["density"]
                        example_support |= footprint["support_mask"]
    if example["xy"]:
        xy = np.concatenate(example["xy"])
        np.savez_compressed(output / "example_npy.npz", x_um=xy[:, 0], y_um=xy[:, 1],
                            is_cfos=np.concatenate(example["labels"]),
                            marker_inside_core=np.concatenate(example["inside"]),
                            edges=np.concatenate(example["edges"]), positive_edges=np.concatenate(example["positive_edges"]),
                            core_x_um=footprint["x_um"], core_y_um=footprint["y_um"],
                            core_mask=example_mask, core_density=example_density, support_mask=example_support,
                            animal_id=np.array(row["animal_id"]), condition=np.array(row["condition"]),
                            pixel_x_um=np.array(sx), pixel_y_um=np.array(sy),
                            native_width_um=np.array(geometry["roi"].shape[1] * sx),
                            native_height_um=np.array(geometry["roi"].shape[0] * sy))
    return dict(strata=strata, memberships=memberships, contours=contours, sensitivity=sensitivity, acquisition_id=acquisition, cohort=cohort)


def aggregate_animals(strata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_columns = ["cohort", "endpoint", "animal_id", "region", "condition", "cage"]
    for key, selected in strata.groupby(group_columns, dropna=False, sort=True):
        row = dict(zip(group_columns, key))
        for column in ("N", "m", "shape_N", "shape_m", "shape_core_N", "core_observed", "core_expected",
                       "graph_N", "graph_m", "graph_edges", "cluster_count", "clustered_cells", "cluster_expected_cells"):
            row[column] = selected[column].sum(min_count=1)
        row["sections"] = selected.acquisition_id.nunique()
        row["strata"] = len(selected)
        row["shape_strata"] = int(selected.shape_available.sum())
        row["shape_sections"] = selected.loc[selected.shape_available, "acquisition_id"].nunique()
        row["core_enrichment"] = row["core_observed"] / row["core_expected"] if row["core_expected"] > 0 else np.nan
        m = row["m"]
        row["cluster_fraction"] = row["clustered_cells"] / m if m > 0 else np.nan
        row["cluster_expected_fraction"] = row["cluster_expected_cells"] / m if m > 0 else np.nan
        row["cluster_excess"] = ((row["clustered_cells"] - row["cluster_expected_cells"]) / m
                                  if m > 0 and row["cluster_expected_cells"] > 1e-12 else np.nan)
        row["largest_cluster"] = int(selected.largest_cluster.max())
        row["core_reason"] = ("estimable" if np.isfinite(row["core_enrichment"]) else
                               "marker_channel_unavailable" if not selected.marker_channel_present.any() else
                               "insufficient_marker_footprint_support" if row["shape_N"] == 0 else
                               "no_positive_cells_in_shape_eligible_strata" if row["shape_m"] == 0 else
                               "zero_expected_positive_cells_inside_footprint")
        row["cluster_reason"] = ("estimable" if np.isfinite(row["cluster_excess"]) else
                                  "marker_channel_unavailable" if np.isnan(m) else
                                  "no_positive_cells" if m == 0 else "zero_expected_clustered_cells")
        rows.append(row)
    return pd.DataFrame(rows)


def unit_metrics(animals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, selected in animals.groupby("cohort"):
        selected = selected.copy()
        selected["unit_id"] = selected.animal_id if cohort == "NPY" else selected.cage
        if selected.groupby("unit_id").condition.nunique().max() > 1:
            raise ValueError("An experimental unit belongs to multiple conditions")
        for key, subset in selected.groupby(["endpoint", "region", "condition", "unit_id"], sort=True):
            row = dict(zip(["endpoint", "region", "condition", "unit_id"], key))
            row.update(cohort=cohort, unit_type="animal" if cohort == "NPY" else "cage",
                       animals=";".join(sorted(subset.animal_id)), animal_count=len(subset))
            for metric in METRICS:
                row[metric] = subset[metric].mean()
                row[metric + "_animals"] = int(subset[metric].notna().sum())
            for count in ("N", "m", "shape_N", "shape_m", "cluster_count", "clustered_cells", "cluster_expected_cells"):
                row[count] = subset[count].sum(min_count=1)
            rows.append(row)
    return pd.DataFrame(rows)


def exact_comparison(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    na, nb = len(a), len(b)
    result = dict(n_a=na, n_b=nb, mean_a=a.mean() if na else np.nan, mean_b=b.mean() if nb else np.nan,
                  effect=a.mean() - b.mean() if na and nb else np.nan, p_raw=np.nan,
                  assignments=math.comb(na + nb, na), minimum_attainable_p=np.nan,
                  reason="fewer_than_two_eligible_units_in_a_group")
    if min(na, nb) < 2:
        return result
    values = np.r_[a, b]
    total = values.sum()
    differences = np.asarray([values[list(index)].sum() / na - (total - values[list(index)].sum()) / nb
                              for index in combinations(range(na + nb), na)])
    result.update(p_raw=float(np.mean(np.abs(differences) >= abs(result["effect"]) - 1e-12)),
                  minimum_attainable_p=float(np.mean(np.abs(differences) >= np.max(np.abs(differences)) - 1e-12)),
                  reason="conditional_exact_two_sided_unit_label_comparison")
    return result


def holm_planned(pvalues: np.ndarray) -> np.ndarray:
    """Unavailable planned tests consume a conservative p=1 multiplicity slot."""
    p = np.asarray(pvalues, dtype=float)
    finite = np.isfinite(p)
    working = np.where(finite, p, 1.0)
    order = np.argsort(working, kind="stable")
    adjusted = np.empty(len(p), dtype=float)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(working[order] * np.arange(len(p), 0, -1)))
    adjusted[~finite] = np.nan
    return adjusted


def treatment_contrasts(units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ("core_enrichment", "cluster_excess"):
        for cohort, regions in REGIONS.items():
            for endpoint in ENDPOINTS:
                for region in regions:
                    selected = units.loc[units.cohort.eq(cohort) & units.endpoint.eq(endpoint) & units.region.eq(region)]
                    for comparison in ("Water", "Sucrose"):
                        result = exact_comparison(selected.loc[selected.condition.eq("Allulose"), metric],
                                                  selected.loc[selected.condition.eq(comparison), metric])
                        result.update(cohort=cohort, endpoint=endpoint, region=region, metric=metric,
                                      group_a="Allulose", group_b=comparison, comparison="Allulose - " + comparison,
                                      family=metric, planned_family_size=12,
                                      unit_type="animal" if cohort == "NPY" else "cage",
                                      interpretation=("Exploratory joint treatment/acquisition comparison; NPY cage independence unresolved"
                                                      if cohort == "NPY" else "Exploratory comparison of equal-weight cage means"))
                        rows.append(result)
    result = pd.DataFrame(rows)
    for _, indices in result.groupby("family").groups.items():
        if len(indices) != 12:
            raise AssertionError("Multiplicity family no longer has twelve planned tests")
        result.loc[indices, "p_holm"] = holm_planned(result.loc[indices, "p_raw"])
    return result


def condition_summary(animals: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, selected in units.groupby(["cohort", "endpoint", "region", "condition"], sort=True):
        row = dict(zip(["cohort", "endpoint", "region", "condition"], key))
        subset = animals.loc[np.logical_and.reduce([animals[column].eq(value) for column, value in row.items()])]
        row.update(animals=subset.animal_id.nunique(), units=len(selected), unit_type=selected.unit_type.iloc[0])
        for metric in METRICS:
            row[metric + "_mean"] = selected[metric].mean()
            row[metric + "_units"] = int(selected[metric].notna().sum())
            row[metric + "_animals"] = int(subset[metric].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


def cluster_counts(animals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, selected in animals.groupby(["cohort", "endpoint", "region", "condition"], sort=True):
        row = dict(zip(["cohort", "endpoint", "region", "condition"], key))
        row.update(total_clusters=selected.cluster_count.sum(min_count=1),
                   clustered_positives=selected.clustered_cells.sum(min_count=1),
                   total_positive_cells=selected.m.sum(min_count=1),
                   animals=selected.animal_id.nunique(), animals_with_clusters=int(selected.cluster_count.gt(0).sum()),
                   marker_available_animals=int(selected.N.notna().sum()),
                   interpretation="Descriptive totals; animals/cages, not clusters, are the comparison units")
        rows.append(row)
    return pd.DataFrame(rows)


def validate_outputs(animals, units, strata, members, contours, contrasts, source_checks):
    checks = dict(source_cell_ids_unique=source_checks["cell_ids_unique"],
                  cell_coordinates_finite=source_checks["coordinates_finite"],
                  condition_assignment_consistent=source_checks["condition_assignment_consistent"],
                  graph_edges_unique=bool(strata.graph_unique_pass.all()),
                  graph_edges_tissue_constrained=bool(strata.graph_boundary_pass.all()),
                  positive_counts_in_range=bool(((strata.m.fillna(0) >= 0) & (strata.m.fillna(0) <= strata.N.fillna(0))).all()),
                  clustered_counts_in_range=bool((strata.clustered_cells.fillna(0) <= strata.m.fillna(0)).all()),
                  cluster_expected_counts_in_range=bool(((strata.cluster_expected_cells.fillna(0) >= -1e-9) & (strata.cluster_expected_cells.fillna(0) <= strata.m.fillna(0) + 1e-8)).all()),
                  shape_positive_counts_in_range=bool((strata.shape_m <= strata.shape_N).all()),
                  shape_counts_exclude_unsupported_strata=bool((strata.loc[~strata.shape_available, "shape_N"] == 0).all()),
                  minimum_component_size_two=bool(members.empty or (members.cluster_size >= 2).all()),
                  component_memberships_unique=bool(members.empty or not members.duplicated(["cohort", "endpoint", "cell_id"]).any()),
                  twelve_planned_contrasts_per_metric=bool(contrasts.groupby("family").size().eq(12).all()),
                  positive_exact_pvalues=bool(contrasts.p_raw.dropna().gt(0).all()),
                  zero_positive_fractions_missing=bool(animals.loc[animals.m.eq(0), "cluster_fraction"].isna().all()),
                  zero_positive_shape_enrichment_missing=bool(animals.loc[animals.shape_m.eq(0), "core_enrichment"].isna().all()),
                  supported_contours_capture_target_mass=bool(contours.loc[contours.available, "achieved_mass"].ge(.8 - 1e-12).all()),
                  cage_condition_assignment_consistent=bool(units.loc[units.cohort.eq("POMC")].groupby("unit_id").condition.nunique().le(1).all()))
    recovered = members.groupby(["cohort", "endpoint", "animal_id", "region"]).size().to_dict() if len(members) else {}
    checks["memberships_reconcile_with_animal_clustered_counts"] = all(
        np.isnan(row.clustered_cells) or int(row.clustered_cells) == recovered.get((row.cohort, row.endpoint, row.animal_id, row.region), 0)
        for row in animals.itertuples())
    checks["all_pass"] = all(checks.values())
    if not checks["all_pass"]:
        raise AssertionError("Validation failed: " + ", ".join(key for key, value in checks.items() if not value))
    return {"schema": SCHEMA, "checks": checks, "animal_rows": len(animals), "unit_rows": len(units),
            "strata": len(strata), "clustered_cell_memberships": len(members), "individual_pvalues_computed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Paper directory")
    parser.add_argument("--output", type=Path, default=Path("analyses/spatial_shape_20260920"))
    parser.add_argument("--bandwidth", type=float, default=40.0, help="Fixed Gaussian bandwidth in µm (primary default 40)")
    parser.add_argument("--bandwidth-sensitivity", action="store_true", help="Also summarize 20/80 µm contours by condition; reuse primary graphs; no new p-values")
    parser.add_argument("--workers", type=int, default=len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1)
    parser.add_argument("--source", type=Path, default=None,
                        help="Cell/manifest/geometry tree; defaults to the 2026-09-19 recovered reanalysis")
    parser.add_argument("--cohorts", default=",".join(REGIONS),
                        help="Comma-separated cohorts to analyze (default: every cohort)")
    args = parser.parse_args()
    if args.bandwidth <= 0:
        parser.error("--bandwidth must be positive")
    cohorts = tuple(name.strip() for name in args.cohorts.split(",") if name.strip())
    unknown = [name for name in cohorts if name not in REGIONS]
    if not cohorts or unknown:
        parser.error(f"--cohorts must name known cohorts {sorted(REGIONS)}; got {unknown or 'nothing'}")
    paper = args.root.resolve()
    output = args.output if args.output.is_absolute() else paper / args.output
    source = (args.source if args.source is None or args.source.is_absolute() else paper / args.source)
    if source is None:
        source = paper / "thesis_and_manuscript/Scientific_Reports_Segura_et_al_2026/spatial_reanalysis_20260919"
    source = source.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "contours").mkdir(exist_ok=True)
    config = dict(k=6, edge_cap_um=75.0, bandwidth_um=float(args.bandwidth), bandwidth_sensitivity=args.bandwidth_sensitivity, grid_um=5.0, contour_mass=.8,
                  minimum_marker_anchors=5, kernel_truncation_sigma=3,
                  identity="accepted_legacy", kernel_line_of_sight=False,
                  graph_line_of_sight_step_um=1.0, adaptive_k="min(6,N-1)",
                  cluster_minimum_cells=2, graph_null="fixed N,m and coordinates within acquisition/region/side",
                  contour_rule="80% of clipped kernel mass, descending raster density; ties included; exterior false padding",
                  shape_null="sum over supported strata of m_i*N_inside_i/N_i",
                  cluster_null="analytic degree-based hypergeometric probability of >=1 positive neighbor",
                  experimental_unit={"NPY": "animal; cage independence unresolved", "POMC": "equal-weight cage mean of animal ratios"},
                  multiplicity="12 planned tests per metric family, unavailable p=1 internally; display unavailable as NA",
                  note="KDE contributions clipped to region/side/tissue support; individual contributions are not line-of-sight filtered",
                  cohorts=list(cohorts))
    tasks = []
    inputs = []
    source_checks = dict(cell_ids_unique=True, coordinates_finite=True, condition_assignment_consistent=True)
    for cohort in cohorts:
        cell_path = source / f"data/{cohort}_cells.csv.gz"
        manifest_path = source / f"audit/{cohort}_acquisition_manifest.csv"
        cells = pd.read_csv(cell_path, low_memory=False)
        manifest = pd.read_csv(manifest_path)
        for column in ("is_cfos", "accepted_legacy", "marker_channel_present", "in_tissue", "include_acquisition"):
            cells[column] = bool_column(cells, column)
        manifest["include_acquisition"] = bool_column(manifest, "include_acquisition")
        if cells.cell_id.duplicated().any():
            raise ValueError(f"Duplicate source cell IDs: {cohort}")
        if not np.isfinite(cells[["x_um", "y_um"]].to_numpy(float)).all():
            raise ValueError(f"Nonfinite coordinates: {cohort}")
        if cells.groupby("animal_id").condition.nunique().max() > 1:
            raise ValueError(f"Inconsistent condition assignment: {cohort}")
        if cohort == "POMC" and cells.groupby("cage").condition.nunique().max() > 1:
            raise ValueError("POMC cage occurs in multiple conditions")
        for path in (cell_path, manifest_path):
            inputs.append(dict(path=str(path), sha256=sha256(path), size_bytes=path.stat().st_size))
        use = cells.loc[cells.include_acquisition & cells.region.isin(REGIONS[cohort])]
        for row in manifest.loc[manifest.include_acquisition].to_dict("records"):
            selected = use.loc[use.acquisition_id.eq(row["acquisition_id"])].copy()
            tasks.append(dict(cells=selected, manifest=row, config=config, output=str(output)))
    for name in ("spatial.py", "extract_cells.py"):
        path = source / "src" / name
        inputs.append(dict(path=str(path), sha256=sha256(path), size_bytes=path.stat().st_size))
    workers = min(max(args.workers, 1), len(tasks))
    print(f"Analyzing {len(tasks)} independent acquisitions with {workers} workers; BLAS threads=1", flush=True)
    start = time.time()
    strata, members, contours, sensitivity = [], [], [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(process_acquisition, task) for task in tasks]
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            strata.extend(result["strata"])
            members.extend(result["memberships"])
            contours.extend(result["contours"])
            sensitivity.extend(result["sensitivity"])
            print(f"{index}/{len(tasks)} {result['cohort']} {result['acquisition_id']}", flush=True)
    strata = pd.DataFrame(strata).sort_values(["cohort", "endpoint", "animal_id", "region", "acquisition_id", "hemifield"])
    members = pd.DataFrame(members, columns=MEMBERSHIP_COLUMNS).sort_values(["cohort", "endpoint", "animal_id", "region", "cluster_id", "cell_id"])
    contours = pd.DataFrame(contours).sort_values(["cohort", "animal_id", "region", "acquisition_id", "hemifield"])
    animals = aggregate_animals(strata)
    units = unit_metrics(animals)
    contrasts = treatment_contrasts(units)
    summaries = condition_summary(animals, units)
    counts = cluster_counts(animals)
    validation = validate_outputs(animals, units, strata, members, contours, contrasts, source_checks)
    if sensitivity:
        sensitivity_summary = []
        for bandwidth, subset in pd.DataFrame(sensitivity).groupby("bandwidth_um"):
            variant_animals = aggregate_animals(subset)
            variant_units = unit_metrics(variant_animals)
            summary = condition_summary(variant_animals, variant_units)
            summary["bandwidth_um"] = bandwidth
            summary["interpretation"] = "Descriptive contour bandwidth sensitivity; same primary graphs; no additional hypothesis tests"
            keep = ["cohort", "endpoint", "region", "condition", "bandwidth_um", "animals", "units", "unit_type",
                    "core_enrichment_mean", "core_enrichment_units", "core_enrichment_animals", "interpretation"]
            sensitivity_summary.append(summary[keep])
        pd.concat(sensitivity_summary, ignore_index=True).to_csv(output / "bandwidth_sensitivity_condition_summary.csv", index=False)
    # Geometry masks fix support and adjacency and are part of source provenance.
    # Hash once after analysis; no large raw image or segmentation files are hashed.
    for geometry_path in sorted({task["manifest"]["geometry_path"] for task in tasks}):
        path = Path(geometry_path)
        inputs.append(dict(path=str(path), sha256=sha256(path), size_bytes=path.stat().st_size, role="fixed_tissue_region_side_geometry"))
    for name, table in (("animal_metrics", animals), ("unit_metrics", units), ("condition_summary", summaries),
                        ("treatment_contrasts", contrasts), ("stratum_metrics", strata),
                        ("cluster_counts_by_condition", counts), ("contour_manifest", contours)):
        table.to_csv(output / f"{name}.csv", index=False)
    members.to_csv(output / "cluster_membership.csv.gz", index=False, compression="gzip")
    validation["elapsed_seconds"] = time.time() - start
    validation["workers"] = workers
    write_json(output / "validation.json", validation)
    metadata = dict(schema=SCHEMA, parameters=config, input_files=inputs,
                    script_sha256=sha256(Path(__file__)), output_directory=str(output),
                    source_directory=str(source), workers=workers, timestamp_utc=pd.Timestamp.now(tz="UTC").isoformat(),
                    statistics_are_exploratory=True, individual_reports_generated=False,
                    NPY_anatomy="Sampled field, not registered ARC", POMC_anatomy="Reviewed ARC and ME")
    write_json(output / "analysis_metadata.json", metadata)
    print(f"Complete: {output}; all validation checks passed; {validation['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
