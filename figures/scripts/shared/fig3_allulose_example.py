#!/usr/bin/env python3
"""Render the reviewed bilateral Allulose field with its accepted anatomy.

The existing native KDE rasters, accepted cell calls and fixed graph definition
are displayed in Panel D's reviewed coordinate frame. No segmentation or
inferential metric is changed.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle
import numpy as np
import pandas as pd
from scipy.ndimage import binary_fill_holes

import shape_spatial_analysis as analysis
import shape_spatial_figures as style

PAPER = Path(__file__).resolve().parents[2]
ANIMAL = "E8_FR6-4"
ACQUISITION = ANIMAL + "__S01"
# Panels H and I count every pair of activated cells closer than r, sweeping r
# from 20 to 150 um. The example draws four of those radii around one cell.
ILLUSTRATED_RADII_UM = (20., 70., 110., 150.)
PAIR_RADIUS_UM = max(ILLUSTRATED_RADII_UM)
# The drawn connections use the innermost radius. At larger r the counted pairs
# fill the field and the picture stops showing structure.
NETWORK_RADIUS_UM = min(ILLUSTRATED_RADII_UM)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mapped_xy(xy, matrix):
    homogeneous = np.column_stack((np.asarray(xy, float), np.ones(len(xy)))) @ matrix.T
    return homogeneous[:, :2] / homogeneous[:, 2:]


def load_cartoon():
    path = PAPER / "Fig3/02_make_cfos_npy_cartoons.py"
    spec = importlib.util.spec_from_file_location("fig3_reviewed_allulose_cartoon", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def prepare(analysis_dir, output_dir):
    analysis_dir, output_dir = Path(analysis_dir).resolve(), Path(output_dir).resolve()
    artifact_dir = output_dir / "examples"
    artifact_dir.mkdir(parents=True, exist_ok=True)
    cartoon = load_cartoon()
    data = cartoon.load_allulose_inputs(PAPER)
    display_dapi, _, _, tissue, lumen = cartoon._matched_layers(data, None)
    tissue_receipt = cartoon._matched_tissue_outline_receipt(data)
    lumen_receipt = cartoon._matched_lumen_receipt(data, None)
    tissue_outline = binary_fill_holes(tissue)
    display_to_native = np.asarray(data.spec.display_to_native_xy, dtype=float)
    native_to_display = np.linalg.inv(display_to_native)
    manifest = pd.read_csv(analysis_dir / "inputs/NPY_acquisition_manifest.csv")
    row = manifest[manifest.acquisition_id.eq(ACQUISITION)]
    if len(row) != 1:
        raise ValueError("Expected exactly one reviewed E8_FR6-4 acquisition")
    row = row.iloc[0]
    scale = (float(row.pixel_x_um), float(row.pixel_y_um))
    geometry_path = Path(row.geometry_path)
    with np.load(geometry_path) as source:
        geometry = {key: source[key] for key in ("tissue", "side", "roi")}
    frame = pd.read_csv(analysis_dir / "inputs/NPY_cells.csv.gz", low_memory=False)
    frame = frame[frame.acquisition_id.eq(ACQUISITION) & frame.region.eq("FIELD")].copy()
    for key in ("accepted_legacy", "is_cfos", "marker_channel_present", "in_tissue"):
        frame[key] = analysis.bool_column(frame, key)
    ix = np.rint(frame.x_px).astype(int).clip(0, geometry["roi"].shape[1]-1)
    iy = np.rint(frame.y_px).astype(int).clip(0, geometry["roi"].shape[0]-1)
    frame["geometry_side"] = geometry["side"][iy, ix]
    marker = frame[frame.accepted_legacy & frame.marker_channel_present & frame.in_tissue].copy()
    rows = pd.read_csv(analysis_dir / "stratum_metrics.csv")
    rows = rows[rows.acquisition_id.eq(ACQUISITION) & rows.endpoint.eq("marker_cfos") & rows.region.eq("FIELD")]
    memberships = pd.read_csv(analysis_dir / "cluster_membership.csv.gz")
    memberships = memberships[memberships.acquisition_id.eq(ACQUISITION) & memberships.endpoint.eq("marker_cfos") & memberships.region.eq("FIELD")]
    cell_ids, points, positives, edges_all, positive_edges_all, cores = [], [], [], [], [], []
    pair_edges_all, pair_counts_all, side_of = [], [], []
    source_contours, graph_checks = [], []
    offset = 0
    for side in sorted(marker.geometry_side.unique()):
        cells = marker[marker.geometry_side.eq(side)].sort_values("cell_id")
        native_um = cells[["x_um", "y_um"]].to_numpy(float)
        positive = cells.is_cfos.to_numpy(bool)
        support = geometry["tissue"] & (geometry["side"] == side) & (geometry["roi"] == 1)
        edges = analysis.graph_edges(native_um, 6, 75.)
        edges = edges[analysis.segment_valid(native_um, edges, support, scale)]
        _, sizes, positive_edges, cluster_sizes = analysis.observed_clusters(len(cells), positive, edges)
        # Every activated pair within the radius, with no neighbour limit and no
        # line-of-sight filter: this is what the scale-resolved analysis counts.
        activated = np.flatnonzero(positive)
        separation = np.linalg.norm(native_um[activated][:, None, :] - native_um[activated][None, :, :], axis=-1)
        close = np.argwhere(np.triu(separation <= NETWORK_RADIUS_UM, k=1))
        pair_edges = activated[close] if len(close) else np.empty((0, 2), dtype=np.int64)
        within = np.triu(np.ones_like(separation, dtype=bool), k=1)
        side_pair_counts = [int(((separation <= radius) & within).sum()) for radius in ILLUSTRATED_RADII_UM]
        hemi = cells.hemifield.iloc[0]
        recorded = rows[rows.hemifield.eq(hemi)]
        if len(recorded) != 1:
            raise ValueError("Missing recorded Allulose stratum")
        recorded = recorded.iloc[0]
        actual = {"N": len(cells), "m": int(positive.sum()), "graph_edges": len(edges),
                  "cluster_count": len(cluster_sizes), "clustered_cells": int((sizes > 0).sum())}
        if any(int(recorded[key]) != value for key, value in actual.items()):
            raise ValueError(f"Reconstructed display graph differs from recorded analysis: {hemi}")
        expected_sizes = memberships[memberships.hemifield.eq(hemi)].set_index("cell_id").cluster_size.to_dict()
        observed_sizes = {cell: int(size) for cell, size in zip(cells.cell_id, sizes) if size > 0}
        if observed_sizes != expected_sizes:
            raise ValueError(f"Display cluster membership changed: {hemi}")
        contour_path = analysis_dir / "contours" / f"NPY_{ACQUISITION}_FIELD_side{int(side)}.npz"
        with np.load(contour_path) as cache:
            if float(cache["bandwidth_um"]) != 40 or float(cache["grid_um"]) != 5 or float(cache["contour_mass"]) != .8:
                raise ValueError("Allulose display requires the current 40 µm / 80% KDE")
            gx, gy = np.meshgrid(cache["x_um"], cache["y_um"])
            native_px = np.column_stack((gx.ravel()/scale[0], gy.ravel()/scale[1]))
            mapped = mapped_xy(native_px, native_to_display)
            cores.append((mapped[:, 0].reshape(gx.shape), mapped[:, 1].reshape(gy.shape), cache["core_mask"].copy()))
        source_contours.append({"path": str(contour_path), "sha256": digest(contour_path)})
        cell_ids.extend(cells.cell_id)
        points.append(cells[["x_px", "y_px"]].to_numpy(float))
        positives.append(positive)
        edges_all.append(edges + offset)
        positive_edges_all.append(positive_edges + offset)
        pair_edges_all.append(pair_edges + offset)
        pair_counts_all.append(side_pair_counts)
        side_of.extend([int(side)] * len(cells))
        graph_checks.append({"hemifield": str(hemi), **actual, "membership_identical": True})
        offset += len(cells)
    native_xy = np.concatenate(points)
    display_xy = mapped_xy(native_xy, native_to_display)
    positive = np.concatenate(positives)
    height, width = data.spec.output_shape_yx
    visible = ((display_xy[:, 0] >= -.5) & (display_xy[:, 0] <= width-.5) &
               (display_xy[:, 1] >= -.5) & (display_xy[:, 1] <= height-.5))
    if (len(positive), int(visible.sum()), int(positive.sum()), int((positive & visible).sum())) != (494, 492, 131, 131):
        raise ValueError("Reviewed Allulose example cell counts changed")
    # A circle only stays a circle if the display transform is a similarity.
    singular = np.linalg.svd(native_to_display[:2, :2], compute_uv=False)
    anisotropy = float(singular[0] / singular[1])
    if anisotropy > 1.01:
        raise ValueError(f"Display transform is not a similarity: anisotropy {anisotropy:.4f}")
    pairs = np.concatenate(pair_edges_all).reshape(-1, 2)
    degree = np.bincount(pairs.ravel(), minlength=len(positive))
    radius_px = PAIR_RADIUS_UM / float(data.spec.um_per_px)
    fits = ((display_xy[:, 0] - radius_px >= -.5) & (display_xy[:, 0] + radius_px <= width - .5) &
            (display_xy[:, 1] - radius_px >= -.5) & (display_xy[:, 1] + radius_px <= height - .5))
    # A ring drawn across the ventricle would encircle cells that are never
    # paired with its centre, so the illustration only uses a cell whose largest
    # ring contains no activated cell from the opposite side.
    sides = np.asarray(side_of)
    scale_um = np.array(scale)
    unpaired_inside = np.zeros(len(positive), dtype=bool)
    for index in np.flatnonzero(positive & visible & fits):
        reach = np.linalg.norm((native_xy - native_xy[index]) * scale_um, axis=1) <= PAIR_RADIUS_UM
        unpaired_inside[index] = bool((reach & positive & (sides != sides[index])).any())
    eligible = np.flatnonzero(positive & visible & fits & (degree >= 4) & ~unpaired_inside)
    if not len(eligible):
        raise ValueError("No activated cell can carry the radius illustration")
    counts = degree[eligible]
    order = np.lexsort((np.asarray(cell_ids, dtype=str)[eligible], counts))
    focal = int(eligible[order[len(order) // 2]])
    # How the focal cell's neighbour count grows with r is the curve in H and I.
    separation_um = np.linalg.norm((native_xy - native_xy[focal]) * scale_um, axis=1)
    # Counted pairs never cross the ventricle, so neither may the illustration.
    activated_all = positive & (sides == sides[focal])
    activated_all[focal] = False
    focal_counts = [int(((separation_um <= radius) & activated_all).sum()) for radius in ILLUSTRATED_RADII_UM]
    artifact = artifact_dir / "example_allulose_npy.npz"
    np.savez_compressed(artifact, cell_id=np.asarray(cell_ids), native_xy_px=native_xy,
                        display_xy_px=display_xy, is_cfos=positive, visible_in_reviewed_frame=visible,
                        edges=np.concatenate(edges_all), positive_edges=np.concatenate(positive_edges_all),
                        pair_edges=pairs, pair_radius_um=np.array(PAIR_RADIUS_UM),
                        network_radius_um=np.array(NETWORK_RADIUS_UM),
                        pair_counts_by_radius=np.asarray(pair_counts_all, dtype=int).sum(axis=0),
                        illustrated_radii_um=np.asarray(ILLUSTRATED_RADII_UM, dtype=float),
                        focal_index=np.array(focal), focal_degree=np.array(int(degree[focal])),
                        focal_counts=np.asarray(focal_counts, dtype=int),
                        display_anisotropy=np.array(anisotropy),
                        core_display_x_px=np.stack([c[0] for c in cores]),
                        core_display_y_px=np.stack([c[1] for c in cores]), core_mask=np.stack([c[2] for c in cores]),
                        tissue_mask=tissue_outline, lumen_mask=lumen,
                        display_to_native_xy=display_to_native, native_to_display_xy=native_to_display,
                        um_per_display_px=np.array(data.spec.um_per_px),
                        animal_id=np.array(ANIMAL), condition=np.array("Allulose"))
    hil_path = PAPER / cartoon.VENTRICLE_HIL_DEFAULT_RELATIVE / "allulose_ventricle_receipt.json"
    provenance = {"schema": "fig3_allulose_anatomical_example_v1", "animal_id": ANIMAL,
                  "condition": "Allulose", "same_field_as_panel": "D", "analysis_changed": False,
                  "new_segmentation": False, "display_frame_yx": [height, width],
                  "orientation": "Unchanged reviewed Panel-D frame; dorsal at top, ventral median eminence at bottom",
                  "display_to_native_xy": display_to_native.tolist(), "native_to_display_xy": native_to_display.tolist(),
                  "native_pixel_um_xy": list(scale), "display_um_per_px": data.spec.um_per_px,
                  "display_rule": "Complete reviewed bilateral window; exact homogeneous point transform; accepted display-only tissue and ventricular outlines",
                  "marker_cells_native": len(positive), "marker_cells_visible": int(visible.sum()),
                  "double_positive_cells_native": int(positive.sum()), "double_positive_cells_visible": int((positive & visible).sum()),
                  "graph_checks": graph_checks, "contours": source_contours,
                  "tissue_outline": tissue_receipt, "ventricle_outline": lumen_receipt,
                  "hil_receipt": {"path": str(hil_path), "sha256": digest(hil_path)},
                  "native_sources": {key: {"path": str(source.path), "sha256": source.sha256} for key, source in data.sources.items()},
                  "input_sha256": {name: digest(analysis_dir/name) for name in
                      ("inputs/NPY_cells.csv.gz", "inputs/NPY_acquisition_manifest.csv", "stratum_metrics.csv", "cluster_membership.csv.gz")},
                  "geometry": {"path": str(geometry_path), "sha256": digest(geometry_path)},
                  "cartoon_script_sha256": digest(PAPER / "Fig3/02_make_cfos_npy_cartoons.py"),
                  "module_sha256": digest(__file__), "example_sha256": digest(artifact),
                  "legacy_example_preserved": {"path": str(analysis_dir / "example_npy.npz"), "sha256": digest(analysis_dir / "example_npy.npz")}}
    receipt = artifact_dir / "example_allulose_npy_provenance.json"
    receipt.write_text(json.dumps(provenance, indent=2, ensure_ascii=False)+"\n")
    return artifact, receipt


def plot_example(example_path, output_dir, language, dpi):
    es = language == "es"
    with np.load(example_path) as source:
        data = {key: source[key] for key in source.files}
    xy = data["display_xy_px"]
    positive = data["is_cfos"].astype(bool)
    visible = data["visible_in_reviewed_frame"].astype(bool)
    height, width = data["tissue_mask"].shape
    # Matching panel E's aspect fills the master row exactly: the cartoon is
    # limited by width, so extra page height would only add white space.
    fig = plt.figure(figsize=(4., 4. * 252.847 / 311.052))
    ax = fig.add_axes([.012, .085, .775, .80])
    ax.set_axis_off()
    field = Rectangle((-.5, -.5), width, height, transform=ax.transData)
    ax.contourf(data["tissue_mask"].astype(float), levels=[.5, 1.5], colors=["#edf0f3"], zorder=0)
    ax.contourf(data["lumen_mask"].astype(float), levels=[.5, 1.5], colors=["white"], zorder=4.5)
    # The accepted outlines are display annotations only; they do not remove cells.
    tissue = np.pad(data["tissue_mask"].astype(float), 1)
    ax.contour(np.arange(-1, width+1), np.arange(-1, height+1), tissue,
               levels=[.5], colors=["#667482"], linewidths=.8, zorder=4)
    ax.contour(data["lumen_mask"].astype(float), levels=[.5], colors=["#1d2631"], linewidths=.85, zorder=5)
    ax.scatter(xy[~positive & visible, 0], xy[~positive & visible, 1], s=2.8,
               color="#828e99", alpha=.75, linewidths=0, zorder=6)
    radii = [float(r) for r in data["illustrated_radii_um"]]
    counts = [int(c) for c in data["focal_counts"]]
    network_um = float(data["network_radius_um"])
    focal = int(data["focal_index"])
    scale_um = float(data["um_per_display_px"])
    network = LineCollection(xy[data["pair_edges"]], colors=style.POSITIVE,
                            linewidths=.5, alpha=.55, zorder=7)
    network.set_clip_path(field)
    ax.add_collection(network)
    ax.scatter(xy[positive & visible, 0], xy[positive & visible, 1], s=7.5,
               color=style.POSITIVE, edgecolors="white", linewidths=.22, zorder=8)
    # One cell carries the distance key: concentric circles at four of the radii
    # the analysis sweeps. Every activated cell inside a circle is one counted
    # pair with the centre cell, and the analysis repeats that for every cell.
    centre = xy[focal]
    for radius in radii:
        ax.add_patch(Circle(tuple(centre), radius / scale_um, facecolor="none", edgecolor="black",
                            linewidth=1.35, linestyle=(0, (4, 2.4)), zorder=8.6))
    ax.scatter([centre[0]], [centre[1]], s=20, facecolor="white", edgecolors="black",
               linewidths=1.1, zorder=9.2)
    for index, radius in enumerate(radii):
        top = centre[1] - radius / scale_um
        ax.text(centre[0], top - 3, f"{radius:g}" + (" µm" if index == len(radii) - 1 else ""),
                ha="center", va="bottom", fontsize=6.6, fontweight="bold", color="black", zorder=9.6,
                bbox={"boxstyle": "round,pad=.08", "facecolor": "white", "edgecolor": "none", "alpha": .92})
    ax.text(684, 152, "3V", ha="center", va="center", fontsize=9.75, fontweight="bold",
            color="black", zorder=9.5,
            bbox={"boxstyle": "round,pad=.14", "facecolor": "white", "edgecolor": "none", "alpha": .82})
    ax.annotate("ME", xy=(755, height-40), xytext=(755, height+94), ha="center", va="center",
                fontsize=9.75, fontweight="bold", color="black", arrowprops={"arrowstyle": "-", "lw": .65, "color": "black"}, zorder=9)
    length = 100 / float(data["um_per_display_px"])
    x0, y0 = width-length-35, height+82
    ax.plot([x0, x0+length], [y0, y0], color="black", linewidth=1.5)
    ax.text(x0+length/2, y0-14, "100 µm", ha="center", va="bottom", fontsize=9.15, color="black")
    ax.set_xlim(-15, width+15)
    ax.set_ylim(height+130, -15)
    ax.set_aspect("equal")
    fig.text(.392, .985, "Pares de células activadas\ndentro de r" if es else "Pairs of activated cells\nwithin r",
             ha="center", va="top", fontsize=11.6, fontweight="bold", color="black")
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor="#828e99", markersize=3, label="NPY⁺"),
               Line2D([], [], marker="o", color="none", markerfacecolor=style.POSITIVE, markersize=3.8, label="c-FOS⁺NPY⁺"),
               Line2D([], [], color="#667482", lw=.9, label="Tejido" if es else "Tissue"),
               Line2D([], [], color="#1d2631", lw=.9, label="Ventrículo" if es else "Ventricle"),
               Line2D([], [], color=style.POSITIVE, lw=.9,
                      label=(f"Par < {network_um:g} µm" if es else f"Pair < {network_um:g} µm")),
               Line2D([], [], color="black", lw=1.2, linestyle=(0, (4, 2.4)),
                      label="Radio r" if es else "Radius r")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(.798, .52), ncol=1, fontsize=6.3,
               frameon=False, handlelength=1.0, handletextpad=.4, labelspacing=.85, labelcolor="black")
    fig.text(.392, .028, (f"E8_FR6-4:Alulosa, r = {radii[0]:g}, {radii[1]:g}, {radii[2]:g}, {radii[3]:g} µm" if es else
                          f"E8_FR6-4:Allulose, r = {radii[0]:g}, {radii[1]:g}, {radii[2]:g}, {radii[3]:g} µm"),
             ha="center", va="center", fontsize=7.6, color="black")
    suffix = "_spanish" if es else ""
    return style.save_figure(fig, Path(output_dir)/"Fig3"/f"Spatial_shape_definition_NPY{suffix}", dpi,
                             "Reviewed bilateral Allulose hypothalamus: fixed tissue and ventricular contours, ventral ME at bottom")


def legend(language, example_path=None):
    es = language == "es"
    radii, counts, network, segments = ILLUSTRATED_RADII_UM, None, NETWORK_RADIUS_UM, None
    if example_path is not None:
        with np.load(example_path) as source:
            radii = [float(r) for r in source["illustrated_radii_um"]]
            counts = [int(c) for c in source["focal_counts"]]
            network = float(source["network_radius_um"])
            segments = int(len(source["pair_edges"]))
    listed = ", ".join(f"{r:g}" for r in radii)
    tally = ", ".join(str(c) for c in counts) if counts else ""
    drawn = f"{segments} " if segments else ""
    if es:
        return (f"Ejemplo del campo hipotalamico bilateral de alulosa E8_FR6-4, tambien mostrado en D. Se reutilizan el "
                f"marco registrado aceptado, el contorno externo del tejido y el poligono ventricular; el tercer "
                f"ventriculo (3V) es dorsal y la eminencia media (ME) es ventral. Estos trazos son anotaciones de "
                f"presentacion y no cambian las regiones ni los denominadores del analisis. Los puntos grises son "
                f"centroides NPY positivos sin c-FOS; los naranjas son centroides dobles positivos c-FOS/NPY. Los "
                f"{drawn}segmentos naranjas unen cada par de celulas activadas separado por menos de {network:g} um, el "
                f"radio mas pequeno que recorre el analisis: son pares contados, no conexiones anatomicas, y la "
                f"proximidad no establece conectividad funcional. A radios mayores los pares contados llenan el campo y "
                f"la imagen dejaria de mostrar estructura. Los cuatro circulos discontinuos negros, centrados en una "
                f"celula activada marcada en blanco, muestran a escala real cuatro de los radios que recorren los "
                f"paneles H e I, rotulados arriba en um: {listed} um. Cada celula activada dentro de un circulo forma un "
                f"par contado con la celula central; para esta celula los recuentos son {tally}, y su crecimiento con r "
                f"es lo que resumen las curvas de H e I. La celula se eligio como la de numero mediano de vecinos al "
                f"radio menor entre las activadas cuyo circulo mayor cabe en el recorte, de modo que es tipica y no "
                f"seleccionada por conveniencia. El analisis repite ese recuento para todas las celulas activadas, no "
                f"solo para la marcada, y dentro de cada lado del tejido. La ventana bilateral aceptada muestra 492 de "
                f"494 centroides NPY y los 131 dobles positivos; dos centroides solo NPY quedan fuera del recorte y "
                f"permanecen en el analisis. Barra de escala, 100 um. Es un ejemplo descriptivo del metodo, no una "
                f"comparacion estadistica individual.")
    return (f"Example from the bilateral Allulose hypothalamic field E8_FR6-4, also shown in D. The accepted registered "
            f"display frame, external tissue outline and ventricular polygon are reused; the third ventricle (3V) is "
            f"dorsal and the median eminence (ME) is ventral. These outlines are display annotations and do not change "
            f"analysis regions or denominators. Gray points are NPY-positive centroids without c-FOS; orange points are "
            f"c-FOS/NPY double-positive centroids. The {drawn}orange segments join every pair of activated cells closer "
            f"than {network:g} um, the smallest radius the analysis sweeps. They are counted pairs, not anatomical "
            f"connections, and proximity does not establish functional connectivity. At larger radii the counted pairs "
            f"fill the field and the picture would stop showing structure. The four black dashed circles, centred on one "
            f"activated cell marked in white, show four of the radii panels H and I sweep, at true scale and labelled "
            f"above in um: {listed} um. Every activated cell inside a circle forms one counted pair with the central "
            f"cell; for this cell the counts are {tally}, and how they grow with r is what the curves in H and I "
            f"summarise. The cell was chosen as the one with the median number of neighbours at the smallest radius "
            f"among activated cells whose largest circle fits inside the displayed frame, so it is typical rather than "
            f"selected for effect. The analysis repeats that count for every activated cell, not only the marked one, "
            f"and within each tissue side. The accepted bilateral window displays 492 of 494 NPY centroids and all 131 "
            f"double positives; two NPY-only centroids lie outside the original display crop and remain in the "
            f"analysis. Scale bar, 100 um. This is a descriptive method example, not an individual statistical "
            f"comparison.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    output = args.output_dir or args.analysis_dir
    style.set_style()
    example, receipt = prepare(args.analysis_dir, output)
    files = [str(example), str(receipt)]
    for language in ("en", "es"):
        files.extend(plot_example(example, output/"spatial_panels", language, args.dpi))
    print(json.dumps({"files": files}, indent=2))


if __name__ == "__main__":
    main()
