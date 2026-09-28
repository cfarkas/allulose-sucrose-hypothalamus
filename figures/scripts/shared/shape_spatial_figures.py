#!/usr/bin/env python3
"""Render the replacement Fig. 3/4 closed-field and positive-cluster panels.

All coordinates and summary values come from shape_spatial_analysis.py outputs.
The example is a method display of one observed NPY field, never a fabricated
anatomical drawing. Statistical units are animals (NPY) and cages (POMC).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


CONDITIONS = ("Water", "Sucrose", "Allulose")
COLORS = {"Water": "#2f91b8", "Sucrose": "#c51b1f", "Allulose": "#159447"}
NAMES_ES = {"Water": "Agua", "Sucrose": "Sacarosa", "Allulose": "Alulosa"}
INK = "#263443"
PURPLE = "#7d459b"
POSITIVE = "#c66a20"


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 8.5,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": "#5e6974", "axes.linewidth": .6,
        "xtick.major.width": .6, "ytick.major.width": .6,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def save_figure(fig, stem: Path, dpi: int, title: str):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), metadata={
        "Title": title, "CreationDate": None, "ModDate": None,
    })
    fig.savefig(stem.with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    return [str(stem.with_suffix(ext)) for ext in (".pdf", ".png")]


def select_rows(frame, cohort, endpoint, region=None):
    selected = frame[(frame.cohort == cohort) & (frame.endpoint == endpoint)]
    if region is not None:
        selected = selected[selected.region == region]
    return selected


def finite_values(frame, column):
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
    return values[np.isfinite(values)]


def p_label(value):
    if not np.isfinite(value):
        return "NE"
    return "<0.001" if value < .001 else f"{value:.3f}"


def adjusted_ps(contrasts, cohort, endpoint, region, metric):
    block = select_rows(contrasts, cohort, endpoint, region)
    values = []
    for other in CONDITIONS[:2]:
        rows = block[(block.metric == metric) & (block.group_a == "Allulose") &
                     (block.group_b == other)]
        values.append(p_label(float(rows.p_holm.iloc[0])) if len(rows) else "NE")
    return " / ".join(values)


def plot_metric(ax, block, metric, language, annotate_n=True):
    """Display all independent units, arithmetic means, and matched null values."""
    es = language == "es"
    is_cluster = metric == "cluster_fraction"
    scale = 100 if is_cluster else 1
    displayed = []
    counts = []
    for index, condition in enumerate(CONDITIONS):
        group = block[block.condition == condition].sort_values("unit_id")
        observed = pd.to_numeric(group[metric], errors="coerce").to_numpy(float) * scale
        valid = np.isfinite(observed)
        count = int(valid.sum())
        counts.append((count, len(group)))
        jitter = np.linspace(-.145, .145, len(group)) if len(group) > 1 else np.zeros(len(group))
        # Units remain in deterministic order so the paired reference is traceable.
        if is_cluster:
            expected = pd.to_numeric(group.cluster_expected_fraction, errors="coerce").to_numpy(float) * 100
            x_obs = index + jitter - .055
            x_exp = index + jitter + .055
            for j in range(len(group)):
                if np.isfinite(observed[j]) and np.isfinite(expected[j]):
                    ax.plot([x_obs[j], x_exp[j]], [observed[j], expected[j]],
                            lw=.55, alpha=.30, color=COLORS[condition], zorder=1)
            exp_valid = np.isfinite(expected)
            ax.scatter(x_exp[exp_valid], expected[exp_valid], marker="D", s=14,
                       facecolors="white", edgecolors=COLORS[condition],
                       linewidths=.75, alpha=.80, zorder=3)
            displayed.extend(expected[exp_valid])
        else:
            x_obs = index + jitter
        ax.scatter(x_obs[valid], observed[valid], s=16, color=COLORS[condition],
                   edgecolors="white", linewidths=.25, alpha=.70, zorder=4)
        if count:
            mean = float(np.mean(observed[valid]))
            ax.plot([index - .24, index + .24], [mean, mean],
                    lw=1.7, color=COLORS[condition], solid_capstyle="round", zorder=5)
            ax.scatter([index], [mean], s=25, color=COLORS[condition],
                       edgecolors="white", linewidths=.4, zorder=6)
            displayed.extend(observed[valid])
        else:
            ax.text(index, .05, "NE", ha="center", va="bottom", fontsize=7,
                    color=COLORS[condition], transform=ax.get_xaxis_transform())
    if is_cluster:
        ax.set_ylim(-3, 104)
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel("En grupos (%)" if es else "Clustered (%)")
    else:
        ax.axhline(1, color="#a0a7ad", linestyle=(0, (3, 2)), lw=.7, zorder=0)
        ymax = max([1.0] + list(displayed))
        ax.set_ylim(0, ymax * 1.20)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        ax.set_ylabel("Enriquecimiento" if es else "Core enrichment")
    labels = []
    for condition, (n, total) in zip(CONDITIONS, counts):
        name = NAMES_ES[condition] if es else condition
        labels.append(name + (f"\nn={n}/{total}" if annotate_n else ""))
    ax.set_xticks(range(3), labels)
    ax.set_xlim(-.48, 2.48)
    ax.spines[["right", "top"]].set_visible(False)
    ax.grid(axis="y", color="#e8eaed", lw=.45, zorder=0)
    ax.set_axisbelow(True)
    return counts


def summary_key(fig, cohort, language, y):
    es = language == "es"
    unit = ("Animal" if cohort == "NPY" else ("Jaula" if es else "Cage"))
    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor="#707d88",
               markeredgecolor="white", markersize=4, label=unit),
        Line2D([], [], marker="o", color="#3d4852", lw=1.5, markersize=4,
               label="Media" if es else "Mean"),
        Line2D([], [], marker="D", color="none", markerfacecolor="white",
               markeredgecolor="#707d88", markersize=4,
               label="Esperado al azar" if es else "Random-label expectation"),
    ]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.53, y),
               ncol=3, frameon=False, fontsize=6.6, handlelength=1.5,
               columnspacing=1.3, handletextpad=.5)


def plot_npy(units, contrasts, endpoint, language, output_dir, dpi):
    es = language == "es"
    fig, axes = plt.subplots(1, 2, figsize=(6.60, 3.12))
    fig.subplots_adjust(left=.095, right=.985, bottom=.31, top=.76, wspace=.40)
    block = select_rows(units, "NPY", endpoint, "FIELD")
    for ax, metric in zip(axes, ("core_enrichment", "cluster_fraction")):
        plot_metric(ax, block, metric, language)
    axes[0].set_title("Contorno de densidad NPY (80 %)" if es else "NPY density contour (80%)", pad=8)
    axes[1].set_title("Grupos espaciales (≥2 células)" if es else "Spatial clusters (≥2 cells)", pad=8)
    title = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺NPY⁺"
    fig.text(.54, .94, title + (": organización espacial" if es else " spatial organization"),
             ha="center", va="center", fontsize=10, fontweight="bold")
    summary_key(fig, "NPY", language, .17)
    core = adjusted_ps(contrasts, "NPY", endpoint, "FIELD", "core_enrichment")
    cluster = adjusted_ps(contrasts, "NPY", endpoint, "FIELD", "cluster_excess")
    lead = "p Holm, alulosa vs. agua / sacarosa" if es else "Holm p, Allulose vs Water / Sucrose"
    body = (f"Contorno: {core}; exceso de agrupación: {cluster}" if es else
            f"Core: {core}; clustering excess: {cluster}")
    fig.text(.54, .092, lead + "\n" + body, ha="center", va="center", fontsize=6.7, linespacing=1.5)
    stem = "Figure3_Spatial_cFOS" + ("_NPY" if endpoint == "marker_cfos" else "") + "_occurrence"
    return save_figure(fig, output_dir / "Fig3" / language / stem, dpi, title + " closed-field and cluster analysis")


def plot_pomc(units, contrasts, endpoint, language, output_dir, dpi, height_factor=1.0):
    """Render existing POMC results with nominal and Holm comparison symbols."""
    es = language == "es"
    fig, axes = plt.subplots(2, 2, figsize=(6.15, 3.12 * height_factor), squeeze=False)
    fig.subplots_adjust(left=.115, right=.985, bottom=.16, top=.83, hspace=.70, wspace=.40)
    annotations = []
    for row, region in enumerate(("ARC", "ME")):
        block = select_rows(units, "POMC", endpoint, region)
        for col, metric in enumerate(("core_enrichment", "cluster_fraction")):
            ax = axes[row, col]
            plot_metric(ax, block, metric, language)
            ax.tick_params(axis="both", labelsize=6.2, pad=1.5)
            ax.yaxis.label.set_size(7)
            tested_metric = "cluster_excess" if metric == "cluster_fraction" else metric
            tested = select_rows(contrasts, "POMC", endpoint, region)
            low, high = ax.get_ylim(); span = high - low
            ax.set_ylim(low, high + .46 * span)
            for level, other in enumerate(("Sucrose", "Water")):
                rows = tested[(tested.metric == tested_metric) & (tested.group_a == "Allulose") & (tested.group_b == other)]
                if len(rows) != 1:
                    raise ValueError(f"Expected one existing POMC comparison: {endpoint} {region} {tested_metric} Allulose/{other}")
                result = rows.iloc[0]; p_raw = float(result.p_raw); p_holm = float(result.p_holm)
                def symbol(value):
                    return "NE" if not np.isfinite(value) else ("****" if value < .0001 else "***" if value < .001 else "**" if value < .01 else "*" if value < .05 else "NS")
                nominal_symbol = symbol(p_raw); holm_symbol = symbol(p_holm)
                label = f"nominal_p: {nominal_symbol}; Holm: {holm_symbol}"
                left = CONDITIONS.index(other); right = CONDITIONS.index("Allulose")
                y = high + (.055 + .23 * level) * span; tick = .022 * span
                ax.plot([left, left, right, right], [y, y + tick, y + tick, y], color=INK, lw=.65, clip_on=False)
                ax.text((left + right) / 2, y + tick + .013 * span, label, ha="center", va="bottom", fontsize=6.2, color="black")
                annotations.append(dict(endpoint=endpoint,region=region,displayed_metric=metric,tested_metric=tested_metric,group_a="Allulose",group_b=other,p_raw=p_raw if np.isfinite(p_raw) else None,p_holm=p_holm if np.isfinite(p_holm) else None,nominal_symbol=nominal_symbol,holm_symbol=holm_symbol,label=label,reason=str(result.reason)))
        axes[row, 0].text(-.26, .5, region, transform=axes[row, 0].transAxes,
                          va="center", ha="center", fontsize=8, fontweight="bold", rotation=90)
    axes[0, 0].set_title("Contorno POMC (80 %)" if es else "POMC contour (80%)", fontsize=8, pad=7)
    axes[0, 1].set_title("Grupos (≥2 células)" if es else "Clusters (≥2 cells)", fontsize=8, pad=7)
    title = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺POMC⁺"
    fig.text(.55, .965, title + (": organización espacial" if es else " spatial organization"),
             ha="center", va="center", fontsize=9, fontweight="bold")
    summary_key(fig, "POMC", language, .047)
    suffix = "_spanish" if es else ""
    stem = "Figure4_Spatial_cFOS" + ("_POMC" if endpoint == "marker_cfos" else "") + "_occurrence" + suffix
    outputs = save_figure(fig, output_dir / "Fig4" / stem, dpi, title + " closed-field and cluster analysis")
    receipt = output_dir / "Fig4" / (stem + "_comparisons.json")
    receipt.write_text(json.dumps({"height_factor":height_factor,"statistical_results_recomputed":False,"annotations":annotations,"symbol_key":{"NS":"p >= 0.05","*":"p < 0.05","**":"p < 0.01","***":"p < 0.001","****":"p < 0.0001","NE":"not estimable"},"comparison_label_format":"nominal_p: [symbol from p_raw]; Holm: [symbol from p_holm]","cluster_annotation_tests":"observed-minus-expected cluster fraction (cluster_excess), while plotted points retain observed and expected fractions"},indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return outputs


def plot_example(example_path, output_dir, language, dpi):
    es = language == "es"
    with np.load(example_path, allow_pickle=False) as example:
        x, y = example["x_um"], example["y_um"]
        positive = example["is_cfos"].astype(bool)
        gx, gy = example["core_x_um"], example["core_y_um"]
        mask = example["core_mask"].astype(float)
        positive_edges = np.asarray(example["positive_edges"], int).reshape(-1, 2)
        animal_id = str(example["animal_id"].item())
    if not (len(x) == len(y) == len(positive)) or len(x) == 0:
        raise ValueError("Example coordinates and positivity are invalid")
    if mask.shape != (len(gy), len(gx)):
        raise ValueError("Example grid must use (y, x) indexing")
    fig = plt.figure(figsize=(4.00, 3.00))
    ax = fig.add_axes([.04, .28, .93, .63])
    ax.set_axis_off()
    # Pad the binary support with an outside ring so every displayed boundary closes.
    dx = float(np.median(np.diff(gx)))
    dy = float(np.median(np.diff(gy)))
    px = np.r_[gx[0] - dx, gx, gx[-1] + dx]
    py = np.r_[gy[0] - dy, gy, gy[-1] + dy]
    padded_mask = np.pad(mask, 1, constant_values=0)
    ax.contourf(px, py, padded_mask, levels=[.5, 1.5], colors=["#f0e7f5"], zorder=0)
    ax.contour(px, py, padded_mask, levels=[.5], colors=[PURPLE], linewidths=1.0, zorder=1)
    ax.scatter(x[~positive], y[~positive], s=3.0, c="#a1a6ae", alpha=.65, lw=0, zorder=2)
    if len(positive_edges):
        xy = np.column_stack((x, y))
        segments = xy[positive_edges]
        ax.add_collection(LineCollection(segments, colors=POSITIVE, linewidths=.8, alpha=.8, zorder=3))
    ax.scatter(x[positive], y[positive], s=9.0, c=POSITIVE, edgecolors="white", linewidths=.25, zorder=4)
    xmin, xmax = float(min(x.min(), gx[mask.any(axis=0)].min())), float(max(x.max(), gx[mask.any(axis=0)].max()))
    ymin, ymax = float(min(y.min(), gy[mask.any(axis=1)].min())), float(max(y.max(), gy[mask.any(axis=1)].max()))
    spanx, spany = xmax - xmin, ymax - ymin
    ax.set_xlim(xmin - .035 * spanx, xmax + .035 * spanx)
    ax.set_ylim(ymax + .065 * spany, ymin - .035 * spany)
    ax.set_aspect("equal", adjustable="box")
    scale_x = xmax - .02 * spanx - 100
    scale_y = ymax + .035 * spany
    ax.plot([scale_x, scale_x + 100], [scale_y, scale_y], lw=1.6, color=INK)
    ax.text(scale_x + 50, scale_y - .012 * spany, "100 µm", ha="center", va="bottom", fontsize=6.4)
    fig.text(.50, .975, "Contorno y grupos espaciales" if es else "Density contour and spatial clusters",
             ha="center", va="top", fontsize=9, fontweight="bold")
    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor="#a1a6ae", markersize=3,
               label="NPY⁺"),
        Line2D([], [], marker="o", color="none", markerfacecolor=POSITIVE, markersize=4,
               label="c-FOS⁺NPY⁺"),
        Patch(facecolor="#f0e7f5", edgecolor=PURPLE,
              label="Contorno NPY 80 %" if es else "80% NPY contour"),
        Line2D([], [], color=POSITIVE, lw=1,
               label="Enlace entre positivas" if es else "Positive–positive link"),
    ]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.5, .20), ncol=2,
               fontsize=6.1, frameon=False, handlelength=1.2, columnspacing=1.2)
    fig.text(.50, .074, (("Ejemplo de campo: " if es else "Field example: ") + animal_id + "\n" +
                         ("Contorno independiente de c-FOS; grupos ≥2 células" if es else
                          "c-FOS-independent contour; clusters ≥2 cells")),
             ha="center", va="center", fontsize=6.0, linespacing=1.4)
    suffix = "_spanish" if es else ""
    return save_figure(fig, output_dir / "Fig3" / f"Spatial_shape_definition_NPY{suffix}", dpi,
                       "Observed NPY field: closed support and positive-cell graph")


def legend_text(cohort, endpoint, language):
    """Complete replacement-panel caption; upstream captions are not touched."""
    marker = cohort
    es = language == "es"
    endpoint_label = "c-FOS⁺" if endpoint == "cfos" else f"c-FOS⁺{marker}⁺"
    eligible_en = "all DAPI-associated nuclear centroids" if endpoint == "cfos" else f"legacy {marker}-associated nuclear centroids"
    eligible_es = "todos los centroides nucleares asociados a DAPI" if endpoint == "cfos" else f"los centroides nucleares asociados a {marker} en el análisis original"
    if es:
        return (f"Organización espacial de {endpoint_label}. El contorno de densidad de {marker} delimita el conjunto de mayor densidad que contiene al menos el 80 % de la masa de una estimación gaussiana de densidad (ancho de banda 40 µm, truncada a tres anchos de banda; malla de 5 µm), calculada a partir de todos los centroides {marker} positivos, independientemente de c-FOS. Puede incluir varias regiones cerradas desconectadas. Se calcula por sección y lado, se restringe al tejido permitido y requiere al menos cinco centroides del marcador. Izquierda: enriquecimiento de células positivas dentro del contorno, observado/esperado bajo etiquetado aleatorio condicionado por los estratos; la línea discontinua indica 1. Derecha: porcentaje de células positivas que pertenece a componentes conectados de dos o más células positivas en un grafo simétrico de hasta seis vecinos (menos en estratos escasos), con distancia máxima de 75 µm y enlaces limitados al tejido, lado y región. El grafo utiliza {eligible_es}. Los diamantes abiertos muestran la fracción esperada, calculada analíticamente bajo etiquetado aleatorio, manteniendo posiciones, grafo y recuentos positivos por estrato; las líneas tenues unen lo observado y esperado de cada unidad. "
                + ("Los puntos pequeños representan animales (tres por condición; un campo por animal), y las barras y puntos grandes sus medias aritméticas. " if cohort == "NPY" else
                   "Las filas muestran ARC y ME. Se agregan recuentos entre secciones de cada animal; los puntos pequeños representan medias de los animales disponibles de cada jaula, y las barras y puntos grandes las medias aritméticas entre jaulas, con igual peso por jaula. ")
                + "n indica unidades estimables/totales para la métrica graficada; NE indica no estimable, nunca ausencia de agrupación. La disponibilidad para las pruebas sobre el exceso de agrupación puede diferir y se detalla en los datos fuente. Las comparaciones exactas entre condiciones utilizan enriquecimiento y exceso de agrupación observado menos esperado; los valores p se ajustan con Holm sobre las 12 comparaciones planificadas de cada familia de métricas, incluyendo comparaciones no estimables como p=1 solo para el ajuste. Los recuentos de componentes por condición y resultados completos se proporcionan como datos fuente. "
                + ("La adquisición/plataforma está confundida con la condición en NPY y no se ha resuelto la independencia entre jaulas; los resultados son exploratorios y no separan ambos efectos. El contorno describe el campo observado, no un límite anatómico. " if cohort == "NPY" else
                   "Los contornos de densidad POMC se definen dentro de ARC o ME y no sustituyen sus límites anatómicos. ")
                + "Los componentes conectados son grupos de proximidad espacial, no circuitos funcionales.")
    return (f"Spatial organization of {endpoint_label}. The {marker} density contour encloses the highest-density set containing at least 80% of a Gaussian kernel density estimate (40 µm bandwidth, truncated at three bandwidths; 5 µm grid), constructed from all {marker}-positive centroids independently of c-FOS. It may include multiple disconnected closed regions. It is estimated by section and side, clipped to allowed tissue, and requires at least five marker centroids. Left: observed/expected positive-cell enrichment within the marker-density contour under stratum-conditioned random labeling; the dashed reference is 1. Right: the percentage of positive cells in connected components of at least two positive cells on a symmetric graph of up to six nearest neighbors (fewer in sparse strata) with a 75 µm edge cap, restricted to supported tissue, side and region. The graph uses {eligible_en}. Open diamonds show the analytically expected fraction under random labeling with fixed positions, graph and within-stratum positive counts; faint connectors pair each unit's observed and expected values. "
            + ("Small points are animals (three per condition; one field per animal); larger points and bars are arithmetic means. " if cohort == "NPY" else
               "Rows show ARC and ME. Counts are aggregated across each animal's sections; small points are means of the available animals in each cage, and larger points and bars are equal-weight arithmetic means across cages. ")
            + "n denotes estimable/total units for the displayed metric; NE denotes non-estimable, never absence of clustering. Availability for the clustering-excess tests can differ and is detailed in the source tables. Exact between-condition tests use core enrichment and observed-minus-expected clustering excess; p values are Holm-adjusted over the 12 planned comparisons within each metric family, counting non-estimable comparisons as p=1 for adjustment only. Component counts by condition and full statistical results are supplied as source data. "
            + ("NPY acquisition/platform is confounded with condition and cage independence remains unresolved; the results are exploratory and cannot separate these effects. The contour describes the observed field, not an anatomical boundary. " if cohort == "NPY" else
               "POMC density contours are defined within ARC or ME and do not replace their anatomical boundaries. ")
            + "Connected components denote spatial proximity groups, not functional circuits.")


def example_legend(language):
    if language == "es":
        return ("Ejemplo metodológico construido con el campo NPY observado WATER_NPY1. Los puntos grises son centroides NPY positivos sin c-FOS; los naranjas son dobles positivos c-FOS/NPY. El área morada delimita el contorno de densidad cerrado que contiene al menos el 80 % de la masa de la densidad NPY, calculado independientemente de la positividad c-FOS y recortado al soporte tisular permitido. Los segmentos naranjas conectan pares positivos del grafo simétrico de hasta seis vecinos (menos en estratos escasos) con distancia máxima de 75 µm; los componentes de dos o más células son los grupos espaciales contabilizados. El dibujo muestra un campo real para explicar las medidas; no es una comparación individual ni un borde anatómico. Barra, 100 µm. Los paneles H–I resumen todos los animales por condición.")
    return ("Method example using the observed WATER_NPY1 NPY field. Gray points are NPY-positive centroids without c-FOS; orange points are c-FOS/NPY double-positive centroids. The purple area delineates the closed density contour containing at least 80% of the NPY density mass, estimated independently of c-FOS positivity and clipped to allowed tissue support. Orange segments connect positive pairs in the symmetric graph of up to six nearest neighbors (fewer in sparse strata) with a 75 µm edge cap; components of two or more cells are counted as spatial clusters. The drawing uses one real field to explain the measures; it is not an individual statistical comparison or an anatomical border. Scale bar, 100 µm. Panels H–I summarize all animals by condition.")


def write_legends(output_dir):
    outputs = []
    for figure, cohort, entries in (("Fig3", "NPY", {"H": "cfos", "I": "marker_cfos"}),
                                    ("Fig4", "POMC", {"I": "cfos", "J": "marker_cfos"})):
        path = output_dir / figure / "legends.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        data = json.loads(path.read_text()) if path.exists() else {}
        for language in ("en", "es"):
            language_legends = data.setdefault(language, {})
            for panel, endpoint in entries.items():
                language_legends[panel] = legend_text(cohort, endpoint, language)
            if figure == "Fig3":
                language_legends["G"] = example_legend(language)
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        outputs.append(str(path))
    return outputs


def generate_figures(analysis_dir: Path, output_dir: Path, dpi=600):
    set_style()
    units = pd.read_csv(analysis_dir / "unit_metrics.csv")
    contrasts = pd.read_csv(analysis_dir / "treatment_contrasts.csv")
    outputs = []
    for language in ("en", "es"):
        for endpoint in ("cfos", "marker_cfos"):
            outputs.extend(plot_npy(units, contrasts, endpoint, language, output_dir, dpi))
            outputs.extend(plot_pomc(units, contrasts, endpoint, language, output_dir, dpi))
        outputs.extend(plot_example(analysis_dir / "example_npy.npz", output_dir, language, dpi))
    outputs.extend(write_legends(output_dir))
    receipt = {"analysis_dir": str(analysis_dir), "dpi": dpi,
               "statistical_units": {"NPY": "animal", "POMC": "cage"},
               "plots": ["closed marker field enrichment", "positive cells in spatial components >=2"],
               "example_is_observed_data": True, "files": outputs}
    receipt_path = output_dir / "figure_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    print(json.dumps(generate_figures(args.analysis_dir, args.output_dir, args.dpi), indent=2))
