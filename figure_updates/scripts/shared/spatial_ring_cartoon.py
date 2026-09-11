#!/usr/bin/env python3
"""Vector explanation of the six DAPI-quantile shells used in Figures 3 and 4.

The deterministic illustrative nuclei are not experimental observations. Their
shells use the same covariance normalization and quantiles as the analyses.
This module renders illustrations only; it never changes analytical inputs.
"""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

COLORS = ("#a36ade", "#c298eb", "#78adef", "#76c7eb", "#8bdfdf", "#b1efe5")
INK = "#263443"
POSITIVE = "#cf7b20"
VENTRICULAR_FLOOR_Y = -1.5
VENTRICLE_SCALE = 2.0
ILLUSTRATIVE_TISSUE_CENTRE_Y = -1.15
ILLUSTRATIVE_TISSUE_HALF_WIDTH = 2.7



def ventricle_half_width(y):
    """Illustrative third-ventricle contour, not an anatomical measurement."""
    y = VENTRICULAR_FLOOR_Y + (np.asarray(y) - VENTRICULAR_FLOOR_Y) / VENTRICLE_SCALE
    # Enlarge the contour about the fixed floor; the analysis never uses it.
    # Narrow superior neck, widening inferiorly to a rounded, nearly flat floor.
    taper = np.clip((1.8 - y) / 3.1, 0, 1)
    width = .055 + .660 * taper**2.2
    corner = .515 + .20 * np.sqrt(np.maximum(0, 1 - ((y + 1.30) / .20)**2))
    return VENTRICLE_SCALE * np.where(y < -1.30, corner, width)


def illustrative_geometry():
    """Return DAPI points around an empty third ventricle and quantile rings."""
    rng = np.random.default_rng(20260910)
    y0 = ILLUSTRATIVE_TISSUE_CENTRE_Y
    w = ILLUSTRATIVE_TISSUE_HALF_WIDTH
    candidates = rng.uniform([-w, y0 - 2.15], [w, y0 + 2.15], (2400, 2))
    x, y = candidates.T
    tissue = (x / w)**2 + ((y - y0) / 2.15)**2 < 1
    lumen = (y > VENTRICULAR_FLOOR_Y) & (np.abs(x) < ventricle_half_width(y))
    points = candidates[tissue & ~lumen][:360]
    centre = points.mean(axis=0)
    centred = points - centre
    eigenvalues, eigenvectors = np.linalg.eigh(np.cov(centred, rowvar=False))
    eigenvalues = np.maximum(eigenvalues, eigenvalues.max() * 1e-8)
    normalized = centred @ eigenvectors @ np.diag(1 / np.sqrt(eigenvalues))
    radius = np.sqrt(np.sum(normalized * normalized, axis=1))
    edges = np.quantile(radius, np.linspace(0, 1, 7))
    shells = np.searchsorted(edges[1:-1], radius, side="right")
    assert np.array_equal(np.bincount(shells), np.repeat(60, 6))
    return points, centre, normalized, radius, edges, shells, eigenvalues, eigenvectors


def caption(marker="NPY", language="en"):
    if language == "es":
        scope = ("Se utiliza un campo por animal." if marker == "NPY" else
                 "Los anillos se calculan por sección; se suman numeradores y denominadores del mismo anillo entre las secciones de cada animal antes de calcular su porcentaje.")
        endpoint_note = ("La única ecuación utiliza todos los núcleos c-FOS positivos para H o los núcleos dobles positivos c-FOS/NPY para I, contando los núcleos del anillo s. Ambos perfiles usan los mismos anillos DAPI y denominador. Las líneas finas de H–I muestran animales individuales; las líneas gruesas y los mapas de calor muestran la media aritmética de los porcentajes por animal. " if marker == "NPY" else
                         "En la ecuación, el numerador selecciona todos los núcleos c-FOS positivos o los dobles positivos del anillo s; ambos perfiles comparten el denominador DAPI. ")
        return ("Esquema ilustrativo de núcleos DAPI alrededor de una representación ampliada del tercer ventrículo, sin datos experimentales ni escala anatómica. Los puntos azules representan núcleos DAPI y el espacio blanco el lumen ventricular. El signo + señala la posición media de los núcleos, no un centro anatómico definido por el ventrículo. Los arcos discontinuos prolongan los límites de distancia sobre el lumen vacío; el contorno ventricular es solo una guía anatómica. "
                "Los seis anillos se construyen normalizando las distancias por la covarianza de las posiciones nucleares, para considerar la orientación y el alargamiento de su distribución. Sus límites corresponden a cuantiles de la distancia normalizada: cada anillo contiene aproximadamente un sexto de los núcleos. "
                "Los anillos internos 1–2, destacados en morado, contienen aproximadamente un tercio; superficie y grosor pueden diferir. La distribución nuclear ilustrativa se sitúa cerca del piso ventricular; el análisis no impone esta posición. Los anillos son posiciones relativas, no capas anatómicas. "
                f"Los límites no dependen de la positividad c-FOS o {marker}. Para cada anillo, el perfil muestra 100 × núcleos c-FOS positivos (o doble positivos c-FOS/{marker}) / todos los núcleos DAPI del mismo anillo. "
                + endpoint_note + scope + " Cada animal aporta un perfil de seis valores. En el dibujo hay 360 núcleos, 60 por anillo.")
    scope = ("One field contributes per animal." if marker == "NPY" else
             "Rings are constructed within each section; numerator and denominator counts for corresponding rings are summed across an animal's sections before calculating its percentages.")
    endpoint_note = ("The single equation uses all c-FOS-positive nuclei for H or c-FOS/NPY double-positive nuclei for I, counting nuclei within ring s. Both endpoints use the same DAPI-defined rings and denominator. Thin traces in H–I show individual animals; thick traces and heat maps show arithmetic means of the animal percentages. " if marker == "NPY" else
                     "In the equation, the numerator selects all c-FOS-positive or double-positive nuclei in ring s; both endpoints share the DAPI denominator. ")
    return ("Illustrative schematic of DAPI nuclei around an enlarged third ventricle, without experimental data or anatomical scale. Blue dots denote DAPI nuclei and the white space the ventricular lumen. The + marks the mean nuclear position, not an anatomical centre defined by the ventricle. Dashed arcs continue the distance boundaries across the empty lumen; the ventricular outline is only an anatomical guide. "
            "Six rings are constructed by covariance-normalizing distances from the mean nuclear position to account for the orientation and elongation of the nuclear distribution. Boundaries follow normalized-distance quantiles, so each ring contains approximately one-sixth of nuclei. "
            "Inner rings 1–2, highlighted in purple, contain approximately one-third; areas and thicknesses can differ. The illustrative nuclear distribution is positioned near the ventricular floor; the analysis does not impose this location. The rings describe relative nuclear positions, not anatomical layers. "
            f"Boundaries are independent of c-FOS or {marker} positivity. For each ring, the plotted profile is 100 × c-FOS-positive (or c-FOS/{marker} double-positive) nuclei / all DAPI nuclei in that ring. "
            + endpoint_note + scope + " Each animal contributes one six-value profile. The drawing contains 360 nuclei, 60 per ring.")


def draw_cartoon(fig, bounds=(0.02, 0.03, 0.96, 0.93), marker="NPY", language="en", font_scale=1.0):
    """Draw one tissue schematic; all methodological detail stays in the legend."""
    left,bottom,width,height = bounds
    ax = fig.add_axes([left,bottom+.28*height,width,.72*height])
    ax.set_axis_off()
    points, centre, normalized, radius, edges, shells, eigenvalues, eigenvectors = illustrative_geometry()
    theta = np.linspace(0, 2*np.pi, 361)
    unit = np.column_stack((np.cos(theta), np.sin(theta)))
    transform = np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
    for index in reversed(range(6)):
        boundary = (unit * edges[index+1]) @ transform + centre
        ax.add_patch(Polygon(boundary, closed=True, facecolor=COLORS[index],
                             edgecolor="#6822a7" if index < 2 else "#277ca9",
                             linewidth=1.5 if index < 2 else 1.05, zorder=1))
    ax.scatter(points[:,0], points[:,1], s=5.5, color="#1454d5", alpha=.98,
               linewidths=0, zorder=2)
    # Keep the ventricular lumen free of nuclei; rings remain defined by DAPI.
    y = np.linspace(VENTRICULAR_FLOOR_Y, VENTRICULAR_FLOOR_Y + VENTRICLE_SCALE * (2.7 - VENTRICULAR_FLOOR_Y), 800)
    half = ventricle_half_width(y)
    contour = np.concatenate((np.column_stack((-half,y)), np.column_stack((half,y))[::-1]))
    lumen_patch = Polygon(contour, closed=True, facecolor="white", edgecolor="#18232f", linewidth=1.8, zorder=3)
    ax.add_patch(lumen_patch)
    # Complete the mathematical ellipses across the empty lumen, using dashed
    # lines to distinguish distance boundaries from ventricular anatomy.
    for index in range(6):
        boundary = (unit * edges[index+1]) @ transform + centre
        line, = ax.plot(boundary[:,0], boundary[:,1],
                        color="#6822a7" if index < 2 else "#277ca9",
                        linewidth=.85, linestyle=(0, (3, 3)), alpha=.60, zorder=3.1)
        line.set_clip_path(lumen_patch)
    ax.plot(*centre, marker="+", color=INK, markersize=7, markeredgewidth=1.1, zorder=4)
    # Spread the labels across the lower arc so the smaller rings stay legible.
    for index, angle in enumerate(np.linspace(-.60, -2.65, 6)):
        direction = np.array([np.cos(angle), np.sin(angle)])
        whitened = direction @ eigenvectors @ np.diag(1 / np.sqrt(eigenvalues))
        scale = np.linalg.norm(whitened)
        low = 0 if index == 0 else edges[index]
        xy = centre + direction * ((low + edges[index+1]) / (2*scale))
        ax.text(*xy, str(index+1), ha="center", va="center", fontsize=10*font_scale,
                fontweight="bold", color=INK,
                bbox=dict(boxstyle="circle,pad=.11",fc=COLORS[index],ec="none"),zorder=5)
    es = language == "es"
    ax.scatter([4.02], [2.5], s=30, color="#1454d5")
    ax.text(4.24,2.5,"DAPI",va="center",fontsize=11*font_scale,color=INK)
    ax.plot(4.02, 1.65, marker="+", color=INK, markersize=7, markeredgewidth=1.1)
    ax.text(4.24,1.65,"Media DAPI" if es else "DAPI mean",va="center",fontsize=10*font_scale,color=INK)
    ax.annotate("Tercer ventrículo" if es else "Third ventricle", xy=(float(ventricle_half_width(1.0)),1.0), xytext=(4.0,.8),
                ha="left",va="center",fontsize=10*font_scale,color=INK,
                arrowprops=dict(arrowstyle="-",color="#18232f",lw=1.05,connectionstyle="angle3,angleA=0,angleB=90"))
    ax.annotate("Anillos internos 1–2" if es else "Inner rings 1–2",xy=(1.0,centre[1]-.5),xytext=(4.0,centre[1]-.5),
                ha="left",va="center",fontsize=10*font_scale,color="#6822a7",fontweight="bold",
                arrowprops=dict(arrowstyle="-",color="#6822a7",lw=1.1))
    ax.set(xlim=(-4.0,6.9),ylim=(-4.1,4.0),aspect="equal")
    formula = fig.add_axes([left+.015*width,bottom,width*.97,height*.23])
    formula.set_axis_off()
    font = dict(fontsize=10*font_scale,color=INK,va="center",transform=formula.transAxes)
    formula.text(.265,.48,"Ocurrencia (%) =" if es else "Occurrence (%) =",ha="right",**font)
    positive = ("c-FOS⁺ (H) o c-FOS⁺"+marker+"⁺ (I)" if es else
                "c-FOS⁺ (H) or c-FOS⁺"+marker+"⁺ (I)") if marker == "NPY" else (
                "c-FOS⁺ o c-FOS⁺"+marker+"⁺" if es else "c-FOS⁺ or c-FOS⁺"+marker+"⁺")
    numerator = positive+"\nNúcleos del anillo s" if es else positive+"\nnuclei in ring s"
    formula.text(.60,.85,numerator,ha="center",linespacing=1.1,**font)
    formula.plot([.285,.915],[.49,.49],color=INK,lw=.8,transform=formula.transAxes)
    formula.text(.60,.13,"Núcleos DAPI del anillo s" if es else "DAPI nuclei in ring s",ha="center",**font)
    formula.text(.937,.48,"×100",ha="left",**font)
    return [ax,formula]


def save_cartoon(outdir, marker="NPY", dpi=600):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for language, suffix in (("en", ""), ("es", "_spanish")):
        fig = plt.figure(figsize=(7.2, 5.4), facecolor="white")
        draw_cartoon(fig, bounds=(.015, .035, .97, .93), marker=marker, language=language, font_scale=1.4 if marker=="NPY" else 1.0)
        stem = outdir/f"Spatial_ring_definition_{marker}{suffix}"
        fig.savefig(stem.with_suffix('.pdf'), metadata={"Title": "DAPI ring definition: illustrative schematic", "CreationDate": None, "ModDate": None})
        fig.savefig(stem.with_suffix('.svg'), metadata={"Date": None})
        svg=stem.with_suffix('.svg')
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
        fig.savefig(stem.with_suffix('.png'), dpi=dpi)
        stem.with_suffix('.txt').write_text(caption(marker, language)+'\n')
        outputs.append(str(stem.with_suffix('.pdf')))
        plt.close(fig)
    geometry = illustrative_geometry()
    receipt = {"schema": "spatial_ring_cartoon_v3", "illustration_only": True,
               "experimental_data_used": False, "marker": marker,
               "illustrative_nuclei": len(geometry[0]),
               "view": "DAPI nuclei surrounding an illustrative third ventricle",
               "ventricle_shape": "narrow superior neck, inferior flare and rounded floor", "shell_counts": np.bincount(geometry[5]).tolist(),
               "inner_shells": [1, 2], "quantile_edges": geometry[4].tolist(),
               "illustrative_nuclear_mean_xy": geometry[1].tolist(),
               "ventricular_floor_y": VENTRICULAR_FLOOR_Y,
               "ventricle_scale_xy": VENTRICLE_SCALE,
               "layout": "enlarged third ventricle; inner rings near its floor; all six rings share the DAPI mean",
               "centre": "mean of eligible DAPI centroids", "distance": "covariance-normalized radius",
               "denominator": "eligible DAPI nuclei in the same shell",
               "formula_count": 1,
               "profile_endpoints": [
                   {"panel": "H" if marker == "NPY" else None, "endpoint": "cfos_occurrence", "numerator": "all c-FOS-positive DAPI nuclei", "denominator": "all DAPI nuclei in the same ring"},
                   {"panel": "I" if marker == "NPY" else None, "endpoint": "double_occurrence", "numerator": f"c-FOS/{marker} double-positive DAPI nuclei", "denominator": "all DAPI nuclei in the same ring"}],
               "plot_values": "100 * positive count / DAPI count, separately for each animal and ring",
               "group_summary": "arithmetic mean of animal percentages",
               "anatomical_outline_used_for_ring_assignment": False,
               "ring_boundaries_in_lumen": "dashed continuation; no nuclei in the illustrative lumen",
               "files": outputs}
    (outdir/f'Spatial_ring_definition_{marker}.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--marker', choices=('NPY', 'POMC'), default='NPY')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    print(json.dumps(save_cartoon(args.output_dir, args.marker, args.dpi), indent=2))
