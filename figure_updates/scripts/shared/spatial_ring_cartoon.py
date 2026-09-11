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
ILLUSTRATIVE_TISSUE_CENTRE_Y = -1.35



def ventricle_half_width(y):
    """Illustrative third-ventricle contour, not an anatomical measurement."""
    y = np.asarray(y)
    # Narrow superior neck, widening inferiorly to a rounded, nearly flat floor.
    taper = np.clip((1.8 - y) / 3.1, 0, 1)
    width = .055 + .660 * taper**2.2
    corner = .515 + .20 * np.sqrt(np.maximum(0, 1 - ((y + 1.30) / .20)**2))
    return np.where(y < -1.30, corner, width)


def illustrative_geometry():
    """Return DAPI points around an empty third ventricle and quantile rings."""
    rng = np.random.default_rng(20260910)
    y0 = ILLUSTRATIVE_TISSUE_CENTRE_Y
    candidates = rng.uniform([-3.0, y0 - 2.15], [3.0, y0 + 2.15], (2400, 2))
    x, y = candidates.T
    tissue = (x / 3.0)**2 + ((y - y0) / 2.15)**2 < 1
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
        return ("Esquema ilustrativo de núcleos DAPI alrededor del tercer ventrículo, sin datos experimentales ni escala anatómica. Los puntos azules representan núcleos DAPI y el espacio blanco el lumen ventricular. El signo + señala la posición media de los núcleos, no un centro anatómico definido por el ventrículo. "
                "Los seis anillos se construyen normalizando las distancias por la covarianza de las posiciones nucleares, para considerar la orientación y el alargamiento de su distribución. Sus límites corresponden a cuantiles de la distancia normalizada: cada anillo contiene aproximadamente un sexto de los núcleos. "
                "Los anillos internos 1–2, destacados en morado, contienen aproximadamente un tercio; superficie y grosor pueden diferir. La distribución nuclear ilustrativa se sitúa cerca del piso ventricular; el análisis no impone esta posición. Los anillos son posiciones relativas, no capas anatómicas. "
                f"Los límites no dependen de la positividad c-FOS o {marker}. Para cada anillo, el perfil muestra 100 × núcleos c-FOS positivos (o doble positivos c-FOS/{marker}) / todos los núcleos DAPI del mismo anillo. "
                + scope + " Cada animal aporta un perfil de seis valores. En el dibujo hay 360 núcleos, 60 por anillo.")
    scope = ("One field contributes per animal." if marker == "NPY" else
             "Rings are constructed within each section; numerator and denominator counts for corresponding rings are summed across an animal's sections before calculating its percentages.")
    return ("Illustrative schematic of DAPI nuclei around the third ventricle, without experimental data or anatomical scale. Blue dots denote DAPI nuclei and the white space the ventricular lumen. The + marks the mean nuclear position, not an anatomical centre defined by the ventricle. "
            "Six rings are constructed by covariance-normalizing distances from the mean nuclear position to account for the orientation and elongation of the nuclear distribution. Boundaries follow normalized-distance quantiles, so each ring contains approximately one-sixth of nuclei. "
            "Inner rings 1–2, highlighted in purple, contain approximately one-third; areas and thicknesses can differ. The illustrative nuclear distribution is positioned near the ventricular floor; the analysis does not impose this location. The rings describe relative nuclear positions, not anatomical layers. "
            f"Boundaries are independent of c-FOS or {marker} positivity. For each ring, the plotted profile is 100 × c-FOS-positive (or c-FOS/{marker} double-positive) nuclei / all DAPI nuclei in that ring. "
            + scope + " Each animal contributes one six-value profile. The drawing contains 360 nuclei, 60 per ring.")


def draw_cartoon(fig, bounds=(0.02, 0.03, 0.96, 0.93), marker="NPY", language="en", font_scale=1.0):
    """Draw one tissue schematic; all methodological detail stays in the legend."""
    left,bottom,width,height = bounds
    ax = fig.add_axes([left,bottom+.20*height,width,.80*height])
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
    ax.scatter(points[:,0], points[:,1], s=9.5, color="#1454d5", alpha=.98,
               linewidths=0, zorder=2)
    # Keep the ventricular lumen free of nuclei; rings remain defined by DAPI.
    y = np.linspace(VENTRICULAR_FLOOR_Y, 2.7, 400)
    half = ventricle_half_width(y)
    contour = np.concatenate((np.column_stack((-half,y)), np.column_stack((half,y))[::-1]))
    ax.add_patch(Polygon(contour, closed=True, facecolor="white", edgecolor="#18232f", linewidth=1.8, zorder=3))
    ax.plot(*centre, marker="+", color=INK, markersize=7, markeredgewidth=1.1, zorder=4)
    # Put the six labels along a right-facing radius in the original tissue view.
    direction = np.array([np.cos(-.38),np.sin(-.38)])
    whitened = direction @ eigenvectors @ np.diag(1 / np.sqrt(eigenvalues))
    scale = np.linalg.norm(whitened)
    for index in range(6):
        low = 0 if index == 0 else edges[index]
        xy = centre + direction * ((low + edges[index+1]) / (2*scale))
        ax.text(*xy, str(index+1), ha="center", va="center", fontsize=10*font_scale,
                fontweight="bold", color=INK,
                bbox=dict(boxstyle="circle,pad=.11",fc=COLORS[index],ec="none"),zorder=5)
    es = language == "es"
    ax.scatter([3.62], [.60], s=30, color="#1454d5")
    ax.text(3.84,.60,"DAPI",va="center",fontsize=11*font_scale,color=INK)
    ax.annotate("Tercer ventrículo" if es else "Third ventricle", xy=(float(ventricle_half_width(1.0)),1.0), xytext=(3.60,-.25),
                ha="left",va="center",fontsize=10*font_scale,color=INK,
                arrowprops=dict(arrowstyle="-",color="#18232f",lw=1.05,connectionstyle="angle3,angleA=0,angleB=90"))
    ax.annotate("Anillos internos 1–2" if es else "Inner rings 1–2",xy=(1.0,centre[1]-.5),xytext=(3.60,centre[1]-.5),
                ha="left",va="center",fontsize=10*font_scale,color="#6822a7",fontweight="bold",
                arrowprops=dict(arrowstyle="-",color="#6822a7",lw=1.1))
    ax.set(xlim=(-3.6,6.3),ylim=(-4.0,1.2),aspect="equal")
    formula = fig.add_axes([left+.07*width,bottom,width*.86,height*.15])
    formula.set_axis_off()
    font = dict(fontsize=10*font_scale,color=INK,va="center",transform=formula.transAxes)
    formula.text(.34,.48,"Ocurrencia (%) =" if es else "Occurrence (%) =",ha="right",**font)
    formula.text(.60,.83,("Núcleos c-FOS⁺"+marker+"⁺ del anillo s") if es else ("c-FOS⁺"+marker+"⁺ nuclei in ring s"),ha="center",**font)
    formula.plot([.36,.84],[.49,.49],color=INK,lw=.8,transform=formula.transAxes)
    formula.text(.60,.14,"Núcleos DAPI del anillo s" if es else "DAPI nuclei in ring s",ha="center",**font)
    formula.text(.87,.48,"×100",ha="left",**font)
    return [ax,formula]


def save_cartoon(outdir, marker="NPY", dpi=600):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for language, suffix in (("en", ""), ("es", "_spanish")):
        fig = plt.figure(figsize=(10, 4.7), facecolor="white")
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
    receipt = {"schema": "spatial_ring_cartoon_v2", "illustration_only": True,
               "experimental_data_used": False, "marker": marker,
               "illustrative_nuclei": len(geometry[0]),
               "view": "DAPI nuclei surrounding an illustrative third ventricle",
               "ventricle_shape": "narrow superior neck, inferior flare and rounded floor", "shell_counts": np.bincount(geometry[5]).tolist(),
               "inner_shells": [1, 2], "quantile_edges": geometry[4].tolist(),
               "illustrative_nuclear_mean_xy": geometry[1].tolist(),
               "ventricular_floor_y": VENTRICULAR_FLOOR_Y,
               "layout": "inner rings near the ventricular floor; all six rings share the DAPI mean",
               "centre": "mean of eligible DAPI centroids", "distance": "covariance-normalized radius",
               "denominator": "eligible DAPI nuclei in the same shell", "files": outputs}
    (outdir/f'Spatial_ring_definition_{marker}.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--marker', choices=('NPY', 'POMC'), default='NPY')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    print(json.dumps(save_cartoon(args.output_dir, args.marker, args.dpi), indent=2))
