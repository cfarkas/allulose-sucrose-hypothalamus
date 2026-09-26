#!/usr/bin/env python3
"""Figure 4 panels I and J: scale-resolved POMC clustering and its global test.

These replicate Figure 3's panels H and I for the POMC cohort in ARC. The
measurement, the exact conditional random-labelling null and the studentized
maximum-deviation global envelope test are identical; only the population and
the experimental unit differ. Figure 4 compares equal-weight cage means, so
each animal's standardized curve is averaged within its cage and the envelope
test enumerates cage-label allocations.

ME is not analysed. It retains 74 marker-associated cells with 28 activated
across 11 animals, and only 6 of those animals have two or more activated
cells, which is too sparse for a curve; forcing one would report noise.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import shape_spatial_figures as style

CONDITIONS = ("Water", "Sucrose", "Allulose")
REGION = "ARC"
# With two cages in one condition the permutation distribution of the
# between-condition F has a very long tail, so the 95% envelope can sit two
# orders of magnitude above the observed curve. On a linear axis the observed
# curve then disappears into the axis line and the panel shows nothing. A
# logarithmic axis is used above this ratio so that the observed curve, and how
# far below the envelope it lies, stay readable.
LOG_AXIS_RATIO = 8.0
PANELS = {"cfos": ("I", "Figure4_Spatial_cFOS_occurrence"),
          "marker_cfos": ("J", "Figure4_Spatial_cFOS_POMC_occurrence")}


def formatted(value) -> str:
    value = float(value)
    if not np.isfinite(value):
        return "NE"
    return f"{value:.2g}" if value < .001 else f"{value:.3f}"


def cage_curves(curves: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    block = curves[curves.endpoint.eq(endpoint) & curves.region.eq(REGION)].dropna(subset=["z"])
    wide = block.pivot_table(index=["condition", "cage", "animal_id"], columns="radius_um", values="z").dropna()
    return wide.groupby(level=["condition", "cage"]).mean()


def plot_panel(curves, envelope, test, endpoint, language, destination, dpi):
    es = language == "es"
    fig, axes = plt.subplots(1, 2, figsize=(6.60, 3.12))
    fig.subplots_adjust(left=.095, right=.985, bottom=.175, top=.785, wspace=.34)
    cages = cage_curves(curves, endpoint)
    band = envelope[envelope.endpoint.eq(endpoint) & envelope.region.eq(REGION)].sort_values("radius_um")
    row = test[test.endpoint.eq(endpoint) & test.region.eq(REGION)].iloc[0]
    radii = np.asarray([float(c) for c in cages.columns], dtype=float)

    left = axes[0]
    left.axhline(0, color="black", linestyle=(0, (3, 2)), lw=.6, zorder=0)
    for condition in CONDITIONS:
        if condition not in cages.index.get_level_values("condition"):
            continue
        group = cages.xs(condition, level="condition")
        for _, values in group.iterrows():
            left.plot(radii, values.to_numpy(float), color=style.COLORS[condition], lw=.7, alpha=.55, zorder=2)
        left.plot(radii, group.to_numpy(float).mean(axis=0), color=style.COLORS[condition], lw=2.1, zorder=3,
                  label=style.NAMES_ES[condition] if es else condition)
    # Headroom for the frameless legend, so no curve runs underneath it.
    drawn = cages.to_numpy(float)
    low, high = float(np.nanmin(drawn)), float(np.nanmax(drawn))
    span = max(high - low, 1e-9)
    left.set_ylim(min(low, 0.) - .06 * span, high + .30 * span)
    left.set_xlabel("Distancia r (µm)" if es else "Distance r (µm)", fontsize=8, fontweight="bold")
    left.set_ylabel("Agrupación frente al azar (DE)" if es else "Clustering vs chance (SD)",
                    fontsize=8, fontweight="bold")
    left.set_title("Curvas por jaula" if es else "Per-cage curves", pad=8, fontsize=8)
    left.legend(frameon=False, fontsize=6.4, loc="upper right", handlelength=1.4)

    right = axes[1]
    observed = band.observed_f.to_numpy(float)
    critical = band.critical_f.to_numpy(float)
    logarithmic = np.nanmax(critical) > LOG_AXIS_RATIO * max(np.nanmax(observed), 1e-9)
    if logarithmic:
        right.set_yscale("log")
        floor = max(min(np.nanmin(observed), np.nanmin(critical)) * .40, 1e-4)
        ceiling = np.nanmax(critical) * 3.6
    else:
        floor, ceiling = 0., max(np.nanmax(observed), np.nanmax(critical)) * 1.30
    right.fill_between(band.radius_um, critical, ceiling, color="#d8dce0", alpha=.55, lw=0, zorder=0,
                       label="Rechazo global" if es else "Global rejection")
    right.plot(band.radius_um, critical, color="black", linestyle=(0, (4, 2)), lw=.9, zorder=2,
               label="Envolvente 95 %" if es else "95% envelope")
    right.plot(band.radius_um, observed, color="black", lw=1.8, zorder=3,
               label="F observado" if es else "Observed F")
    exits = band.exits_envelope.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    if exits.any():
        right.scatter(band.radius_um.to_numpy()[exits], observed[exits], s=13, color="black", zorder=4)
    right.set_ylim(floor, ceiling)
    right.set_xlabel("Distancia r (µm)" if es else "Distance r (µm)", fontsize=8, fontweight="bold")
    axis_label = "F entre condiciones" if es else "Between-condition F"
    if logarithmic:
        axis_label += " (log)"
    right.set_ylabel(axis_label, fontsize=8, fontweight="bold")
    right.set_title("Prueba de envolvente global" if es else "Global envelope test", pad=8, fontsize=8)
    # A white-backed frame, because the envelope can pass through the upper
    # left corner at these sample sizes and the key must stay readable.
    right.legend(loc="upper left", fontsize=6.1, handlelength=1.4, frameon=True,
                 facecolor="white", edgecolor="none", framealpha=.88, borderpad=.3)

    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺POMC⁺"
    fig.text(.54, .962, label + (": organización espacial en ARC" if es else " spatial organization in ARC"),
             ha="center", va="center", fontsize=10, fontweight="bold")
    counts = "/".join(str(int(row[f"n_{condition.lower()}"])) for condition in CONDITIONS)
    span = f"{radii.min():g}–{radii.max():g} µm"
    if es:
        headline = (f"Envolvente global, {span} ({counts} jaulas): p={formatted(row.p_value)}; "
                    f"p Holm={formatted(row.p_holm_two_endpoint_family)}")
    else:
        headline = (f"Global envelope, {span} ({counts} cages): p={formatted(row.p_value)}; "
                    f"Holm p={formatted(row.p_holm_two_endpoint_family)}")
    fig.text(.54, .885, headline, ha="center", va="center", fontsize=7.4, fontweight="bold")
    suffix = "_spanish" if es else ""
    stem = PANELS[endpoint][1] + suffix
    return style.save_figure(fig, Path(destination) / stem, dpi,
                             label + " scale-resolved clustering and global envelope test in ARC")


def legend_text(curves, envelope, test, endpoint, language):
    es = language == "es"
    row = test[test.endpoint.eq(endpoint) & test.region.eq(REGION)].iloc[0]
    band = envelope[envelope.endpoint.eq(endpoint) & envelope.region.eq(REGION)].sort_values("radius_um")
    radii = band.radius_um.to_numpy(float)
    exits = band.exits_envelope.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    cages = cage_curves(curves, endpoint)
    animals = curves[curves.endpoint.eq(endpoint) & curves.region.eq(REGION)].dropna(subset=["z"]).animal_id.nunique()
    counts = "/".join(str(int(row[f"n_{condition.lower()}"])) for condition in CONDITIONS)
    span = f"{radii.min():g}–{radii.max():g} µm"
    population = ("todos los centroides asociados a DAPI" if endpoint == "cfos" else "centroides POMC positivos") if es else (
        "all DAPI-associated centroids" if endpoint == "cfos" else "POMC-positive centroids")
    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺POMC⁺"
    crossing = (f"{radii[exits].min():g}–{radii[exits].max():g} µm" if exits.any() else None)
    observed, critical = band.observed_f.to_numpy(float), band.critical_f.to_numpy(float)
    logarithmic = np.nanmax(critical) > LOG_AXIS_RATIO * max(np.nanmax(observed), 1e-9)
    axis_note = ("" if not logarithmic else
                 ("El eje vertical derecho es logarítmico: con dos jaulas en una condición la envolvente queda "
                  "hasta dos órdenes de magnitud por encima de la curva observada, que en escala lineal sería "
                  "indistinguible del eje. " if es else
                  "The right-hand vertical axis is logarithmic: with two cages in one condition the envelope sits "
                  "up to two orders of magnitude above the observed curve, which on a linear axis would be "
                  "indistinguishable from the axis line. "))
    if es:
        text = (f"Organización espacial de {label} en ARC, resuelta por distancia. Para cada animal se cuenta el número "
                f"de pares de células activadas separadas por menos de r dentro de {population}, por adquisición y lado "
                f"del tejido, y se expresa en desviaciones estándar respecto de su propia hipótesis nula de etiquetado "
                f"aleatorio; la media y la varianza nulas son exactas, de modo que la geometría del campo se cancela. "
                f"La unidad independiente es la jaula: las curvas de los animales se promedian con igual peso dentro de "
                f"cada jaula. Izquierda: líneas finas, jaulas; líneas gruesas, medias por condición; referencia 0. "
                f"Derecha: F entre condiciones en cada radio con la envolvente global del 95 % obtenida por enumeración "
                f"exhaustiva de {int(row.allocations)} asignaciones de etiquetas de jaula. La prueba es la desviación "
                f"máxima estudentizada sobre todo el intervalo {span}, por lo que un único valor p cubre todas las "
                f"distancias: p={formatted(row.p_value)}, p Holm={formatted(row.p_holm_two_endpoint_family)} con "
                f"{counts} jaulas ({animals} animales). Con estos tamaños ninguna prueba exacta puede dar "
                f"p<{formatted(row.minimum_attainable_p)}. ")
        text += (f"La curva observada supera la envolvente entre {crossing}. " if crossing else
                 "La curva observada no supera la envolvente en ningún radio. ")
        text += axis_note
        text += ("La eminencia media no se analiza: conserva 74 células asociadas al marcador con 28 activadas en 11 "
                 "animales, y solo 6 tienen dos o más activadas, insuficiente para una curva. El método es idéntico al "
                 "de los paneles H e I de la Figura 3; solo cambian la población y la unidad experimental. La "
                 "proximidad no establece conectividad funcional.")
    else:
        text = (f"Scale-resolved spatial organization of {label} in ARC. For each animal the number of activated-cell "
                f"pairs closer than r is counted within {population}, by acquisition and tissue side, and expressed in "
                f"standard deviations from its own random-labelling null; the null mean and variance are exact, so "
                f"field geometry cancels. The independent unit is the cage: animal curves are averaged with equal "
                f"weight within each cage. Left: thin lines are cages, thick lines are condition means, and the "
                f"reference is 0. Right: the between-condition F at each radius with its global 95% envelope, from an "
                f"exhaustive enumeration of {int(row.allocations)} cage-label allocations. The test is the studentized "
                f"maximum deviation over the whole {span} range, so one p-value covers every distance: "
                f"p={formatted(row.p_value)}, Holm p={formatted(row.p_holm_two_endpoint_family)} with {counts} cages "
                f"({animals} animals). At these sizes no exact test can return p<{formatted(row.minimum_attainable_p)}. ")
        text += (f"The observed curve exits the envelope between {crossing}. " if crossing else
                 "The observed curve does not exit the envelope at any radius. ")
        text += axis_note
        text += ("The median eminence is not analysed: it retains 74 marker-associated cells with 28 activated across "
                 "11 animals, and only 6 of those animals have two or more activated cells, too sparse for a curve. "
                 "The method is identical to Figure 3 panels H and I; only the population and the experimental unit "
                 "differ. Proximity does not establish functional connectivity.")
    return text


def generate(analysis_dir: Path, output_dir: Path, dpi: int = 600) -> dict:
    analysis_dir, output_dir = Path(analysis_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    curves = pd.read_csv(analysis_dir / "scale_resolved_curves.csv")
    envelope = pd.read_csv(analysis_dir / "scale_resolved_envelope.csv")
    test = pd.read_csv(analysis_dir / "scale_resolved_test.csv")
    if set(test.region) != {REGION}:
        raise ValueError(f"Figure 4's scale-resolved analysis must cover {REGION} only; found {sorted(set(test.region))}")
    style.set_style()
    plt.rcParams.update({"text.color": "black", "axes.labelcolor": "black", "axes.titlecolor": "black",
                         "xtick.color": "black", "ytick.color": "black", "xtick.labelcolor": "black",
                         "ytick.labelcolor": "black", "axes.edgecolor": "black"})
    files, legends = [], {"en": {}, "es": {}}
    for language in ("en", "es"):
        for endpoint, (letter, _) in PANELS.items():
            files.extend(plot_panel(curves, envelope, test, endpoint, language, output_dir, dpi))
            legends[language][letter] = legend_text(curves, envelope, test, endpoint, language)
    path = output_dir / "legends.json"
    path.write_text(json.dumps(legends, indent=2, ensure_ascii=False) + "\n")
    files.append(str(path))
    receipt = dict(schema="fig4-scale-resolved-panels-v1", analysis_dir=str(analysis_dir.resolve()),
                   output_dir=str(output_dir.resolve()), dpi=dpi, region=REGION,
                   unit="biological cage", me_analysed=False,
                   endpoints={endpoint: letter for endpoint, (letter, _) in PANELS.items()},
                   files=[str(Path(f).relative_to(output_dir)) for f in files])
    (output_dir / "panel_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    print(json.dumps(generate(args.analysis_dir, args.output_dir, args.dpi), indent=2))


if __name__ == "__main__":
    main()
