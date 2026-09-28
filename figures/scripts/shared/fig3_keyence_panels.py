#!/usr/bin/env python3
"""Render Figure 3 animal summaries from the documented Keyence extension.

All inferential results are read from supplied tables. This renderer computes
only descriptive means/SDs and writes the existing E/F/G/H/I panel contracts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

import shape_spatial_figures as spatial
import fig3_allulose_example as allulose_example


CONDITIONS = spatial.CONDITIONS
COLORS = spatial.COLORS
BAR_COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
PAIRS = (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose"))
QUANT_SIZE = (311.052 / 72, 252.847 / 72)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pformat(value):
    value = float(value)
    if not np.isfinite(value):
        return "NE"
    return f"{value:.2g}" if value < .001 else f"{value:.3f}"


def require(frame, columns, name):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{name} missing columns: {sorted(missing)}")


def one_row(frame, endpoint, name):
    rows = frame[frame.endpoint.eq(endpoint)]
    if len(rows) != 1:
        raise ValueError(f"Expected one {name} row for {endpoint}; found {len(rows)}")
    return rows.iloc[0]


def keyence(value):
    return "keyence" in str(value).lower()


def platform_handles(language, expected=False):
    es = language == "es"
    handles = [Line2D([], [], marker="o", linestyle="none", markerfacecolor="#9aa3ab",
                      markeredgecolor="black", markersize=4, label="Animal")]
    if expected:
        handles += [
            Line2D([], [], color="black", lw=1.7, label="Media" if es else "Mean"),
            Line2D([], [], marker="D", linestyle="none", markerfacecolor="white",
                   markeredgecolor="black", markersize=4,
                   label="Esperado al azar" if es else "Random-label expectation"),
        ]
    return handles


def included_rows(frame, endpoint=None):
    name = endpoint + "_include_in_main" if endpoint else "include_in_main"
    if name not in frame:
        return frame
    flag = frame[name].astype(str).str.lower()
    if not flag.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Main-analysis eligibility must be explicitly boolean")
    return frame.loc[flag.isin(["true", "1"])].copy()


def comparison_symbol(value):
    value = float(value)
    if not np.isfinite(value): return "NE"
    if value < .0001: return "****"
    if value < .001: return "***"
    if value < .01: return "**"
    if value < .05: return "*"
    return "NS"


BRACKETS = (("Water", "Sucrose", 0, 1, 1.025), ("Allulose", "Sucrose", 1, 2, 1.025),
            ("Allulose", "Water", 0, 2, 1.14))


def comparison_brackets(ax, contrasts, endpoint, metric):
    block = contrasts[contrasts.endpoint.eq(endpoint) & contrasts.metric.eq(metric)]
    transform = ax.get_xaxis_transform()
    for first, second, left, right, y in BRACKETS:
        row = block[(block.group_a.eq(first) & block.group_b.eq(second)) |
                    (block.group_a.eq(second) & block.group_b.eq(first))]
        if len(row) != 1:
            raise ValueError(f"Expected one Figure 3 comparison: {endpoint}/{metric}/{first}-{second}")
        label = f"p={pformat(row.iloc[0].p_raw)}"
        ax.plot([left, left, right, right], [y, y+.025, y+.025, y],
                transform=transform, clip_on=False, color="black", lw=.7)
        ax.text((left+right)/2, y+.034, label, transform=transform,
                ha="center", va="bottom", fontsize=6.5, color="black", clip_on=False)


def abundance_rows(stats, endpoint):
    block = stats[stats.endpoint.eq(endpoint)]
    omnibus = block[block.test.eq("one_way_ANOVA")]
    if len(omnibus) != 1:
        raise ValueError(f"Expected one ANOVA row for {endpoint}")
    pairs = {}
    for first, second in PAIRS:
        row = block[block.test.eq("exact_two_sided_Mann_Whitney") &
                    ((block.group_a.eq(first) & block.group_b.eq(second)) |
                     (block.group_a.eq(second) & block.group_b.eq(first)))]
        if len(row) != 1:
            raise ValueError(f"Expected one exact MWU row for {endpoint}: {first}/{second}")
        pairs[(first, second)] = float(row.iloc[0].p_value)
    return omnibus.iloc[0], pairs


def robust_rows(robust, endpoint):
    """The unequal-variance omnibus and pairwise rows for one abundance endpoint."""
    if robust is None:
        return None, None
    block = robust.loc[robust.domain.eq("abundance") & robust.endpoint.eq(endpoint)]
    if block.empty:
        return None, None
    omnibus = block.loc[block.test.eq("exact_permutation_welch_anova")]
    pairs = block.loc[block.test.eq("exact_studentized_permutation")]
    omnibus = omnibus.iloc[0] if len(omnibus) and np.isfinite(float(omnibus.iloc[0].p_value)) else None
    return omnibus, (pairs if len(pairs) else None)


def plot_abundance(values, stats, endpoint, language, destination, dpi, robust=None):
    es = language == "es"
    fig, ax = plt.subplots(figsize=QUANT_SIZE)
    anova, pairs = abundance_rows(stats, endpoint)
    # The unequal-variance block adds two lines above the axes, so it needs more headroom.
    robust_omnibus, robust_pairs = robust_rows(robust, endpoint)
    fig.subplots_adjust(left=.185, right=.985, bottom=.24, top=.77)
    heights = []
    for index, condition in enumerate(CONDITIONS):
        rows = included_rows(values[values.condition.eq(condition)], endpoint).sort_values("animal")
        observed = pd.to_numeric(rows[endpoint], errors="coerce").to_numpy(float) * 100
        valid = np.isfinite(observed)
        y = observed[valid]
        mean = float(y.mean()) if len(y) else np.nan
        sd = float(y.std(ddof=1)) if len(y) > 1 else 0
        if len(y):
            ax.bar(index, mean, yerr=sd, width=.52, color=BAR_COLORS[condition],
                   edgecolor="black", linewidth=1.1, capsize=3,
                   error_kw={"elinewidth": 1.1, "capthick": 1.1}, zorder=2)
        jitter = np.linspace(-.13, .13, len(rows)) if len(rows) > 1 else np.zeros(len(rows))
        for j, (_, row) in enumerate(rows.iterrows()):
            if valid[j]:
                ax.scatter(index + jitter[j], observed[j], marker="o",
                           facecolor=BAR_COLORS[condition], edgecolor="black", linewidth=.8,
                           s=25, zorder=4, clip_on=False)
        top = max([mean + sd] + list(y)) if len(y) else 0
        heights.append(top)
    ymax = max(max(heights), .01) * 1.27
    ax.set_ylim(0, ymax)
    for index, condition in enumerate(CONDITIONS):
        n = int(pd.to_numeric(included_rows(values[values.condition.eq(condition)], endpoint)[endpoint], errors="coerce").notna().sum())
        ax.text(index, heights[index] + .035 * ymax, f"n={n}", ha="center", va="bottom",
                fontsize=8, fontweight="bold")
    names = [spatial.NAMES_ES[c] if es else c for c in CONDITIONS]
    ax.set_xticks(range(3), names, rotation=24, ha="right", fontsize=9, fontweight="bold")
    ax.set_xlim(-.5, 2.5)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ylabel = ("Núcleos c-FOS⁺ (% DAPI)" if es else "c-FOS⁺ nuclei (% DAPI)") if endpoint == "cfos_over_dapi" else (
        "Núcleos c-FOS⁺NPY⁺ (% NPY)" if es else "c-FOS⁺NPY⁺ nuclei (% NPY)")
    ax.set_ylabel(ylabel, fontsize=9, fontweight="bold")
    ax.tick_params(width=1, length=3, labelsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
    title = "Activación c-FOS total" if es else "Total c-FOS activation"
    if endpoint == "double_over_npy":
        title = "Activación c-FOS en núcleos NPY" if es else "c-FOS activation of NPY nuclei"
    fig.text(.56, .966, title, ha="center", va="top", fontsize=10, fontweight="bold")
    df = int(sum(int(anova[f"n_{c.lower()}"]) for c in CONDITIONS) - 3)
    initials = {"Water": "A" if es else "W", "Sucrose": "S", "Allulose": "Al" if es else "A"}
    if robust_omnibus is not None:
        label = "p de permutación exacta" if es else "Exact permutation p"
        fig.text(.56, .871, f"{label}={pformat(robust_omnibus.p_value)}",
                 ha="center", fontsize=7.8, fontweight="bold")
        pairtext = "; ".join(f"{initials[a]}–{initials[b]} {pformat(pairs[(a,b)])}" for a,b in PAIRS)
        fig.text(.56, .811, ("p MWU exacta: " if es else "Exact MWU p: ") + pairtext, ha="center", fontsize=6.7)
    else:
        fig.text(.56, .871, f"ANOVA: F(2,{df})={float(anova.statistic):.2f}; p={pformat(anova.p_value)}",
                 ha="center", fontsize=7.6, fontweight="bold")
        pairtext = "; ".join(f"{initials[a]}–{initials[b]} {pformat(pairs[(a,b)])}" for a,b in PAIRS)
        fig.text(.56, .811, ("p MWU exacta (nominal): " if es else "Exact MWU p (nominal): ") + pairtext,
                 ha="center", fontsize=6.7)
    fig.legend(handles=platform_handles(language), loc="center", bbox_to_anchor=(.55, .045),
               ncol=2, frameon=False, fontsize=6.6, handletextpad=.4, columnspacing=1.2)
    old_letter = "D" if endpoint == "cfos_over_dapi" else "E"
    suffix = "_spanish" if es else ""
    return spatial.save_figure(fig, destination / f"Figure_3_cFos_NPY_Panel_{old_letter}{suffix}",
                               dpi, title)


def plot_metric(ax, block, metric, language):
    """Existing closed-contour display, with acquisition platform point shapes."""
    es = language == "es"
    cluster = metric == "cluster_fraction"
    displayed, counts = [], []
    for index, condition in enumerate(CONDITIONS):
        recorded = block[block.condition.eq(condition)]
        group = included_rows(recorded).sort_values("unit_id")
        y = pd.to_numeric(group[metric], errors="coerce").to_numpy(float) * (100 if cluster else 1)
        valid = np.isfinite(y)
        counts.append((int(valid.sum()), len(recorded)))
        jitter = np.linspace(-.145, .145, len(group)) if len(group) > 1 else np.zeros(len(group))
        x = index + jitter - (.055 if cluster else 0)
        if cluster:
            expected = pd.to_numeric(group.cluster_expected_fraction, errors="coerce").to_numpy(float) * 100
            xe = index + jitter + .055
            for j in range(len(group)):
                if np.isfinite(y[j]) and np.isfinite(expected[j]):
                    ax.plot([x[j], xe[j]], [y[j], expected[j]], lw=.55, alpha=.30,
                            color=COLORS[condition], zorder=1)
            good = np.isfinite(expected)
            ax.scatter(xe[good], expected[good], marker="D", s=14, facecolors="white",
                       edgecolors=COLORS[condition], linewidths=.75, alpha=.8, zorder=3)
            displayed.extend(expected[good])
        for j, (_, row) in enumerate(group.iterrows()):
            if valid[j]:
                ax.scatter(x[j], y[j], marker="o",
                           s=16, color=COLORS[condition],
                           edgecolors="white", linewidths=.25, alpha=.8, zorder=4)
        if valid.any():
            mean = float(y[valid].mean())
            ax.plot([index-.24, index+.24], [mean, mean], lw=1.7, color=COLORS[condition],
                    solid_capstyle="round", zorder=5)
            displayed.extend(y[valid])
        else:
            ax.text(index, .05, "NE", ha="center", fontsize=7, transform=ax.get_xaxis_transform())
    if cluster:
        ax.set_ylim(-3, 104)
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel("En grupos (%)" if es else "Clustered (%)")
    else:
        ax.axhline(1, color="black", linestyle=(0, (3, 2)), lw=.6, zorder=0)
        ax.set_ylim(0, max([1.] + list(displayed)) * 1.20)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        ax.set_ylabel("Enriquecimiento" if es else "Core enrichment")
    ax.set_xticks(range(3), [(spatial.NAMES_ES[c] if es else c) + f"\nn={n}/{total}"
                            for c, (n, total) in zip(CONDITIONS, counts)])
    ax.set_xlim(-.48, 2.48)
    ax.spines[["right", "top"]].set_visible(False)
    ax.grid(axis="y", color="#e8eaed", lw=.45, zorder=0)
    ax.set_axisbelow(True)


def per_measure_p(per_feature, endpoint, feature, language):
    """The single-measure exact permutation p shown above its own sub-plot."""
    if per_feature is None or "feature" not in getattr(per_feature, "columns", []):
        return ""
    row = per_feature.loc[per_feature.endpoint.eq(endpoint) & per_feature.feature.eq(feature)]
    if len(row) != 1 or not np.isfinite(float(row.iloc[0].p_value_exact)):
        return ""
    label = "p permutación exacta" if language == "es" else "exact permutation p"
    return f"\n{label}={pformat(row.iloc[0].p_value_exact)}"


def plot_scale_resolved(curves, envelope, test, endpoint, language, destination, dpi):
    """Panels H and I: clustering against chance at every distance, and its global test."""
    es = language == "es"
    fig, axes = plt.subplots(1, 2, figsize=(6.60, 3.12))
    # The condition keys sit inside the axes, so the space a bottom legend used
    # to occupy belongs to the plots.
    fig.subplots_adjust(left=.095, right=.985, bottom=.175, top=.785, wspace=.34)
    block = curves[curves.endpoint.eq(endpoint)].dropna(subset=["z"])
    band = envelope[envelope.endpoint.eq(endpoint)].sort_values("radius_um")
    row = test[test.endpoint.eq(endpoint)].iloc[0]

    left = axes[0]
    left.axhline(0, color="black", linestyle=(0, (3, 2)), lw=.6, zorder=0)
    for condition in CONDITIONS:
        group = block[block.condition.eq(condition)]
        for _, animal in group.groupby("animal_id"):
            animal = animal.sort_values("radius_um")
            left.plot(animal.radius_um, animal.z, color=COLORS[condition], lw=.7, alpha=.55, zorder=2)
        mean = group.groupby("radius_um").z.mean().sort_index()
        left.plot(mean.index, mean.to_numpy(), color=COLORS[condition], lw=2.1, zorder=3,
                  label=spatial.NAMES_ES[condition] if es else condition)
    left.set_xlabel("Distancia r (µm)" if es else "Distance r (µm)", fontsize=8, fontweight="bold")
    left.set_ylabel("Agrupación frente al azar (DE)" if es else "Clustering vs chance (SD)",
                    fontsize=8, fontweight="bold")
    left.set_title("Curvas por animal" if es else "Per-animal curves", pad=8, fontsize=8)
    left.legend(frameon=False, fontsize=6.4, loc="upper right", handlelength=1.4)

    right = axes[1]
    radii = band.radius_um.to_numpy(float)
    observed = band.observed_f.to_numpy(float)
    critical = band.critical_f.to_numpy(float)
    ceiling = max(np.nanmax(observed), np.nanmax(critical)) * 1.14
    right.fill_between(radii, critical, ceiling, color="#d8dce0", alpha=.55, lw=0, zorder=0,
                       label="Rechazo global" if es else "Global rejection")
    right.plot(radii, critical, color="black", linestyle=(0, (4, 2)), lw=.9, zorder=2,
               label="Envolvente 95 %" if es else "95% envelope")
    right.plot(radii, observed, color="black", lw=1.8, zorder=3,
               label="F observado" if es else "Observed F")
    exits = band.exits_envelope.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    if exits.any():
        right.scatter(radii[exits], observed[exits], s=13, color="black", zorder=4)
    right.set_ylim(0, ceiling)
    right.set_xlabel("Distancia r (µm)" if es else "Distance r (µm)", fontsize=8, fontweight="bold")
    right.set_ylabel("F entre condiciones" if es else "Between-condition F", fontsize=8, fontweight="bold")
    right.set_title("Prueba de envolvente global" if es else "Global envelope test", pad=8, fontsize=8)
    right.legend(frameon=False, fontsize=6.1, loc="upper left", handlelength=1.4)

    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺NPY⁺"
    fig.text(.54, .962, label + (": organización espacial" if es else " spatial organization"),
             ha="center", va="center", fontsize=10, fontweight="bold")
    counts = "/".join(str(int(row[f"n_{condition.lower()}"])) for condition in CONDITIONS)
    span = f"{radii.min():g}–{radii.max():g} µm"
    if es:
        headline = (f"Envolvente global, {span} (n={counts}): p={pformat(row.p_value)}; "
                    f"p Holm={pformat(row.p_holm_two_endpoint_family)}")
    else:
        headline = (f"Global envelope, {span} (n={counts}): p={pformat(row.p_value)}; "
                    f"Holm p={pformat(row.p_holm_two_endpoint_family)}")
    fig.text(.54, .885, headline, ha="center", va="center", fontsize=7.4, fontweight="bold")
    stem = "Figure3_Spatial_cFOS" + ("_NPY" if endpoint == "marker_cfos" else "") + "_occurrence"
    outputs = spatial.save_figure(fig, destination / "Fig3" / language / stem, dpi,
                                  label + " scale-resolved clustering and global envelope test")
    receipt = dict(statistics_recomputed=False, endpoint=endpoint,
                   radii_um=[float(r) for r in radii], p_value=float(row.p_value),
                   p_holm=float(row.p_holm_two_endpoint_family),
                   exits_envelope_at=[float(r) for r in radii[exits]])
    path = destination / "Fig3" / language / (stem + "_envelope.json")
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    outputs.append(str(path))
    return outputs


def plot_measure_scatter(units, omnibus, contrasts, endpoint, language, destination, dpi):
    """Fallback for bundles without scale-resolved curves: one test per measure."""
    es = language == "es"
    fig, axes = plt.subplots(1, 2, figsize=(6.60, 3.12))
    fig.subplots_adjust(left=.095, right=.985, bottom=.27, top=.67, wspace=.40)
    block = spatial.select_rows(units, "NPY", endpoint, "FIELD")
    for ax, metric, test_metric in zip(axes, ("core_enrichment", "cluster_fraction"), ("core_enrichment", "cluster_excess")):
        plot_metric(ax, block, metric, language)
        comparison_brackets(ax, contrasts, endpoint, test_metric)
    axes[0].set_title(("Contorno de densidad NPY (80 %)" if es else "NPY density contour (80%)")
                      + per_measure_p(omnibus, endpoint, "core_enrichment", language), pad=27, fontsize=8)
    axes[1].set_title(("Grupos espaciales (≥2 células)" if es else "Spatial clusters (≥2 cells)")
                      + per_measure_p(omnibus, endpoint, "cluster_excess", language), pad=27, fontsize=8)
    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺NPY⁺"
    fig.text(.54, .955, label + (": organización espacial" if es else " spatial organization"),
             ha="center", va="center", fontsize=10, fontweight="bold")
    fig.legend(handles=platform_handles(language, expected=True), loc="center", bbox_to_anchor=(.53, .09),
               ncol=3, frameon=False, fontsize=6.4, handlelength=1.3, columnspacing=1.1, handletextpad=.4)
    stem = "Figure3_Spatial_cFOS" + ("_NPY" if endpoint == "marker_cfos" else "") + "_occurrence"
    outputs = spatial.save_figure(fig, destination / "Fig3" / language / stem, dpi,
                                  label + " closed-density and positive-cluster organization")
    rows = contrasts[contrasts.endpoint.eq(endpoint)]
    annotations = [dict(endpoint=endpoint, metric=row.metric, group_a=row.group_a, group_b=row.group_b,
                        p_nominal=float(row.p_raw), p_holm=float(row.p_holm),
                        label=f"p={pformat(row.p_raw)}", displayed="nominal p only; Holm retained in source tables")
                   for row in rows.itertuples()]
    path = destination / "Fig3" / language / (stem + "_comparisons.json")
    path.write_text(json.dumps(dict(statistics_recomputed=False, annotations=annotations), indent=2)+"\n")
    outputs.append(str(path))
    return outputs


def robust_text(robust, endpoint, language):
    """Unequal-variance omnibus, variance check and studentized pairwise p values."""
    if robust is None:
        return "", ""
    block = robust.loc[robust.domain.eq("abundance") & robust.endpoint.eq(endpoint)]
    if block.empty:
        return "", ""
    def one(name):
        found = block.loc[block.test.eq(name)]
        return found.iloc[0] if len(found) else None
    welch, variance = one("welch_one_way_anova"), one("brown_forsythe_equal_variance")
    enumerated = one("exact_permutation_welch_anova")
    if enumerated is None or not np.isfinite(float(enumerated.p_value)):
        return "", ""
    pairs = block.loc[block.test.eq("exact_studentized_permutation")]
    detail = "; ".join(
        f"{spatial.NAMES_ES[row.group_a] if language == 'es' else row.group_a}–"
        f"{spatial.NAMES_ES[row.group_b] if language == 'es' else row.group_b}: p={pformat(row.p_value)}"
        for row in pairs.itertuples())
    variance_p = pformat(variance.p_value) if variance is not None and np.isfinite(float(variance.p_value)) else "NE"
    floor = float(pairs.minimum_attainable_p.max()) if len(pairs) else float("nan")
    welch_p = pformat(welch.p_value) if welch is not None and np.isfinite(float(welch.p_value)) else "NE"
    if language == "es":
        head = (f" El grupo agua es mucho más disperso que los grupos de azúcar y su dispersión coincide con el lote de "
                f"adquisición, por lo que el contraste global no supone varianzas iguales. F robusto de Welch con nulo "
                f"enumerado sobre las {int(enumerated.enumerated_labelings)} asignaciones distintas de animales: "
                f"F={float(enumerated.statistic):.4g}, p={pformat(enumerated.p_value)}. Prueba de Brown–Forsythe de igualdad "
                f"de varianzas: p={variance_p}. La aproximación F de Welch da p={welch_p}; sus grados de libertad se estiman "
                f"a partir de varianzas de tres animales y es anticonservadora con estos tamaños.")
        tail = (f" Permutación estudentizada exacta (t de Welch sobre todas las asignaciones): {detail}."
                + (f" Con estos tamaños de grupo ninguna prueba exacta puede dar p<{floor:.3g}." if np.isfinite(floor) else "")
                if detail else "")
    else:
        head = (f" The Water group is far more dispersed than either sugar group and its spread coincides with acquisition "
                f"batch, so the omnibus does not assume a common variance. Welch's robust F with its null enumerated over "
                f"all {int(enumerated.enumerated_labelings)} distinct animal allocations: F={float(enumerated.statistic):.4g}, "
                f"p={pformat(enumerated.p_value)}. Brown–Forsythe test of equal variance: p={variance_p}. The Welch F "
                f"approximation gives p={welch_p}; its degrees of freedom are estimated from three-animal variances and it "
                f"is anti-conservative at these group sizes.")
        tail = (f" Exact studentized permutation (Welch t over every animal-label allocation): {detail}."
                + (f" At these group sizes no exact test can return p<{floor:.3g}." if np.isfinite(floor) else "")
                if detail else "")
    return head, tail


def abundance_legend(values, stats, endpoint, language, robust=None):
    es = language == "es"
    omnibus, pairs = abundance_rows(stats, endpoint)
    n = [int(pd.to_numeric(included_rows(values[values.condition.eq(c)], endpoint)[endpoint], errors="coerce").notna().sum()) for c in CONDITIONS]
    total = int(len(values))
    flagged = int(values[endpoint + "_iqr_excluded"].astype(str).str.lower().isin(["true", "1"]).sum()) if endpoint + "_iqr_excluded" in values else 0
    iqrtext_en = ("The within-condition 1.5\u00d7IQR rule was applied separately by endpoint to this experiment; "
                  + (f"no outliers were flagged and all {total} animals were retained."
                     if not flagged else f"{flagged} of {total} animals were flagged and excluded from this endpoint."))
    iqrtext_es = ("Se aplic\u00f3 la regla 1,5\u00d7RIC dentro de cada condici\u00f3n y desenlace de este experimento; "
                  + (f"no se detectaron valores at\u00edpicos y se conservaron los {total} animales."
                     if not flagged else f"se marcaron y excluyeron {flagged} de {total} animales en este desenlace."))
    pairtext = "; ".join(f"{spatial.NAMES_ES[a] if es else a}–{spatial.NAMES_ES[b] if es else b}: p={pformat(pairs[(a,b)])}"
                         for a,b in PAIRS)
    robust_head, robust_tail = robust_text(robust, endpoint, language)
    anova_label_es = "ANOVA de una vía con varianzas iguales (sensibilidad)" if robust_head else "ANOVA de una vía"
    anova_label_en = "Equal-variance one-way ANOVA (sensitivity)" if robust_head else "One-way ANOVA"
    if es:
        measure = "Núcleos c-FOS positivos como porcentaje de DAPI" if endpoint == "cfos_over_dapi" else "Núcleos dobles positivos c-FOS/NPY como porcentaje de núcleos NPY positivos"
        return (f"{measure}. Los recuentos se suman entre campos/secciones de cada animal antes de calcular el porcentaje; cada animal recibe el mismo peso. Agua/sacarosa/alulosa: n={n[0]}/{n[1]}/{n[2]}. Los círculos representan animales de la cohorte original; barras, media ± DE. {iqrtext_es}{robust_head}{robust_tail} {anova_label_es}: F={float(omnibus.statistic):.4g}, p={pformat(omnibus.p_value)}. Mann–Whitney bilateral con enumeración exacta de etiquetas de animales y rangos medios para empates: {pairtext}; valores p nominales, sin ajuste, conforme al análisis de abundancia original. La combinación de cohortes y plataformas es exploratoria; la condición y la adquisición no están completamente equilibradas, por lo que estas comparaciones no aíslan un efecto de tratamiento.")
    measure = "c-FOS-positive nuclei as a percentage of DAPI nuclei" if endpoint == "cfos_over_dapi" else "c-FOS/NPY double-positive nuclei as a percentage of NPY-positive nuclei"
    return (f"{measure}. Counts are summed across each animal's fields/sections before computing its percentage; animals have equal weight. Water/Sucrose/Allulose: n={n[0]}/{n[1]}/{n[2]}. Circles denote animals from the original experiment; bars show mean ± SD. {iqrtext_en}{robust_head}{robust_tail} {anova_label_en}: F={float(omnibus.statistic):.4g}, p={pformat(omnibus.p_value)}. Two-sided Mann–Whitney tests use exact animal-label enumeration with average ranks for ties: {pairtext}; p values are nominal, without adjustment, preserving the original abundance analysis. The pooled cohort/platform comparison is exploratory; condition and acquisition are incompletely balanced, so these comparisons do not isolate a treatment effect.")


def measure_sentences(omnibus, adjusted, endpoint, language):
    """One clause per separately tested measure; there is no joint test to report."""
    es = language == "es"
    names = {"core_enrichment": ("enriquecimiento en el contorno de densidad" if es else "density-contour enrichment"),
             "cluster_excess": ("exceso de agrupación" if es else "clustering excess")}
    parts = []
    for _, row in omnibus[omnibus.endpoint.eq(endpoint)].iterrows():
        label = names.get(row.feature, row.feature)
        sugar = adjusted[adjusted.endpoint.eq(endpoint) & adjusted.feature.eq(row.feature)]
        if es:
            clause = (f"{label}: prueba exacta de tres condiciones con {int(row.enumerated_labelings)} asignaciones, "
                      f"p={pformat(row.p_value_exact)}, p Holm={pformat(row.p_holm_four_measure_family)}")
            if len(sugar):
                clause += (f"; alulosa–sacarosa dentro de la misma plataforma, "
                           f"p={pformat(sugar.iloc[0].p_value_exact)}, p Holm={pformat(sugar.iloc[0].p_holm_four_measure_family)}")
        else:
            clause = (f"{label}: exact three-condition test over {int(row.enumerated_labelings)} allocations, "
                      f"p={pformat(row.p_value_exact)}, Holm p={pformat(row.p_holm_four_measure_family)}")
            if len(sugar):
                clause += (f"; Allulose–Sucrose within the shared platform, "
                           f"p={pformat(sugar.iloc[0].p_value_exact)}, Holm p={pformat(sugar.iloc[0].p_holm_four_measure_family)}")
        parts.append(clause)
    return ". ".join(parts)


def scale_resolved_legend(curves, envelope, test, omnibus, adjusted, endpoint, language):
    es = language == "es"
    row = test[test.endpoint.eq(endpoint)].iloc[0]
    band = envelope[envelope.endpoint.eq(endpoint)].sort_values("radius_um")
    radii = band.radius_um.to_numpy(float)
    exits = band.exits_envelope.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    span = f"{radii.min():g}–{radii.max():g} µm"
    if exits.any():
        crossing = f"{radii[exits].min():g}–{radii[exits].max():g} µm"
    else:
        crossing = None
    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺NPY⁺"
    graph = ("todos los centroides asociados a DAPI" if endpoint == "cfos" else "centroides NPY positivos") if es else (
        "all DAPI-associated centroids" if endpoint == "cfos" else "NPY-positive centroids")
    counts = "/".join(str(int(row[f"n_{condition.lower()}"])) for condition in CONDITIONS)
    if es:
        text = (f"Organización espacial de {label} resuelta por distancia. Izquierda: para cada animal, el número "
                f"observado de pares positivo–positivo separados por menos de r, expresado en desviaciones estándar "
                f"respecto de su propia hipótesis nula de etiquetado aleatorio. El etiquetado aleatorio mantiene fijas "
                f"todas las posiciones celulares, el estrato y el número de positivos, y solo permuta qué células son "
                f"positivas; la media y la varianza nulas son exactas, de modo que la geometría del campo se cancela y "
                f"no interviene ninguna corrección de borde. Los pares se cuentan dentro de {graph}, por adquisición y "
                f"lado del tejido. Líneas finas, animales; líneas gruesas, medias por condición; referencia 0. "
                f"Derecha: F entre condiciones en cada radio, con la envolvente global del 95 % obtenida por "
                f"enumeración exhaustiva de {int(row.allocations)} asignaciones de etiquetas de animales. "
                f"La prueba es la desviación máxima estudentizada sobre todo el intervalo {span}, de modo que un único "
                f"valor p cubre todas las distancias: p={pformat(row.p_value)}, p Holm={pformat(row.p_holm_two_endpoint_family)} "
                f"con n={counts} animales estimables. Con estos tamaños de grupo ninguna prueba exacta puede dar "
                f"p<{pformat(row.minimum_attainable_p)}. ")
        text += (f"La curva observada supera la envolvente entre {crossing}. " if crossing else
                 "La curva observada no supera la envolvente en ningún radio. ")
        text += (f"Las medidas escalares previas se conservan y se prueban por separado, sin prueba conjunta: "
                 f"{measure_sentences(omnibus, adjusted, endpoint, language)}. "
                 "La proximidad no establece conectividad funcional, y la adquisición del agua sigue confundida con la condición.")
    else:
        text = (f"Scale-resolved spatial organization of {label}. Left: for each animal, the observed number of "
                f"positive–positive pairs closer than r, expressed in standard deviations from its own "
                f"random-labelling null. Random labelling holds every cell position, stratum and positive count "
                f"fixed and permutes only which cells are positive; the null mean and variance are exact, so field "
                f"geometry cancels and no edge correction applies. Pairs are counted within {graph}, by acquisition "
                f"and tissue side. Thin lines are animals, thick lines are condition means, and the reference is 0. "
                f"Right: the between-condition F at each radius with its global 95% envelope, obtained by exhaustive "
                f"enumeration of {int(row.allocations)} animal-label allocations. The test is the studentized maximum "
                f"deviation over the whole {span} range, so one p-value covers every distance: p={pformat(row.p_value)}, "
                f"Holm p={pformat(row.p_holm_two_endpoint_family)} with n={counts} estimable animals. At these group "
                f"sizes no exact test can return p<{pformat(row.minimum_attainable_p)}. ")
        text += (f"The observed curve exits the envelope between {crossing}. " if crossing else
                 "The observed curve does not exit the envelope at any radius. ")
        text += (f"The earlier scalar measures are retained and tested separately, with no joint test: "
                 f"{measure_sentences(omnibus, adjusted, endpoint, language)}. "
                 "Proximity does not establish functional connectivity, and Water acquisition remains confounded with condition.")
    return text


def spatial_legend(omnibus, adjusted, endpoint, language):
    """Legend for the fallback measure-scatter rendering."""
    es = language == "es"
    label = "c-FOS⁺" if endpoint == "cfos" else "c-FOS⁺NPY⁺"
    head = (f"Organización espacial de {label}. " if es else f"Spatial organization of {label}. ")
    body = ("El contorno cerrado encierra al menos el 80 % de la masa KDE de los centroides NPY positivos, "
            "estimada independientemente de c-FOS con corrección de borde. Izquierda: enriquecimiento "
            "observado/esperado dentro del contorno. Derecha: porcentaje observado de positivas en componentes "
            "conectados de ≥2 positivas, con la esperanza bajo etiquetado aleatorio. Cada medida se prueba por "
            "separado; no se realiza ninguna prueba conjunta. " if es else
            "The closed contour encloses at least 80% of the KDE mass of NPY-positive centroids, estimated "
            "independently of c-FOS and edge corrected. Left: observed/expected enrichment inside the contour. "
            "Right: observed percentage of positives in connected components of at least two, with the "
            "random-labelling expectation. Each measure is tested on its own; no joint test is performed. ")
    return head + body + measure_sentences(omnibus, adjusted, endpoint, language) + "."


def comparison_legend(contrasts, endpoint, language, excluded=0):
    block = contrasts[contrasts.endpoint.eq(endpoint)]
    details = "; ".join(f"{r.metric}, {r.group_a}–{r.group_b}: p={pformat(r.p_raw)}, Holm p={pformat(r.p_holm)}" for r in block.itertuples())
    excluded = int(excluded)
    iqr_en = ("The within-condition 1.5\u00d7IQR rule excluded no animals."
              if not excluded else f"The within-condition 1.5\u00d7IQR rule excluded {excluded} animal(s) from this endpoint.")
    iqr_es = ("La regla 1,5\u00d7RIC por condici\u00f3n no excluy\u00f3 ning\u00fan animal."
              if not excluded else f"La regla 1,5\u00d7RIC por condici\u00f3n excluy\u00f3 {excluded} animal(es) en este desenlace.")
    if language == "es":
        return (" Los corchetes muestran las tres comparaciones entre condiciones mediante pruebas exactas bilaterales de diferencia de medias entre animales. Cada corchete muestra únicamente su p nominal sin ajustar; los valores Holm figuran en esta leyenda y en las tablas fuente. La corrección Holm abarca seis comparaciones por métrica dentro de Figura 3 (dos desenlaces × tres contrastes). Los corchetes de agrupación prueban el exceso observado menos esperado, aunque el eje muestra la fracción observada. Para cada valor: NS, p≥0,05; *, p<0,05; **, p<0,01; ***, p<0,001; ****, p<0,0001; NE, no estimable. Las pruebas PERMANOVA conjuntas se informan aquí y en las tablas fuente, separadas de los contrastes por métrica. " + iqr_es + " " + details + ".")
    return (" Brackets show all three between-condition comparisons using exact two-sided animal-label mean-difference tests. Each bracket displays only its nominal, unadjusted p; the Holm-adjusted values appear in this legend and the source tables. Holm correction covers six comparisons per metric within Figure 3 (two endpoints × three contrasts). Clustering brackets test observed-minus-expected excess, while the axis displays the observed fraction. For either value: NS, p≥0.05; *, p<0.05; **, p<0.01; ***, p<0.001; ****, p<0.0001; NE, not estimable. Joint PERMANOVA results are reported in this legend and source tables separately from the metric-specific contrasts. " + iqr_en + " " + details + ".")


def generate(analysis_dir, output_dir, dpi=600, example_path=None):
    analysis_dir, output_dir = Path(analysis_dir), Path(output_dir)
    names = ["animal_plot_values.csv", "abundance_statistics.csv", "unit_metrics.csv",
             "spatial_permanova.csv", "spatial_sugar_platform_sensitivity.csv", "treatment_contrasts.csv"]
    scale_names = ["scale_resolved_curves.csv", "scale_resolved_envelope.csv", "scale_resolved_test.csv"]
    scale_resolved = all((analysis_dir / name).is_file() for name in scale_names)
    for optional in ["abundance_robust_statistics.csv"] + (scale_names if scale_resolved else []):
        if (analysis_dir / optional).is_file():
            names = names + [optional]
    tables = {name: pd.read_csv(analysis_dir / name) for name in names}
    robust = tables.get("abundance_robust_statistics.csv")
    values, stats, units = (tables[n] for n in names[:3])
    omnibus, adjusted = tables["spatial_permanova.csv"], tables["spatial_sugar_platform_sensitivity.csv"]
    contrasts = tables["treatment_contrasts.csv"]
    if values.animal.str.contains("Keyence", case=False).any():
        raise ValueError("Figure 3 must use only its own experiment's animals")
    recorded = values.groupby("condition").size().to_dict()
    if set(recorded) != set(CONDITIONS) or min(recorded.values()) < 3:
        raise ValueError(f"Figure 3 needs at least three animals in each condition; recorded {recorded}")
    require(values, ["animal", "condition", "cfos_over_dapi", "double_over_npy", "platform"], names[0])
    require(stats, ["endpoint", "test", "group_a", "group_b", "statistic", "p_value",
                    "n_water", "n_sucrose", "n_allulose"], names[1])
    require(units, ["unit_id", "cohort", "endpoint", "region", "condition", "core_enrichment",
                    "cluster_fraction", "cluster_expected_fraction", "cluster_excess"], names[2])
    if values.animal.duplicated().any():
        raise ValueError("Abundance values must contain one row per animal")
    if set(values.condition) != set(CONDITIONS):
        raise ValueError("Expected Water, Sucrose and Allulose conditions")
    for metric in ("cfos_over_dapi", "double_over_npy"):
        data = pd.to_numeric(values[metric], errors="coerce").dropna()
        if ((data < 0) | (data > 1)).any():
            raise ValueError(f"{metric} must contain fractions, not percentages")
    platform_map = values.set_index("animal").platform.to_dict()
    units = units[units.cohort.eq("NPY") & units.region.eq("FIELD")].copy()
    units["platform"] = units.unit_id.map(platform_map)
    if units.platform.isna().any():
        raise ValueError("Spatial animals lack platform mapping in animal_plot_values.csv")
    if units.duplicated(["endpoint", "unit_id"]).any():
        raise ValueError("Spatial unit table must have one NPY FIELD row per endpoint/animal")
    for metric in ("cfos_over_dapi", "double_over_npy"):
        row, _ = abundance_rows(stats, metric)
        for condition in CONDITIONS:
            count = int(pd.to_numeric(included_rows(values[values.condition.eq(condition)], metric)[metric], errors="coerce").notna().sum())
            if int(row[f"n_{condition.lower()}"]) != count:
                raise ValueError(f"Abundance statistics and plotted n differ: {metric}/{condition}")
    for endpoint in ("cfos", "marker_cfos"):
        block = units[units.endpoint.eq(endpoint)]
        complete = np.isfinite(block[["core_enrichment", "cluster_excess"]].to_numpy(float)).all(axis=1)
        if "include_in_main" in block:
            complete &= block.include_in_main.astype(str).str.lower().isin(["true", "1"]).to_numpy()
        for _, row in omnibus[omnibus.endpoint.eq(endpoint)].iterrows():
            for condition in CONDITIONS:
                group = block.condition.eq(condition).to_numpy()
                if int(row[f"n_{condition.lower()}"]) != int((complete & group).sum()):
                    raise ValueError(f"Per-measure n differs from source metrics: {endpoint}/{row.feature}/{condition}")
                if int(row[f"n_recorded_{condition.lower()}"]) != int(group.sum()):
                    raise ValueError(f"Recorded n differs from source metrics: {endpoint}/{row.feature}/{condition}")
    if example_path is not None:
        raise ValueError("Panel G uses the fixed reviewed Allulose field; --example-path is no longer supported")
    example_path, example_receipt = allulose_example.prepare(analysis_dir, output_dir)
    spatial.set_style()
    # Figure 3 writes in black; the shared style's slate ink stays with the other figures.
    plt.rcParams.update({"text.color": "black", "axes.labelcolor": "black",
                         "axes.titlecolor": "black", "xtick.color": "black", "ytick.color": "black",
                         "xtick.labelcolor": "black", "ytick.labelcolor": "black",
                         "axes.edgecolor": "black"})
    quantitative_dir = output_dir / "quantitative_inputs"
    spatial_dir = output_dir / "spatial_panels"
    files = [str(output_dir / "examples" / path.name) for path in (example_path, example_receipt)]
    quantitative_legends, shape_legends = {}, {}
    for language in ("en", "es"):
        quantitative_legends[language], shape_legends[language] = {}, {"G": allulose_example.legend(language, example_path)}
        for letter, endpoint in (("E", "cfos_over_dapi"), ("F", "double_over_npy")):
            files.extend(plot_abundance(values, stats, endpoint, language, quantitative_dir, dpi, robust))
            quantitative_legends[language][letter] = abundance_legend(values, stats, endpoint, language, robust)
            quantitative_legends[language][letter] += (
                " Se mantienen las regiones celulares aceptadas y las reglas de asociación de Figura 3, sin nueva segmentación ni redibujo de regiones."
                if language == "es" else
                " Accepted cell ROIs and the Figure 3 association rules are retained; no new segmentation or ROI redrawing is performed.")
        for letter, endpoint in (("H", "cfos"), ("I", "marker_cfos")):
            if scale_resolved:
                files.extend(plot_scale_resolved(tables["scale_resolved_curves.csv"],
                                                 tables["scale_resolved_envelope.csv"],
                                                 tables["scale_resolved_test.csv"],
                                                 endpoint, language, spatial_dir, dpi))
                shape_legends[language][letter] = scale_resolved_legend(
                    tables["scale_resolved_curves.csv"], tables["scale_resolved_envelope.csv"],
                    tables["scale_resolved_test.csv"], omnibus, adjusted, endpoint, language)
            else:
                files.extend(plot_measure_scatter(units, omnibus, contrasts, endpoint, language, spatial_dir, dpi))
                shape_legends[language][letter] = spatial_legend(omnibus, adjusted, endpoint, language)
                excluded_here = int(units.loc[units.endpoint.eq(endpoint), "iqr_excluded"].astype(str).str.lower().isin(["true", "1"]).sum()) if "iqr_excluded" in units else 0
                shape_legends[language][letter] += comparison_legend(contrasts, endpoint, language, excluded_here)
        files.extend(allulose_example.plot_example(example_path, spatial_dir, language, dpi))
    for path, data in ((output_dir / "quantitative_legends.json", quantitative_legends),
                       (spatial_dir / "Fig3" / "legends.json", shape_legends)):
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
        files.append(str(path))
    receipt = {"analysis_dir": str(analysis_dir.resolve()), "output_dir": str(output_dir.resolve()),
               "dpi": dpi, "statistics_recomputed": False, "statistical_unit": "animal",
               "shape_metrics_unchanged": ["core_enrichment", "cluster_fraction", "cluster_expected_fraction", "cluster_excess"],
               "point_shapes": {"original_animals": "circle"},
               "input_sha256": {name: digest(analysis_dir / name) for name in names},
               "example_sha256": digest(example_path), "source_sha256": digest(__file__),
               "files": {str(Path(path).relative_to(output_dir)): digest(path) for path in files}}
    (output_dir / "panel_render_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--example-path", type=Path)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    if args.dpi < 300:
        parser.error("Publication exports require at least 300 dpi")
    print(json.dumps(generate(args.analysis_dir, args.output_dir or args.analysis_dir,
                              args.dpi, args.example_path), indent=2))
