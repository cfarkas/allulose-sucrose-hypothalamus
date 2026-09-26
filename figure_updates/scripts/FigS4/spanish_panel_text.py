#!/usr/bin/env python3
"""Translate live Matplotlib panel labels without changing plotted geometry or data.

Combined publication figures remain English-only. Renderers call this helper only
after saving the English master and English panel crops, then export additional
``_spanish`` crops from the same live axes.
"""

from __future__ import annotations

from collections.abc import Mapping

from matplotlib.figure import Figure
from matplotlib.text import Text


GENERIC_REPLACEMENTS: tuple[tuple[str, str], ...] = (
    ("Primary Day-6 consumption endpoint", "Consumo · día 6"),
    ("Primary Day-6 body-weight endpoint", "Peso corporal · día 6"),
    ("Day-6 consumption endpoint", "Consumo · día 6"),
    ("Day-6 body-weight endpoint", "Peso corporal · día 6"),
    ("Exploratory within-sex cage analyses · no treatment-by-sex interaction test", "Exploratorio por sexo · sin prueba sexo×tratamiento"),
    ("Cumulative consumption", "Consumo acumulado"),
    ("Cage-balanced body-weight change", "Cambio de peso por jaula"),
    ("Consumption from Day 1 (mL/cage)", "Consumo desde el día 1 (mL/jaula)"),
    ("Consumption (mL/cage)", "Consumo (mL/jaula)"),
    ("Cage-mean change from baseline (%)", "Cambio medio por jaula desde basal (%)"),
    ("Change from Day 1 (%)", "Cambio desde el día 1 (%)"),
    ("Change from baseline (g)", "Cambio desde basal (g)"),
    ("Body-weight change", "Cambio de peso corporal"),
    ("body-weight change", "cambio de peso corporal"),
    ("Body weight", "Peso corporal"),
    ("Solution fraction", "Fracción de solución"),
    ("Total fluid", "Líquido total"),
    ("Solution", "Solución"),
    ("Food", "Alimento"),
    ("Study day", "Día del estudio"),
    ("thin: cages", "fino: jaulas"),
    ("thick/band: equal-cage mean/95% CI", "grueso/banda: media por jaula/IC del 95 %"),
    ("pale: mice", "pálido: ratones"),
    ("thin: cage-day", "fino: jaula-día"),
    ("One-way ANOVA", "ANOVA de una vía"),
    ("Holm exact", "Holm exacta"),
    ("95% CI not estimable", "IC del 95 % no estimable"),
    ("95% CI\nnot estimable", "IC del 95 %\nno estimable"),
    ("(one cage)", "(una jaula)"),
    ("Water", "Agua"),
    ("Sucrose", "Sacarosa"),
    ("Allulose", "Alulosa"),
    ("Female", "Hembra"),
    ("Male", "Macho"),
    (" mice", " ratones"),
    (" mouse", " ratón"),
    ("W-S", "Ag-Sac"),
    ("W-A", "Ag-Alu"),
    ("S-A", "Sac-Alu"),
    (" cages", " jaulas"),
    (" cage", " jaula"),
    ("Day 6", "Día 6"),
    ("Day 3", "Día 3"),
    ("Day 1", "Día 1"),
)


def translate_text(value: str, extra: Mapping[str, str] | None = None) -> str:
    """Return a deterministic Spanish label while preserving numbers and symbols."""
    replacements = list(GENERIC_REPLACEMENTS)
    if extra:
        replacements = list(extra.items()) + replacements
    result = value
    for source, target in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        result = result.replace(source, target)
    return result


def translate_figure_texts_to_spanish(
    figure: Figure, extra: Mapping[str, str] | None = None
) -> dict[str, str]:
    """Translate every nonempty live text artist and return an English→Spanish receipt."""
    figure.canvas.draw()
    receipt: dict[str, str] = {}
    for artist in figure.findobj(match=Text):
        source = artist.get_text()
        if not source:
            continue
        translated = translate_text(source, extra=extra)
        if translated != source:
            artist.set_text(translated)
            receipt[source] = translated
    # Axis formatters regenerate tick-label Text objects during draw. Pin only
    # the axes whose categorical labels changed, preserving the exact positions.
    for axis in figure.axes:
        x_source = [label.get_text() for label in axis.get_xticklabels()]
        x_translated = [translate_text(value, extra=extra) for value in x_source]
        if x_translated != x_source:
            axis.set_xticks(axis.get_xticks(), labels=x_translated)
            receipt.update({a: b for a, b in zip(x_source, x_translated) if a != b})
        y_source = [label.get_text() for label in axis.get_yticklabels()]
        y_translated = [translate_text(value, extra=extra) for value in y_source]
        if y_translated != y_source:
            axis.set_yticks(axis.get_yticks(), labels=y_translated)
            receipt.update({a: b for a, b in zip(y_source, y_translated) if a != b})
    figure.canvas.draw()
    return dict(sorted(receipt.items()))
