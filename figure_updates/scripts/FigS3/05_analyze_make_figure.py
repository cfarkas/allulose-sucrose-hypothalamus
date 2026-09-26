#!/usr/bin/env python3
"""Analyze animal-level Figure S3 outputs and render the publication figure."""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shapely
import sklearn
import statsmodels.api as sm
import pyvips
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from scipy import stats
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

from figs3_common import (
    FIGS3_DIR,
    ORGANS,
    TREATMENTS,
    FigS3Error,
    atomic_json,
    read_csv,
    require_regular_file,
    sha256_file,
)


def shared_helper_directory() -> Path:
    """Locate scripts/shared without assuming a fixed depth or a mount point.

    The paper root is recognised by its marker files, so the tree still works
    when a reader unpacks the archive under another name, and the search also
    succeeds from a scripts/<figure> mirror, where climbing one directory would
    land in Paper/scripts instead of Paper. A sibling copy of the helper, the
    arrangement the supplementary figures already use, is accepted last.
    """
    for parent in FIGS3_DIR.parents:
        marker = (parent / "scripts" / "setup").is_dir() and (
            parent / "README.txt"
        ).is_file()
        if parent.name == "Paper" or marker:
            candidate = parent / "scripts" / "shared"
            if (candidate / "spanish_panel_text.py").is_file():
                return candidate
    if (FIGS3_DIR / "spanish_panel_text.py").is_file():
        return FIGS3_DIR
    raise FigS3Error(
        "Cannot find scripts/shared/spanish_panel_text.py for the Spanish panels"
    )


# Isolated subpanels ship in English and Spanish; the multipanel master stays
# English-only. The translation runs on the live text artists of an already
# saved English panel, so geometry and data are identical in both languages.
sys.path.insert(0, str(shared_helper_directory()))
from spanish_panel_text import translate_figure_texts_to_spanish, translate_text


ANALYSIS_VERSION = "figs3-animal-composition-1.8"
RAW_CLASS_ORDER = [
    "Cancer cell",
    "Lymphocytes",
    "Fibroblasts",
    "Plasmocytes",
    "Eosinophils",
    "Neutrophils",
    "Macrophages",
    "Muscle Cell",
    "Endothelial Cell",
    "Red blood cell",
    "Epithelial",
    "Apoptotic Body",
    "Mitotic Figures",
    "Minor Stromal Cell",
]
CLASS_ORDER = [name for name in RAW_CLASS_ORDER if name != "Cancer cell"]
CLASS_SHORT = {
    "Lymphocytes": "Lymphocyte-like",
    "Fibroblasts": "Fibroblast-like",
    "Plasmocytes": "Plasmocyte-like",
    "Eosinophils": "Eosinophil-like",
    "Neutrophils": "Neutrophil-like",
    "Macrophages": "Macrophage-like",
    "Muscle Cell": "Muscle-like",
    "Endothelial Cell": "Endothelial-like",
    "Red blood cell": "RBC-like",
    "Epithelial": "Epithelial-like",
    "Apoptotic Body": "Apoptotic-like",
    "Mitotic Figures": "Mitotic-like",
    "Minor Stromal Cell": "Minor stromal-like",
}
TREATMENT_COLORS = {"Water": "#2166AC", "Sucrose": "#D73027", "Allulose": "#008F86"}
ORGAN_COLORS = {"Kidney": "#0097ff", "Liver": "#22b14c", "Spleen": "#b446c8"}
ORGAN_LABEL_VALUE = {"Kidney": 1, "Liver": 2, "Spleen": 3}
CELL_COLORS = {
    name: color_value
    for name, color_value in zip(
        CLASS_ORDER,
        [
            "#2166ac",
            "#f46d43",
            "#762a83",
            "#e6ab02",
            "#67a9cf",
            "#1b7837",
            "#8c510a",
            "#008571",
            "#b2182b",
            "#4d9221",
            "#7b3294",
            "#e08214",
            "#555555",
        ],
    )
}
CONTRASTS = ("Sucrose vs Water", "Allulose vs Water")
# Deterministic Spanish for the isolated subpanels. scripts/shared applies the
# longest key first, so composite labels are translated before their parts.
SPANISH_TERMS = {
    "Multi-organ WSI workflow and within-slide organ segmentation": "Flujo WSI multiorgánico y segmentación de órganos en el portaobjetos",
    "Deviation from Water composition": "Desviación respecto a la composición de Agua",
    "Three-organ histology and representative phenotype examples": "Histología de tres órganos y ejemplos fenotípicos representativos",
    "Compositional ordination (CLR-PCA)": "Ordenación composicional (ACP sobre CLR)",
    "Classical H&E effects vs Water (95% CI; permutation-p stars)": "Efectos clásicos de H&E frente a Agua (IC 95 %; asteriscos de p por permutación)",
    "Model-predicted cell-type composition effects vs Water": "Efectos de composición celular predichos por el modelo frente a Agua",
    "Workflow": "Flujo de trabajo",
    "24 animals": "24 animales",
    "One H&E WSI": "Una WSI de H&E",
    "Liver + kidney + spleen": "Hígado + riñón + bazo",
    "L0/L2 + CV masks": "L0/L2 + máscaras de visión por computador",
    "within-slide organ labels": "rótulos de órganos dentro del portaobjetos",
    "24 L0 tiles/organ": "24 teselas L0 por órgano",
    "spatially balanced": "equilibradas espacialmente",
    "HistoPLUS + classical CV": "HistoPLUS + visión por computador clásica",
    "animal-level composition": "composición a nivel de animal",
    "Organ segmentation": "Segmentación de órganos",
    "Kidney": "Riñón",
    "Liver": "Hígado",
    "Spleen": "Bazo",
    "cross-validated CLR\nMahalanobis distance": "distancia de Mahalanobis CLR\ncon validación cruzada",
    " vs Water": " frente a Agua",
    "instances": "instancias",
    "BH-significant:": "Significativo tras BH:",
    "of displayed outputs": "de las salidas mostradas",
    "standardized effect vs Water": "efecto estandarizado frente a Agua",
    "CLR effect vs Water": "efecto CLR frente a Agua",
    "No phenotype-composition contrast survived global BH q<0.05": "Ningún contraste de composición fenotípica superó la corrección global BH q<0,05",
}
# "Eosin" is a substring of "Eosinophil", so the two live in different panels'
# glossaries instead of one shared map that could corrupt the longer word.
SPANISH_TERMS_BY_PANEL = {
    "E": {
        "Nuclear density": "Densidad nuclear",
        "Nuclear area": "Área nuclear",
        "Eccentricity": "Excentricidad",
        "Solidity": "Solidez",
        "Hematoxylin": "Hematoxilina",
        "Eosin": "Eosina",
    },
    "F": {
        "Lymphocyte": "Linfocito",
        "Fibroblast": "Fibroblasto",
        "Plasmocyte": "Plasmocito",
        "Eosinophil": "Eosinófilo",
        "Neutrophil": "Neutrófilo",
        "Macrophage": "Macrófago",
        "Muscle": "Muscular",
        "Endothelial": "Endotelial",
        "RBC": "Eritrocito",
        "Epithelial": "Epitelial",
        "Apoptotic": "Apoptótico",
        "Mitotic": "Mitótico",
        "Minor stromal": "Estromal menor",
    },
    "C": {
        "lymphocyte-like": "tipo linfocito",
        "fibroblast-like": "tipo fibroblasto",
        "plasmocyte-like": "tipo plasmocito",
        "eosinophil-like": "tipo eosinófilo",
        "neutrophil-like": "tipo neutrófilo",
        "macrophage-like": "tipo macrófago",
        "muscle-like": "tipo muscular",
        "endothelial-like": "tipo endotelial",
        "rbc-like": "tipo eritrocito",
        "epithelial-like": "tipo epitelial",
        "apoptotic-like": "tipo apoptótico",
        "mitotic-like": "tipo mitótico",
        "minor stromal-like": "tipo estromal menor",
    },
}
MIN_DETECTED_ANIMALS = 18
MIN_TOTAL_CLASS_COUNT = 24
# Significance asterisks are drawn at twice their previous size, so the axis
# padding below is widened to keep a doubled star inside its own axes.
STAR_FONTSIZE_BARS = 18.0
STAR_FONTSIZE_FOREST = 15.0
PLOT_JITTER_SEED = 1707


def deterministic_stripplot(*args, seed: int = PLOT_JITTER_SEED, **kwargs):
    """Draw a Seaborn stripplot without leaking or consuming global RNG state.

    Seaborn 0.13 implements categorical jitter with ``np.random.uniform``.
    Without a local seed, Panel B moves thousands of pixels between otherwise
    identical reproductions. Saving and restoring the legacy NumPy RNG state
    makes the raster stable while leaving every caller's random stream intact.
    """
    random_state = np.random.get_state()
    try:
        np.random.seed(seed)
        return sns.stripplot(*args, **kwargs)
    finally:
        np.random.set_state(random_state)


# Master geometry. The canvas was tightened from 16.0 x 13.6 inches so the
# panels sit closer together; the inter-panel gaps shrink far more than the
# panels themselves, and every long label was shortened to keep the tighter
# layout free of collisions.
FIGURE_WIDTH = 15.0
FIGURE_HEIGHT = 12.3
FIGURE_LEFT = 0.033
FIGURE_RIGHT = 0.992
FIGURE_BOTTOM = 0.038
FIGURE_TOP = 0.982
FIGURE_HSPACE = 0.16
FIGURE_WSPACE = 0.20
ROW_RATIOS = (0.74, 1.10, 1.04, 1.10)
_COLUMN_UNITS = 12 + 11 * FIGURE_WSPACE
_COLUMN = (FIGURE_RIGHT - FIGURE_LEFT) * FIGURE_WIDTH / _COLUMN_UNITS
_ROW_UNITS = sum(ROW_RATIOS) + 3 * FIGURE_HSPACE * sum(ROW_RATIOS) / len(ROW_RATIOS)
_ROW = (FIGURE_TOP - FIGURE_BOTTOM) * FIGURE_HEIGHT / _ROW_UNITS
PANEL_WIDTH_LEFT = 5 * _COLUMN + 4 * FIGURE_WSPACE * _COLUMN
PANEL_WIDTH_RIGHT = 7 * _COLUMN + 6 * FIGURE_WSPACE * _COLUMN
PANEL_WIDTH_FULL = 12 * _COLUMN + 11 * FIGURE_WSPACE * _COLUMN
PANEL_HEIGHT_TOP = ROW_RATIOS[0] * _ROW
PANEL_HEIGHT_WIDE = ROW_RATIOS[1] * _ROW
PANEL_HEIGHT_MIDDLE = ROW_RATIOS[2] * _ROW


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figs3-root", type=Path, default=FIGS3_DIR)
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Fresh output root for analysis tables, bilingual panels, legends and masters.",
    )
    parser.add_argument("--permutations", type=int, default=99999)
    parser.add_argument("--seed", type=int, default=1707)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_complete_inputs(
    root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inference_summary = json.loads(
        require_regular_file(
            root / "histoplus" / "inference_summary.json", "inference summary"
        ).read_text()
    )
    if (
        inference_summary.get("status") != "complete"
        or inference_summary.get("complete_slides") != 24
    ):
        raise FigS3Error(
            "Figure analysis requires 24/24 completed HistoPLUS slide receipts"
        )
    segmentation_summary = json.loads(
        require_regular_file(
            root / "organ_segmentation" / "segmentation_summary.json",
            "organ-segmentation summary",
        ).read_text()
    )
    if (
        segmentation_summary.get("status") != "complete"
        or segmentation_summary.get("slides") != 24
        or segmentation_summary.get("image_registration_performed") is not False
    ):
        raise FigS3Error(
            "Figure analysis requires 24/24 segmentation-only slide receipts"
        )
    if (root / "registration").exists():
        raise FigS3Error(
            "Obsolete Figure S3 image-registration artifacts must be absent"
        )
    slide = pd.read_csv(root / "source_data" / "Figure_S3_slide_manifest.csv")
    phenotype = pd.read_csv(
        root / "histoplus" / "Figure_S3_HistoPLUS_phenotype_summary.csv"
    )
    tile_cv = pd.read_csv(
        root / "histoplus" / "Figure_S3_classical_tile_cv_features.csv"
    )
    if (
        slide.sample_id.nunique() != 24
        or phenotype.sample_id.nunique() != 24
        or tile_cv.sample_id.nunique() != 24
    ):
        raise FigS3Error("An aggregate input does not contain all 24 animals")
    return slide, phenotype, tile_cv


def composition_eligibility(
    organ_counts: pd.DataFrame,
) -> tuple[pd.Series, pd.Series, pd.Series]:
    detected = (organ_counts[CLASS_ORDER] > 0).sum(axis=0)
    total = organ_counts[CLASS_ORDER].sum(axis=0)
    eligible = (detected >= MIN_DETECTED_ANIMALS) & (total >= MIN_TOTAL_CLASS_COUNT)
    if int(eligible.sum()) < 2:
        raise FigS3Error(
            "Fewer than two phenotype components pass the prevalence filter"
        )
    return eligible, detected, total


def composition_values(phenotype: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    index = ["sample_id", "treatment", "cohort", "organ"]
    counts = phenotype.pivot_table(
        index=index,
        columns="model_class_raw",
        values="cell_count_in_sampled_tiles",
        aggfunc="sum",
        fill_value=0,
    ).reindex(columns=RAW_CLASS_ORDER, fill_value=0)
    if len(counts) != 72:
        raise FigS3Error(
            f"Expected 72 animal-organ composition rows, found {len(counts)}"
        )
    counts_frame = counts.reset_index()
    counts_frame["total_cells"] = counts[RAW_CLASS_ORDER].sum(axis=1).to_numpy()
    clr_parts = []
    for organ in ORGANS:
        organ_counts = counts_frame[counts_frame.organ == organ].reset_index(drop=True)
        eligible, detected, total = composition_eligibility(organ_counts)
        used = eligible.index[eligible].tolist()
        adjusted = organ_counts[used].astype(float) + 0.5
        log_values = np.log(adjusted.div(adjusted.sum(axis=1), axis=0))
        transformed = log_values.sub(log_values.mean(axis=1), axis=0)
        for class_name in CLASS_ORDER:
            part = organ_counts[index].copy()
            part["model_class_raw"] = class_name
            part["clr_value"] = (
                transformed[class_name].to_numpy()
                if bool(eligible[class_name])
                else np.nan
            )
            part["analysis_eligible"] = bool(eligible[class_name])
            part["detected_animals"] = int(detected[class_name])
            part["total_class_count"] = int(total[class_name])
            part["eligibility_rule"] = (
                f"detected_animals>={MIN_DETECTED_ANIMALS} and "
                f"total_class_count>={MIN_TOTAL_CLASS_COUNT} within organ"
            )
            part["pseudocount_applied"] = bool(eligible[class_name])
            clr_parts.append(part)
    clr_frame = pd.concat(clr_parts, ignore_index=True)
    clr_frame["phenotype_display"] = clr_frame["model_class_raw"].map(CLASS_SHORT)
    clr_frame["pseudocount"] = 0.5
    return counts_frame, clr_frame


def regression_one(
    frame: pd.DataFrame,
    value_column: str,
    permutations: int,
    seed: int,
) -> list[dict[str, float | str | int]]:
    data = frame[[value_column, "treatment", "cohort"]].dropna().copy()
    if set(data["treatment"]) != set(TREATMENTS):
        raise FigS3Error(f"Missing treatment level while fitting {value_column}")
    y = data[value_column].to_numpy(float)
    x = np.column_stack(
        [
            np.ones(len(data)),
            (data["treatment"] == "Sucrose").astype(float),
            (data["treatment"] == "Allulose").astype(float),
            (data["cohort"] == "POMC").astype(float),
        ]
    )
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise FigS3Error("Treatment/cohort design matrix is rank deficient")
    model = sm.OLS(y, x).fit(cov_type="HC3")
    if permutations < 1:
        raise FigS3Error("At least one residual permutation is required")
    generator = np.random.default_rng(seed)
    permutations_index = np.vstack(
        [generator.permutation(len(y)) for _ in range(permutations)]
    )
    inverse = np.linalg.inv(x.T @ x)
    projection = x @ inverse
    leverage = np.sum(projection * x, axis=1)
    if np.any(leverage >= 1.0 - 1e-12):
        raise FigS3Error("HC3 is undefined for a unit-leverage observation")
    output = []
    for contrast, coefficient_index in zip(CONTRASTS, (1, 2)):
        # Freedman--Lane tests this coefficient conditional on all other
        # columns, including the other treatment contrast as a nuisance term.
        nuisance_columns = [j for j in range(x.shape[1]) if j != coefficient_index]
        x0 = x[:, nuisance_columns]
        fitted = x0 @ np.linalg.lstsq(x0, y, rcond=None)[0]
        residual = y - fitted
        y_permuted = fitted[None, :] + residual[permutations_index]
        beta_permuted = y_permuted @ projection
        residual_permuted = y_permuted - beta_permuted @ x.T
        # Use the SAME HC3-studentized statistic in the observed and every
        # permuted fit. An ordinary-OLS permutation SE is not interchangeable.
        se_permuted = np.sqrt(
            (residual_permuted**2 / (1.0 - leverage)[None, :]**2)
            @ (projection[:, coefficient_index]**2)
        )
        t_permuted = beta_permuted[:, coefficient_index] / np.maximum(se_permuted, 1e-15)
        observed_t = float(model.params[coefficient_index] / model.bse[coefficient_index])
        permutation_p = float(
            (1 + np.sum(np.abs(t_permuted) >= abs(observed_t))) / (permutations + 1)
        )
        effect = float(model.params[coefficient_index])
        output.append(
            {
                "contrast": contrast,
                "effect": effect,
                "robust_se": float(model.bse[coefficient_index]),
                "ci95_low": float(model.conf_int()[coefficient_index, 0]),
                "ci95_high": float(model.conf_int()[coefficient_index, 1]),
                "hc3_p_value": float(model.pvalues[coefficient_index]),
                "freedman_lane_permutation_p_value": permutation_p,
                "n_animals": len(y),
                "n_water": int((data["treatment"] == "Water").sum()),
                "n_sucrose": int((data["treatment"] == "Sucrose").sum()),
                "n_allulose": int((data["treatment"] == "Allulose").sum()),
                "cohort_adjusted": "true",
            }
        )
    return output


def add_fdr(frame: pd.DataFrame, p_column: str) -> pd.DataFrame:
    output = frame.copy()
    output["q_value_global"] = multipletests(output[p_column], method="fdr_bh")[1]
    output["q_value_within_organ"] = np.nan
    for organ, indices in output.groupby("organ").groups.items():
        output.loc[indices, "q_value_within_organ"] = multipletests(
            output.loc[indices, p_column], method="fdr_bh"
        )[1]
    return output


def composition_contrasts(
    clr: pd.DataFrame, permutations: int, seed: int
) -> pd.DataFrame:
    rows = []
    for organ_index, organ in enumerate(ORGANS):
        for class_index, class_name in enumerate(CLASS_ORDER):
            subset = clr[(clr.organ == organ) & (clr.model_class_raw == class_name)]
            eligible = bool(subset.analysis_eligible.iloc[0])
            eligibility_fields = {
                "analysis_eligible": eligible,
                "detected_animals": int(subset.detected_animals.iloc[0]),
                "total_class_count": int(subset.total_class_count.iloc[0]),
                "eligibility_rule": str(subset.eligibility_rule.iloc[0]),
                "exclusion_reason": (
                    ""
                    if eligible
                    else "below prespecified organ-level prevalence/count threshold"
                ),
            }
            if not eligible:
                for contrast in CONTRASTS:
                    rows.append(
                        {
                            "organ": organ,
                            "model_class_raw": class_name,
                            "phenotype_display": CLASS_SHORT[class_name],
                            "outcome": "CLR model-predicted phenotype composition",
                            "effect_scale": "CLR log-ratio units",
                            "multiplicative_relative_abundance": np.nan,
                            "contrast": contrast,
                            "effect": np.nan,
                            "robust_se": np.nan,
                            "ci95_low": np.nan,
                            "ci95_high": np.nan,
                            "hc3_p_value": 1.0,
                            "freedman_lane_permutation_p_value": 1.0,
                            "n_animals": len(subset),
                            "n_water": int((subset.treatment == "Water").sum()),
                            "n_sucrose": int((subset.treatment == "Sucrose").sum()),
                            "n_allulose": int((subset.treatment == "Allulose").sum()),
                            "cohort_adjusted": "not_tested",
                            **eligibility_fields,
                        }
                    )
                continue
            result = regression_one(
                subset,
                "clr_value",
                permutations,
                seed + organ_index * 1000 + class_index,
            )
            for row in result:
                rows.append(
                    {
                        "organ": organ,
                        "model_class_raw": class_name,
                        "phenotype_display": CLASS_SHORT[class_name],
                        "outcome": "CLR model-predicted phenotype composition",
                        "effect_scale": "CLR log-ratio units",
                        "multiplicative_relative_abundance": math.exp(
                            float(row["effect"])
                        ),
                        **eligibility_fields,
                        **row,
                    }
                )
    return add_fdr(pd.DataFrame(rows), "freedman_lane_permutation_p_value")


def deviation_values(counts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for organ in ORGANS:
        organ_counts = counts[counts.organ == organ].reset_index(drop=True)
        eligible, _, _ = composition_eligibility(organ_counts)
        used = eligible.index[eligible].tolist()
        matrix = organ_counts[used].to_numpy(float) + 0.5
        proportions = matrix / matrix.sum(axis=1, keepdims=True)
        log_values = np.log(proportions)
        clr = log_values - log_values.mean(axis=1, keepdims=True)
        water_indices = np.flatnonzero(organ_counts.treatment.to_numpy() == "Water")
        for index, sample in organ_counts.iterrows():
            reference_indices = (
                water_indices[water_indices != index]
                if sample.treatment == "Water"
                else water_indices
            )
            if len(reference_indices) < 6:
                raise FigS3Error(
                    "Water-reference deviation requires at least six reference animals"
                )
            estimator = LedoitWolf().fit(clr[reference_indices])
            distance = float(np.sqrt(estimator.mahalanobis(clr[index : index + 1])[0]))
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "treatment": sample.treatment,
                    "cohort": sample.cohort,
                    "organ": organ,
                    "water_reference_deviation": distance,
                    "log1p_water_reference_deviation": math.log1p(distance),
                    "reference_water_n": len(reference_indices),
                    "reference_rule": "leave-one-out for Water; all Water controls for sugar groups",
                    "composition_classes_used": "|".join(used),
                    "composition_class_count": len(used),
                }
            )
    return pd.DataFrame(rows)


def deviation_contrasts(
    values: pd.DataFrame, permutations: int, seed: int
) -> pd.DataFrame:
    rows = []
    for organ_index, organ in enumerate(ORGANS):
        results = regression_one(
            values[values.organ == organ],
            "log1p_water_reference_deviation",
            permutations,
            seed + 9000 + organ_index,
        )
        for row in results:
            rows.append(
                {
                    "organ": organ,
                    "outcome": "log1p Water-reference CLR Mahalanobis deviation",
                    **row,
                }
            )
    return add_fdr(pd.DataFrame(rows), "freedman_lane_permutation_p_value")


def deviation_anova(values: pd.DataFrame) -> pd.DataFrame:
    """One-way treatment ANOVA requested for the Water-deviation display."""
    rows = []
    for organ in ORGANS:
        subset = values[values.organ == organ]
        groups = [
            subset.loc[
                subset.treatment == treatment, "water_reference_deviation"
            ].to_numpy(float)
            for treatment in TREATMENTS
        ]
        result = stats.f_oneway(*groups)
        rows.append(
            {
                "organ": organ,
                "test": "one-way ANOVA",
                "outcome": "cross-validated CLR Mahalanobis distance from Water",
                "factor": "treatment",
                "f_statistic": float(result.statistic),
                "p_value": float(result.pvalue),
                "df_between": len(TREATMENTS) - 1,
                "df_within": len(subset) - len(TREATMENTS),
                "n_animals": len(subset),
            }
        )
    return pd.DataFrame(rows)


def pca_scores(counts: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[float]]]:
    rows = []
    variance = {}
    for organ in ORGANS:
        subset = counts[counts.organ == organ].reset_index(drop=True)
        eligible, _, _ = composition_eligibility(subset)
        used = eligible.index[eligible].tolist()
        matrix = subset[used].to_numpy(float) + 0.5
        composition = matrix / matrix.sum(axis=1, keepdims=True)
        clr = np.log(composition)
        clr -= clr.mean(axis=1, keepdims=True)
        pca = PCA(n_components=2, svd_solver="full").fit(clr)
        scores = pca.transform(clr)
        variance[organ] = [float(value) for value in pca.explained_variance_ratio_]
        for index, sample in subset.iterrows():
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "treatment": sample.treatment,
                    "cohort": sample.cohort,
                    "organ": organ,
                    "PC1": scores[index, 0],
                    "PC2": scores[index, 1],
                    "PC1_variance_fraction": variance[organ][0],
                    "PC2_variance_fraction": variance[organ][1],
                    "composition_classes_used": "|".join(used),
                    "composition_class_count": len(used),
                }
            )
    return pd.DataFrame(rows), variance


def classical_animal_values(tile_cv: pd.DataFrame) -> pd.DataFrame:
    weighted = ("hematoxylin_mean", "eosin_mean")
    rows = []
    for keys, frame in tile_cv.groupby(
        ["sample_id", "treatment", "cohort", "organ"], sort=True
    ):
        area = frame.analyzed_tissue_area_mm2.to_numpy(float)
        row = dict(zip(("sample_id", "treatment", "cohort", "organ"), keys))
        row["selected_tile_count"] = len(frame)
        row["analyzed_tissue_area_mm2"] = float(area.sum())
        row["classical_nuclei_count"] = int(frame.classical_nuclei_count.sum())
        row["classical_nuclei_density_mm2"] = (
            row["classical_nuclei_count"] / row["analyzed_tissue_area_mm2"]
        )
        for name in weighted:
            row[name] = float(np.average(frame[name], weights=area))
        for name in (
            "tissue_fraction",
            "laplacian_variance",
            "classical_nucleus_area_median_um2",
            "classical_nucleus_eccentricity_median",
            "classical_nucleus_solidity_median",
        ):
            row[name] = float(frame[name].median())
        rows.append(row)
    result = pd.DataFrame(rows)
    if len(result) != 72 or set(result.selected_tile_count) != {24}:
        raise FigS3Error(
            "Classical CV aggregation is not balanced at 24 tiles per animal-organ"
        )
    return result


CV_FEATURES = {
    "classical_nuclei_density_mm2": "Nuclear density",
    "classical_nucleus_area_median_um2": "Nuclear area",
    "classical_nucleus_eccentricity_median": "Eccentricity",
    "classical_nucleus_solidity_median": "Solidity",
    "hematoxylin_mean": "Hematoxylin",
    "eosin_mean": "Eosin",
}


def classical_contrasts(
    values: pd.DataFrame, permutations: int, seed: int
) -> pd.DataFrame:
    rows = []
    for organ_index, organ in enumerate(ORGANS):
        subset = values[values.organ == organ].copy()
        for feature_index, (feature_name, display) in enumerate(CV_FEATURES.items()):
            mean = subset[feature_name].mean()
            std = subset[feature_name].std(ddof=1)
            subset["standardized_value"] = (subset[feature_name] - mean) / std
            results = regression_one(
                subset,
                "standardized_value",
                permutations,
                seed + 12000 + organ_index * 100 + feature_index,
            )
            for row in results:
                rows.append(
                    {
                        "organ": organ,
                        "feature": feature_name,
                        "feature_display": display,
                        "outcome": "animal-level classical H&E computer-vision feature",
                        "effect_scale": "within-organ standard deviations",
                        **row,
                    }
                )
    return add_fdr(pd.DataFrame(rows), "freedman_lane_permutation_p_value")


def save_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def displayed_cell_bar_values(phenotype: pd.DataFrame) -> pd.DataFrame:
    """Animal-level observed counts/fractions for the 13 displayed model outputs."""
    keys = ["sample_id", "treatment", "cohort", "organ"]
    displayed = phenotype[phenotype.model_class_raw.isin(CLASS_ORDER)].copy()
    displayed = displayed[
        keys + ["model_class_raw", "cell_count_in_sampled_tiles"]
    ].rename(columns={"cell_count_in_sampled_tiles": "displayed_class_count"})
    totals = (
        displayed.groupby(keys)
        .displayed_class_count.sum()
        .rename("displayed_total_count")
    )
    displayed = displayed.merge(totals, on=keys, how="left")
    displayed["displayed_class_fraction"] = (
        displayed.displayed_class_count / displayed.displayed_total_count
    )
    displayed["displayed_class_percent"] = displayed.displayed_class_fraction * 100
    displayed["log10_displayed_class_count_plus_1"] = np.log10(
        displayed.displayed_class_count + 1
    )
    displayed["phenotype_display"] = displayed.model_class_raw.map(CLASS_SHORT)
    expected = 24 * 3 * len(CLASS_ORDER)
    if len(displayed) != expected:
        raise FigS3Error(
            f"Expected {expected} animal-organ displayed-cell rows, found {len(displayed)}"
        )
    return displayed


def translate_panel_to_spanish(figure, extra: dict[str, str]) -> dict[str, str]:
    """Translate a saved English panel in place, including inset tick labels.

    Every panel here draws into inset axes, and an inset is a child axes that
    never appears in ``figure.axes``. The shared helper therefore re-pins the
    categorical tick labels of the outer axes only, and a redraw regenerates the
    inset's labels from its fixed formatter in English. Walking the child axes
    keeps organ, cell-type and feature tick labels in Spanish.
    """
    receipt = dict(translate_figure_texts_to_spanish(figure, extra=extra))
    pending = list(figure.axes)
    seen: set[int] = set()
    while pending:
        axis = pending.pop()
        if id(axis) in seen:
            continue
        seen.add(id(axis))
        pending.extend(getattr(axis, "child_axes", []))
        for get_ticks, get_labels, set_ticks in (
            (axis.get_xticks, axis.get_xticklabels, axis.set_xticks),
            (axis.get_yticks, axis.get_yticklabels, axis.set_yticks),
        ):
            source = [label.get_text() for label in get_labels()]
            translated = [translate_text(value, extra=extra) for value in source]
            if translated == source:
                continue
            set_ticks(get_ticks(), labels=translated)
            receipt.update({a: b for a, b in zip(source, translated) if a != b})
    figure.canvas.draw()
    return dict(sorted(receipt.items()))


def significance_stars(value: float) -> str:
    if not np.isfinite(value):
        return ""
    if value < 0.0001:
        return "****"
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""


def style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_letter(ax, letter: str, title: str) -> None:
    """Place the bold letter in a fixed points-sized gutter, never in axes units.

    An axes-fraction offset scales with the panel, so the wide right-hand panels
    pushed their letter far enough left to collide with the left neighbour's
    text. A points offset keeps every letter the same short distance outside its
    own panel and inside the column gutter.
    """
    for text, offset, size in ((letter, -10.0, 15.0), (title, 7.0, 9.0)):
        ax.annotate(
            text,
            xy=(0.0, 1.0),
            xycoords="axes fraction",
            xytext=(offset, 5.0),
            textcoords="offset points",
            fontsize=size,
            fontweight="bold",
            ha="left",
            va="baseline",
            annotation_clip=False,
        )


def component_title(ax, title: str) -> None:
    ax.set_title(title, fontsize=7.5, fontweight="bold", pad=3)


def panel_side_legend(ax, handles, ncol: int, fontsize: float = 6.3) -> None:
    """Park a legend on the panel's title row so it cannot cover plotted data."""
    ax.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.055),
        ncol=ncol,
        fontsize=fontsize,
        frameon=False,
        handletextpad=0.45,
        columnspacing=1.1,
        borderaxespad=0.0,
    )


def plot_workflow(ax, standalone: bool = True) -> None:
    ax.set_axis_off()
    if standalone:
        panel_letter(ax, "A", "Multi-organ WSI workflow")
    else:
        component_title(ax, "Workflow")
    boxes = [
        ("24 animals\nWater 7 | Sucrose 9 | Allulose 8", "#e8f1fa"),
        ("One H&E WSI\nLiver + kidney + spleen", "#f4e8fa"),
        ("L0/L2 + CV masks\nwithin-slide organ labels", "#e9f6ec"),
        ("24 L0 tiles/organ\n0.5 µm/px; spatially balanced", "#fff3dd"),
        ("HistoPLUS + classical CV\nanimal-level composition", "#fce8e8"),
    ]
    y_positions = np.linspace(0.845, 0.055, len(boxes))
    for index, ((text, face), y) in enumerate(zip(boxes, y_positions)):
        box = FancyBboxPatch(
            (0.05, y),
            0.90,
            0.145,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor=face,
            edgecolor="#555555",
            linewidth=0.8,
            transform=ax.transAxes,
        )
        ax.add_patch(box)
        ax.text(
            0.5,
            y + 0.0725,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7.0,
        )
        if index < len(boxes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (0.5, y - 0.004),
                    (0.5, y_positions[index + 1] + 0.149),
                    arrowstyle="-|>",
                    mutation_scale=8,
                    color="#555555",
                    transform=ax.transAxes,
                )
            )


def representative_sample(
    root: Path, slide: pd.DataFrame
) -> tuple[str, dict[str, object]]:
    """Select a Water example without using cross-slide image alignment.

    Tissue fraction is an intrinsic per-slide segmentation measurement. The
    Water animal nearest the Water median is deterministic, transparent, and
    independent of every treatment-comparison endpoint displayed later.
    """
    water_ids = sorted(slide.loc[slide.treatment == "Water", "sample_id"].tolist())
    if len(water_ids) != 7:
        raise FigS3Error(
            f"Representative selection requires seven Water animals, found {len(water_ids)}"
        )
    fractions: dict[str, float] = {}
    for sample_id in water_ids:
        receipt = json.loads(
            require_regular_file(
                root / "organ_segmentation" / sample_id / "segmentation_receipt.json",
                f"organ-segmentation receipt for {sample_id}",
            ).read_text(encoding="utf-8")
        )
        fraction = float(receipt.get("tissue_fraction", math.nan))
        if not math.isfinite(fraction) or not 0.0 < fraction < 1.0:
            raise FigS3Error(f"Invalid tissue fraction for {sample_id}: {fraction}")
        fractions[sample_id] = fraction
    median_fraction = float(np.median(list(fractions.values())))
    selected = min(
        water_ids,
        key=lambda sample_id: (abs(fractions[sample_id] - median_fraction), sample_id),
    )
    selection = {
        "method": "Water animal closest to median whole-slide tissue fraction; sample_id breaks exact ties",
        "metric_source": "organ_segmentation/<sample_id>/segmentation_receipt.json:tissue_fraction",
        "eligible_treatment": "Water",
        "eligible_sample_count": len(water_ids),
        "candidate_tissue_fractions": fractions,
        "median_tissue_fraction": median_fraction,
        "selected_sample_id": selected,
        "selected_tissue_fraction": fractions[selected],
    }
    return selected, selection


def plot_organ_overview(
    ax,
    root: Path,
    sample_id: str,
    standalone: bool = True,
) -> None:
    image = np.asarray(
        Image.open(root / "organ_segmentation" / sample_id / "organ_overlay.png")
    )
    overview_width_px = int(image.shape[1])
    conversion_receipt = json.loads(
        require_regular_file(
            root / "exports" / sample_id / "conversion_receipt.json",
            "WSI conversion receipt",
        ).read_text(encoding="utf-8")
    )
    source_mpp = float(conversion_receipt["source_mpp"])
    true_width_l0 = float(conversion_receipt["true_width"])
    overview_mpp = source_mpp * true_width_l0 / overview_width_px
    if not math.isfinite(overview_mpp) or overview_mpp <= 0:
        raise FigS3Error(
            f"Invalid organ-overview calibration for {sample_id}: {overview_mpp}"
        )
    labels = np.asarray(
        Image.open(root / "organ_segmentation" / sample_id / "organ_label_mask.png")
    )
    # Trim the empty slide background so the three organs, not the scanner
    # margin, use the panel area. Image and labels are cropped identically, so
    # the organ label centroids below stay correct.
    label_rows, label_columns = np.nonzero(labels > 0)
    if label_rows.size == 0:
        raise FigS3Error(f"Empty organ label mask for {sample_id}")
    margin = max(4, int(round(0.02 * max(labels.shape))))
    row_0 = max(0, int(label_rows.min()) - margin)
    row_1 = min(labels.shape[0], int(label_rows.max()) + 1 + margin)
    column_0 = max(0, int(label_columns.min()) - margin)
    column_1 = min(labels.shape[1], int(label_columns.max()) + 1 + margin)
    image = image[row_0:row_1, column_0:column_1]
    labels = labels[row_0:row_1, column_0:column_1]
    ax.imshow(image)
    ax.set_axis_off()
    if standalone:
        panel_letter(ax, "B", f"Three organs segmented in one WSI ({sample_id})")
    else:
        component_title(ax, f"Organ segmentation ({sample_id})")
    for organ in ORGANS:
        y, x = np.nonzero(labels == ORGAN_LABEL_VALUE[organ])
        if len(x) == 0:
            raise FigS3Error(f"Missing {organ} pixels in representative organ mask")
        ax.text(
            float(np.median(x)),
            float(np.median(y)),
            organ,
            ha="center",
            va="center",
            color="white",
            fontsize=8,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.22",
                "facecolor": "#111111",
                "edgecolor": "white",
                "alpha": 0.82,
            },
        )

    scale_um = 2000.0
    scale_px = scale_um / overview_mpp
    if scale_px >= 0.48 * image.shape[1]:
        raise FigS3Error(
            f"Organ-overview 2-mm scale bar is too wide for {sample_id}: {scale_px:.1f}px"
        )
    # The representative spleen and its label occupy the lower-right field.
    # Keep the calibrated bar in the empty lower-left scanner background.
    x_0 = max(18.0, 0.045 * image.shape[1])
    x_1 = x_0 + scale_px
    y_bar = image.shape[0] - max(18.0, 0.040 * image.shape[0])
    ax.plot(
        [x_0, x_1],
        [y_bar, y_bar],
        color="black",
        linewidth=6.0,
        solid_capstyle="butt",
        zorder=20,
    )
    ax.plot(
        [x_0, x_1],
        [y_bar, y_bar],
        color="white",
        linewidth=3.3,
        solid_capstyle="butt",
        zorder=21,
    )
    label = ax.text(
        (x_0 + x_1) / 2.0,
        y_bar - max(7.0, 0.012 * image.shape[0]),
        "2 mm",
        color="white",
        ha="center",
        va="bottom",
        fontsize=7.0,
        fontweight="bold",
        zorder=22,
    )
    label.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="black")])


def plot_combined_panel_a(
    ax,
    root: Path,
    sample_id: str,
) -> None:
    """Compact Panel A containing the workflow and within-slide organ map."""
    ax.set_axis_off()
    panel_letter(
        ax,
        "A",
        "Multi-organ WSI workflow and within-slide organ segmentation",
    )
    workflow = ax.inset_axes([0.005, 0.02, 0.615, 0.895])
    organs = ax.inset_axes([0.645, 0.02, 0.335, 0.895])
    plot_workflow(workflow, standalone=False)
    plot_organ_overview(organs, root, sample_id, standalone=False)


def vips_crop_rgb(
    path: Path, x: int, y: int, width: int, height: int, target: int = 840
) -> Image.Image:
    slide = pyvips.Image.new_from_file(str(path), access="random")
    crop = slide.crop(x, y, width, height)
    crop = crop.resize(target / width, vscale=target / height)
    if crop.format != "uchar":
        crop = crop.cast("uchar")
    array = np.frombuffer(crop.write_to_memory(), dtype=np.uint8).reshape(
        crop.height, crop.width, crop.bands
    )
    image = Image.fromarray(array[..., :3], mode="RGB")
    return image.resize((target, target), Image.Resampling.LANCZOS)


def draw_scale_bar(
    image: Image.Image,
    crop_width_l0: float,
    length_um: int,
    source_mpp: float = 0.26178,
) -> None:
    """Add a calibrated, high-contrast scale bar to an L0-derived crop."""
    pixel_length = max(
        8,
        int(round(image.width * length_um / (crop_width_l0 * source_mpp))),
    )
    x1 = image.width - 22
    x0 = x1 - pixel_length
    y = image.height - 24
    draw = ImageDraw.Draw(image)
    draw.line((x0, y, x1, y), fill="#111111", width=10)
    draw.line((x0, y, x1, y), fill="white", width=6)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 22)
    except OSError:
        font = ImageFont.load_default()
    draw.text(
        (x0, y - 30),
        f"{length_um} µm",
        fill="white",
        font=font,
        stroke_width=3,
        stroke_fill="#111111",
    )


def organ_histology_overlays(
    root: Path, sample_id: str
) -> dict[str, tuple[Image.Image, str, int]]:
    manifest = pd.read_csv(
        root / "histoplus" / sample_id / "selected_tile_manifest.csv"
    )
    cells = pd.read_parquet(
        root / "histoplus" / sample_id / "cells.parquet",
        columns=[
            "tile_uid",
            "organ",
            "model_class_raw",
            "centroid_x_l0",
            "centroid_y_l0",
        ],
    )
    l0_path = root / "exports" / sample_id / "1_L0_rgb.tif"
    outputs: dict[str, tuple[Image.Image, str, int]] = {}
    for organ in ORGANS:
        candidates = manifest[manifest.organ == organ].copy()
        counts = cells[cells.organ == organ].groupby("tile_uid").size()
        candidates["cell_count"] = candidates.tile_uid.map(counts).fillna(0).astype(int)
        median_count = float(candidates.cell_count.median())
        candidates["distance_to_median"] = (candidates.cell_count - median_count).abs()
        chosen = candidates.sort_values(["distance_to_median", "sampling_rank"]).iloc[0]
        image = vips_crop_rgb(
            l0_path,
            int(chosen.x_l0),
            int(chosen.y_l0),
            int(chosen.width_l0),
            int(chosen.height_l0),
        )
        draw = ImageDraw.Draw(image)
        tile_cells = cells[
            (cells.tile_uid == chosen.tile_uid)
            & (cells.model_class_raw != "Cancer cell")
        ]
        x_scale = image.width / float(chosen.width_l0)
        y_scale = image.height / float(chosen.height_l0)
        for cell in tile_cells.itertuples():
            x = (float(cell.centroid_x_l0) - float(chosen.x_l0)) * x_scale
            y = (float(cell.centroid_y_l0) - float(chosen.y_l0)) * y_scale
            color_value = CELL_COLORS.get(str(cell.model_class_raw), "#333333")
            draw.ellipse(
                (x - 2.5, y - 2.5, x + 2.5, y + 2.5), outline=color_value, width=2
            )
        draw_scale_bar(image, float(chosen.width_l0), 100)
        outputs[organ] = (image, str(chosen.tile_uid), int(len(tile_cells)))
    return outputs


def _draw_target_polygons(
    image: Image.Image,
    cells: pd.DataFrame,
    x_l0: float,
    y_l0: float,
    width_l0: float,
    height_l0: float,
) -> None:
    """Draw true model-output nuclear contours on an L0 crop."""
    draw = ImageDraw.Draw(image)
    x_scale = image.width / width_l0
    y_scale = image.height / height_l0
    for cell in cells.itertuples():
        geometry = shapely.from_wkb(cell.geometry)
        polygons = (
            list(geometry.geoms) if geometry.geom_type == "MultiPolygon" else [geometry]
        )
        for polygon in polygons:
            if polygon.is_empty or polygon.geom_type != "Polygon":
                continue
            points = [
                (
                    (float(x) - x_l0) * x_scale,
                    (float(y) - y_l0) * y_scale,
                )
                for x, y in polygon.exterior.coords
            ]
            if len(points) < 3:
                continue
            draw.line(points, fill="#111111", width=5, joint="curve")
            draw.line(points, fill="#FFD92F", width=3, joint="curve")


def significant_cell_examples(
    root: Path,
    phenotype: pd.DataFrame,
    composition_stats: pd.DataFrame,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    """Build paired L0 insets only for globally BH-significant phenotype effects."""
    significant = (
        composition_stats[
            composition_stats.analysis_eligible
            & (composition_stats.q_value_global < 0.05)
        ]
        .sort_values(["q_value_global", "organ", "model_class_raw", "contrast"])
        .head(1)
    )
    examples: list[dict[str, object]] = []
    records: list[dict[str, object]] = []
    for result in significant.itertuples():
        sugar = str(result.contrast).split(" vs ")[0]
        displayed = phenotype[
            (phenotype.organ == result.organ)
            & phenotype.model_class_raw.isin(CLASS_ORDER)
        ].copy()
        totals = (
            displayed.groupby(["sample_id", "treatment"])["cell_count_in_sampled_tiles"]
            .sum()
            .rename("displayed_total")
        )
        target = (
            displayed[displayed.model_class_raw == result.model_class_raw]
            .set_index(["sample_id", "treatment"])["cell_count_in_sampled_tiles"]
            .rename("target_count")
        )
        fractions = pd.concat([totals, target], axis=1).fillna(0).reset_index()
        fractions["displayed_class_fraction"] = (
            fractions.target_count / fractions.displayed_total
        )
        pair: dict[str, dict[str, object]] = {}
        for treatment in ("Water", sugar):
            group = fractions[fractions.treatment == treatment].copy()
            group_median = float(group.displayed_class_fraction.median())
            group["median_distance"] = (
                group.displayed_class_fraction - group_median
            ).abs()
            chosen_animal = group.sort_values(["median_distance", "sample_id"]).iloc[0]
            sample_id = str(chosen_animal.sample_id)
            cells = pd.read_parquet(
                root / "histoplus" / sample_id / "cells.parquet",
                columns=["geometry", "tile_uid", "organ", "model_class_raw"],
            )
            cells = cells[
                (cells.organ == result.organ) & cells.model_class_raw.isin(CLASS_ORDER)
            ].copy()
            total_by_tile = cells.groupby("tile_uid").size().rename("displayed_total")
            target_cells = cells[cells.model_class_raw == result.model_class_raw]
            target_by_tile = (
                target_cells.groupby("tile_uid").size().rename("target_count")
            )
            manifest = pd.read_csv(
                root / "histoplus" / sample_id / "selected_tile_manifest.csv"
            )
            manifest = manifest[manifest.organ == result.organ].copy()
            manifest = manifest.merge(total_by_tile, on="tile_uid", how="left")
            manifest = manifest.merge(target_by_tile, on="tile_uid", how="left")
            manifest[["displayed_total", "target_count"]] = manifest[
                ["displayed_total", "target_count"]
            ].fillna(0)
            manifest["tile_fraction"] = (
                manifest.target_count / manifest.displayed_total.replace(0, np.nan)
            ).fillna(0)
            manifest["animal_distance"] = (
                manifest.tile_fraction - float(chosen_animal.displayed_class_fraction)
            ).abs()
            positive = manifest[manifest.target_count > 0]
            tile_pool = positive if not positive.empty else manifest
            chosen_tile = tile_pool.sort_values(
                ["animal_distance", "sampling_rank", "tile_uid"]
            ).iloc[0]
            image = vips_crop_rgb(
                root / "exports" / sample_id / "1_L0_rgb.tif",
                int(chosen_tile.x_l0),
                int(chosen_tile.y_l0),
                int(chosen_tile.width_l0),
                int(chosen_tile.height_l0),
                target=680,
            )
            tile_targets = target_cells[target_cells.tile_uid == chosen_tile.tile_uid]
            _draw_target_polygons(
                image,
                tile_targets,
                float(chosen_tile.x_l0),
                float(chosen_tile.y_l0),
                float(chosen_tile.width_l0),
                float(chosen_tile.height_l0),
            )
            draw_scale_bar(image, float(chosen_tile.width_l0), 50)
            pair[treatment] = {
                "image": image,
                "sample_id": sample_id,
                "tile_uid": str(chosen_tile.tile_uid),
                "sample_fraction": float(chosen_animal.displayed_class_fraction),
                "group_median_fraction": group_median,
                "tile_fraction": float(chosen_tile.tile_fraction),
                "target_count": int(chosen_tile.target_count),
                "displayed_total": int(chosen_tile.displayed_total),
            }
            records.append(
                {
                    "organ": result.organ,
                    "model_class_raw": result.model_class_raw,
                    "phenotype_display": CLASS_SHORT[str(result.model_class_raw)],
                    "contrast": result.contrast,
                    "effect_clr": float(result.effect),
                    "multiplicative_relative_abundance": float(
                        result.multiplicative_relative_abundance
                    ),
                    "q_value_global": float(result.q_value_global),
                    "treatment_displayed": treatment,
                    "sample_id": sample_id,
                    "group_median_displayed_class_fraction": group_median,
                    "sample_displayed_class_fraction": float(
                        chosen_animal.displayed_class_fraction
                    ),
                    "tile_uid": str(chosen_tile.tile_uid),
                    "tile_displayed_class_fraction": float(chosen_tile.tile_fraction),
                    "target_cell_count_in_tile": int(chosen_tile.target_count),
                    "displayed_cell_count_in_tile": int(chosen_tile.displayed_total),
                    "selection_rule": (
                        "animal closest to treatment median displayed-class fraction; "
                        "positive tile closest to selected-animal displayed-class fraction"
                    ),
                    "outline_definition": "true HistoPLUS/LazySlide model-output nuclear polygon",
                    "source_l0_path": str(
                        (root / "exports" / sample_id / "1_L0_rgb.tif").resolve()
                    ),
                }
            )
        examples.append(
            {
                "organ": str(result.organ),
                "model_class_raw": str(result.model_class_raw),
                "phenotype_display": CLASS_SHORT[str(result.model_class_raw)],
                "contrast": str(result.contrast),
                "sugar": sugar,
                "relative_abundance": float(result.multiplicative_relative_abundance),
                "q_value_global": float(result.q_value_global),
                "pair": pair,
            }
        )
    return examples, pd.DataFrame.from_records(records)


def plot_histoplus(
    ax,
    overlays: dict[str, tuple[Image.Image, str, int]],
    significant_examples: list[dict[str, object]],
    sample_id: str,
) -> None:
    ax.set_axis_off()
    panel_letter(
        ax,
        "C",
        "Three-organ histology and representative phenotype examples",
    )
    # The crops are square, so the inset boxes are square in figure inches too;
    # otherwise each image letterboxes and the strip reads as five widely
    # separated tiles instead of one continuous row.
    for x, organ in zip((0.014, 0.2006, 0.3872), ORGANS):
        inset = ax.inset_axes([x, 0.030, 0.1686, 0.845])
        image, tile_uid, cell_count = overlays[organ]
        inset.imshow(image)
        inset.set_axis_off()
        inset.set_title(
            f"{organ} | Water ({sample_id})\n{cell_count:,} instances",
            fontsize=7,
            fontweight="bold",
            pad=2,
        )
    if significant_examples:
        example = significant_examples[0]
        pair = example["pair"]
        ax.text(
            0.8084,
            0.998,
            (
                f"BH-significant: {example['organ']} "
                f"{str(example['phenotype_display']).lower()} · "
                f"{example['sugar']} vs Water "
                f"{example['relative_abundance']:.2f}×; q={example['q_value_global']:.3f}"
            ),
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.5,
            fontweight="bold",
        )
        for x, treatment in ((0.6308, "Water"), (0.8174, str(example["sugar"]))):
            data = pair[treatment]
            inset = ax.inset_axes([x, 0.030, 0.1686, 0.845])
            inset.imshow(data["image"])
            inset.set_xticks([])
            inset.set_yticks([])
            for spine in inset.spines.values():
                spine.set_color(TREATMENT_COLORS[treatment])
                spine.set_linewidth(2.2)
            inset.set_title(
                f"{treatment}\n{data['sample_fraction'] * 100:.1f}% of displayed outputs",
                fontsize=5.9,
                color=TREATMENT_COLORS[treatment],
                fontweight="bold",
                pad=2,
            )
    else:
        ax.text(
            0.5,
            0.20,
            "No phenotype-composition contrast survived global BH q<0.05",
            transform=ax.transAxes,
            ha="center",
            fontsize=7,
            fontweight="bold",
        )


def effect_matrix(
    frame: pd.DataFrame, row_name: str, row_order: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ordered_columns = [(organ, contrast) for organ in ORGANS for contrast in CONTRASTS]
    effects = frame.pivot(
        index=row_name, columns=["organ", "contrast"], values="effect"
    ).reindex(index=row_order, columns=pd.MultiIndex.from_tuples(ordered_columns))
    qvalues = frame.pivot(
        index=row_name, columns=["organ", "contrast"], values="q_value_global"
    ).reindex(index=row_order, columns=pd.MultiIndex.from_tuples(ordered_columns))
    effects.columns = [
        f"{organ}\n{contrast.split()[0]}" for organ, contrast in ordered_columns
    ]
    qvalues.columns = effects.columns
    return effects, qvalues


def plot_effect_heatmap(
    ax,
    contrasts: pd.DataFrame,
    standalone: bool = True,
) -> None:
    effects, qvalues = effect_matrix(contrasts, "model_class_raw", CLASS_ORDER)
    vmax = max(0.35, float(np.nanquantile(np.abs(effects.to_numpy()), 0.85)))
    display = effects.copy()
    display.index = [CLASS_SHORT[name] for name in display.index]
    sns.heatmap(
        display,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        cbar_kws={"label": "CLR effect vs Water", "shrink": 0.62, "pad": 0.02},
        linewidths=0.55,
        linecolor="#f7f7f7",
    )
    for row_index in range(qvalues.shape[0]):
        for column_index in range(qvalues.shape[1]):
            stars = significance_stars(float(qvalues.iloc[row_index, column_index]))
            if not stars:
                continue
            label = ax.text(
                column_index + 0.5,
                row_index + 0.5,
                stars,
                ha="center",
                va="center",
                color="white",
                fontsize=11,
                fontweight="bold",
                zorder=8,
            )
            label.set_path_effects(
                [path_effects.withStroke(linewidth=2.2, foreground="#111111")]
            )
    ax.set_facecolor("#d9d9d9")
    if standalone:
        panel_letter(ax, "E", "Sugar-associated phenotype-composition effects")
    else:
        component_title(ax, "Animal-level CLR effects vs Water")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    ax.text(
        1.0,
        -0.125,
        "global BH q: *<0.05, **<0.01, ***<0.001, ****<0.0001; gray = filtered",
        transform=ax.transAxes,
        ha="right",
        fontsize=6.2,
    )


def plot_cell_barplots(
    fraction_ax,
    count_ax,
    contrasts: pd.DataFrame,
    cell_values: pd.DataFrame,
) -> None:
    significant = contrasts[
        contrasts.analysis_eligible & contrasts.q_value_global.lt(0.05)
    ].sort_values(["q_value_global", "organ", "model_class_raw", "contrast"])
    if significant.empty:
        fraction_ax.set_axis_off()
        count_ax.set_axis_off()
        fraction_ax.text(
            0.5,
            0.5,
            "No global BH-significant cell output",
            ha="center",
            va="center",
            transform=fraction_ax.transAxes,
        )
        return
    result = significant.iloc[0]
    data = cell_values[
        (cell_values.organ == result.organ)
        & (cell_values.model_class_raw == result.model_class_raw)
    ].copy()
    for axis, y_name, y_label in (
        (fraction_ax, "displayed_class_percent", "% of displayed outputs"),
        (
            count_ax,
            "log10_displayed_class_count_plus_1",
            "log10(objects / 24 tiles + 1)",
        ),
    ):
        sns.barplot(
            data=data,
            x="treatment",
            y=y_name,
            order=TREATMENTS,
            hue="treatment",
            hue_order=TREATMENTS,
            palette=TREATMENT_COLORS,
            errorbar="se",
            capsize=0.12,
            width=0.70,
            edgecolor="#222222",
            linewidth=0.9,
            legend=False,
            ax=axis,
        )
        deterministic_stripplot(
            data=data,
            x="treatment",
            y=y_name,
            order=TREATMENTS,
            hue="treatment",
            hue_order=TREATMENTS,
            palette=TREATMENT_COLORS,
            dodge=False,
            size=3.6,
            edgecolor="#111111",
            linewidth=0.35,
            alpha=0.95,
            legend=False,
            ax=axis,
        )
        axis.set_xlabel("")
        axis.set_ylabel(y_label, fontsize=6.3)
        axis.tick_params(labelsize=5.8)
        axis.grid(axis="x", visible=False)
    fraction_ax.set_title(
        f"Observed {result.organ} {CLASS_SHORT[result.model_class_raw].lower()}",
        fontsize=7,
        fontweight="bold",
        pad=2,
    )
    fraction_ax.set_xticks(range(len(TREATMENTS)))
    fraction_ax.set_xticklabels([])
    sugar = str(result.contrast).split(" vs ")[0]
    sugar_index = TREATMENTS.index(sugar)
    maximum = float(data.displayed_class_percent.max())
    bracket_y = maximum * 1.12
    bracket_step = max(maximum * 0.035, 0.8)
    fraction_ax.plot(
        [0, 0, sugar_index, sugar_index],
        [bracket_y, bracket_y + bracket_step, bracket_y + bracket_step, bracket_y],
        color="#111111",
        linewidth=1.0,
        clip_on=False,
    )
    star_label = fraction_ax.text(
        sugar_index / 2,
        bracket_y + bracket_step * 1.02,
        significance_stars(float(result.q_value_global)),
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color="white",
    )
    star_label.set_path_effects(
        [path_effects.withStroke(linewidth=2.2, foreground="#111111")]
    )
    fraction_ax.text(
        0.98,
        0.97,
        f"global BH q={result.q_value_global:.3f}",
        transform=fraction_ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
    )
    fraction_ax.set_ylim(top=bracket_y + bracket_step * 5.0)
    count_ax.set_xticks(range(len(TREATMENTS)))
    count_ax.set_xticklabels(
        [
            f"{treatment}\nn={int((data.treatment == treatment).sum())}"
            for treatment in TREATMENTS
        ]
    )
    count_ax.text(
        1.0,
        -0.38,
        "bars = mean ± SE; points = animals; barplots are descriptive",
        transform=count_ax.transAxes,
        ha="right",
        fontsize=5.4,
    )


def plot_effects_with_cell_bars(
    ax,
    contrasts: pd.DataFrame,
    cell_values: pd.DataFrame,
) -> None:
    ax.set_axis_off()
    panel_letter(
        ax,
        "E",
        "Phenotype-composition effects with observed cell-output barplots",
    )
    heatmap = ax.inset_axes([0.015, 0.08, 0.62, 0.84])
    fraction_bar = ax.inset_axes([0.70, 0.55, 0.27, 0.34])
    count_bar = ax.inset_axes([0.70, 0.11, 0.27, 0.34])
    plot_effect_heatmap(heatmap, contrasts, standalone=False)
    plot_cell_barplots(fraction_bar, count_bar, contrasts, cell_values)


def plot_cell_effect_barplots(ax, contrasts: pd.DataFrame) -> None:
    """Full-width replacement for the heatmap: organ-wise cell-type effect bars."""
    ax.set_axis_off()
    panel_letter(
        ax,
        "F",
        "Model-predicted cell-type composition effects vs Water",
    )
    finite_ci = contrasts[["ci95_low", "ci95_high"]].to_numpy(float)
    limit = max(1.0, float(np.nanmax(np.abs(finite_ci))) * 1.26)
    colors = {
        "Sucrose vs Water": TREATMENT_COLORS["Sucrose"],
        "Allulose vs Water": TREATMENT_COLORS["Allulose"],
    }
    offsets = {"Sucrose vs Water": -0.19, "Allulose vs Water": 0.19}
    for organ_index, organ in enumerate(ORGANS):
        inset = ax.inset_axes([0.035 + organ_index * 0.322, 0.175, 0.300, 0.735])
        organ_rows = contrasts[(contrasts.organ == organ) & contrasts.analysis_eligible]
        classes = [
            class_name
            for class_name in CLASS_ORDER
            if class_name in set(organ_rows.model_class_raw)
        ]
        x = np.arange(len(classes), dtype=float)
        for contrast in CONTRASTS:
            rows = (
                organ_rows[organ_rows.contrast == contrast]
                .set_index("model_class_raw")
                .reindex(classes)
            )
            effect = rows.effect.to_numpy(float)
            low = rows.ci95_low.to_numpy(float)
            high = rows.ci95_high.to_numpy(float)
            positions = x + offsets[contrast]
            inset.bar(
                positions,
                effect,
                width=0.36,
                color=colors[contrast],
                edgecolor="#222222",
                linewidth=0.55,
                alpha=0.92,
                zorder=2,
            )
            inset.errorbar(
                positions,
                effect,
                yerr=np.vstack([effect - low, high - effect]),
                fmt="none",
                ecolor="#222222",
                elinewidth=0.8,
                capsize=1.8,
                zorder=3,
            )
            for class_index, row in enumerate(rows.itertuples()):
                stars = significance_stars(float(row.q_value_global))
                if not stars:
                    continue
                positive = float(row.effect) >= 0
                y_position = (
                    float(row.ci95_high) + limit * 0.035
                    if positive
                    else float(row.ci95_low) - limit * 0.035
                )
                star = inset.text(
                    positions[class_index],
                    y_position,
                    stars,
                    ha="center",
                    va="bottom" if positive else "top",
                    fontsize=STAR_FONTSIZE_BARS,
                    fontweight="bold",
                    color="#111111",
                    zorder=5,
                )
                star.set_path_effects(
                    [path_effects.withStroke(linewidth=2.6, foreground="white")]
                )
        inset.axhline(0, color="#222222", linewidth=0.8, zorder=1)
        inset.set_ylim(-limit, limit)
        inset.set_xlim(-0.65, len(classes) - 0.35)
        inset.set_title(organ, fontsize=8, fontweight="bold", pad=3)
        inset.set_xticks(x)
        inset.set_xticklabels(
            [CLASS_SHORT[name].replace("-like", "") for name in classes],
            rotation=58,
            ha="right",
            fontsize=5.6,
        )
        inset.tick_params(axis="y", labelsize=6)
        inset.grid(axis="x", visible=False)
        inset.grid(axis="y", color="#dedede", linewidth=0.55)
        inset.set_xlabel("")
        inset.set_ylabel("CLR effect vs Water" if organ_index == 0 else "", fontsize=7)
    panel_side_legend(
        ax,
        [
            Patch(
                facecolor=TREATMENT_COLORS["Sucrose"],
                edgecolor="#222222",
                label="Sucrose",
            ),
            Patch(
                facecolor=TREATMENT_COLORS["Allulose"],
                edgecolor="#222222",
                label="Allulose",
            ),
        ],
        ncol=2,
        fontsize=6.5,
    )


def plot_pca(ax, scores: pd.DataFrame) -> None:
    ax.set_axis_off()
    panel_letter(ax, "D", "Compositional ordination (CLR-PCA)")
    for index, organ in enumerate(ORGANS):
        inset = ax.inset_axes([0.045, 0.075 + (2 - index) * 0.300, 0.945, 0.240])
        data = scores[scores.organ == organ]
        for treatment in TREATMENTS:
            subset = data[data.treatment == treatment]
            inset.scatter(
                subset.PC1,
                subset.PC2,
                color=TREATMENT_COLORS[treatment],
                marker="o",
                s=27,
                edgecolor="#222222",
                linewidth=0.35,
                alpha=0.95,
            )
        inset.axhline(0, color="#dddddd", lw=0.5)
        inset.axvline(0, color="#dddddd", lw=0.5)
        # Each organ is its own ordination, so the variance explained belongs in
        # that row's own header rather than in one shared axis label that would
        # read as if it applied to all three. The header goes inside the axes on
        # an opaque patch, because an axes title would land on the x tick labels
        # of the row above once the rows are packed this close, and the y limit
        # is opened at the top so the header never covers an animal's point.
        low, high = inset.get_ylim()
        inset.set_ylim(low, high + 0.30 * (high - low))
        inset.text(
            0.010,
            0.955,
            f"{organ} · PC1 {data.PC1_variance_fraction.iloc[0] * 100:.0f}%"
            f" · PC2 {data.PC2_variance_fraction.iloc[0] * 100:.0f}%",
            transform=inset.transAxes,
            ha="left",
            va="top",
            fontsize=6.2,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.2},
            zorder=6,
        )
        inset.set_xlabel(
            "PC1" if index == len(ORGANS) - 1 else "", fontsize=6, labelpad=1.5
        )
        inset.set_ylabel("PC2", fontsize=6, labelpad=1.5)
        inset.tick_params(labelsize=5, pad=1.2)
    panel_side_legend(
        ax,
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=TREATMENT_COLORS[t],
                markeredgecolor="#222222",
                label=t,
                markersize=5,
            )
            for t in TREATMENTS
        ],
        ncol=3,
        fontsize=6.0,
    )


def plot_deviation(ax, values: pd.DataFrame, anova: pd.DataFrame) -> None:
    # Drawn into an inset for the same reason as every other panel here: the
    # rotated y-axis label of a directly plotted axes reaches into the gutter and
    # collides with the panel letter and with the left neighbour's text.
    ax.set_axis_off()
    panel_letter(ax, "B", "Deviation from Water composition")
    inset = ax.inset_axes([0.075, 0.150, 0.918, 0.760])
    sns.boxplot(
        data=values,
        x="organ",
        y="water_reference_deviation",
        hue="treatment",
        hue_order=TREATMENTS,
        palette=TREATMENT_COLORS,
        showfliers=False,
        width=0.72,
        linewidth=1.15,
        saturation=1,
        boxprops={"alpha": 0.88, "edgecolor": "#222222"},
        medianprops={"color": "#111111", "linewidth": 1.25},
        whiskerprops={"color": "#333333", "linewidth": 1.0},
        capprops={"color": "#333333", "linewidth": 1.0},
        ax=inset,
    )
    deterministic_stripplot(
        data=values,
        x="organ",
        y="water_reference_deviation",
        hue="treatment",
        hue_order=TREATMENTS,
        palette=TREATMENT_COLORS,
        dodge=True,
        size=3.6,
        edgecolor="#111111",
        linewidth=0.35,
        alpha=0.92,
        ax=inset,
        legend=False,
    )
    handles, labels = inset.get_legend_handles_labels()
    legend = inset.get_legend()
    if legend is not None:
        legend.remove()
    # The legend used to sit inside the axes at upper left, directly under the
    # Kidney ANOVA annotation. On the title row it can never overlap either.
    panel_side_legend(ax, handles[:3], ncol=3, fontsize=6.3)
    inset.set_ylabel(
        "cross-validated CLR\nMahalanobis distance", fontsize=6.6, labelpad=1.5
    )
    inset.set_xlabel("")
    inset.tick_params(axis="y", labelsize=6)
    inset.tick_params(axis="x", labelsize=7.5)
    maximum = values.water_reference_deviation.max()
    for index, organ in enumerate(ORGANS):
        row = anova[anova.organ == organ].iloc[0]
        p_text = "p<0.001" if row.p_value < 0.001 else f"p={row.p_value:.3f}"
        inset.text(
            index,
            maximum * 1.035,
            f"ANOVA {p_text}",
            ha="center",
            va="bottom",
            fontsize=6.4,
            fontweight="bold" if row.p_value < 0.05 else "normal",
            color="#111111",
        )
    inset.set_ylim(top=maximum * 1.16)


def plot_cv_forest(ax, contrasts: pd.DataFrame) -> None:
    order = list(CV_FEATURES)
    ax.set_axis_off()
    panel_letter(
        ax,
        "E",
        "Classical H&E effects vs Water (95% CI; permutation-p stars)",
    )
    data_limit = max(
        1.0,
        float(np.nanmax(np.abs(contrasts[["ci95_low", "ci95_high"]].to_numpy(float))))
        * 1.06,
    )
    limit = data_limit * 1.34
    y = np.arange(len(order))[::-1]
    # The two contrasts are pushed further apart than the marker size needs,
    # because the doubled significance asterisk is what sets the spacing here.
    offsets = {"Sucrose vs Water": 0.24, "Allulose vs Water": -0.24}
    colors = {
        "Sucrose vs Water": TREATMENT_COLORS["Sucrose"],
        "Allulose vs Water": TREATMENT_COLORS["Allulose"],
    }
    for index, organ in enumerate(ORGANS):
        inset = ax.inset_axes([0.115 + index * 0.298, 0.115, 0.272, 0.800])
        organ_rows = contrasts[contrasts.organ == organ]
        for contrast in CONTRASTS:
            rows = (
                organ_rows[organ_rows.contrast == contrast]
                .set_index("feature")
                .reindex(order)
            )
            effect = rows.effect.to_numpy(float)
            low = rows.ci95_low.to_numpy(float)
            high = rows.ci95_high.to_numpy(float)
            inset.errorbar(
                effect,
                y + offsets[contrast],
                xerr=np.vstack([effect - low, high - effect]),
                fmt="o",
                color=colors[contrast],
                ecolor=colors[contrast],
                markersize=3.7,
                elinewidth=1.0,
                capsize=1.8,
                markeredgecolor="#222222",
                markeredgewidth=0.3,
                zorder=3,
            )
            permutation_p = rows.freedman_lane_permutation_p_value.to_numpy(float)
            for feature_index, p_value in enumerate(permutation_p):
                stars = significance_stars(float(p_value))
                if not stars:
                    continue
                positive = effect[feature_index] >= 0
                x_position = (
                    high[feature_index] + limit * 0.035
                    if positive
                    else low[feature_index] - limit * 0.035
                )
                star = inset.text(
                    x_position,
                    y[feature_index] + offsets[contrast],
                    stars,
                    ha="left" if positive else "right",
                    va="center",
                    color=colors[contrast],
                    fontsize=STAR_FONTSIZE_FOREST,
                    fontweight="bold",
                    zorder=5,
                )
                star.set_path_effects(
                    [path_effects.withStroke(linewidth=2.2, foreground="white")]
                )
        inset.axvline(0, color="#222222", linewidth=0.8, linestyle="--", zorder=1)
        inset.set_xlim(-limit, limit)
        inset.set_ylim(-0.6, len(order) - 0.4)
        inset.set_title(organ, fontsize=7, fontweight="bold", pad=2)
        inset.set_yticks(y)
        if index == 0:
            inset.set_yticklabels([CV_FEATURES[name] for name in order], fontsize=5.8)
        else:
            inset.set_yticklabels([])
        inset.tick_params(axis="x", labelsize=5.5)
        inset.grid(axis="x", color="#e1e1e1", linewidth=0.55)
        inset.grid(axis="y", visible=False)
        if index == 1:
            inset.set_xlabel("standardized effect vs Water", fontsize=6)
    panel_side_legend(
        ax,
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color=TREATMENT_COLORS["Sucrose"],
                label="Sucrose",
                markersize=4.5,
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color=TREATMENT_COLORS["Allulose"],
                label="Allulose",
                markersize=4.5,
            ),
        ],
        ncol=2,
        fontsize=6.3,
    )


def render_all(
    input_root: Path,
    output_root: Path,
    slide: pd.DataFrame,
    composition_stats: pd.DataFrame,
    significant_examples: list[dict[str, object]],
    cell_bar_values: pd.DataFrame,
    scores: pd.DataFrame,
    deviations: pd.DataFrame,
    deviation_anova_values: pd.DataFrame,
    cv_stats: pd.DataFrame,
    dpi: int,
) -> tuple[Path, Path, str, dict[str, object], Path]:
    style()
    panels_dir = output_root / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    sample_id, representative_selection = representative_sample(input_root, slide)
    histology_overlays = organ_histology_overlays(input_root, sample_id)
    # Panels are lettered in reading order: A/B share the top row, C is the
    # full-width histology strip, D/E share the third row and F is the bottom
    # row. The letters used to run A, G, D, F, H, E because two early panels
    # were merged away, which made the figure unreadable as a sequence.
    plotters = {
        "A": lambda ax: plot_combined_panel_a(ax, input_root, sample_id),
        "B": lambda ax: plot_deviation(ax, deviations, deviation_anova_values),
        "C": lambda ax: plot_histoplus(
            ax, histology_overlays, significant_examples, sample_id
        ),
        "D": lambda ax: plot_pca(ax, scores),
        "E": lambda ax: plot_cv_forest(ax, cv_stats),
        "F": lambda ax: plot_cell_effect_barplots(ax, composition_stats),
    }
    # Isolated panels reuse the master's own panel rectangle, so a reader sees
    # exactly the composition the multipanel shows, and the declashing checked
    # on the master also holds for every standalone export.
    panel_sizes = {
        "A": (PANEL_WIDTH_LEFT, PANEL_HEIGHT_TOP),
        "B": (PANEL_WIDTH_RIGHT, PANEL_HEIGHT_TOP),
        "C": (PANEL_WIDTH_FULL, PANEL_HEIGHT_WIDE),
        "D": (PANEL_WIDTH_LEFT, PANEL_HEIGHT_MIDDLE),
        "E": (PANEL_WIDTH_RIGHT, PANEL_HEIGHT_MIDDLE),
        "F": (PANEL_WIDTH_FULL, PANEL_HEIGHT_WIDE),
    }
    translations: dict[str, str] = {}
    for letter, plotter in plotters.items():
        figure = plt.figure(figsize=panel_sizes[letter])
        axis = figure.add_axes([0.0, 0.0, 1.0, 1.0])
        plotter(axis)
        stem = panels_dir / f"Figure_S3_Panel_{letter}"
        figure.savefig(
            stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white"
        )
        figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        translations.update(
            translate_panel_to_spanish(
                figure,
                extra={**SPANISH_TERMS, **SPANISH_TERMS_BY_PANEL.get(letter, {})},
            )
        )
        spanish = panels_dir / f"Figure_S3_Panel_{letter}_spanish"
        figure.savefig(
            spanish.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white"
        )
        figure.savefig(
            spanish.with_suffix(".pdf"), bbox_inches="tight", facecolor="white"
        )
        plt.close(figure)
    translation_path = (
        output_root / "provenance" / "Figure_S3_spanish_translation_receipt.json"
    )
    atomic_json(
        translation_path,
        {
            "status": "complete",
            "scope": "isolated subpanels only; the multipanel master is English-only",
            "helper": "scripts/shared/spanish_panel_text.py",
            "panels": sorted(plotters),
            "translations": dict(sorted(translations.items())),
        },
    )

    figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    grid = figure.add_gridspec(
        4,
        12,
        height_ratios=list(ROW_RATIOS),
        left=FIGURE_LEFT,
        right=FIGURE_RIGHT,
        bottom=FIGURE_BOTTOM,
        top=FIGURE_TOP,
        hspace=FIGURE_HSPACE,
        wspace=FIGURE_WSPACE,
    )
    axes = {
        "A": figure.add_subplot(grid[0, 0:5]),
        "B": figure.add_subplot(grid[0, 5:12]),
        "C": figure.add_subplot(grid[1, 0:12]),
        "D": figure.add_subplot(grid[2, 0:5]),
        "E": figure.add_subplot(grid[2, 5:12]),
        "F": figure.add_subplot(grid[3, 0:12]),
    }
    for letter in ("A", "B", "C", "D", "E", "F"):
        plotters[letter](axes[letter])
    png = output_root / "Figure_S3.png"
    pdf = output_root / "Figure_S3.pdf"
    figure.savefig(png, dpi=dpi, facecolor="white")
    figure.savefig(pdf, facecolor="white")
    plt.close(figure)
    return png, pdf, sample_id, representative_selection, translation_path


def write_legend(
    root: Path,
    representative: str,
    example_manifest: pd.DataFrame,
    permutations: int,
) -> list[Path]:
    """Write the bilingual figure legend and the bilingual per-panel legends."""
    if example_manifest.empty:
        example_sentence = (
            "No phenotype-composition contrast survived global BH q<0.05, so no "
            "cell-example inset was displayed."
        )
        spanish_phenotype = translate_text(
            str(example.phenotype_display).lower(),
            extra=SPANISH_TERMS_BY_PANEL["C"],
        )
        spanish_organ = translate_text(str(example.organ), extra=SPANISH_TERMS)
        spanish_contrast = translate_text(str(example.contrast), extra=SPANISH_TERMS)
        example_sentence_spanish = (
            "Ningún contraste de composición fenotípica superó la corrección "
            "global BH q<0,05, por lo que no se muestra ningún recuadro de "
            "células de ejemplo."
        )
    else:
        example = example_manifest.iloc[0]
        example_sentence = (
            "The paired mini-insets show the globally BH-significant "
            f"{example.phenotype_display} result in {example.organ}: "
            f"{example.contrast} ({example.multiplicative_relative_abundance:.2f}-fold "
            f"relative abundance; q={example.q_value_global:.3f}). Yellow contours are "
            "the true target-class nuclear polygons. Animals were selected closest to "
            "the treatment-median displayed-class fraction, then a positive tile was "
            "selected closest to that animal's fraction."
        )
        spanish_phenotype = translate_text(
            str(example.phenotype_display).lower(),
            extra=SPANISH_TERMS_BY_PANEL["C"],
        )
        spanish_organ = translate_text(str(example.organ), extra=SPANISH_TERMS)
        spanish_contrast = translate_text(str(example.contrast), extra=SPANISH_TERMS)
        example_sentence_spanish = (
            "Los recuadros emparejados muestran el resultado significativo tras "
            f"corrección global BH de {spanish_phenotype} en {spanish_organ}: "
            f"{spanish_contrast} (abundancia relativa "
            f"{example.multiplicative_relative_abundance:.2f}×; q={example.q_value_global:.3f}). "
            "Los contornos amarillos son los polígonos nucleares reales de la clase "
            "diana. Se eligió el animal más próximo a la mediana de la fracción de la "
            "clase mostrada en su tratamiento y, dentro de él, la tesela positiva más "
            "próxima a esa fracción."
        )

    legend_dir = root / "legends"
    legend_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    def write(name: str, body: str) -> None:
        path = legend_dir / name
        path.write_text(body, encoding="utf-8")
        written.append(path)

    text = f"""FIGURE S3 — MULTI-ORGAN H&E COMPUTER VISION

(A) Compact fused panel containing the reproducible workflow and representative within-slide segmentation ({representative}). Liver, Kidney, and Spleen are labeled directly in their verified top-to-bottom order. Each of 24 animal-level Motic whole-slide images contains all three organs. Scanner-native L0 and L2 images were exported with physical-resolution and SHA-256 receipts. Within-slide organ masks defined spatially balanced sampling, and exactly 24 direct-pixel-QC tiles per organ were analyzed at 0.5 µm/px. The displayed Water animal is the one closest to the median whole-slide tissue fraction among Water animals. (B) Cross-validated Mahalanobis deviation from the Water composition, with organ-wise one-way ANOVA p values. (C) Full-width horizontal strip of matched median-density Kidney, Liver, and Spleen H&E tiles from Water animal {representative}, followed by paired Water/Allulose significant-cell examples; the tumor-specific model class is omitted from downstream composition analysis. {example_sentence} (D) Organ-specific CLR-PCA colored by treatment only; each row reports its own PC1 and PC2 variance fractions. (E) Model-independent classical H&E standardized effects and HC3 95% confidence intervals shown as forest plots. Panel E asterisks encode unadjusted Freedman–Lane permutation p values; no classical feature survives global BH q<0.05. (F, bottom) Full-width organ-specific barplots of cohort-adjusted centered-log-ratio cell-type effects for Sucrose and Allulose relative to Water. Whiskers are HC3 95% confidence intervals; high-contrast asterisks encode global BH q-value thresholds.

Statistics: the animal is the biological unit (Water n=7, Sucrose n=9, Allulose n=8). The tumor-specific model class was excluded from classification/composition displays but remains in raw model-output tables. Within each organ, remaining phenotype components were eligible for CLR, PCA, and Water-reference deviation only when detected in at least 18/24 animals and represented by at least 24 total objects; filtered raw counts remain available and contrasts are marked not tested. A 0.5 pseudocount was applied only to eligible components. OLS models included treatment and a prespecified source-cohort covariate; HC3 confidence intervals and {permutations:,} Freedman–Lane residual permutations were computed. For each contrast, the reduced model retained the other treatment contrast and cohort; the observed and permuted statistics were both studentized using HC3 standard errors. Global Benjamini–Hochberg correction covers displayed contrasts within each endpoint family. Panel F global-BH q stars use *<0.05, **<0.01, ***<0.001, and ****<0.0001. Panel E uses the same thresholds for unadjusted permutation p values and explicitly distinguishes them from FDR significance. Panel B additionally reports requested one-way treatment ANOVA p values. The panels themselves print no footnote: the asterisk thresholds, the bar and whisker definitions and the FDR caveat are defined here and nowhere else, so the graphics stay uncluttered and the enlarged asterisks remain readable at print size. Sampled counts were not extrapolated to whole organs, and failures were never encoded as zero.

Interpretation boundary: HistoPLUS was trained on human tumor H&E and has not been validated for normal mouse kidney, liver, or spleen. Its labels are therefore exploratory model-predicted nuclear phenotypes (the panel uses “-like” display labels), not validated mouse cell identities. Associations are described as sugar-associated deviations from Water, not causal effects. Classical H&E features provide a model-independent sensitivity analysis. AGUA_Sin_cerebro (m26-016) is included as Water; “brain not collected” is retained only as a provenance note.

Language: this multipanel master is English-only. The isolated subpanels are published in English and Spanish, and every legend in this directory has a Spanish counterpart.
"""
    write("Figure_S3_LEGEND.txt", text)

    text_spanish = f"""FIGURA S3 — VISIÓN POR COMPUTADOR SOBRE H&E MULTIORGÁNICA

(A) Panel compacto que fusiona el flujo de trabajo reproducible y la segmentación representativa dentro del portaobjetos ({representative}). Hígado, riñón y bazo se rotulan directamente en su orden verificado de arriba abajo. Cada una de las 24 imágenes Motic de portaobjetos completo, una por animal, contiene los tres órganos. Los niveles L0 y L2 nativos del escáner se exportaron con recibos de resolución física y SHA-256. Las máscaras de órganos dentro de cada portaobjetos definieron un muestreo equilibrado espacialmente, y se analizaron exactamente 24 teselas por órgano, con control de calidad sobre los píxeles reales, a 0,5 µm/px. El animal de Agua mostrado es el más próximo a la mediana de la fracción de tejido del portaobjetos completo entre los animales de Agua. (B) Desviación de Mahalanobis con validación cruzada respecto a la composición de Agua, con valores p de ANOVA de una vía por órgano. (C) Banda horizontal de ancho completo con teselas de H&E de riñón, hígado y bazo de densidad mediana del animal de Agua {representative}, seguida de ejemplos emparejados Agua/Alulosa de células significativas; la clase específica de tumor se omite del análisis de composición. {example_sentence_spanish} (D) ACP sobre CLR por órgano, coloreado solo por tratamiento; cada fila indica sus propias fracciones de varianza de PC1 y PC2. (E) Efectos estandarizados clásicos de H&E, independientes del modelo, con intervalos de confianza HC3 del 95 % en gráficos de bosque. Los asteriscos del panel E corresponden a valores p de permutación de Freedman–Lane sin ajustar; ninguna variable clásica supera la corrección global BH q<0,05. (F, abajo) Barras de ancho completo por órgano con los efectos de tipo celular en razón logarítmica centrada, ajustados por cohorte, para Sacarosa y Alulosa frente a Agua. Los bigotes son intervalos de confianza HC3 del 95 % y los asteriscos de alto contraste codifican los umbrales de q tras corrección global BH.

Estadística: la unidad biológica es el animal (Agua n=7, Sacarosa n=9, Alulosa n=8). La clase específica de tumor se excluyó de las representaciones de clasificación y composición, pero permanece en las tablas de salida cruda del modelo. Dentro de cada órgano, los componentes fenotípicos restantes solo fueron elegibles para CLR, ACP y desviación respecto a la referencia de Agua cuando se detectaron en al menos 18 de 24 animales y con al menos 24 objetos en total; los recuentos crudos filtrados siguen disponibles y sus contrastes se marcan como no evaluados. Se aplicó un pseudoconteo de 0,5 solo a los componentes elegibles. Los modelos de mínimos cuadrados incluyeron el tratamiento y una covariable prespecificada de cohorte de origen; se calcularon intervalos de confianza HC3 y {permutations:,} permutaciones de residuos de Freedman–Lane. Para cada contraste, el modelo reducido conservó el otro contraste de tratamiento y la cohorte; tanto el estadístico observado como el de cada permutación utilizaron errores estándar HC3. La corrección global de Benjamini–Hochberg cubre los contrastes mostrados dentro de cada familia de variables. Los asteriscos de q global del panel F usan *<0,05, **<0,01, ***<0,001 y ****<0,0001. El panel E emplea los mismos umbrales para valores p de permutación sin ajustar y los distingue explícitamente de la significación tras control de la tasa de falsos descubrimientos. El panel B añade los valores p de ANOVA de una vía solicitados. Los paneles no llevan nota al pie: los umbrales de los asteriscos, la definición de barras y bigotes y la advertencia sobre la tasa de falsos descubrimientos se definen aquí y en ningún otro sitio, de modo que las figuras quedan despejadas y los asteriscos, ahora del doble de tamaño, se leen bien al tamaño de impresión. Los recuentos muestreados no se extrapolaron a órganos completos y ningún fallo se codificó como cero.

Límite de interpretación: HistoPLUS se entrenó con H&E de tumores humanos y no está validado para riñón, hígado ni bazo de ratón normales. Sus etiquetas son, por tanto, fenotipos nucleares exploratorios predichos por el modelo (los paneles usan etiquetas de “tipo”), no identidades celulares validadas en ratón. Las asociaciones se describen como desviaciones asociadas al azúcar respecto a Agua, no como efectos causales. Las variables clásicas de H&E aportan un análisis de sensibilidad independiente del modelo. AGUA_Sin_cerebro (m26-016) se incluye como Agua; “cerebro no recolectado” se conserva solo como nota de procedencia.

Idioma: este máster multipanel es solo en inglés. Los subpaneles aislados se publican en inglés y español, y cada leyenda de este directorio tiene su equivalente en español.
"""
    write("Figure_S3_LEGEND_spanish.txt", text_spanish)

    panels_english = {
        "A": (
            "Reproducible workflow for the 24 Motic H&E whole-slide images and a "
            f"representative within-slide organ segmentation of Water animal {representative}. "
            "Liver, Kidney and Spleen are labeled in their verified top-to-bottom order; "
            "the map is cropped to the tissue so the three organs, not the slide "
            "background, fill the panel."
        ),
        "B": (
            "Cross-validated centered-log-ratio Mahalanobis distance of each animal from "
            "the Water composition, by organ. Boxes are the interquartile range, points "
            "are individual animals, and Water controls are evaluated leave-one-out. "
            "One-way treatment ANOVA p values are printed above each organ."
        ),
        "C": (
            f"Matched median-density Kidney, Liver and Spleen H&E tiles from Water animal {representative}, "
            "with 100 µm scale bars and per-tile retained instance counts, followed by the "
            "paired Water/sugar example tiles for the strongest globally BH-significant "
            f"phenotype contrast with 50 µm scale bars. {example_sentence}"
        ),
        "D": (
            "Organ-specific principal-component analysis of the centered-log-ratio "
            "composition, colored by treatment only. Each organ is an independent "
            "ordination, so each row reports its own PC1 and PC2 variance fractions."
        ),
        "E": (
            "Model-independent classical H&E and nuclear features as forest plots of the "
            "standardized effect versus Water with HC3 95% confidence intervals. Points "
            "are the estimate and bars the HC3 95% confidence interval. Asterisks mark "
            "unadjusted Freedman–Lane permutation p values, *<0.05, **<0.01, ***<0.001 "
            "and ****<0.0001; they are not FDR-corrected, and no feature in this panel "
            "survives global BH q<0.05."
        ),
        "F": (
            "Cohort-adjusted centered-log-ratio effects of Sucrose and Allulose relative "
            "to Water for every eligible model-predicted cell type, within each organ. "
            "Bars are the cohort-adjusted CLR effect and whiskers are HC3 95% confidence "
            "intervals. Asterisks encode global Benjamini–Hochberg q values, *<0.05, "
            "**<0.01, ***<0.001 and ****<0.0001."
        ),
    }
    panels_spanish = {
        "A": (
            "Flujo de trabajo reproducible para las 24 imágenes Motic de H&E de "
            "portaobjetos completo y segmentación representativa de los órganos dentro "
            f"del portaobjetos del animal de Agua {representative}. Hígado, riñón y bazo se "
            "rotulan en su orden verificado de arriba abajo; el mapa se recorta al tejido "
            "para que los tres órganos, y no el fondo del portaobjetos, ocupen el panel."
        ),
        "B": (
            "Distancia de Mahalanobis con validación cruzada, sobre la razón logarítmica "
            "centrada, de cada animal respecto a la composición de Agua, por órgano. Las "
            "cajas son el rango intercuartílico, los puntos son animales individuales y "
            "los controles de Agua se evalúan dejando uno fuera. Sobre cada órgano se "
            "indica el valor p del ANOVA de una vía por tratamiento."
        ),
        "C": (
            "Teselas de H&E de densidad mediana de riñón, hígado y bazo del animal de Agua "
            f"{representative}, con barras de escala de 100 µm y el recuento de instancias "
            "retenidas por tesela, seguidas de las teselas de ejemplo emparejadas "
            "Agua/azúcar del contraste fenotípico más fuerte significativo tras corrección "
            f"global BH, con barras de escala de 50 µm. {example_sentence_spanish}"
        ),
        "D": (
            "Análisis de componentes principales por órgano de la composición en razón "
            "logarítmica centrada, coloreado solo por tratamiento. Cada órgano es una "
            "ordenación independiente, por lo que cada fila indica sus propias fracciones "
            "de varianza de PC1 y PC2."
        ),
        "E": (
            "Variables clásicas de H&E y nucleares, independientes del modelo, como "
            "gráficos de bosque del efecto estandarizado frente a Agua con intervalos de "
            "confianza HC3 del 95 %. Los puntos son la estimación y las barras el "
            "intervalo de confianza HC3 del 95 %. Los asteriscos indican valores p de "
            "permutación de Freedman–Lane sin ajustar, *<0,05, **<0,01, ***<0,001 y "
            "****<0,0001; no están corregidos por tasa de falsos descubrimientos y "
            "ninguna variable de este panel supera la corrección global BH q<0,05."
        ),
        "F": (
            "Efectos en razón logarítmica centrada, ajustados por cohorte, de Sacarosa y "
            "Alulosa respecto a Agua para cada tipo celular predicho por el modelo que "
            "resultó elegible, dentro de cada órgano. Las barras son el efecto CLR "
            "ajustado por cohorte y los bigotes el intervalo de confianza HC3 del 95 %. "
            "Los asteriscos codifican los valores q tras corrección global de "
            "Benjamini–Hochberg, *<0,05, **<0,01, ***<0,001 y ****<0,0001."
        ),
    }
    for letter in ("A", "B", "C", "D", "E", "F"):
        write(
            f"Figure_S3_Panel_{letter}_LEGEND.txt",
            f"Figure S3, panel {letter}. {panels_english[letter]}\n",
        )
        write(
            f"Figure_S3_Panel_{letter}_LEGEND_spanish.txt",
            f"Figura S3, panel {letter}. {panels_spanish[letter]}\n",
        )
    return written


def main() -> int:
    args = parse_args()
    if args.permutations < 999:
        raise FigS3Error("Use at least 999 permutations")
    input_root = args.figs3_root.resolve()
    if args.output_root is None:
        output_root = input_root
    else:
        output_root = args.output_root.expanduser().resolve()
        if output_root.exists() or output_root.is_symlink():
            raise FigS3Error(f"Fresh --output-root already exists: {output_root}")
        output_root.mkdir(parents=True, mode=0o700)
    slide, phenotype, tile_cv = read_complete_inputs(input_root)
    counts, clr = composition_values(phenotype)
    cell_bar_values = displayed_cell_bar_values(phenotype)
    composition_stats = composition_contrasts(clr, args.permutations, args.seed)
    significant_examples, example_manifest = significant_cell_examples(
        input_root, phenotype, composition_stats
    )
    deviations = deviation_values(counts)
    deviation_stats = deviation_contrasts(deviations, args.permutations, args.seed)
    deviation_anova_values = deviation_anova(deviations)
    scores, variance = pca_scores(counts)
    classical_values = classical_animal_values(tile_cv)
    classical_stats = classical_contrasts(
        classical_values, args.permutations, args.seed
    )
    source = output_root / "source_data"
    source.mkdir(parents=True, exist_ok=True)
    save_table(counts, source / "Figure_S3_animal_phenotype_counts.csv")
    save_table(clr, source / "Figure_S3_animal_composition_clr.csv")
    save_table(composition_stats, source / "Figure_S3_composition_contrasts.csv")
    save_table(
        cell_bar_values,
        source / "Figure_S3_displayed_cell_barplot_values.csv",
    )
    save_table(
        example_manifest,
        source / "Figure_S3_significant_cell_example_manifest.csv",
    )
    save_table(deviations, source / "Figure_S3_Water_reference_deviation_values.csv")
    save_table(
        deviation_stats, source / "Figure_S3_Water_reference_deviation_contrasts.csv"
    )
    save_table(
        deviation_anova_values, source / "Figure_S3_Water_reference_deviation_ANOVA.csv"
    )
    save_table(scores, source / "Figure_S3_composition_pca_scores.csv")
    save_table(classical_values, source / "Figure_S3_classical_cv_animal_values.csv")
    save_table(classical_stats, source / "Figure_S3_classical_cv_contrasts.csv")
    png, pdf, representative, representative_selection, translation_receipt = (
        render_all(
            input_root,
            output_root,
            slide,
            composition_stats,
            significant_examples,
            cell_bar_values,
            scores,
            deviations,
            deviation_anova_values,
            classical_stats,
            args.dpi,
        )
    )
    legends = write_legend(output_root, representative, example_manifest, args.permutations)
    output_files = [
        png,
        pdf,
        translation_receipt,
        *legends,
        source / "Figure_S3_animal_phenotype_counts.csv",
        source / "Figure_S3_animal_composition_clr.csv",
        source / "Figure_S3_composition_contrasts.csv",
        source / "Figure_S3_displayed_cell_barplot_values.csv",
        source / "Figure_S3_significant_cell_example_manifest.csv",
        source / "Figure_S3_Water_reference_deviation_values.csv",
        source / "Figure_S3_Water_reference_deviation_contrasts.csv",
        source / "Figure_S3_Water_reference_deviation_ANOVA.csv",
        source / "Figure_S3_composition_pca_scores.csv",
        source / "Figure_S3_classical_cv_animal_values.csv",
        source / "Figure_S3_classical_cv_contrasts.csv",
    ]
    significant_composition = int((composition_stats.q_value_global < 0.05).sum())
    significant_deviation = int((deviation_stats.q_value_global < 0.05).sum())
    significant_classical = int((classical_stats.q_value_global < 0.05).sum())
    receipt = {
        "status": "complete",
        "analysis_version": ANALYSIS_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "biological_unit": "animal",
        "animal_count": 24,
        "treatment_counts": slide.treatment.value_counts().sort_index().to_dict(),
        "agual_sin_cerebro_handling": "m26-016 included in Water; note retained only as provenance",
        "organs_per_wsi": list(ORGANS),
        "image_registration_performed": False,
        "tiles_per_animal_organ": 24,
        "composition_pseudocount": 0.5,
        "excluded_model_classes_from_composition_and_figure": ["Cancer cell"],
        "panel_letters_in_reading_order": ["A", "B", "C", "D", "E", "F"],
        "master_figure_size_inches": [FIGURE_WIDTH, FIGURE_HEIGHT],
        "language_policy": {
            "multipanel_master": "English only",
            "isolated_panels": "English and Spanish",
            "legends": "English and Spanish, figure-level and per panel",
        },
        "significant_cell_example_insets": {
            "display_rule": "global BH q<0.05 only; strongest contrast displayed",
            "displayed_contrasts": (
                int(example_manifest.contrast.nunique())
                if not example_manifest.empty
                else 0
            ),
            "manifest": example_manifest.to_dict(orient="records"),
        },
        "cell_barplots": {
            "rows": int(len(cell_bar_values)),
            "displayed_classes": len(CLASS_ORDER),
            "fraction_denominator": "all 13 retained displayed model outputs within animal-organ",
            "count_transform": "log10(displayed class objects in 24 sampled tiles + 1)",
            "inferential_role": "descriptive; inference is the CLR contrast",
        },
        "panel_F_effect_barplots": {
            "source": "Figure_S3_composition_contrasts.csv",
            "panels": "Kidney, Liver, Spleen",
            "bars": "cohort-adjusted CLR treatment effect vs Water",
            "whiskers": "HC3 95% confidence interval",
            "stars": "global BH q-value thresholds",
            "heatmap_displayed": False,
        },
        "composition_prevalence_filter": {
            "minimum_detected_animals_within_organ": MIN_DETECTED_ANIMALS,
            "minimum_total_class_count_within_organ": MIN_TOTAL_CLASS_COUNT,
            "eligible_classes_by_organ": {
                organ: sorted(
                    clr.loc[
                        (clr.organ == organ) & clr.analysis_eligible, "model_class_raw"
                    ]
                    .unique()
                    .tolist()
                )
                for organ in ORGANS
            },
            "filtered_classes_retained_as_raw_counts_and_marked_not_tested": True,
        },
        "permutations": args.permutations,
        "permutation_method": "Freedman-Lane residual permutation; contrast-specific nuisance model; matched observed/permuted HC3 t statistics",
        "permutation_seed": args.seed,
        "permutation_nuisance_columns": "intercept, other treatment contrast, source cohort",
        "permutation_statistic": "HC3-studentized treatment coefficient",
        "permutation_p_estimator": "(1 + extreme permutations) / (1 + permutations)",
        "multiple_testing": "Benjamini-Hochberg globally within endpoint family and within organ",
        "pca_explained_variance": variance,
        "Water_reference_deviation_one_way_ANOVA": deviation_anova_values.to_dict(
            orient="records"
        ),
        "significant_global_q_lt_0_05": {
            "composition_contrasts": significant_composition,
            "Water_reference_deviation_contrasts": significant_deviation,
            "classical_cv_contrasts": significant_classical,
        },
        "representative_sample_id": representative,
        "representative_sample_selection": representative_selection,
        "histoPLUS_cross_species_validated": False,
        "interpretation": "exploratory model-predicted nuclear phenotypes; sugar-associated, not causal",
        "sampled_counts_extrapolated_to_whole_organ": False,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": mpl.__version__,
            "seaborn": sns.__version__,
            "statsmodels": sm.__version__,
            "sklearn": sklearn.__version__,
            "shapely": shapely.__version__,
        },
        "outputs": [
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in output_files
        ],
    }
    atomic_json(output_root / "provenance" / "Figure_S3_analysis_receipt.json", receipt)
    print(
        json.dumps(
            {
                "status": "complete",
                "figure_png": str(png),
                "figure_pdf": str(pdf),
                "representative_sample_id": representative,
                "significant_global_q_lt_0_05": receipt["significant_global_q_lt_0_05"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
