FIGURE S2 — FEMALE SINGLE-BOTTLE BEHAVIOR

Status: publication figure ready and reproducible from figure-local raw_data.
The combined master is English only. Each A-D panel is supplied in English and
as an identically cropped *_spanish PDF/PNG with a Spanish panel legend.

Contents
  Figure_S2_single_bottle_female.pdf/png English 600-dpi master
  panels/                                English and Spanish A-D panel pairs
  legends/                               manuscript legend and panel legends
  source_data/, provenance/              plotted values, statistics, receipts
  raw_data/                               independent input copies + manifest
  scripts/                                analyzer, renderer, translation helper

Reproduce safely into fresh paths (do not render over this directory):

  python scripts/01_analyze_single_bottle_female.py \
    --consumption raw_data/consumption.csv --weights raw_data/weights.csv \
    --outdir /tmp/FigS2_analysis_rebuild
  python scripts/02_make_figure_s2_female.py \
    --analysis-dir /tmp/FigS2_analysis_rebuild \
    --output-dir /tmp/FigS2_figure_rebuild --dpi 600

Inference uses cage as the experimental unit. Day-6 omnibus values are ordinary
one-way ANOVA; exhaustive cage-label omnibus results are sensitivity analyses,
and pairwise comparisons are exact with Holm adjustment. Water has one female
cage, so ANOVA assumptions and confidence intervals are especially limited;
the analysis is exploratory, not confirmatory.
