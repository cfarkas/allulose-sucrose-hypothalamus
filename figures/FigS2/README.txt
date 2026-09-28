FIGURE S2 — FEMALE SINGLE-BOTTLE FLUID CONSUMPTION AND BODY WEIGHT

This figure presents the female animals from the single-bottle experiment.
The independent experimental unit is the cage; animal weight changes are
averaged within cages and each cage receives equal weight.

Panels
  A  Cumulative bottle-volume removal across the observation days.
  B  Day-6 bottle-volume endpoint.
  C  Percentage body-weight change from each animal's baseline.
  D  Day-6 cage-mean percentage body-weight change.

Day-6 omnibus tests use ordinary one-way ANOVA. Exhaustive cage-label tests
provide sensitivity analyses, and exact pairwise comparisons use Holm
adjustment within sex and endpoint. Small cage counts limit precision.
Sex-stratified summaries do not establish a treatment-by-sex interaction.
The Water group has one female cage, which particularly limits variance estimation.

Files
  Figure_S2_single_bottle_female.pdf/.png  Complete English figure.
  panels/       Individual A–D panels in English and Spanish, PNG and PDF.
  legends/      Figure and panel captions in both languages.
  source_data/  Plotted values, tests, intervals and source audits.
  raw_data/     Consumption and weight tables with their input manifest.
  provenance/   Input hashes and figure validation records.
All publication PNGs are exported at 600 dpi.

Reproduce from Paper/ with the paper environment and fresh output paths:
  python FigS2/01_analyze_single_bottle_female.py \
    --consumption FigS2/raw_data/consumption.csv \
    --weights FigS2/raw_data/weights.csv --outdir /tmp/FigS2_analysis
  python FigS2/02_make_figure_s2_female.py \
    --analysis-dir /tmp/FigS2_analysis \
    --output-dir /tmp/FigS2_render --dpi 600
