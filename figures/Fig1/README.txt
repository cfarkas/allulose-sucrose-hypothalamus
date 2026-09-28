FIGURE 1 — SINGLE-BOTTLE FLUID CONSUMPTION AND BODY WEIGHT

Figure 1 compares water, sucrose and allulose during single-bottle exposure.
The independent experimental unit is the cage. Weight measurements include
24 mice in 12 cages: Water/Sucrose/Allulose 5/9/10 mice and 3/5/4 cages.

Panels
  A  Cumulative bottle-volume removal on days 1, 3 and 6.
  B  Day-6 bottle-volume endpoint, with one point per cage.
  C  Percentage body-weight change from each animal's baseline.
  D  Day-6 cage-mean percentage body-weight change.

Animal weight changes are averaged within each cage, and cages receive equal
weight. Day-6 omnibus tests use ordinary one-way ANOVA; exact cage-label tests
provide sensitivity analyses. Pairwise exact comparisons use Holm adjustment
within outcome. Individual mouse trajectories are descriptive. Bottle-volume
removal is measured in mL/cage, without an independent correction for losses.

The treatment assignment of FR5-4 is Allulose through the day-6 endpoint;
its transfer to Water occurred after that measurement. E10/FR6-2 is Water.
Nonpositive weights are invalid measurements. Source labels and audit tables
document these rules.

Files
  Figure_1_behavior_experiment_1.pdf/.png  Complete English figure, 600-dpi PNG.
  panels/       Individual A–D panels in English and Spanish, PNG and PDF.
  legends/      Figure and panel captions in both languages.
  source_data/  Plotted values, tests, intervals and sample-size simulations.
  raw_data/     Consumption and weight tables with their input manifest.
  provenance/   Source hashes and validation records.

Reproduce from Paper/ with Python from the paper environment, using fresh paths:
  python Fig1/01_analyze_behavior_experiment1.py \
    --consumption Fig1/raw_data/consumption.csv \
    --weights Fig1/raw_data/weights.csv --outdir /tmp/figure1_analysis
  python Fig1/02_make_figure_1_behavior.py \
    --analysis-dir /tmp/figure1_analysis \
    --output-dir /tmp/figure1_render --dpi 600

The pilot-based planning simulations used by Figure S8 can be recomputed with:
  python Fig1/03_estimate_sample_size_monte_carlo.py --replicates 500000
These simulations describe assumptions for future cage-level studies; they do
not estimate achieved power for the neuronal or glial experiments.
