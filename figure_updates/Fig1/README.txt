FIGURE 1 — FINAL REPRODUCIBLE PACKAGE
=====================================

Status
------
This directory contains the current Figure 1 publication package. Superseded
promoted files are retained under legacy/ by the promotion utility.

Contents
--------
- Figure_1_behavior_experiment_1.pdf: English combined vector PDF.
- Figure_1_behavior_experiment_1.png: English combined 600-dpi PNG.
- panels/: exact A-D crops as English PDF/PNG and matching _spanish PDF/PNG.
- legends/: manuscript-style English master legend, English panel legends,
  concise Spanish standalone-panel legends, and the caption text.
- source_data/: plotted values, statistics, confidence intervals, and source
  hashes used by the renderer.
- provenance/: statistical audit, deterministic figure manifest, Spanish
  translation receipt, and renderer notes.
- raw_data/: independent physical copies of consumption.csv and weights.csv,
  plus their copy/hash manifest.
- 01_analyze_behavior_experiment1.py / 02_make_figure_1_behavior.py: current
  analysis/renderer scripts; exact mirrors are in scripts/Fig1/.

Statistical unit
----------------
The independent experimental unit is the cage. The displayed omnibus tests are
ordinary equal-variance one-way ANOVAs on one Day-6 endpoint per cage. Exact
cage-label tests are retained as sensitivity analyses, and the three pairwise
comparisons for each outcome are adjusted using Holm. Mouse trajectories are descriptive;
they are not treated as independent inferential replicates.

Reproduce safely
----------------
From the directory holding README.txt and Fig1/ (Paper or the GitHub
figure_updates directory), use Python from the paper/figure environment:

  python Fig1/01_analyze_behavior_experiment1.py \
    --consumption Fig1/raw_data/consumption.csv \
    --weights Fig1/raw_data/weights.csv --outdir /tmp/figure1_rebuild/analysis
  python Fig1/02_make_figure_1_behavior.py \
    --analysis-dir /tmp/figure1_rebuild/analysis \
    --output-dir /tmp/figure1_rebuild/figure --dpi 600

These commands read only the small raw tables. Use fresh output directories;
do not point the renderer at this mixed final package. The 2026-09-07 terminology
correction labels the exact comparison as "exact, Holm-adjusted". It changes
statistical wording, with the comparisons and p-values retained.

Canonical analysis outputs are also present in analyses/Fig1/results. The
figure-local source_data directory contains the subset directly consumed by or
exported for the publication figure.

Acceptance checks completed on 2026-08-23
-----------------------------------------
- Fresh analysis and rendering from the figure-local raw files: PASS.
- 18 expected graphics (English master plus A-D bilingual pairs): PASS.
- All PDFs are one page; fonts are embedded/subset and no Type 3 fonts: PASS.
- English/Spanish PNG crop dimensions match exactly for every panel: PASS.
- Combined figure is English-only; Spanish is limited to standalone panels:
  PASS.
- PNG master is 6720 x 7320 pixels with 600-dpi metadata: PASS.
- Raw source/destination SHA-256 values match and inodes are distinct: PASS.
- No symlinks: PASS.

Current graphic hashes — 8 September 2026
-----------------------------------------
Combined PDF:
  06ff14fd5b4f73ba6df15161071debdddbc1a57575b8b11559e5ccf38a2e295d
Combined PNG:
  4f382e286730e3b786990b55d99e7ad8a5c9cf9b44fa891037e84dbc273509c5

All panel hashes are recorded in
provenance/Figure_1_behavior_experiment_1_manifest.json.

Author reconciliation and sample-size reconstruction — 8 September 2026
---------------------------------------------------------------------
FR5-4 moved from E2 Allulose to E10 Water for dehydration immediately after
the day-6 measurement. Its earlier measurements retain their original group.
E10/FR6-2, originally labelled Control, is included as Water. The current
weight summary contains 24 mice (Water/Sucrose/Allulose 5/9/10) in 12 cages
(3/5/4). Original source labels and invalid zero entries remain available.

03_estimate_sample_size_monte_carlo.py reconstructs planning scenarios from
the day-6 cage endpoints. It uses 500,000 simulations per candidate size,
three balanced groups, alpha 0.05 and target power 0.80. Group-specific
variances with Welch ANOVA give minima of 5 cages/group for volume removed
and 14 for cage-mean weight change. Common-variance ordinary ANOVA gives
4 and 9, respectively. These pilot-based assumptions do not establish
achieved power of completed cohorts or neuronal/glial endpoints.

Run from Paper/:
  python Fig1/03_estimate_sample_size_monte_carlo.py --replicates 500000

The portable source tables are in source_data/sample_size_pilot/. Outputs
are in analyses/Fig1/results/sample_size_monte_carlo_20260908/. The dated
reconstruction is distinguished from the historical planning simulation.
The verified 500,000-replicate output tables are copied into
source_data/sample_size_monte_carlo/; their execution settings and input/script
hashes are in provenance/SAMPLE_SIZE_MONTE_CARLO_RECEIPT.json.
