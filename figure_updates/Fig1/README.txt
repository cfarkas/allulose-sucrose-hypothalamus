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

Authoritative graphic hashes
----------------------------
Combined PDF:
  96a22118214642f1b9779e76592a326fcdae0c966d337eaa2ff7a56fe5a6506a
Combined PNG:
  b90958f3103ea314b71618b4df7d640438695c09e9e745d00427127d871c4d76

All panel hashes are recorded in
provenance/Figure_1_behavior_experiment_1_manifest.json.
