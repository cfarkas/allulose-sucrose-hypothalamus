SUPPLEMENTARY FIGURE S6 — MICROGLIAL CLASSIFICATION BENCHMARK
===========================================================

Created 2026-09-07 at the user's request to accompany Figure 5.

The same 203 human-reviewed cells from 15 animals are used for both methods.
The existing morphometric proposals reach 73.9% balanced accuracy; frozen
DINOv2 features plus a linear classifier reach 57.0%. Balanced accuracy is the
unweighted mean recall of the four human-reviewed morphology classes, not
ordinary cell-level accuracy. Ordinary accuracy is 74.4% versus 60.6%.

The English master has four panels:
  A  Balanced accuracy with whole-animal bootstrap 95% intervals.
  B  Recall for Ramified, Rod-like, Activated and Amoeboid reference classes.
  C  Confusion matrix for the existing morphometric proposals.
  D  Confusion matrix for animal-held-out DINOv2 predictions.

The paired baseline-minus-DINOv2 difference is 16.9 percentage points, with a
95% animal-bootstrap interval of 7.8 to 25.1 points. Intervals use 2,000 paired
resamples of animals and fixed predictions; they do not include retraining
variability. The DINOv2 classifier, normalization, PCA and regularization
selection exclude each entire outer held-out animal. The morphometric
comparator is a previously computed unsupervised proposal set, whose clustering
used the full cell collection without fitting human class labels.

This is an internal benchmark on a selected human-review sample, not an
independent biological validation cohort or a population prevalence estimate.
The morphological labels do not independently establish functional activation.
Figure 5's current classifier and publication figures have not been replaced.

Delivered files
---------------
  Figure_S6.pdf / Figure_S6.png    English master; PNG is 600 dpi.
  panels/                         Exact English/Spanish A-D PDF/PNG panel pairs.
  legends/                        Full and per-panel captions in both languages.
  source_data/                    Plotted metrics, paired predictions and draws.
  provenance/                     Source-freeze, build, validation and install receipts.
  raw_data/                       Frozen prediction/reference tables and manifest.
  01_make_figure_s6_microglia_classifier.py
                                  Recomputes every plotted metric and figure.

Local raw_data contains frozen derived inputs for this supplement. Original
microscopy and the human-review campaign remain in Figure 5. The full trained
candidate comparison is at:
  analyses/Fig5/results/transfer_reviewed_20260907_v1/

Rebuild in a fresh output directory, from the Paper root
-------------------------------------------------------
  PAPER_PYTHON=python  # Python from the installed paper or figure environment
  "$PAPER_PYTHON" FigS6/01_make_figure_s6_microglia_classifier.py \
    --output-dir /tmp/figs6_benchmark_fresh --dpi 600
  "$PAPER_PYTHON" scripts/utilities/07_validate_figure_outputs.py \
    --figure FigS6 --root /tmp/figs6_benchmark_fresh

The script is mirrored byte-for-byte under scripts/FigS6/. It uses only numpy,
pandas and matplotlib; no GPU, encoder download or retraining is required to
reproduce this supplement. The same benchmark is included in the eleven-figure
root reproduction workflow.

Check the frozen inputs and independently recomputed metrics without writing:
  "$PAPER_PYTHON" FigS6/01_make_figure_s6_microglia_classifier.py --validate-only

Install a verified rebuild through the normal canonical figure installer:
  "$PAPER_PYTHON" scripts/utilities/04_promote_figure.py \
    --figure FigS6 --from /tmp/figs6_benchmark_fresh

Each input is bound to RAW_DATA_MANIFEST.csv by size and SHA-256. The renderer
requires identical cell/animal/reference-label identities for both methods,
checks every outer and inner animal partition, and independently recomputes
scores, confusion matrices and paired bootstrap intervals. All must match the
frozen Figure 5 experiment before any output directory is created.
