FIGURE S3 — KIDNEY, LIVER AND SPLEEN HISTOLOGY

Figure S3 analyzes 24 H&E whole-slide images, each containing kidney, liver
and spleen. The cohort comprises 7 Water, 9 Sucrose and 8 Allulose animals.
AGUA and AGUA_Sin_cerebro both denote Water; the latter records a sample note.

Panels
  A  Analysis workflow and within-slide organ segmentation.
  B  Cross-validated deviation from Water cell-type composition.
  C  Representative three-organ histology and segmented-cell examples.
  D  Organ-specific centered-log-ratio principal-component analysis.
  E  Model-independent H&E and nuclear-feature effects.
  F  Model-predicted nuclear-phenotype composition effects.

The workflow exports scanner-native L0/L2 images, segments the organs, selects
24 spatially balanced L0 tiles per organ, runs HistoPLUS, and extracts classical
H&E and nuclear features. Counts describe the sampled tissue area and are not
extrapolated to whole organs. Representative overview selection uses proximity
to the median whole-slide tissue fraction; the selected Water animal is m26-015.

Animals are the biological replicates. Eligible composition components occur
in at least 18 of 24 animals and contain at least 24 objects. Counts receive
a 0.5 pseudocount before centered-log-ratio transformation. Treatment effects
adjust for source cohort. Coefficient-specific Freedman–Lane nulls retain the
other treatment contrast and source cohort; observed and permuted statistics
both use HC3 standard errors. Tests use 99,999 permutations and seed 1707, with
Benjamini–Hochberg correction within endpoint families. Water-reference
distances use shrinkage Mahalanobis distance and leave-one-out controls.

HistoPLUS labels are exploratory model-predicted nuclear phenotypes, displayed
as '-like', because the model was trained on human tumor H&E. They are not
validated mouse cell identities or diagnostic findings. The tumor-specific
class is excluded from composition displays and inference. Classical image
features provide a model-independent sensitivity analysis.

Files
  Figure_S3.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/            Individual A–F panels in English and Spanish, PNG and PDF.
  legends/           Figure and panel captions in both languages.
  source_data/       Metadata, animal values, ordinations and statistical tests.
  exports/           Full-resolution and working-resolution slide exports.
  organ_segmentation/ and histoplus/  Masks, tiles, predictions and receipts.
  provenance/        Source hashes, method references and validation records.

Reproduce all figures using python3 reproduce.py from the GitHub root.
From a complete Paper/ tree, the Figure S3 workflow is:
  FigS3/06_run_figure_s3.sh
It is resumable and requires the slide inputs, its scientific environments,
CUDA for inference, and the locally authorized HistoPLUS model weights.
FIGS3_SYSTEM_PYTHON, FIGS3_LAZYSLIDE_PYTHON and HISTOPLUS_WEIGHT_FILE select
the interpreters and weight file. Completed, verified inference can be reused.

The numbered stages validate inputs, prepare metadata, export slide levels,
segment organs, perform inference, and compute/render the figure. The source
images are read-only. Recompute the figure from completed inference with:
  python FigS3/05_analyze_make_figure.py --permutations 99999 --seed 1707 --dpi 600

HistoPLUS weights carry their own CC-BY-NC-ND-4.0 terms. Method and software
references are listed in provenance/ and the manuscript.
