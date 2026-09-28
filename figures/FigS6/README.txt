FIGURE S6 — MICROGLIAL MORPHOLOGY CLASSIFICATION BENCHMARK

Figure S6 compares morphometric proposals with frozen DINOv2 image features
and a linear classifier on the same 203 human-reviewed cells from 15 animals.
Balanced accuracy is the mean recall across Ramified, Rod-like, Activated and
Amoeboid reference morphology classes.

Panels
  A  Balanced accuracy with whole-animal bootstrap 95% intervals.
  B  Recall for each of the four morphology classes.
  C  Confusion matrix for morphometric proposals.
  D  Confusion matrix for animal-held-out DINOv2 predictions.

Intervals use 2,000 paired animal resamples and fixed predictions. They do not
include retraining variation. DINOv2 normalization, PCA, classifier fitting
and regularization selection exclude each outer held-out animal. The unsupervised
morphometric comparator clusters the full cell collection without fitting
human labels. This is an internal benchmark on a selected review sample;
shape classification does not establish functional microglial activation.

Files
  Figure_S6.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/            Individual A–D panels in English and Spanish, PNG and PDF.
  legends/           Figure and panel captions in both languages.
  source_data/       Metrics, paired predictions and bootstrap draws.
  raw_data/          Prediction/reference tables and their input manifest.
  provenance/        Source hashes, input partitions and validation records.

Reproduce from Paper/ with the paper or figure environment and a fresh path:
  python FigS6/01_make_figure_s6_microglia_classifier.py \
    --output-dir /tmp/figs6_render --dpi 600
  python scripts/utilities/07_validate_figure_outputs.py \
    --figure FigS6 --root /tmp/figs6_render

The same commands work from the GitHub root with figures/ prepended
to script paths. Only numpy, pandas and matplotlib are needed; the renderer
uses saved predictions and does not require microscopy, a GPU or retraining.
It verifies matching cell identities, reference labels and animal partitions,
then recomputes all displayed metrics. Check inputs and metrics without writing:
  python FigS6/01_make_figure_s6_microglia_classifier.py --validate-only
