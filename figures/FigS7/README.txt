FIGURE S7 — HUMAN-SUPERVISED CELLPOSE TRAINING

Panels
  A  Segmentation, human mask correction and model-training workflow.
  B  DAPI, c-FOS and POMC training-image crops and reviewed mask overlays.
  C  Saved training-loss curves for five checkpoints.

The 512 × 512-pixel crops and masks are bound to their source images by
SHA-256. Signal/overlay pairs use identical display contrast. POMC overlays
apply a local DAPI-area size filter at rendering time to the ten stored fields.
The saved model losses correspond to the training annotations used for each
checkpoint, without that display-time filtering. They describe training fit,
not independent test performance. Per-object display decisions are in
source_data/pomc_training_size_qc.csv.

Files
  Figure_S7.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/            Individual A–C panels in English and Spanish, PNG and PDF.
  legends/           Figure and panel captions in both languages.
  raw_data/          Verified training crops, masks and loss values.
  provenance/        Input manifest and crop-source records.

Reproduce from Paper/ using the figure environment and a fresh output path:
  python FigS7/01_make_figure_s7_training.py \
    --output-dir /tmp/figure_s7_render --dpi 600
From the GitHub root, prepend figures/ to the script path.
The renderer verifies provenance/input_manifest.json and reads the supplied
inputs. Model training is separate from figure rendering. Shared plotting code
is in scripts/shared/render_method_figures.py.
