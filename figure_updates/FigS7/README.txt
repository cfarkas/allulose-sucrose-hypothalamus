Figure S7: human-supervised Cellpose training

Regenerate from the small verified inputs, into a fresh output directory:
  python FigS7/01_make_figure_s7_training.py --output-dir /tmp/figure_s7_new --dpi 600

The complete master is English. Isolated panels and legends are supplied in
English and Spanish. The renderer verifies provenance/input_manifest.json
before reading source data. Shared plotting code is in
scripts/shared/render_method_figures.py. No model retraining or new Monte
Carlo simulations occur during rendering. Original images, masks, loss values
and pilot-derived power values are not modified.

POMC overlays now apply a local DAPI-area size heuristic to all ten stored
training fields. Historical masks remain in the verified NPZ inputs; filtered
overlays are derived at rendering time. The models were not retrained, and
the displayed loss curves still describe their historical training. Per-object
decisions are in source_data/pomc_training_size_qc.csv.
