Figure S8: pilot-based Monte Carlo sample-size estimates

Regenerate from the small verified inputs, into a fresh output directory:
  python FigS8/01_make_figure_s8_power.py --output-dir /tmp/figure_s8_new --dpi 600

The complete master is English. Isolated panels and legends are supplied in
English and Spanish. The renderer verifies provenance/input_manifest.json
before reading source data. Shared plotting code is in
scripts/shared/render_method_figures.py. No model retraining or new Monte
Carlo simulations occur during rendering. Original images, masks, loss values
and pilot-derived power values are not modified.
