FIGURE S8 — PILOT-BASED MONTE CARLO SAMPLE-SIZE ESTIMATES

Panels
  A  Power curves for day-6 cumulative bottle-volume removal.
  B  Power curves for day-6 percentage body-weight change.

The pilot contains 12 independent cages: 3 Water, 5 Sucrose and 4 Allulose.
Solid curves use group-specific pilot SDs and Welch ANOVA; dashed curves use
a pooled common variance and ordinary ANOVA. Each group size uses 500,000
normally distributed simulations, alpha 0.05 and PCG64 seed 20260908. The
horizontal line marks 80% power. Curves are displayed through 24 cages/group.
These are conditional planning estimates for cage-level behavioral outcomes;
they do not establish achieved power or sample sizes for neuronal/glial assays.

Files
  Figure_S8.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/            Individual A–B panels in English and Spanish, PNG and PDF.
  legends/           Figure and panel captions in both languages.
  source_data/       Saved simulation inputs and plotted values.
  provenance/        Input hashes and rendering records.

Render the supplied curves from Paper/ into a fresh output directory:
  python FigS8/01_make_figure_s8_power.py \
    --output-dir /tmp/figure_s8_render --dpi 600
From the GitHub root, prepend figures/ to the script path.
The renderer verifies provenance/input_manifest.json before reading data.
It does not run new simulations. To recompute the planning simulations:
  python Fig1/03_estimate_sample_size_monte_carlo.py --replicates 500000
The assumptions and full simulation tables are supplied with Figure 1.
