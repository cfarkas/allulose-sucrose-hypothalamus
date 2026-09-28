FIGURE S4 — TWO-BOTTLE CONSUMPTION, FOOD REMOVAL AND BODY WEIGHT

Figure S4 presents a two-bottle assay with eight cages and two mice per cage:
2 Water, 3 Sucrose and 3 Allulose cages. Control cages receive two water bottles;
the other cages receive water and the assigned sugar solution.

Panels
  A  Consumption from the bottle on the assigned-solution side.
  B  Water-bottle consumption.
  C  Assigned-solution fraction of total fluid consumption.
  D  Food removal.
  E  Percentage body-weight change.

Bottle mass changes and food removal are reported in g/cage/day. Summaries
exclude initialization records and use seven daily observations per cage.
The solution fraction in controls describes the two water-bottle positions.
Weights are expressed relative to the pre-exposure baseline. Cage is the
independent unit. Displayed endpoints use ordinary one-way ANOVA and raw exact
two-sided Mann–Whitney comparisons; those pairwise values are unadjusted.
Small cage counts and unmeasured fluid losses constrain interpretation.

Files
  Figure_S4_behavior_preference.pdf/.png  Complete English figure.
  panels/       Individual A–E panels in English and Spanish, PNG and PDF.
  legends/      Figure and panel captions in both languages.
  source_data/  Plotted values, endpoint tests and source-audit tables.
  raw_data/UCSC_MICE.xlsx  Input workbook, supplied through Zenodo.
  provenance/   Input and output hashes and validation records.
All publication PNGs are exported at 600 dpi.

Reproduce from Paper/ with the paper environment and fresh output paths:
  python FigS4/01_analyze_behavior_experiment2.py \
    --workbook FigS4/raw_data/UCSC_MICE.xlsx --output-dir /tmp/FigS4_analysis
  python FigS4/02_make_figure_s4_behavior_preference.py \
    --workbook FigS4/raw_data/UCSC_MICE.xlsx \
    --analysis-dir /tmp/FigS4_analysis --output-dir /tmp/FigS4_render --dpi 600

RAW_DATA_MANIFEST.csv identifies the packaged workbook and its scientific
source. The analysis reads its cells, formulas and recorded measurements.
