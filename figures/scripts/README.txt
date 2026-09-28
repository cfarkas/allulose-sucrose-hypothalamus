PAPER SCRIPTS

Run from the Paper/ root with the configured scientific environments:
  ./reproduce_all_figures.sh --check-only
  ./reproduce_all_figures.sh --plan
  ./reproduce_all_figures.sh --output "$PWD" --force

The launcher stages, validates and installs thirteen figure packages. Complete
figures are English; individual panels and captions are English and Spanish.
PNGs use 600 dpi. PAPER_PYTHON and PAPER_S3_PYTHON select the interpreters.

Figure scripts
  Fig1/   Single-bottle behavior.
  Fig2/   Oral-challenge microscopy and regional c-FOS labeling.
  Fig3/   NPY-associated labeling and spatial clustering.
  Fig4/   POMC-associated labeling and ARC spatial clustering.
  Fig5/   GFAP/Iba1 signals and microglial morphology.
  FigS1/  Male single-bottle behavior.
  FigS2/  Female single-bottle behavior.
  FigS3/  Kidney, liver and spleen histology.
  FigS4/  Two-bottle behavior.
  FigS5/  First-exposure ACTH/CLIP-associated labeling.
  FigS6/  Human-reference morphology benchmark.
  FigS7/  Cellpose training documentation.
  FigS8/  Pilot-based sample-size estimates.

The figure-root scripts are mirrored here. Each Paper/FigN/README.txt describes
its panels, methods, input requirements and reproduction command.

Shared tools
  shared/recheck_current_statistics.py   Portable NPY/POMC numerical checks.
  shared/01_annotate_regions.py          ARC/ME/VMN anatomical review.
  shared/02_validate_raw_cellpose_masks.py  Raw/mask coverage validation.
  shared/03_safe_stage_fig45.py          Staging for Figures 4 and 5.
  shared/05_annotate_microglia.py        Microglial morphology review.
  shared/spanish_panel_text.py          Spanish panel text.

Validation and installation
  utilities/01_validate_bundle.py         Reconstructed-tree validation.
  utilities/03_check_reproducibility.py    Syntax, dependencies and mirrors.
  utilities/04_promote_figure.py          Validated figure installation.
  utilities/07_validate_figure_outputs.py  English masters and bilingual panels.
  utilities/09_audit_canonical_inputs.py   Raw-data and active-script audit.

The behavioral experiments use cages for shared measurements. NPY inference
uses animals with undocumented cage membership; POMC spatial inference uses
biological cages. Figure S3 phenotypes and Figure 5 morphology labels have the
interpretation limits described in their figure guides.
