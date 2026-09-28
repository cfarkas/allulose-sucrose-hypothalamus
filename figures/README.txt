REPRODUCE THE ALLULOSE–SUCROSE HYPOTHALAMUS FIGURES

This Paper tree contains five main figures and eight supplementary figures.
Each figure has an English PDF and 600-dpi PNG master, individual English and
Spanish panels, captions, source tables, and a README describing its methods.

Figure map
  Fig1   Single-bottle fluid consumption and body-weight change.
  Fig2   Oral-challenge design, microscopy and regional c-FOS labeling.
  Fig3   NPY-associated c-FOS labeling and spatial clustering.
  Fig4   POMC-associated c-FOS labeling in ARC/ME and ARC spatial clustering.
  Fig5   Regional glial signals and predicted microglial morphology.
  FigS1  Male single-bottle behavioral measurements.
  FigS2  Female single-bottle behavioral measurements.
  FigS3  Kidney, liver and spleen histology.
  FigS4  Two-bottle consumption, food removal and body weight.
  FigS5  ACTH/CLIP-associated c-FOS after first oral exposure.
  FigS6  Microglial morphology classification benchmark.
  FigS7  Human-supervised Cellpose training.
  FigS8  Pilot-based Monte Carlo sample-size estimates.

Start from GitHub
  git clone https://github.com/cfarkas/allulose-sucrose-hypothalamus.git
  cd allulose-sucrose-hypothalamus
  python3 reproduce.py --check-only
  python3 reproduce.py

The repository README and DETAILED_GUIDE.md describe computer requirements,
downloads, environments and logs. The launcher verifies the Zenodo data,
reconstructs Paper/, installs the repository's figure files, and runs the
thirteen-figure workflow. Data availability for the NPY raw acquisitions is
described in Fig3/README.txt; their reviewed coordinate inputs support the
portable numerical check below.

Run from a prepared Paper/ tree
  ./reproduce_all_figures.sh --check-only
  ./reproduce_all_figures.sh --plan
  ./reproduce_all_figures.sh --output "$PWD" --force

PAPER_PYTHON and PAPER_S3_PYTHON select the configured scientific interpreters.
The launcher renders into fresh staging directories, validates the figure
outputs, installs them in FigN/, and verifies installed hashes. It checks raw
inputs and scripts before and after the run. Use each figure's README for an
individual build and its required inputs.

Portable neuronal statistical check, from the GitHub root
  python3 -m venv .figures-venv
  .figures-venv/bin/python -m pip install -r requirements-figures.txt
  .figures-venv/bin/python figures/scripts/shared/recheck_current_statistics.py

This recomputes NPY/POMC abundance and spatial statistics from reviewed cell
coordinates without downloading microscopy. Animals are the NPY units;
biological cages are the POMC spatial units. Acquisition and sampling differences
remain limitations of the measured populations.

Microglial review
Normal reproduction reuses 203 saved human choices and the matching candidate.
To make independent choices, add --do_microglial_choices to reproduce.py or
reproduce_all_figures.sh. The terminal supplies a local review URL; Fig5/README.txt
explains the review and training workflow. Predicted morphology classes are
shape assignments and do not independently establish functional activation.

Validation
  python scripts/utilities/07_validate_figure_outputs.py --all-live
This checks the English masters, matching bilingual panel sets, PNG resolution,
PDF signatures and figure legends. Full thirteen-figure end-to-end reproduction
has not been benchmarked; individual numerical, rendering and installation
checks define the validation scope.
