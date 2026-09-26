# Figure guide

The reproduction workflow covers five main figures and eight supplementary
figures. Complete figures are in English; individual panels and captions are
available in English and Spanish.

[Reproduction instructions](README.md) · [Manual setup and archive checks](DETAILED_GUIDE.md)

## Code and source inputs

Analysis scripts, figure assets and small source inputs are in
[`figure_updates/`](figure_updates/), which follows the layout of the
reconstructed `Paper/` directory. Large microscopy datasets are retrieved from
the four [Zenodo records](README.md#paper-and-data-on-zenodo). Model binaries,
training-loss arrays and their roles are described in the [model guide](models/README.md).

The normal command is:

```bash
python3 reproduce.py
```

It verifies the Zenodo archives, reconstructs `Paper/`, checks the repository's
figure files against `figure-updates.json`, and installs those files before
running the thirteen-figure launcher. The installer checks SHA-256 hashes and
file modes and stops on unrecognized local edits. Installation backups are
saved under `.figure-update-backups/`. Repeat the command to resume a run.
`--figure-updates` is an explicit alias for this workflow.

The `--archived-release` option reproduces the ten-figure `v1.0.0` package from
an unmodified reconstruction. Its archive identities and manual replay commands
are documented in the [detailed guide](DETAILED_GUIDE.md).

## Figures 1, S1 and S2: single-bottle behavior

[Figure 1](figure_updates/Fig1/) presents fluid consumption and body-weight
change. [S1](figure_updates/FigS1/) and [S2](figure_updates/FigS2/) show the
sex-specific results. Pairwise comparisons use exact tests and Holm adjustment.
The weight analysis includes 24 mice in 12 cages.

The humane transfer of FR5-4 from allulose cage E2 to water cage E10 occurred
after the day-6 measurement. Linked records identify one animal and retain its
allulose assignment through that endpoint; E10/FR6-2 belongs to the water group.
The package includes the source consumption and weight tables, analysis scripts,
figures and captions. See the [Figure 1 instructions](figure_updates/Fig1/README.txt)
for reproduction from the small data package.

## Figure 2: hypothalamic labeling

Figure 2 describes the oral-challenge protocol and regional hypothalamic c-FOS
labeling. Its analysis scripts and inputs come from the Zenodo archives;
reproduced figures, source values and legends are saved in `Paper/Fig2/`.

## Figure 3: NPY-associated labeling and spatial clustering

[Figure 3](figure_updates/Fig3/) contains microscopy and nuclear maps for water,
sucrose and allulose, with investigator-reviewed anatomy on native DAPI images.
The complete figure has panels A–I. Panels E/F/G occupy the middle row and
spatial curves H/I appear below. The expanded cohort has 5/3/3
water/sucrose/allulose animals, each contributing one accepted field.

Figure 3G illustrates positive-cell pair counting in reviewed allulose field
E8_FR6-4. Curves H/I compare clustering over 20–150 µm using exact
random-labeling moments conditional on cell positions, tissue sides and
positive counts. The global envelope controls across radii; Holm correction
covers the two endpoints. These distances describe proximity, not anatomical
subregions or connectivity.

The total c-FOS comparison gives p=0.220130. The NPY-associated comparison gives
p=0.0285714 and Holm p=0.0571429. WATER_NPY3 and WATER_NPY4 have zero and one
double-positive nuclei and no estimable marker-associated curve; WATER_NPY5
has two and is included. Thus panel I uses 3/3/3 animals. Both added water
animals remain in the whole-field abundance analysis, whose enumerated Welch
omnibus p values are 0.0746753 (c-FOS/DAPI) and 0.0166667 (c-FOS/NPY); these
are unadjusted, and allulose–sucrose pairwise p values are 0.100.

The [portable numerical bundle](figure_updates/Fig3/source_data/current_statistics/)
recomputes these values with `scripts/shared/recheck_current_statistics.py`.
The former six-shell PERMANOVA and dispersion outputs remain historical.
Two added Zeiss water acquisitions improve instrument overlap, but mixed
acquisition/staining batches and unknown NPY cage membership still constrain
inference. Exact enumeration removes simulation error; finite-sample validity
requires exchangeability.

## Figure 4: POMC labeling and microscopy

[Figure 4](figure_updates/Fig4/) has panels A–J and examines the familiarized
cohort after a four-hour fast. A uniform local-DAPI area threshold is applied
across conditions in 76 reconstructed sections: 780 of 3,403 processed POMC
objects are excluded. ARC c-FOS/POMC labeling is higher under allulose than
water in the nominal exact comparison (p=0.047619), with global p=0.176454 and
cage-mean p=0.200. The direction is consistent across 0.5/1.0/1.5 size thresholds;
nominal significance depends on the criterion. See the
[size-filter methods and sensitivity analysis](figure_updates/Fig4/provenance/POMC_SIZE_QC_20260909.md).

Current panels I/J use ARC clustering curves and global envelopes across
20–150 µm. Animal curves are averaged equally within biological cages.
The c-FOS endpoint includes 14 animals in 2/3/3 water/sucrose/allulose cages
(p=0.592857); the POMC-associated endpoint includes 13 animals in 2/2/3 cages
(p=0.714286). Both two-endpoint Holm values are 1.000. Sparse positive labeling
precludes the corresponding ME analysis. The
[portable numerical bundle](figure_updates/Fig4/source_data/current_statistics/)
contains the current inputs and results; the former shell tests are historical.

Panel F and the POMC/NPY-GFP inset display measured channel intensities within
accepted marker regions. POMC uses the full-field normalization and amber tint
of panel C. Screen compositing retains channel-intensity variation; region
membership defines where a signal is displayed. Display checks cover dark
source pixels, intensity variation, label numbering, region membership,
bright-channel compositing and registered geometry.

## Figures 5 and S6: glial signals and microglial morphology

[Figure 5](figure_updates/Fig5/) provides glial analyses, mask refinement,
human review, supervised CNN training, frozen DINOv2 features and evaluation
that holds out whole animals. Classifier-derived microglial labels are
provisional morphological assignments.

[S6](figure_updates/FigS6/) compares morphometric proposals with frozen DINOv2
image features on the same **203 human-reviewed cells from 15 animals**.
Balanced accuracy is **73.9% versus 57.0%**, respectively. The figure includes
paired animal-bootstrap intervals, class recall and confusion matrices.
DINOv2 preprocessing and classifier selection exclude each held-out animal.
The unsupervised morphometric proposals use the full cell collection without
fitting human class labels. This is an internal benchmark on a selected review
sample; morphology alone does not establish functional microglial activation.

### Render S6 without the large data download

Create an environment in the Git repository:

```bash
python3 -m venv .figures-venv
.figures-venv/bin/python -m pip install -r requirements-figures.txt
```

Render into a fresh output directory:

```bash
.figures-venv/bin/python figure_updates/FigS6/01_make_figure_s6_microglia_classifier.py --validate-only
.figures-venv/bin/python figure_updates/FigS6/01_make_figure_s6_microglia_classifier.py --output-dir figure-renders/S6 --dpi 600
.figures-venv/bin/python figure_updates/scripts/utilities/07_validate_figure_outputs.py --figure FigS6 --root figure-renders/S6
```

The renderer recomputes scores and checks frozen input hashes and animal
partitions. It needs no GPU, microscopy download, pretrained weights or
retraining. Use a different output directory for a subsequent render.

### Reuse saved choices or make your own

Normal runs load the 203 saved choices, their masks and the matching trained
candidate without opening a reviewer. The inputs are in
[Fig5/microglial_review_data](figure_updates/Fig5/microglial_review_data/): an open
CSV and ten ZIPs, about 74 MB total, containing 300 review crops/masks, 14,415
classifier input crops and soma masks, source tables and the saved CNN candidate.

To make independent choices during Figure 5, run:

```bash
python3 reproduce.py --do_microglial_choices
```

Or, from a reconstructed `Paper/` tree:

```bash
./reproduce_all_figures.sh --do_microglial_choices
```

At the Figure 5 step, the terminal displays `MICROGLIAL_REVIEW_URL`. Open that
address in a browser. The session has **300 targets and zero prefilled human
choices**. Review masks, paint corrections if needed and select classes. Finish
after 200–300 choices with class/animal coverage; at 300, training starts
automatically. The workflow continues after training succeeds. On a remote
server, forward the printed localhost port through SSH to your computer.

Independent choices are saved separately from the default 203-label set. The
printed review directory can be resumed with `--review-dir`. Candidate results
and the label receipt are saved in
`Paper/analyses/Fig5/results/current_microglial_choices/`. This branch is separate
from the fixed Figure 5 panels and S6 benchmark.

For review without the microscopy download, see the
[Figure 5 commands](figure_updates/Fig5/README.txt). DINOv2 training uses the
verified upstream encoder obtained by `11_prepare_dinov2_assets.py`.

## Figure S3: peripheral histology

[S3](figure_updates/FigS3/) analyzes liver, kidney and spleen histology. Observed
and permuted treatment coefficients use matched HC3 standard errors. Each
coefficient-specific Freedman–Lane null retains the other treatment contrast
and source cohort.

With 99,999 permutations, the hepatic epithelial-like Allulose/Water CLR ratio
is 4.06 (95% CI 2.22–7.42), p=0.00031 and global BH q=0.02418. Six classical
contrasts are nominal, including lower hepatic nuclear density; none passes BH.
The phenotype labels are exploratory model predictions. See the
[statistical methods and validation](figure_updates/FigS3/provenance/STATISTICAL_CORRECTION_20260909.md).

## Figure S4: two-bottle behavior

S4 contains the two-bottle consumption and food-intake analyses. Its scripts
and inputs come from the Zenodo archives; reproduced source data, figures and
captions are saved in `Paper/FigS4/`.

## Figure S5: first exposure and ACTH/CLIP labeling

[S5](figure_updates/FigS5/) examines first oral exposure after a 16-hour fast,
using ACTH/CLIP as a POMC-related signal. It includes three animals per
condition, including NPY-M, a male NPY-GFP water animal. Four-channel assignments
and DAPI-based anatomy are recorded for all nine animals.

Spatial c-FOS and double-positive profiles give exact PERMANOVA p=0.400 and
p=0.464286, respectively, with BH q=0.464286 for both. The eight-WT subset is
supplied as a sensitivity comparison.

## Figure S7: Cellpose training

[S7](figure_updates/FigS7/) contains training-image crops, reviewed masks and
recorded loss curves, linked to their sources by SHA-256. The displayed POMC
training masks include the size filter; the archived models were not retrained
on those filtered masks. Human correction of training masks is distinct from
independent evaluation.

## Figure S8: sample-size planning

[S8](figure_updates/FigS8/) presents Monte Carlo power curves from the
[Figure 1 estimator](figure_updates/Fig1/03_estimate_sample_size_monte_carlo.py).
The simulations use day-6 outcomes from 12 pilot cages, 500,000 simulations per
group size, group-specific SDs and Welch ANOVA. Five cages per group reach
80.16% power for bottle-volume disappearance; fourteen reach 80.77% for
cage-mean weight change. These estimates are conditional on the pilot data
and do not establish achieved power or neuronal/glial sample sizes.
[Source values and curves](figure_updates/Fig1/source_data/sample_size_monte_carlo/)
are included.

## Verification

Verification covers the following scopes:

- An end-to-end replay of the archived ten-figure package reproduced all 136
  reference PNGs byte for byte.
- Separate Figure 1/S6 checks reproduced 36 PDF/PNG outputs byte for byte.
- Figure 3 checks cover ring geometry, stored counts, spatial calculations and
  agreement with numerical reference tables.
- Figure 4 checks cover POMC size filtering, threshold sensitivity and native
  microscopy display.
- S3 checks include an independent per-permutation HC3 reference calculation
  and full data/figure validation.
- S5 reconstruction from native channels, saved segmentations and accepted
  anatomy agrees with the nine-animal regional and spatial results.
- S7/S8 checks reproduce the reference PNG figures and panels byte for byte.

A complete end-to-end replay of all thirteen figures has not been performed.
The launcher, file installation and individual reconstruction checks are
validated separately. Repository checks can be run with:

```bash
python3 validate_repository.py --root .
python3 -m unittest discover -s tests -v
python3 validate_figure_updates.py
```

Full-data pipeline tests require the reconstructed `Paper/` tree and its
microscopy exports. Code is MIT-licensed; original data and figures are
CC BY 4.0. See the [license](LICENSE) for scope and third-party notices.
