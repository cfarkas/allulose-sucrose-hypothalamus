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

## Figure 3: NPY-associated labeling and spatial profiles

[Figure 3](figure_updates/Fig3/) contains microscopy and nuclear maps for water,
sucrose and allulose, with investigator-reviewed anatomy on native DAPI images.
The complete figure has panels A–I. Panels E/F/G occupy the middle row and
spatial profiles H/I appear below.

Figure 3G illustrates six rings around a schematic third ventricle. The ring
origin is the mean position of the field's eligible DAPI nuclei. Ring boundaries
place approximately equal numbers of eligible nuclei in each ring, so inner
rings 1–2 contain approximately one-third of the nuclei. Equal nuclear counts
do not imply equal areas, thicknesses or atlas-registered anatomy.

The schematic uses illustrative DAPI positions, with 60 nuclei per ring,
concentric boundaries and inner rings near the floor of the ventricle. A
DAPI-mean key identifies the origin. One equation defines the within-ring
percentages: H counts c-FOS-positive nuclei, I counts c-FOS/NPY double-positive
nuclei, and both use all eligible DAPI nuclei in that ring as the denominator.
The NPY analysis uses one field per animal. All 108 animal-by-ring percentages
agree with their stored counts; the six nuclear counts for each animal differ
by at most one.

The [shared renderer](figure_updates/scripts/shared/spatial_ring_cartoon.py)
reproduces the schematic. Figure 4 refers to this same illustration. The NPY
profile PERMANOVA gives p=0.010714 and q=0.021429; the dispersion test gives
p=0.014286 and q=0.028571.

Radial shell methods have published precedents. Mehta et al. (2013),
[IMACULAT](https://doi.org/10.1371/journal.pone.0061386), used equal-area
elliptical shells within individual nuclei. This exploratory tissue-field
implementation uses approximately equal DAPI nuclear counts and does not
establish anatomical or functional validation of those bands.

## Figure 4: POMC labeling and microscopy

[Figure 4](figure_updates/Fig4/) has panels A–J and examines the familiarized
cohort after a four-hour fast. A uniform local-DAPI area threshold is applied
across conditions in 76 reconstructed sections: 780 of 3,403 processed POMC
objects are excluded. ARC c-FOS/POMC labeling is higher under allulose than
water in the nominal exact comparison (p=0.047619), with global p=0.176454 and
cage-mean p=0.200. The direction is consistent across 0.5/1.0/1.5 size thresholds;
nominal significance depends on the criterion. See the
[size-filter methods and sensitivity analysis](figure_updates/Fig4/provenance/POMC_SIZE_QC_20260909.md).

POMC rings are constructed within each reconstructed section. Counts from
corresponding rings are summed across sections before calculating each animal's
six percentages. The spatial double-positive comparison gives p=q=0.920080.

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
