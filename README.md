# Reproduce the allulose–sucrose paper

This study compares allulose, a sweetener with few calories, with sucrose
(table sugar) and water in mice. It examines feeding behavior and markers of
brain activity after prior exposure, focusing on the hypothalamus, which helps
regulate appetite.

This repository contains code and source tables for **all 13 figures
(1–5 and S1–S8)**. It includes analysis scripts, figure source tables, saved
human microglial labels and trained-model downloads. The workflow retrieves
the archived microscopy datasets from Zenodo. Later figure updates are tracked
separately from that fixed archive.

## Current NPY and POMC statistics — 26 September 2026

The current Figure 3 includes **five water, three sucrose and three allulose
animals**. Its NPY-associated spatial comparison uses three estimable water
animals and gives nominal p=0.0285714, **Holm p=0.0571429**. The Figure 4
spatial comparisons use biological cages and are nonsignificant. These are
exploratory analyses with acquisition/batch and independence limitations.

The [27 September interpretation notes](MICROSCOPY_AND_INTERPRETATION.md) explain
what ROI normalization controls, why magnification alone does not invalidate
these measurements, the remaining acquisition differences, behavioral units
and the limits of neuronal/glial marker interpretation. Numerical results are
unchanged.

Recompute all current cell-pair moments, spatial curves, envelope tests and
NPY abundance statistics without downloading microscopy:

```bash
python3 -m venv .figures-venv
.figures-venv/bin/python -m pip install -r requirements-figures.txt
.figures-venv/bin/python figure_updates/scripts/shared/recheck_current_statistics.py
```

The [NPY inputs](figure_updates/Fig3/source_data/current_statistics/) and
[POMC inputs](figure_updates/Fig4/source_data/current_statistics/) contain
reviewed cell coordinates, classifications, geometry maps and reference
results. This check begins after segmentation. The two September water
acquisitions postdate the Zenodo archive; the archived download alone cannot
recreate their raw-image segmentation. Old six-shell Figure 3–4 tests are
historical and should not be substituted for these results.

## Figures and analyses

| Figures | Contents |
| --- | --- |
| 1, S1 and S2 | Single-bottle fluid consumption and body-weight change |
| 2–4 | Hypothalamic c-FOS labeling, NPY/POMC co-labeling and spatial profiles |
| 5 and S6 | Glial signals, microglial morphology and comparison with human review |
| S3 | Liver, kidney and spleen histology |
| S4 | Two-bottle consumption and food intake |
| S5 | ACTH/CLIP-associated c-FOS labeling after first oral exposure |
| S7 and S8 | Cellpose training and pilot-based sample-size estimates |

The [figure guide](FIGURE_GUIDE.md) explains the analyses, source inputs and
validation for each figure, including the cell-pair illustration in
[Figure 3G](figure_updates/Fig3/Figure_3_cFos_NPY.pdf).

**Normal runs reuse the 203 saved human microglial choices without prompting.**
To classify the cells yourself and retrain, add the optional flag:
`python3 reproduce.py --do_microglial_choices`.

## Trained models

[Download the 21 project-trained model binaries and 13 training-loss arrays](models/README.md).
The model guide provides SHA-256 checksums, model roles and a download helper,
and identifies the manuscript checkpoints and development candidates.
The [model implementation notes](figure_updates/Fig5/source_data/deep_learning_implementation.md)
describe TinyMorphCNN's architecture and 79,780 parameters, frozen DINOv2
features, downstream fitting, training settings and computational hardware.

## Paper and data on Zenodo

Four open-access records provide the software archive and source data. The
reproduction command downloads these records and combines them with the figure
scripts and inputs supplied in this repository.

| Contents | Zenodo record |
| --- | --- |
| Software, manuscript, and figure files | [22265239](https://zenodo.org/records/22265239) |
| Main-figure data (1–5) | [22265241](https://zenodo.org/records/22265241) |
| Supplementary-figure data (S1–S5) | [22265243](https://zenodo.org/records/22265243) |
| Full-resolution tissue scans for Figure S3 | [22284125](https://zenodo.org/records/22284125) |

## 1. Prepare your computer

Use a computer with:

- **A 64-bit Linux operating system and an Intel or AMD processor.** A GPU is not needed.
- **About 800 GB of free disk space** and **at least 32 GiB of available RAM**.
- **Python 3.10 or newer, Git, and Conda.** [Install these first if needed →](SETUP.md)

**The download is about 240 GB.** Keeping the download and unpacked data takes
about 560 GB, with extra space needed while the program runs. If the data and
`/tmp` are on different drives, allow 700 GB on the data drive and 100 GB in
`/tmp`. The program checks these space and memory allowances before downloading.

Native Windows and macOS are not supported by this workflow. Use a suitable
Linux computer or server.

## 2. Download this repository

Open **Terminal** in the folder where you want to keep the project. Copy and
run these two lines:

```bash
git clone https://github.com/cfarkas/allulose-sucrose-hypothalamus.git
cd allulose-sucrose-hypothalamus
```

## 3. Reproduce the paper

Run:

```bash
python3 reproduce.py
```

The program checks your computer, downloads and unpacks the data from Zenodo,
installs the required software, and recreates all thirteen figures. It reuses
the saved human choices and candidate classifier for the microglial review
branch, applies the tested settings and checks the figure outputs.

**Keep the terminal open and the computer awake.** Downloading may take several
hours. Total runtime for the thirteen-figure workflow has not been benchmarked.

When everything has passed, you will see:

```text
SUCCESS: All 13 figures reproduced and verified.
```

## Find your results

Inside this repository, open:

- `Paper/Fig1/` through `Paper/Fig5/` for the main figures.
- `Paper/FigS1/` through `Paper/FigS8/` for the supplementary figures.

Each folder contains PNG/PDF figures and a `source_data/` folder with plotted
values and statistics. Complete figures are in English; individual panels and
legends are available in English and Spanish. Run logs are saved in `logs/`.

## If something stops

Read the message, fix the reported problem, and run `python3 reproduce.py`
again. Completed downloads are reused. The program stops if a check fails.

To check your computer without downloading or installing anything:

```bash
python3 reproduce.py --check-only
```

For manual commands, checksums, troubleshooting and software requirements,
see the [detailed guide](DETAILED_GUIDE.md).

## Verification and interpretation

Validation includes an end-to-end replay of the archived ten-figure package,
with all 136 PNGs matching its reference files, and separate reconstruction,
numerical and display checks for the thirteen-figure workflow. A complete
end-to-end replay of all thirteen figures has not been performed. The
[figure guide](FIGURE_GUIDE.md#verification) describes the scope of those checks.

Figure 5's classifier-derived labels are provisional morphological assignments.
The S6 benchmark uses an internal human-review sample; morphology alone does
not establish functional microglial activation.

To reproduce the ten-figure package associated with `v1.0.0`, use
`python3 reproduce.py --archived-release` in an unmodified reconstruction.

**Reusing this work:** cite the relevant Zenodo records using their DOIs.
Code is MIT-licensed; original documents and data are CC BY 4.0
([license details](LICENSE)).
