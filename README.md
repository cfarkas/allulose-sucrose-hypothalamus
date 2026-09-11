# Reproduce the allulose–sucrose paper
The 10 September update introduces the DAPI ring method once in [Figure 3G](figure_updates/Fig3/Figure_3_cFos_NPY.pdf); [Figure 4](figure_updates/Fig4/Figure_4.pdf) retains its A–J data panels, with unchanged experimental statistics. See [update details](FIGURE_UPDATES.md).


This study compares allulose, a sweetener with few calories, with sucrose
(table sugar) and water in mice. It examines feeding behavior and markers of
brain activity after prior exposure, focusing on the hypothalamus, which helps
regulate appetite.

Recreate **all 13 figures (1–5 and S1–S8)**. Figure 1 includes the corrected
Holm terminology and reconciled cage transfer. Figure 3 includes sucrose and
its reviewed boundaries; S5 includes NPY-M (three animals per condition).
S6 compares 73.9% versus 57.0% balanced accuracy. S7 documents Cellpose training; S8 shows pilot-based sample-size estimates. S3 uses corrected HC3-studentized residual permutations (hepatic epithelial-like p=0.00031, BH q=0.02418).

Figure 4 now applies uniform local-DAPI size QC to POMC masks and reports threshold sensitivity. Its familiarized cohort had a 4-hour fast; the independent ACTH/CLIP experiment in S5 had a 16-hour fast. Original images and masks remain preserved.

**Normal runs reuse the 203 saved human microglial choices without prompting.**
To classify the cells yourself and retrain, add the optional flag:
`python3 reproduce.py --do_microglial_choices`.

The [figure scripts, new review data, and S6 inputs](FIGURE_UPDATES.md) are included
directly in GitHub. The original microscopy is downloaded from the existing
Zenodo records.

## Paper and data on Zenodo

All four records are open access. The command below downloads them automatically.
They contain the original ten-figure release; the verified GitHub update adds the
current scripts, saved reviews, corrected figures, pilot-based sample-size
calculations and S6–S8.

| Contents | Zenodo record |
| --- | --- |
| Software, manuscript, and figure files | [22265239](https://zenodo.org/records/22265239) |
| Main-figure data (1–5) | [22265241](https://zenodo.org/records/22265241) |
| Supplementary-figure data (S1–S5) | [22265243](https://zenodo.org/records/22265243) |
| Full-resolution tissue scans for Figure S3 | [22284125](https://zenodo.org/records/22284125) |

## 1. Prepare your computer

Use a computer with:

- **A 64-bit Linux computer with an Intel or AMD processor.** A GPU is not needed.
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

The program checks your computer, downloads and unpacks all the data from
Zenodo, installs the required software for this project, and recreates all thirteen
figures. It reuses the saved human choices and candidate classifier for the
additional microglial review branch. It applies the tested settings automatically and checks the results
against the published figures.

**Keep the terminal open and the computer awake.** Downloading may take several
hours. The archived ten-figure calculations took about 37 minutes on our test
server. The current thirteen-figure workflow includes additional reconstruction
steps; its total runtime has not been benchmarked.

When everything has passed, you will see:

```text
SUCCESS: All 13 figures reproduced and verified.
```

## Find your results

Inside this repository, open:

- `Paper/Fig1/` through `Paper/Fig5/` for the main figures.
- `Paper/FigS1/` through `Paper/FigS8/` for the supplementary figures.

Each folder contains PNG/PDF figures and a `source_data/` folder with plotted
values and statistics. Main figures are in English; individual panels and
legends are available in English and Spanish.

The archived Zenodo payload also contains its original documents. The current
thesis and manuscript are maintained locally and are not included in this
GitHub figure/code update.
Run logs are saved in `logs/`.

## If something stops

Read the message, fix the reported problem, and run `python3 reproduce.py`
again. Completed downloads are reused. The program stops if a check fails.

To check your computer without downloading or installing anything:

```bash
python3 reproduce.py --check-only
```

For manual commands, checksums, troubleshooting, and scientific limitations,
see the [detailed guide](DETAILED_GUIDE.md).

The original ten-figure replay reproduced all 136 PNGs exactly. Separate checks
reproduced all 36 corrected Figure 1/S6 PDF/PNG files byte for byte. A complete
new thirteen-figure replay has not been rerun. The complete Figure 3 and
nine-animal S5 have separate reconstruction checks; the update guide describes
their scope. Figure 5's candidate class labels
remain provisional, as explained in the [update guide](FIGURE_UPDATES.md).

For the original ten-figure release in an unmodified reconstruction, use
`python3 reproduce.py --archived-release`.

**Reusing this work:** cite the relevant Zenodo records using their DOIs.
Code is MIT-licensed; original documents and data are CC BY 4.0
([license details](LICENSE)).

Figure 4F now displays microscopy intensities inside the accepted POMC ROIs,
including the POMC/NPY-GFP inset. No uniform marker fills are used in those
microscopy views. The matching scripts, legends, and display regression checks
are included in [the figure updates](FIGURE_UPDATES.md).

Figure 3G places the illustrative inner rings near the third-ventricle floor.
All six rings still follow the DAPI mean and count quantiles; the updated shared
renderer and Figure 3 scripts reproduce this placement. The cartoon appears once.
