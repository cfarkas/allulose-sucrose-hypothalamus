# Reproduce the allulose–sucrose paper

This study compares allulose, a sweetener with few calories, with sucrose
(table sugar) and water in mice. It examines feeding behavior and markers of
brain activity after prior exposure, focusing on the hypothalamus, which helps
regulate appetite.

Recreate **all 11 figures (1–5 and S1–S6)**. Figure 1 includes the corrected
Holm terminology, and S6 compares 73.9% versus 57.0% balanced accuracy.

**Normal runs reuse the 203 saved human microglial choices without prompting.**
To classify the cells yourself and retrain, add the optional flag:
`python3 reproduce.py --do_microglial_choices`.

The [figure scripts, new review data, and S6 inputs](FIGURE_UPDATES.md) are included
directly in GitHub. The original microscopy is downloaded from the existing
Zenodo records.

## Paper and data on Zenodo

All four records are open access. The command below downloads them automatically.
They contain the original ten-figure release; the verified GitHub update adds the
current scripts, saved microglial review data and S6.

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
Zenodo, installs the required software for this project, and recreates all eleven
figures. It reuses the saved human choices and candidate classifier for the
additional microglial review branch. It applies the tested settings automatically and checks the results
against the published figures.

**Keep the terminal open and the computer awake.** Downloading may take several
hours. The figure calculations took about 37 minutes on our test server;
allow longer on a slower computer.

When everything has passed, you will see:

```text
SUCCESS: All 11 figures reproduced and verified.
```

## Find your results

Inside this repository, open:

- `Paper/Fig1/` through `Paper/Fig5/` for the main figures.
- `Paper/FigS1/` through `Paper/FigS6/` for the supplementary figures.

Each folder contains PNG/PDF figures and a `source_data/` folder with plotted
values and statistics. Main figures are in English; individual panels and
legends are available in English and Spanish.

The manuscript/thesis PDF and Word files are directly inside `Paper/`.
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
new eleven-figure replay has not been rerun. Figure 5's candidate class labels
remain provisional, as explained in the [update guide](FIGURE_UPDATES.md).

For the original ten-figure release in an unmodified reconstruction, use
`python3 reproduce.py --archived-release`.

**Reusing this work:** cite the relevant Zenodo records using their DOIs.
Code is MIT-licensed; original documents and data are CC BY 4.0
([license details](LICENSE)).
