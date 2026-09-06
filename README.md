# Reproduce the allulose–sucrose paper

Recreate **all 10 figures (1–5 and S1–S5)** from the published data.
You do not need to edit code or choose analysis settings.

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
Zenodo, installs the required software for this project, and recreates all ten
figures. It applies the tested settings automatically and checks the results
against the published figures.

**Keep the terminal open and the computer awake.** Downloading may take several
hours. The figure calculations took about 37 minutes on our test server;
allow longer on a slower computer.

When everything has passed, you will see:

```text
SUCCESS: All 10 figures reproduced and verified.
```

## Find your results

Inside this repository, open:

- `Paper/Fig1/` through `Paper/Fig5/` for the main figures.
- `Paper/FigS1/` through `Paper/FigS5/` for the supplementary figures.

Each folder contains the figure images, PDFs, and related outputs. Run logs
are saved in `logs/`.

## If something stops

Read the message, fix the reported problem, and run `python3 reproduce.py`
again. Completed downloads are reused. The program stops if a check fails.

To check your computer without downloading or installing anything:

```bash
python3 reproduce.py --check-only
```

For manual commands, Zenodo records, checksums, and troubleshooting, see the
[detailed guide](DETAILED_GUIDE.md).

The tested full replay reproduced all 136 figure PNGs exactly.
[License and data-use terms](LICENSE).
