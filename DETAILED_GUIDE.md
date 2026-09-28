# Detailed reproduction guide

Use this guide to prepare the data and software for Figures 1–5 and S1–S8.
The [figure guide](FIGURE_GUIDE.md) links each complete figure, its individual
panels, source tables and analysis instructions.

## Requirements and setup

Use 64-bit Linux, Python 3.10 or newer, Git and Conda. See [setup](SETUP.md).
Plan for at least 32 GiB available RAM and about 800 GB free disk space.
If the data and temporary directories are on different filesystems, allow
700 GB on the data filesystem and 100 GB in `/tmp`. These are planning allowances.
The data download is about 240 GB; downloaded and extracted files use about
560 GB before environments and working space.

```bash
git clone https://github.com/cfarkas/allulose-sucrose-hypothalamus.git
cd allulose-sucrose-hypothalamus
python3 reproduce.py --check-only
python3 reproduce.py
```

The program checks the computer, verifies downloads, reconstructs `Paper/`,
installs the repository figure package, prepares scientific environments,
checks readiness, and runs the thirteen-figure launcher. A failed check stops
the run. Repeat the command to resume verified downloads and setup stages.

## Options and logs

- `--check-only`: check the computer and storage without downloading or installing.
- `--prepare-only`: download, reconstruct, install software and check readiness;
  stop before rendering.
- `--archive-root /absolute/path/to/archives`: use a complete local set of
  canonical ZIP archives; all hashes are checked again.
- `--do_microglial_choices`: open an independent microglial review and train
  a candidate. The default reuses 203 saved choices.

Logs under `logs/` record computer checks, download URLs and byte counts,
software setup, readiness checks and figure execution. `Paper/FigN/` holds
the resulting figures, panels, captions and source tables. `Paper/analyses/`
holds generated analyses. Large data and environment directories are ignored by Git.

## Scientific environments

Conda may be on `PATH` or in a standard installation such as `~/miniforge3`.
The launcher creates `.paper-envs/main` and `.paper-envs/s3`, keeping package
caches in the project. It applies the environment locks, OpenCV 4.10.0.84
compatibility requirements and the numerical thread settings in
`reproduction-runtime.env`. Those settings are needed for the reference
microscopy masks and raster comparisons.

To use prepared environments, set `PAPER_PYTHON` and `PAPER_S3_PYTHON` to their
absolute Python executable paths. The launcher checks their versions and
scientific dependencies. It does not install packages into those custom
environments. Inspect `requirements-replay-compatibility.txt` and the setup
files in `figures/scripts/setup/` for the precise package requirements.

## Data and numerical scope

The GitHub `figures/` tree contains all 13 figure graphics, bilingual
panels, captions, scripts and small source inputs. Large microscopy datasets
come from the four records below. `figure-updates.json` verifies the repository
files installed into the reconstructed tree.

Raw images for WATER_NPY4–5 are not in the Zenodo datasets. Their reviewed
coordinates and classifications are included in the GitHub numerical bundle.
Full raw-image Figure 3 segmentation requires those images and the complete
local analysis inputs described in [its README](figures/Fig3/README.txt).

From the GitHub root, recompute NPY/POMC statistics without microscopy:

```bash
python3 -m venv .figures-venv
.figures-venv/bin/python -m pip install -r requirements-figures.txt
.figures-venv/bin/python figures/scripts/shared/recheck_current_statistics.py
```

Full Figure S3 inference requires the slide inputs, an inference environment,
CUDA and locally authorized HistoPLUS weights. Rendering can reuse verified
inference outputs. Figure S6 uses supplied prediction tables, while S7/S8 use
saved crops, loss values and simulation curves.

## Public Zenodo records

| Package | Zenodo record ID | DOI | Public record |
| --- | --- | --- | --- |
| Software and figure-reproduction package | `22265239` | [`10.5281/zenodo.22265239`](https://doi.org/10.5281/zenodo.22265239) | [https://zenodo.org/records/22265239](https://zenodo.org/records/22265239) |
| Main-figure source data (Figures 1–5) | `22265241` | [`10.5281/zenodo.22265241`](https://doi.org/10.5281/zenodo.22265241) | [https://zenodo.org/records/22265241](https://zenodo.org/records/22265241) |
| Supplementary-figure source data (Figures S1–S5) | `22265243` | [`10.5281/zenodo.22265243`](https://doi.org/10.5281/zenodo.22265243) | [https://zenodo.org/records/22265243](https://zenodo.org/records/22265243) |
| Figure S3 full-resolution H&E WSI exports (24 L0 BigTIFFs) | `22284125` | [`10.5281/zenodo.22284125`](https://doi.org/10.5281/zenodo.22284125) | [https://zenodo.org/records/22284125](https://zenodo.org/records/22284125) |

All four records and all files are public/open. All four records are required
to reconstruct the complete tree. The independently downloadable chunk files
are transport units, not separate datasets.

The logical Figure S3 WSI package contains the 24 full-resolution L0 BigTIFF
exports. Its chunks are hosted in the dedicated WSI record except for one final
quota-overflow transport chunk, which is hosted in the CC BY 4.0
supplementary-figure source-data record. The URL map records the physical host
of every chunk, so reconstruction is automatic and checksum-identical. The
matched L2 working-resolution exports, per-slide conversion receipts, aggregate
SHA-256 manifests, and downstream analysis/QC outputs also remain in the
supplementary-figure source-data record. Figure S3 therefore requires both
records.

The manifest excludes machine-local working state, the unrelated literature
scrape and its downloaded third-party article copies, and unused exploratory
inputs. It retains the manifest-bound acquisition receipts and analysis outputs
that the active reproduction reads.

Figure S3 uses the within-slide organ segmentation masks directly
throughout its active reproduction.

## Verified archives

| Logical record | Canonical archive | Compressed size | SHA-256 | Zenodo transport |
| --- | --- | ---: | --- | --- |
| `software_core` | `paper-software-core.part001.zip` | 1,912,685,403 bytes (1.91 GB; 1.78 GiB) | `e8b240899f3e3d19e50f17a2e6ebbd039184cbeaa2e5f0839d8085eb5ea694b2` | 15 independently verified chunks |
| `main_figure_data` | `paper-main-figure-data.part001.zip` | 21,914,164,144 bytes (21.91 GB; 20.41 GiB) | `e34fd3172a5216d73518b4044b979db498bd2770cf2f41b8305bc00ff0234ba2` | 28 independently verified chunks |
| `main_figure_data` | `paper-main-figure-data.part002.zip` | 15,246,200,736 bytes (15.25 GB; 14.20 GiB) | `4c50829eb6a2c71ef01615fe7b59fad8d93ae50cf15c0073ddcafc71c8dbcabf` | 20 independently verified chunks |
| `main_figure_data` | `paper-main-figure-data.part003.zip` | 21,935,080,784 bytes (21.94 GB; 20.43 GiB) | `f2fa219419729795b03e893a82ee1094deae70f467e9bfb50ca4ed456166ad4a` | 28 independently verified chunks |
| `main_figure_data` | `paper-main-figure-data.part004.zip` | 5,460,029,662 bytes (5.46 GB; 5.09 GiB) | `072981ddc9f36c66edf84a623eb621cf3ab14e70ebf64fac0759bb667dbf4ef4` | 7 independently verified chunks |
| `supplementary_figure_data` | `paper-supplementary-figure-data.part001.zip` | 24,396,735,136 bytes (24.40 GB; 22.72 GiB) | `f5e6d7dc44e4e616ec6feeae3db4421e67bee7c7fe9e85de57a33937f10edeac` | 48 independently verified chunks |
| `figs3_wsi` | `paper-figs3-wsi.part001.zip` | 37,963,064,474 bytes (37.96 GB; 35.36 GiB) | `af9e309f2f7f3b17b889db63e63220b4e69784c087dc39df370cd060f997ece1` | 25 independently verified chunks |
| `figs3_wsi` | `paper-figs3-wsi.part002.zip` | 37,187,908,584 bytes (37.19 GB; 34.63 GiB) | `f9359794aa65068a8847a6648aada0877da3f811960949ed14a6ad8a4d9a2aba` | 24 independently verified chunks |
| `figs3_wsi` | `paper-figs3-wsi.part003.zip` | 42,302,285,852 bytes (42.30 GB; 39.40 GiB) | `30ce3b2432d3b642c1436eceefb53710fc8551f12e077128bf19dee572361468` | 28 independently verified chunks |
| `figs3_wsi` | `paper-figs3-wsi.part004.zip` | 31,935,357,591 bytes (31.94 GB; 29.74 GiB) | `02e07dbb26dbbe1edaf4e187f7b3e61a025d94887ec4aa9a89f75d98c881c792` | 21 independently verified chunks |

The authoritative manifest has SHA-256 `df2ebb51339fdc3cc84dd8bbfb232647570d7e2cd0f12c19a5538ff7b7dd071d`; the sealed
transport manifest has SHA-256 `1ef5cfb1511cd4836d78208e51c7b61e9ba3cd9e5c0097c8d9ab9e0bae828055`. The four Zenodo
records contain 246 transport files in total: the
244 chunks plus the canonical and transport manifests.
The reconstructed tree SHA-256 is `3d9b410828413cd8a4a69b4347b70707d443d628e8491ee3f140ec7875d8f4df`. The reconstructor verifies
each chunk before concatenation, each canonical archive after concatenation,
every listed ZIP member, file size, file SHA-256, POSIX mode, path, and final
tree inventory. It rejects missing or extra URL-map entries, redirects with
credentials, unsafe paths, unlisted ZIP members, and an existing destination.


Data-package identifier in `release-metadata.json`: `v1.0.0`.

## Run a prepared Paper tree

After the repository launcher has installed the inputs and figure files:

```bash
cd Paper
./reproduce_all_figures.sh --check-only
./reproduce_all_figures.sh --plan
./reproduce_all_figures.sh --output "$PWD" --force
```

Figure computation uses fresh staging directories. The launcher validates
the English masters, bilingual panels, numerical inputs and output hashes
before installation. Raw-data and active-script audits bracket the run.
Individual figure READMEs provide scoped commands and their input requirements.

The complete thirteen-figure workflow has not been replayed end to end or
benchmarked. Source verification, individual figure reconstruction, numerical
tests and installation checks define the validation coverage. Check logs and
final validator results before treating a run as successful.

## Troubleshooting

- For insufficient disk or memory, provide the required resources and repeat
  the command. A preflight pass does not guarantee the scientific environment.
- For interrupted downloads, repeat the command; completed verified chunks
  are reused. Do not rename or combine transport chunks manually.
- For a checksum failure, identify the named input from the log and restore
  that file from its documented source.
- For environment failures, inspect the software log and required package
  versions. Dependency checks still apply when a setup receipt exists.
- For figure-input failures, consult the corresponding README and supply its
  raw, reviewed or derived inputs. Do not substitute missing observations with zeros.

## Interpretation and reuse

Experimental units are endpoint-specific; cells and sections are subsamples.
Figure 5 morphology classes are provisional assignments. Figure S3 labels are
exploratory model-predicted nuclear phenotypes. Acquisition and sample-size
limits are documented with the figures and in the [interpretation notes](MICROSCOPY_AND_INTERPRETATION.md).

Code and software-support files use MIT; authors' original documents, figures
and data use CC BY 4.0. See `LICENSE` and `LICENSES/`. File-specific and
third-party notices take precedence. Cite the applicable Zenodo DOIs when
reusing their datasets.
