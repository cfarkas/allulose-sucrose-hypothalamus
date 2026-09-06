# Allulose–sucrose hypothalamus figure reproduction

This repository is the small, public reconstruction companion for:

> Oral allulose after prior familiarization is associated with a distinct
> spatial c-FOS/NPY profile without detectable POMC or glial activation in mice

Repository: [https://github.com/cfarkas/allulose-sucrose-hypothalamus](https://github.com/cfarkas/allulose-sucrose-hypothalamus)  
Release: `v1.0.0`

The large inputs are public on Zenodo; they are not committed to GitHub. The
checksum-bound `manifest.json`, `transport_manifest.json`, and
`zenodo-urls.json` in this repository fetch 244
independently verified chunks, reassemble the 10 canonical ZIP
archives, verify their sizes and SHA-256 digests, and reconstruct the exact
`Paper/` tree used to reproduce Figures
1–5 and S1–S5.

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

Figure S3 uses the original within-slide organ segmentation masks directly
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

## 1. Clone the GitHub companion and validate the release

Use a **64-bit Linux machine**, preferably Linux x86-64, and open a Bash
terminal on the large filesystem where you want to keep the reproduction.
Python 3.10 or newer and Git must already be available. These first commands
download the small companion repository; they do not download the paper data.

```bash
git clone https://github.com/cfarkas/allulose-sucrose-hypothalamus.git
cd allulose-sucrose-hypothalamus
python3 validate_repository.py --root .
sha256sum manifest.json
```

The validator must report `[PASS]`. The manifest command must print:

```text
df2ebb51339fdc3cc84dd8bbfb232647570d7e2cd0f12c19a5538ff7b7dd071d  manifest.json
```

**STORAGE WARNING — read before downloading anything from Zenodo.** The Git
clone is small, but the complete release requires **240,253,512,366 bytes
(240.25 GB; 223.75 GiB) of archives** and **319,524,773,963 bytes (319.52 GB;
297.58 GiB) of extracted files**. Keeping both uses **559,778,286,329 bytes
(559.78 GB; 521.33 GiB)** before filesystem overhead, Conda environments,
package caches, logs, generated analyses and temporary figure builds.

Plan for **at least 700 GB free on the filesystem holding this clone**, plus
**100 GB free on the filesystem holding `/tmp`**. If the clone and `/tmp`
share one filesystem, plan for **800 GB free on that filesystem** so that the
same free bytes are not counted twice. These are planning allowances, not a
measured peak or a guarantee. Conda's environments/cache may be on a third
filesystem; that filesystem needs its own free space. Account for user quotas,
other running jobs and storage used after this check.

The following instructions keep the reconstructed paper **inside your clone**:

```text
allulose-sucrose-hypothalamus/
  README.md
  manifest.json
  check_reproduction.py
  reconstruct_paper.py
  zenodo-archives/              completed canonical ZIP archive cache
  Paper/                       reconstructed paper; all final figures live here
  logs/                        your download, preflight and reproduction logs
```

During downloads, `.zenodo-archives.downloading/` is present instead of the
completed `zenodo-archives/`. During extraction, `.Paper.reconstructing.PID/`
is present instead of `Paper/`. These temporary directories are siblings
inside the clone. `.gitignore` excludes the large data, caches and local logs
from Git. Do not create an empty `Paper/` or `zenodo-archives/` yourself.

## 2. Check this machine before the 240 GB download

Run the portable checker from the clone:

```bash
python3 check_reproduction.py --root . --stage download --check-zenodo
```

It uses the Python standard library and reports:

| Check | Requirement or interpretation |
| --- | --- |
| Release metadata | Manifest, transport map, checksums, record IDs, README and license validation must pass. |
| Operating system | Full reproduction requires 64-bit Linux. Linux x86-64 is the reference platform for the explicit Figure S3 locks. Other Linux architectures get a warning and have not been certified by this check. |
| Bootstrap Python | Python 3.10 or newer runs validation, downloads and this checker. The scientific environments have their own pinned Python versions. |
| Shell/tools | Bash 4.1 or newer, GNU coreutils/findutils/diffutils and the logging command `tee` must be available on `PATH`. The checker reports missing commands and verifies GNU `cp` support for the installation flags. |
| Filesystem | The clone, reconstructed `Paper/` and `/tmp` must be writable, case-sensitive, preserve POSIX modes and allow executable scripts. Small temporary probe files are created and removed. Avoid FAT/exFAT or mounts with `noexec`. |
| Storage | Exact archive/tree sizes come from the validated manifest. Free bytes are checked on the actual destination and `/tmp` filesystems (including a separately mounted existing `Paper/`); shared storage is counted once with the reserves added together. Free inodes are checked where available. |
| RAM | Available RAM and exposed container memory limits are reported. **32 GiB available** is a planning recommendation, not a demonstrated minimum. Scheduler limits and other processes can reduce usable RAM. |
| CPU | Available CPU affinity and exposed container CPU quotas are reported. The reference run used roughly 16 cores; fewer cores will generally take longer. Fewer than four available cores triggers a warning. |
| Figure 3 raster compatibility | Apply `requirements-replay-compatibility.txt` to the main environment: OpenCV 4.10.0.84 is required for exact agreement with the frozen HIL polygons. The preflight checks actual polygon/mask equality. |
| Numerical threads | Source `reproduction-runtime.env` before reproduction. The frozen Figure 5 CPU CNN requires the tested 64-thread library profile; this is separate from the number of physical CPU cores. |
| GPU | A GPU is not required for the portable paper replay, which supports the CPU fallback and uses released WSI exports/receipts. |
| Conda | Discovers Conda from `CONDA_EXE`, `PATH` or common installation directories under your home directory. Reports the main and Figure S3 interpreters if already installed. |
| Public Zenodo access | With `--check-zenodo`, makes one HEAD request to a sample file in each of the four physical records and compares the advertised size. This checks connectivity, not every file or its contents. Full verification happens during download. |

`[FAIL]` means the requested stage is blocked; the checker exits with code 2.
`[WARN]` describes a recommendation or limitation; warnings are visible but
normally do not block the stage. Add `--strict` to make warnings block too:

```bash
python3 check_reproduction.py --root . --stage download --check-zenodo --strict
```

A fresh clone can report missing scientific environments as `[INFO]`: their
pinned setup files arrive with `Paper/`. This does not claim that figure
reproduction is ready. A successful first check ends with
`[PASS] READY TO DOWNLOAD/RECONSTRUCT`. After extraction and environment setup,
you must run the separate reproduction check in step 5.

For manual inspection of the same disks:

```bash
df -h . /tmp
df -i . /tmp
python3 --version
uname -srm
```

If the clone is on a small filesystem, place a new clone on a larger Linux
filesystem and continue there. Setting `TMPDIR` alone does **not** relocate the
published launcher's scratch space: it explicitly creates stages under
`/tmp`. Arrange sufficient space on that mount before reproduction.

## 3. Download all four Zenodo records and reconstruct `Paper/`

### What “cloning Zenodo” means here

Git clones only the small GitHub companion. Zenodo hosts ordinary public
files, not another Git repository. `reconstruct_paper.py` uses
`zenodo-urls.json` to fetch the exact 244 chunk files across the four public
records, verifies them, assembles ten ZIP archives and extracts the complete
manifest-bound paper. No Zenodo account, personal token or authentication
header is required for these public downloads.

Do not try `git clone` on a Zenodo DOI, manually concatenate chunks or unzip
each `.chunk...` file. A chunk is a byte range of one canonical ZIP, not an
independent ZIP. The URL map also handles the Figure S3 overflow chunk hosted
in the supplementary-data record. All four records are necessary; downloading
only the software record cannot reproduce the paper.

### First invocation, from the clone root

Remain in `allulose-sucrose-hypothalamus/`, the directory containing
`manifest.json`. Define the clone location for the commands below. If you open
a new terminal later, return here and define it again.

```bash
REPRO_REPO="$PWD"
mkdir -p "$REPRO_REPO/logs"
set -o pipefail

python3 check_reproduction.py --root "$REPRO_REPO" --stage download --check-zenodo \
  2>&1 | tee "$REPRO_REPO/logs/01-download-preflight.log"
```

Proceed only if that check succeeds. The download and paper destination paths
must be absent on the first invocation. The reconstructor enforces this and
refuses to merge with an existing directory.

```bash
python3 -u reconstruct_paper.py \
  --manifest "$REPRO_REPO/manifest.json" \
  --url-map "$REPRO_REPO/zenodo-urls.json" \
  --download-dir "$REPRO_REPO/zenodo-archives" \
  --output "$REPRO_REPO/Paper" \
  --verbose \
  2>&1 | tee "$REPRO_REPO/logs/02-reconstruction.log"
```

`set -o pipefail` preserves a failing Python exit status when output is piped
through `tee`. Check the command's exit status immediately if in doubt:

```bash
echo "$?"
```

Zero indicates success. Do not continue to environment creation after a
nonzero status or a `[FAIL-CLOSED]` message. The log shows the failing path,
network error or checksum mismatch.

### What the verbose command does, in order

1. Reads the local manifest and URL map, checks their structure and binding,
   validates archive paths and checks that the output does not exist.
2. Creates a private `.zenodo-archives.downloading/` stage with a receipt bound
   to this manifest and this exact URL map. On a resumed run it validates that
   receipt and the retained inventory first.
3. Processes canonical archives in the manifest's record/shard order. For each
   archive it logs the archive number and name, the chunk number/name, public
   URL, expected byte count and expected SHA-256.
4. Downloads each chunk over HTTPS. During longer transfers it prints byte
   progress about every ten seconds when bytes are arriving. A blocked network
   read may remain quiet until the socket timeout. The timeout defaults to 120
   seconds; it is not a limit on total download time.
5. Checks the chunk's exact size and SHA-256 before appending verified bytes to
   the canonical ZIP partial. It reports how many verified bytes have been
   assembled. Transient failures are retried up to five attempts per chunk
   with increasing delays, and each retry is logged.
6. Verifies the complete canonical ZIP's size and SHA-256 after concatenation.
   The stage retains completed ZIPs and at most one canonical partial and one
   transient chunk. It does not retain an extra permanent copy of all chunks.
7. After all ten ZIPs pass, installs `zenodo-archives/` and removes the download
   receipt. The final cache contains the logical record subdirectories
   `software_core/`, `main_figure_data/`, `supplementary_figure_data/` and
   `figs3_wsi/` with their canonical `.zip` files.
8. Creates a private `.Paper.reconstructing.PID/` directory. Rechecks each
   canonical archive before extraction, including when using a local cache.
   Large archive hashing also prints periodic byte progress.
9. Validates every ZIP member's inventory, order, name, size, file type and
   POSIX mode. Extracts each file while hashing its bytes, checks its SHA-256,
   and logs `VERIFIED MEMBER` with the file's path and size. Large members
   print intermediate extraction progress. The thousands of member messages
   are intentional and are saved in the reconstruction log.
10. Checks the final directory/file inventory and modes against the manifest.
    Only after all checks pass does it install the staged directory as
    `Paper/`. The stage is renamed; a second complete extracted copy is not
    needed for this installation.

Successful completion ends with two lines identifying the **verified archive
root** and the **exact reconstructed Paper tree**. At this point, the files
have been reconstructed and verified; the scientific analyses have not yet
been rerun.

### Download duration and interruptions

The release transfers about 240.25 GB. Ideal network-only times are about 53.4
hours at 10 Mbit/s, 5.34 hours at 100 Mbit/s and 32 minutes at 1 Gbit/s.
Zenodo throttling, retries, hashing, disk writes, assembly and extraction add
time. Keep the machine awake and use a persistent terminal such as `tmux` if
working over SSH.

If the network or process stops **before `zenodo-archives/` exists**, leave
`.zenodo-archives.downloading/` untouched. Keep `Paper/` and `zenodo-archives/`
absent and rerun the same reconstruction command. Use `tee -a` on the resumed
invocation if you want to append to the previous log. The downloader rehashes
completed ZIPs and committed chunk boundaries, discards only an uncommitted
partial tail, and resumes without downloading already verified chunks.

Do not edit the manifest or URL map between attempts, and do not rename,
modify or add files inside the hidden stage. Changed bindings, symlinks,
unexpected files and corrupted retained bytes cause a closed failure.

If **`zenodo-archives/` exists but `Paper/` does not**, downloading already
finished. Resume at extraction using the local-archive mode; rerunning the
URL-map command would correctly refuse the existing download directory:

```bash
python3 -u reconstruct_paper.py \
  --manifest "$REPRO_REPO/manifest.json" \
  --archive-root "$REPRO_REPO/zenodo-archives" \
  --output "$REPRO_REPO/Paper" \
  --verbose \
  2>&1 | tee -a "$REPRO_REPO/logs/02-reconstruction.log"
```

Extraction restarts in a fresh stage and checks all archives again. It does
not merge into an unfinished or existing `Paper/`. A force-killed or interrupted extraction
can leave `.Paper.reconstructing.PID/`; that directory is not a successful
output and is not reused. Account for its occupied space before restarting.
Never delete it while its process is still running.

If **`Paper/` already exists**, continue with environment setup and the
reproduction check below. Do not rerun reconstruction over it. Preserve an
existing reproduction and use a separate clone if you need another pristine
manifest reconstruction.

### Reusing a canonical archive cache already downloaded elsewhere

This optional path avoids another 240 GB transfer. Supply a directory with the
four logical record subdirectories and their canonical ZIPs. It must be
separate from `Paper/`; a folder of unassembled chunk files is not sufficient.
The following variable is the only location you need to adapt for this route:

```bash
ARCHIVE_CACHE=/absolute/path/to/existing/canonical-archives
python3 check_reproduction.py --root "$REPRO_REPO" --stage download \
  --archive-root "$ARCHIVE_CACHE"
python3 -u reconstruct_paper.py \
  --manifest "$REPRO_REPO/manifest.json" \
  --archive-root "$ARCHIVE_CACHE" \
  --output "$REPRO_REPO/Paper" \
  --verbose \
  2>&1 | tee "$REPRO_REPO/logs/02-reconstruction.log"
```

The space check credits existing archives and checks their sizes. The
reconstructor independently verifies every archive hash and extracted file;
it does not trust the cache's name or an earlier receipt. The resulting paper
still lives inside this clone.

## 4. Create the two scientific environments

Conda/Miniforge/Miniconda must be installed first. The bootstrap Python used to
download files is not the scientific environment. Make your Conda executable
available to both setup scripts and the launcher:

```bash
conda --version
export CONDA_EXE="${CONDA_EXE:-$(type -P conda)}"
test -x "$CONDA_EXE"
cd "$REPRO_REPO/Paper"

./Fig1/00_create_conda_envs.sh \
  2>&1 | tee "$REPRO_REPO/logs/03-main-environment.log"
./scripts/setup/01_create_s3_environment.sh \
  2>&1 | tee "$REPRO_REPO/logs/04-s3-environment.log"
```

This preserves the absolute `CONDA_EXE` set by Conda shell initialization.
Otherwise `type -P conda` locates its executable on `PATH`, including when
`conda` is also a shell function. `test -x` must succeed. If it does not, set
`CONDA_EXE` to the absolute path of your installation's `bin/conda`. Both
setup scripts accept that executable path. Stop and resolve an
environment command's nonzero status before continuing.

The main environment, `paper_apotome_repro`, is defined in
`Paper/scripts/setup/environment.yml` (Python 3.10 and pinned scientific,
microscopy, Cellpose/PyTorch and figure-rendering dependencies). The separate
`paper_apotome_s3_repro` uses Python 3.12.7, libvips, Poppler, PyArrow and the
pinned Figure S3/review packages. Its Linux x86-64 setup uses
`environment_s3_conda_explicit_linux-64.txt` plus
`requirements_s3_pip_lock.txt`; other platforms use `environment_s3.yml` and
are less extensively tested.

**Required v1.0.0 replay correction: OpenCV 4.10.0.84.** The immutable main
archive's YAML pins `opencv-python-headless==4.12.0.88`. An end-to-end test from
this clone found that OpenCV 4.12 rasterizes the accepted Figure 3 HIL polygons
differently: 12 Water pixels and 31 Allulose pixels differ from their recorded
masks, stopping the renderer with `HIL mask does not match its recorded polygon`.
OpenCV 4.10.0.84 reproduces both accepted masks with zero differing pixels using
the same NumPy 2.2.6. The separate S3 environment already uses OpenCV 4.10.0.84.

After creating the main environment, apply the explicit correction provided by
this GitHub companion:

```bash
conda run -n paper_apotome_repro python -m pip install --no-deps \
  -r "$REPRO_REPO/requirements-replay-compatibility.txt"
```

This changes only the main environment's OpenCV package. Keep the archived
YAML, manifest, HIL masks and HIL receipts intact; the requirements file records
the runtime correction separately. The checker expects this corrected pin,
checks the other package pins against the archived YAML, and actually
rasterizes both polygons before declaring the machine ready. A version/import
check alone did not catch this release incompatibility. If you use a custom
main interpreter, run the same pip command through that interpreter.

The scripts create an absent named environment or check an existing one.
They preserve an existing environment unless explicitly given their
`--force` option. Reusing an existing environment is not proof that all its
versions match the effective replay pins; use the supplied specifications for a fresh
installation and retain its setup logs.

Record portable interpreter paths without assuming the authors' username or
Conda installation directory:

```bash
export PAPER_PYTHON="$(conda run -n paper_apotome_repro python -c 'import sys; print(sys.executable)')"
export PAPER_PYTHON_CZI="$PAPER_PYTHON"
export PAPER_S3_PYTHON="$(conda run -n paper_apotome_s3_repro python -c 'import sys; print(sys.executable)')"
export PAPER_REVIEW_PYTHON="$PAPER_S3_PYTHON"
```

The setup scripts also print the interpreter paths. You may set these variables
to those absolute paths directly. For custom environment names, use the
matching names in the commands and set `PAPER_ENV` / `PAPER_S3_ENV` as
appropriate. Do not prepend the S3 environment to `PATH` globally; the launcher
scopes it to the figure steps that need it.

**Required numerical runtime profile: 64 library threads.** The frozen Figure 5
CNN also depends on the CPU numerical execution settings. Merely matching the
Python/package versions and random seed does not ensure the same training
trajectory. A complete test with 16 library threads rendered all ten figures,
but the final reference-image gate correctly rejected seven Figure 5 PNGs
(the master and the English/Spanish C, D and H panels). The original cell
counts, QC results and morphology clusters matched; CNN-derived values differed.

A controlled first-epoch comparison using the same inputs and PyTorch 2.12.1
identified the missing constraint:

| Numerical library threads | Training loss | Validation loss | Comparison with reference training |
| ---: | ---: | ---: | --- |
| 8 | 0.7644083776797677 | 0.6287438972839606 | Different |
| 16 | 0.7640612057135058 | 0.6216919878405228 | Different |
| 64 | 0.7645012029834801 | 0.6550887733401917 | Exact first-epoch agreement |

Load the companion's explicit profile in the **same Bash terminal** used for
the check and full reproduction:

```bash
source "$REPRO_REPO/reproduction-runtime.env"
```

This sets `OMP_NUM_THREADS=64`, `MKL_NUM_THREADS=64` and
`OPENBLAS_NUM_THREADS=64`. The checker requires these settings and verifies that
the actual main interpreter reports 64 PyTorch threads. Do not lower the
thread count to speed up this frozen replay. **64 library threads is not a
64-core hardware minimum**; a machine with fewer available cores can schedule
those threads but may take longer. Different CPU implementations can still
change floating-point execution, so the final reference-image comparison
remains the proof of a matching run. The immutable launcher, CNN code, accepted
inputs and reference images are preserved.

The portable Figure S3 replay uses the released L0/L2 exports, original organ
segmentation masks and completed inference/QC receipts. It recomputes the
paper statistics and figures from those released inputs. It does not repeat
external MDS conversion or retrain/rerun the external HistoPLUS GPU pipeline;
those source stores, model weights and CUDA installation are not required for
this replay. This is the supported reproduction route documented here.

## 5. Check that this clone is ready to reproduce

Return to the clone root for the new checker:

```bash
cd "$REPRO_REPO"
python3 check_reproduction.py --root . --stage reproduce --strict \
  2>&1 | tee "$REPRO_REPO/logs/05-reproduction-preflight.log"
```

This recalculates storage for an already extracted tree, checks the system,
explicit numerical thread profile and both interpreter paths, imports the scientific modules in each environment,
reports Python/pip version drift against the released specifications plus the
documented OpenCV replay correction, checks Figure 3 polygon-to-mask agreement
with the actual rendering interpreter, and verifies
the manifest SHA-256 and executable modes of
the shipped Python/shell scripts, and runs the published launcher with
`--check-only` **inside this clone's `Paper/`**. The launcher checks scientific
imports, bundle structure, script mirrors, raw inputs, the 24 L0/L2 WSI pairs,
Figure S3 replay/Parquet dependencies and accepted Figure 3/S5 HIL inputs.
Version differences produce warnings (blocking with `--strict`); missing or
broken imports fail. The verbose output is streamed into the log. Checks can take several minutes
and may prepare review/cache files; they do not render the ten figures.

Success must end with `[PASS] READY TO REPRODUCE`. This is still a readiness
result, not proof of a completed analysis. If it fails, fix the reported
missing input, environment, command or space condition and run the check
again. Do not suppress a failure and proceed as if the paper had passed.

For the launcher's direct readiness command and a readable execution plan:

```bash
cd "$REPRO_REPO/Paper"
./reproduce_all_figures.sh --check-only --output "$PWD"
./reproduce_all_figures.sh --plan --output "$PWD" --force
```

## 6. Reproduce all ten figures inside the cloned repository

From the reconstructed `Paper/` directory, run:

```bash
cd "$REPRO_REPO/Paper"
set -o pipefail
./reproduce_all_figures.sh --output "$PWD" --force \
  2>&1 | tee "$REPRO_REPO/logs/06-reproduction.log"
```

Here `--output "$PWD"` refers to
`allulose-sucrose-hypothalamus/Paper/`, not the small clone root. The launcher
requires its own reconstructed Paper directory as the canonical output. It
will refuse another destination. `--force` here authorizes installation of
the rebuilt, frozen-reference-validated Figure 3 as well as the other nine
figures. It does not mean “ignore validation failures.”

The launcher performs the following sequence:

1. Rechecks dependencies, repository structure, raw inputs and accepted HIL
   receipts, then snapshots raw-input metadata and active script hashes.
2. Creates a fresh `/tmp/apotome_full_rebuild_...` stage. Figures 2 and 3 also
   use their own isolated internal computation stages under `/tmp`.
3. Rebuilds Figure 1 (single-bottle behavior), Figure 2
   (design/anatomy/microscopy/c-FOS), Figure 3 (c-FOS/NPY), Figure 4
   (POMC/c-FOS), Figure S5 (ACTH/CLIP), Figure 5 (GFAP/Iba1), Figure S1,
   Figure S2, Figure S3 (organ WSI replay) and Figure S4, in that order.
4. Validates each staged figure, including the frozen Figure 3 reference,
   English multipanel masters, English/Spanish isolated panels and legends.
   Rasters are exported at 600 dpi with PDF companions. The progress counter
   advances only after a figure's staged output validates.
5. Compares every generated publication PNG byte for byte with the shipped
   reference images. A completed 10/10 rendering counter alone is not success;
   this final comparison can still fail, as the CPU-thread test demonstrated.
6. Only after all ten staged figures and their reference images pass, installs the generated analyses,
   figures, panels, legends, source data and provenance into the canonical
   `Paper/analyses/` and `Paper/Fig*/` locations.
7. Verifies installed artifact hashes against the staged outputs, audits raw
   data and active scripts again, and runs the final repository/bundle checks.
   Successful runs clean up their temporary build stage. Previous generated
   analyses, when present, remain recoverable at the `/tmp` location printed
   by the launcher.

Some figure analyses emit `[WARNING]` messages about using an accepted HIL
receipt from a relocated stage, or about stitch-registration JSON byte-hash
drift. A clone necessarily has different absolute paths from the original
review workstation. The released code accepts these cases only after checking
the recorded image SHA-256, dimensions, section identity and other stated
registration/review invariants. The messages explicitly list the matching
invariants; they are not permission to ignore an image/mask checksum failure.
Keep the original receipts intact. A traceback, `[FAIL]`, `[FAIL-CLOSED]`,
nonzero exit status or stopped figure progress still requires diagnosis.

### Follow detailed progress during slower figure steps

The launcher prints `Temporary validation root: /tmp/apotome_full_rebuild_...`.
Some analyses also write a more detailed log inside that stage. In particular,
Figure 5's CPU CNN prints epoch progress to
`Fig5/analysis/keyence20x_pipeline.log`; the main terminal can be quiet while
training is still progressing. Once Figure 5 has started, open a **second
terminal**, set `REPRO_REPO` to your clone's absolute path, and follow that log:

```bash
REPRO_STAGE="$(sed -n 's/^Temporary validation root: //p' "$REPRO_REPO/logs/06-reproduction.log" | tail -n 1)"
test -n "$REPRO_STAGE"
tail -f "$REPRO_STAGE/Fig5/analysis/keyence20x_pipeline.log"
```

The file appears when the Figure 5 analyzer starts. `Ctrl-C` in this second
terminal stops only `tail`; keep the original reproduction terminal running.
After a successful installation, the same analysis log is preserved at
`Paper/analyses/Fig5/results/human_final_run/keyence20x_pipeline.log`.

The Figure 5 launcher explicitly selects `--cnn-device cpu`. On this test
server, the CUDA-enabled PyTorch wheel printed a warning about an older
installed NVIDIA driver while CPU training continued. That warning alone does
not require a driver upgrade for this replay; check the recorded execution
device, analysis progress and final validation status.

The released README's reference clean-room run on the authors' server took
44 minutes 4 seconds with the CPU-fallback CNN and roughly 16 cores. This is
an observed reference, not a runtime promise or a minimum hardware claim.
Slower CPUs, constrained RAM, busy disks or network storage can make it take
several hours; the Zenodo download is additional. Do not interrupt the final
installation stage. If a figure fails, read its error and the printed stage
location in `logs/06-reproduction.log`, resolve the cause, and rerun the
launcher; it rebuilds stages rather than treating an incomplete run as valid.

## 7. Recognize and inspect a successful reproduction

The final launcher output must include:

```text
[PASS] Reproduction workflow completed.
Raw data and active scripts are present and unchanged.
All ten figures were installed, including certified Figure 3 (--force).
```

It also prints the canonical output path inside your clone. A successful
`git clone`, metadata validation, download, environment creation or readiness
check alone is not this result.

Inspect `Paper/Fig1/` through `Paper/Fig5/` and `Paper/FigS1/` through
`Paper/FigS5/` for the publication outputs. Each figure's `README.txt` gives
its specific filenames, panels, source-data tables, legends and provenance;
`Paper/README.txt` describes the complete delivered layout. Working analyses
are under `Paper/analyses/`, and this guide's logs stay under the clone's
`logs/`. Keep the archive cache if you want to reconstruct another pristine
tree without another download.

For an additional post-run check, from `Paper/`:

```bash
"$PAPER_PYTHON" scripts/utilities/03_check_reproducibility.py
"$PAPER_PYTHON" scripts/utilities/01_validate_bundle.py
```

### Completed cloned-directory validation

A full local integration run on **6 September 2026** rebuilt and installed all
ten figures in this clone using the OpenCV correction and 64-thread runtime
profile above. **All 136 publication PNGs matched the shipped references byte
for byte.** The final installed-artifact checks, raw-input/script audit and
bundle checks passed. The figure replay took **36 minutes 58 seconds** on the
test server; reconstruction and downloads are additional.

That run reconstructed all 6,350 manifest files from an existing canonical
archive cache, checking every archive and extracted member again. Separate
live tests downloaded and SHA-256 verified one chunk from each physical Zenodo
record (951,907,127 bytes total). The existing scientific environments were
reused with a local main-environment overlay; a fresh Conda installation on an
empty machine was not part of that test. These results establish the tested
replay, while the system checks and final reference-image gate remain required
on another machine.

Both must report `OVERALL: PASS`. The original manifest identifies the
pristine release reconstruction. Generated outputs and provenance can change
when the scientific workflow reruns; use the launcher's staged/installed
artifact checks and final validators to assess that run rather than demanding
that all regenerated files still match the original archive manifest.

## Notes

- Figure 5 has complete HIL-reviewed ARC/ME/VMN anatomy, but the current
  GFAP/Iba1 microglial-state result remains classifier-provisional. Cells and
  sections are not biological replicates.
- Figure S3 HistoPLUS phenotype labels are exploratory model-predicted outputs,
  not validated mouse cell identities.

These qualifications are part of the package and must be retained when its
figures or outputs are reused.

## License

This is a mixed-license release:

| Material | License |
| --- | --- |
| Executable source code and software-support files | MIT |
| Authors' original documentation, manuscript/thesis materials, figures, legends, masks, provenance, and data | Creative Commons Attribution 4.0 International (CC BY 4.0) |

See `LICENSE`, `LICENSES/MIT.txt`, and `LICENSES/CC-BY-4.0.txt` for the
scope and terms. A file-specific or third-party notice takes precedence.
Third-party works are not relicensed. Cite the applicable Zenodo DOI(s) above
when reusing the package.
