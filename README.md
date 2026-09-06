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

Clone this lightweight GitHub package, then validate it before downloading:

```bash
git clone https://github.com/cfarkas/allulose-sucrose-hypothalamus.git
cd allulose-sucrose-hypothalamus
```

```bash
python3 validate_repository.py --root .
sha256sum manifest.json
```

The first command must report `[PASS]`; the second must print
`df2ebb51339fdc3cc84dd8bbfb232647570d7e2cd0f12c19a5538ff7b7dd071d`.

## Download and reconstruct

Requirements are Python 3.10 or newer, HTTPS access to Zenodo, and a
case-sensitive filesystem with enough free space. Linux is the reference and
supported full-reproduction platform.

Choose a large-volume parent directory. On the first invocation, the download
and output paths below must not already exist; this is a deliberate overwrite
guard.

```bash
REPRO_PARENT=/absolute/path/on/a/large-volume
mkdir -p "$REPRO_PARENT"
test ! -e "$REPRO_PARENT/zenodo-archives"
test ! -e "$REPRO_PARENT/Paper"

python3 reconstruct_paper.py \
  --manifest manifest.json \
  --url-map zenodo-urls.json \
  --download-dir "$REPRO_PARENT/zenodo-archives" \
  --output "$REPRO_PARENT/Paper"
```

This single command downloads every public Zenodo chunk in URL-map order,
verifies each chunk, reassembles and verifies all ten canonical ZIP archives,
reconstructs into a private staging directory, verifies the complete tree, and
only then installs `Paper/` atomically. Expect:

- archive downloads: 240,253,512,366 bytes (240.25 GB; 223.75 GiB);
- reconstructed `Paper/` payload: 319,524,773,963 bytes (319.52 GB; 297.58 GiB);
- downloads plus final tree: 559,778,286,329 bytes (559.78 GB; 521.33 GiB).

Those figures exclude Conda environments and temporary figure-build stages.
Budget substantially more space than the combined minimum; about 700 GB free
is a practical starting point, not a guarantee. The figure launcher stages
work under `/tmp`, so that filesystem must also have ample free space.

The URL-map-v2 downloader is resumable and fail-closed. Until every canonical
ZIP is complete, the requested download directory remains absent and verified
progress is retained in its sibling
`$REPRO_PARENT/.zenodo-archives.downloading`. After a network interruption or
stopped process, leave that hidden directory in place, keep both requested
target paths absent, and rerun the exact same command. The rerun checks a
receipt bound to the exact manifest and URL map, rehashes completed archives
and segment boundaries, discards only an uncommitted mid-segment tail, and
continues without fetching already verified bytes.

Do not edit, rename, or add files to the hidden staging directory. Any changed
binding, unexpected file, symlink, corrupt completed segment, or corrupt
archive is rejected instead of trusted. At most one canonical-ZIP partial and
one transport chunk are retained beyond the already completed archives, so an
interruption does not accumulate multiple transient copies of the same archive.

## Create the pinned environments

Install Miniforge, Miniconda, or Conda first. From the reconstructed Paper
root, create the main environment and the separately pinned Figure S3
replay/review environment:

```bash
cd "$REPRO_PARENT/Paper"
./Fig1/00_create_conda_envs.sh
./scripts/setup/01_create_s3_environment.sh
```

On Linux x86-64, Figure S3 uses its explicit Conda and pip locks. Other
platforms fall back to the pinned YAML and are less extensively tested. The
portable reproduction uses the checksum-bound Figure S3 L0 exports from the
logical WSI package together with the matched L2 exports and receipts from the
supplementary record. It does not require the external MDS source store,
HistoPLUS checkpoint, or GPU pipeline.

## Reproduce every figure

First run the non-writing readiness check:

```bash
./reproduce_all_figures.sh --check-only --output "$PWD"
```

Then rebuild and validate Figures 1–5 and S1–S5:

```bash
./reproduce_all_figures.sh --output "$PWD" --force
```

`--force` is required to install the newly rendered, frozen-reference-
validated Figure 3. Every figure is rendered in isolated temporary storage and
validated before the generated analyses and publication artifacts are
installed into their canonical locations. The launcher audits raw inputs and
active scripts before and after the run.

The measured clean-room replay on the authors' server took 44 minutes 4
seconds with the CPU-fallback CNN and roughly 16 CPU cores. That is a reference
measurement, not a runtime guarantee. Downloading 240+ GB may dominate, and
other hardware may require several hours. Do not interrupt the final
transactional installation stage.

Successful completion ends with all ten figure packages passing. Detailed
per-figure inputs, outputs, validation, and provenance are documented in
`Paper/README.txt` and each `Paper/Fig*/README.txt`.

## Scientific-status notes

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
