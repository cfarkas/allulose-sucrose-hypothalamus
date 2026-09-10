# Ring-definition schematic — 10 September 2026

Figure **3G** alone shows DAPI nuclei around a schematic third ventricle, six rings with the inner two highlighted, and a single within-ring percentage formula. The third ventricle has a narrow superior neck and a flared, rounded base; DAPI points and ring colors have stronger contrast. Methods and caveats are in the legends. Inner rings 1–2 contain approximately one-third of eligible nuclei. The schematic uses illustrative DAPI nuclei; experimental measurements and inference are unchanged.

The NPY profile uses one field per animal. POMC rings are constructed within each reconstructed section, and corresponding-ring numerators and denominators are summed across sections before calculating each animal's six percentages. Equal nuclear counts do not imply equal ring areas, ring thicknesses or atlas-registered anatomy.

Figure 3 uses the available middle-row space: E/F at left, ring schematic G at right, and spatial panels H/I below. The English master figures contain 9 and 10 panels, respectively; Figure 4 refers readers to Figure 3G and has no duplicate ring cartoon; all individual panels and legends are available in English and Spanish. The shared [cartoon renderer](figure_updates/scripts/shared/spatial_ring_cartoon.py) recreates the tissue schematic and uses the same ring assignments as the analysis. The Figure 3 finishing step also renders current spatial-panel layouts from certified numerical tables, preserving the historical sealed analysis.

Checks include both affected spatial reruns, exact agreement with ten prior numerical tables, within-ring denominators and transforms, mathematical checks against both production shell functions, and complete figure output validation. NPY PERMANOVA p=0.010714, q=0.021429 and dispersion p=0.014286, q=0.028571 are unchanged. POMC spatial double-positive p=q=0.920080 is unchanged. A complete 13-figure raw-data replay was not required or performed for this presentation update. Current Word documents remain local; the submitted thesis was not modified.

Radial shell methods have published precedents. Mehta et al. (2013), [IMACULAT](https://doi.org/10.1371/journal.pone.0061386), used equal-area elliptical shells within individual nuclei. Our custom exploratory implementation uses approximately equal numbers of DAPI nuclei across a tissue field. No claim of first use or validation of the exact hypothalamic implementation is made.

# Current figure scripts and data — 9 September 2026

The current update contains thirteen figures: five main figures and S1–S8.
The complete masters remain English; isolated panels and captions are bilingual.
The current thesis and manuscripts remain local and are not uploaded.

- **Figure 4 POMC size QC:** all 76 reconstructed sections were reanalysed using the same local-DAPI area threshold across conditions. Of 3,403 processed POMC objects, 780 were excluded. In the familiarized, 4-hour-fast cohort, ARC c-FOS/POMC was higher under allulose than water (nominal exact p=0.047619; global p=0.176454; cage-mean p=0.200). Direction was retained across 0.5/1.0/1.5 thresholds, while nominal significance depended on the criterion. Original masks and images were retained; [methods and sensitivity](figure_updates/Fig4/provenance/POMC_SIZE_QC_20260909.md) are provided. S7 displays derived training masks after size QC, without retraining archived models.
- **S3 statistical correction:** observed and permuted treatment coefficients now use the same HC3 standard errors. Each coefficient-specific Freedman–Lane null retains the other treatment contrast and source cohort. The fixed 99,999-permutation analysis gives hepatic epithelial-like Allulose/Water CLR ratio 4.06 (95% CI 2.22–7.42), p=0.00031 and global BH q=0.02418. Six classical contrasts are nominal, including lower hepatic nuclear density; none passes BH. Images, sampling, effects and confidence intervals are unchanged. The [correction report](figure_updates/FigS3/provenance/STATISTICAL_CORRECTION_20260909.md) and complete previous/corrected comparison are included.
- **S7 training documentation:** [code and inputs](figure_updates/FigS7/) include original training-image crops, reviewed masks and recorded loss curves, bound to their sources by SHA-256. Human correction is distinguished from independent evaluation.
- **S8 sample-size planning:** [code and inputs](figure_updates/FigS8/) render the recorded pilot-based Monte Carlo power curves. The underlying simulation is [Figure 1's estimator](figure_updates/Fig1/03_estimate_sample_size_monte_carlo.py). These are planning estimates for the pilot outcomes.

Verification covers the independent per-permutation HC3 reference test, the full S3 data/figure gate, and byte-identical PNG replays for all S7/S8 masters and panels. The full thirteen-figure raw-data replay has not been rerun; the launcher and authenticated update installation are checked separately. Existing archived Zenodo hashes and records remain unchanged.

Repository checks: `python -m pytest tests -q` (59 passed). Targeted numerical checks: `python -m pytest figure_updates/Fig4/tests/test_pomc_size_qc.py figure_updates/FigS3/tests/test_hc3_permutation.py -q` (5 passed). The bundled full-data pipeline tests require the reconstructed Paper tree and its raw-data exports; they are run in that tree.

## 8 September update, retained


The 8 September update adds three completed changes:

- **Figure 1 and sex-specific S1/S2:** the humane transfer of FR5-4 from E2
  allulose to E10 water occurred after the day-6 measurement. Linked records
  identify one animal and retain its original group through that endpoint;
  E10/FR6-2 is water. The weight analysis includes 24 mice in 12 cages.
- **Figure 3:** the complete English A–H master now includes the sucrose
  microscopy and nuclear map. Its third ventricle was delineated by the
  investigator on native DAPI, and the tissue outline follows the procedure
  used for water and allulose. Quantified counts and statistics are unchanged.
- **Figure S5:** male NPY-M is included as an NPY-GFP transgenic water animal,
  giving three animals per condition. Its four-channel assignment and DAPI-only
  anatomy/third-ventricle contours are saved with the other eight reviews.
  Spatial c-FOS and double-positive profiles give exact PERMANOVA p=0.400
  and p=0.464286, respectively (BH q=0.464286 for both). The original eight-WT
  results are supplied as a sensitivity comparison. S5 tests first oral exposure
  after a 16-hour fast using ACTH/CLIP as a POMC-related signal.

The [Monte Carlo script](figure_updates/Fig1/03_estimate_sample_size_monte_carlo.py)
reconstructs sample-size scenarios from the day-6 values of 12 pilot cages.
With 500,000 simulations per size, group-specific SDs and Welch ANOVA, five
cages per group reach 80.16% power for bottle-volume disappearance and fourteen
reach 80.77% for cage-mean weight change. These pilot-conditional omnibus
estimates do not establish achieved power or neuronal/glial sample sizes.
[Source values and curves](figure_updates/Fig1/source_data/sample_size_monte_carlo/)
and an [execution receipt](figure_updates/Fig1/provenance/SAMPLE_SIZE_MONTE_CARLO_RECEIPT.json)
are included. See [Figure 1 instructions](figure_updates/Fig1/README.txt).

The update includes figure assets and analysis code. The current thesis and
manuscript remain local and are not distributed in this GitHub update.

The current scripts and small figure inputs are directly available in
[`figure_updates/`](figure_updates/). This directory has the same layout as the
reconstructed `Paper/` directory. Its current launcher includes thirteen figures (1–5 and S1–S8).

- [Figure 1](figure_updates/Fig1/): the renderer now says **“exact, Holm-adjusted”**;
  the exact comparisons and Holm adjustment are unchanged. Corrected reference
  graphics, analysis scripts, raw consumption/weight tables and captions are included.
- [Figure 5](figure_updates/Fig5/): cell/soma mask refinement, painting and human
  review, supervised CNN retraining, frozen DINOv2 features, nested animal-held-out
  evaluation, audit/next-review-queue scripts and verified encoder setup are included.
  These are candidate workflows; the published Figure 5 panels are unchanged.
- [Figure S6](figure_updates/FigS6/): renderer, frozen human references and predictions,
  validation partitions, source tables, English/Spanish panels, captions and provenance.
  [View the PDF](figure_updates/FigS6/Figure_S6.pdf).
- [Shared scripts](figure_updates/scripts/), the root reproduction launcher,
  [archive builder](figure_updates/release/build_archives.py), and source-export
  script include the new supplement. Figure-local scripts and their `scripts/`
  mirrors are identical.

S6 compares **73.9% versus 57.0% balanced accuracy** on the same 203 human-reviewed
cells from 15 animals. It includes paired animal-bootstrap intervals, class
recall and confusion matrices. DINOv2 preprocessing and classifier selection
exclude each held-out animal. The existing unsupervised morphometric proposals
used the full cell collection without fitting human class labels. This is a
selected internal review sample, not an independent biological validation cohort;
morphology alone does not establish functional microglial activation.

## Reproduce S6 without the large data download

From the Git repository, create an environment once:

```bash
python3 -m venv .figures-venv
.figures-venv/bin/python -m pip install -r requirements-figures.txt
```

Then render into a fresh output directory:

```bash
.figures-venv/bin/python figure_updates/FigS6/01_make_figure_s6_microglia_classifier.py --validate-only
.figures-venv/bin/python figure_updates/FigS6/01_make_figure_s6_microglia_classifier.py --output-dir figure-renders/S6 --dpi 600
.figures-venv/bin/python figure_updates/scripts/utilities/07_validate_figure_outputs.py --figure FigS6 --root figure-renders/S6
```

The renderer recomputes all scores and checks the frozen input hashes and animal
partitions. It needs no GPU, microscopy download, pretrained weights or retraining.
Use a different output directory for a subsequent render. Figure 1 can also be
rerun directly from the small package; see its [instructions](figure_updates/Fig1/README.txt).

## Reproduce all thirteen figures

With the same computer requirements as the [main guide](README.md), run:

```bash
python3 reproduce.py
```

This downloads/reuses the original Zenodo release, authenticates the source
update, installs its files with backups, and invokes the thirteen-figure launcher.
The original `manifest.json`, archive chunks, hashes and DOIs remain fixed.
`figure-updates.json` separately binds each added/replaced file by hash and mode.
Unrecognized local edits stop installation before any files are replaced; use a
fresh clone for the update if you have edited the reconstructed publication files.
Backups are saved under `.figure-update-backups/`. Repeat the same command to
resume. An updated `Paper/` uses the normal command; `--figure-updates` is retained as an
explicit alias for the current workflow.

Use `python3 reproduce.py --archived-release` to reproduce the archived
10-figure release in an unmodified reconstruction. The previously reported
136-PNG full replay applies to that archived release. The corrected Figure 1, complete Figure 3, nine-animal S5 and S6 have separate
checks. The S5 update was rebuilt from native channels, saved segmentations and
accepted anatomy/3V masks in a fresh analysis directory; all nine per-animal
regional values and new spatial statistics agree with the reviewed analysis.
This update does not claim a newly completed full thirteen-figure replay.

## Reuse the current choices or make your own

The normal command loads the **203 saved human choices** and their matching masks
and trained candidate, without opening a reviewer. These original choices retain
their exact SHA-256 and are never replaced by model predictions. The new data are
in [Fig5/microglial_review_data](figure_updates/Fig5/microglial_review_data/): an
open CSV of the 203 choices plus ten ZIPs (about 74 MB total) containing the 300
review crops/masks, 14,415 classifier input crops and soma masks, source tables,
and the saved CNN candidate. The original microscopy acquisitions have not changed.

To make independent choices during Figure 5, use the optional flag:

```bash
python3 reproduce.py --do_microglial_choices
```

Or, from a reconstructed Paper tree:

```bash
./reproduce_all_figures.sh --do_microglial_choices
```

At the Figure 5 step, the terminal displays `MICROGLIAL_REVIEW_URL`. Open that
address in a browser. The new session has **300 targets and zero prefilled human
choices**. Review the cell/soma masks, paint corrections if needed, then select
classes. Finish after 200–300 choices with class/animal coverage; at 300 it starts
training automatically. The workflow continues after training succeeds. On a
remote server, forward the printed localhost port through SSH to your computer.

New choices are saved separately and do not alter the default 203-label set.
The review directory is printed and can be resumed with the helper's
`--review-dir` option. Candidate results and the label receipt are retained in
`Paper/analyses/Fig5/results/current_microglial_choices/`. This review/training
branch is separate from the fixed published Figure 5 panels and S6 benchmark;
fresh labels do not silently redefine a published accuracy value.

The review inputs can also be used without the microscopy download, with the
scientific paper environment and `Fig5/12_microglial_choices.py`; see the
[Figure 5 commands](figure_updates/Fig5/README.txt). DINOv2 retraining additionally
uses the verified upstream encoder obtained by `11_prepare_dinov2_assets.py`.

## Zenodo release impact

The live records checked on 7 September still describe the archived release.
The new review data and S6 tables are available directly on GitHub, so updating
Zenodo is not required to run this workflow now. For an updated DOI archive of
these new files, the affected records would be:

| Record | Changed content |
| --- | --- |
| [Software/core 22265239](https://zenodo.org/records/22265239) | Corrected Figure 1, Figure 5 code and reviewed campaign, S6 code/graphics, shared scripts and provenance |
| [Supplementary data 22265243](https://zenodo.org/records/22265243) | Frozen S6 reference/prediction tables and validation design under `FigS6/raw_data/` |

The main raw-data and Figure S3 whole-slide-image payloads do not change for this
update. New versions should retain links to those existing records and include a
fresh archive manifest and download map. The portable microglial training bundle is now included in the GitHub update
and would belong in the next core archive. Local `analyses/` working artifacts
remain excluded from the archive contract. The historical release plan is retained as a record of v1.0.0.

Zenodo [recommends versioning for significant file changes](https://help.zenodo.org/docs/deposit/manage-files/#modify-files-after-publication).
Each [new version has its own persistent identifier](https://help.zenodo.org/docs/deposit/manage-versions/)
linked to earlier versions. **No new Zenodo version was published by this GitHub update.**

Code remains MIT-licensed; original data/figures remain CC BY 4.0. The separate
DINOv2 asset downloader retrieves upstream Apache-2.0 source and weights directly
from their public hosts and verifies their frozen hashes.

Local document work is now stored in `thesis_and_manuscript/`. Input resolvers accept this name and the earlier dated folder; both document-workspace names remain excluded from public bundles.
