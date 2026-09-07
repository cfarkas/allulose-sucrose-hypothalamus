# Figure scripts and supplementary Figure S6 — 7 September 2026

The current scripts and small figure inputs are directly available in
[`figure_updates/`](figure_updates/). This directory has the same layout as the
reconstructed `Paper/` directory. Its 11-figure launcher includes S6.

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

## Reproduce all eleven figures

With the same computer requirements as the [main guide](README.md), run:

```bash
python3 reproduce.py
```

This downloads/reuses the original Zenodo release, authenticates the source
update, installs its files with backups, and invokes the eleven-figure launcher.
The original `manifest.json`, archive chunks, hashes and DOIs remain fixed.
`figure-updates.json` separately binds each added/replaced file by hash and mode.
Unrecognized local edits stop installation before any files are replaced; use a
fresh clone for the update if you have edited the reconstructed publication files.
Backups are saved under `.figure-update-backups/`. Repeat the same command to
resume. An updated `Paper/` uses the normal command; `--figure-updates` is retained as an
explicit alias for the current workflow.

Use `python3 reproduce.py --archived-release` to reproduce the archived
10-figure release in an unmodified reconstruction. The previously reported
136-PNG full replay applies to that archived release. The new S6 and corrected
Figure 1 have separate render checks; this update does not claim a newly completed
full eleven-figure replay.

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
