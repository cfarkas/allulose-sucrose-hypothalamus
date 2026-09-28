PAPER REPRODUCIBILITY UTILITIES
===============================

Run the supported read-only checks from the Paper root:

  ./reproduce_all_figures.sh --check-only

Individual checks:

  "$PAPER_PYTHON" scripts/utilities/03_check_reproducibility.py
  "$PAPER_PYTHON" scripts/utilities/01_validate_bundle.py
  "$PAPER_PYTHON" scripts/utilities/06_import_figure_s5_hil.py --source-only

`01_validate_bundle.py` validates the active publishable figure set: Fig1-Fig5
and FigS1-FigS6. `03_check_reproducibility.py` checks numbered scripts and their
scripts/<figure> mirrors. Figure S3, the multi-organ WSI package renumbered from
Figure 7, has no scripts/FigS3 mirror and is validated directly from its active
front scripts and completed receipts.

Promotion is explicit and non-destructive:

  "$PAPER_PYTHON" scripts/utilities/04_promote_figure.py \
    --figure FigS5 --from /absolute/path/to/reviewed/render

Any differing promoted file is copied first to the figure's timestamped legacy
folder. The promoter never writes raw_data, raw, scripts, hil_review, legacy or
README.txt. Certified Figure 3 is rejected unless --force is explicit; the root
launcher exposes that override as --force.

In root --replace-current mode, every staged artifact is atomically rewritten,
including byte-identical files. The promoter verifies all source/destination
SHA-256 pairs, and the root launcher repeats a ten-package byte-for-byte audit
before deleting the stage.

`08_figure_progress.py` provides the launcher's tqdm counter. Each of its ten
units is credited only after one figure package renders and validates.

`05_refresh_promoted_timestamps.py` is a compatibility atomic refresh helper.
`06_import_figure_s5_hil.py` understands the documented historical Fig5-to-FigS5
relocation, verifies the frozen HIL receipt and every hash, and writes only into
an absent/empty fresh work directory.

`19_update_single_bottle_consumption_wording.py` performs the narrow final-thesis
update associated with Figure 1, S1 and S2 terminology. It rewrites only the
relevant document text and nine matching English-master/Spanish-panel image
members, then proves every other DOCX ZIP member stayed byte-identical.
