PAPER SCRIPTS — CURRENT MAP
===========================

Updated: 2026-08-30

Start here
----------
From the Paper root:

  export PAPER_PYTHON=/home/server/anaconda3/envs/paper_apotome_repro/bin/python
  ./reproduce_all_figures.sh --check-only
  ./reproduce_all_figures.sh --plan

A normal rebuild renders only into fresh /tmp stages until all thirteen figures pass.
It then installs current generated analyses and publication artifacts into the
canonical Paper/analyses and Paper/FigN locations. Paper itself is the delivered
complete package, so active scripts and raw/raw_data remain present and are
audited before and after installation. No nested apotome_rebuild_* output is
created. Figures 2 and 3 retain their additional isolated temporary computation
stages because their certified runners require them. Pass --force to install the
validated certified Figure 3 rebuild.
Every generated English master, bilingual isolated panel and bilingual legend
is atomically rewritten from the successful stage, including byte-identical
outputs. A staged-to-canonical SHA-256 audit must pass for every artifact before
the launcher can report success and remove the stage.


Active figure map
-----------------
  scripts/Fig1/    Figure 1 behavior
  scripts/Fig2/    Figure 2 microscopy/HIL
  scripts/Fig3/    Figure 3 c-FOS/NPY
  scripts/Fig4/    Figure 4 POMC/c-FOS
  scripts/Fig5/    Figure 5 GFAP/Iba1 (classifier gate still pending)
  scripts/FigS1/   Supplementary Figure S1
  scripts/FigS2/   Supplementary Figure S2
  scripts/FigS4/   Supplementary Figure S4 behavior preference
  scripts/FigS5/   Supplementary Figure S5 ACTH/CLIP, bilingual 600 dpi
  scripts/FigS6/   Supplementary Figure S6 microglial classifier benchmark, bilingual 600 dpi

There is no active Figure 6 or Figure 7 package. The GFAP/Iba1 package is
Figure 5 and was called Figure 6 only while the July ACTH/CLIP experiment held
the main Figure 5 slot. The multi-organ WSI package became Figure S3
on 2026-08-29, which shifted behavior preference to S4 and ACTH/CLIP to S5.
Historical frozen receipts may retain the old numbers by design. Figure S3 is
complete and has no scripts/FigS3 mirror; its front scripts, completed receipts,
downstream analysis and bilingual figure are included in reproduction checks.

Shared entry points
-------------------
  scripts/shared/01_annotate_regions.py
      HIL ARC/ME/VMN reviewer for Figure 4 and Figure 5.
  scripts/shared/02_validate_raw_cellpose_masks.py
      Read-only raw/Cellpose coverage validation.
  scripts/shared/03_safe_stage_fig45.py
      Fail-closed no-delete staging for Figures 4 and 5.
  scripts/shared/05_annotate_microglia.py
      Independent Figure 5 microglial-state HIL reviewer.
  scripts/shared/spanish_panel_text.py
      Shared English-to-Spanish panel text support.

Utilities
---------
  scripts/utilities/01_validate_bundle.py
      Read-only publishable-tree validation.
  scripts/utilities/03_check_reproducibility.py
      Syntax, dependency, sibling and front/mirror checks.
  scripts/utilities/04_promote_figure.py
      Atomic all-file install plus staged-to-canonical byte-for-byte verification.
  scripts/utilities/05_refresh_promoted_timestamps.py
      Compatibility refresh using the same atomic installer and SHA-256 check.
  scripts/utilities/06_import_figure_s5_hil.py
      Hash-validating import of the accepted S5 HIL set into fresh work.
  scripts/utilities/08_figure_progress.py
      Tqdm progress display for the thirteen validated figure packages.
  scripts/utilities/09_audit_canonical_inputs.py
      Pre/post raw-data metadata and active-script SHA-256 audit.

Scientific readiness
--------------------
- Figure 4: accepted HIL anatomy complete; see Fig4/README.txt.
- Figure 5: 61/61 HIL anatomy reviews complete; classifier provenance remains the
  publication-final gate; see Fig5/README.txt.
- Figure S5: eight accepted DAPI-only fields, Water n=2, Sucrose n=3,
  Allulose n=3, NPY-M excluded; see FigS5/README.txt.
- Figure 3 remains frozen/certified and is installed only by the explicit
  root-launcher option --force after frozen-reference validation.
- Figure S4 still requires its documented workbook for a fresh analysis.

All publication raster figures are produced at 600 dpi. Superseded figure
graphics are not active outputs.

Figure S7 documents human-supervised Cellpose training; Figure S8 presents
pilot-based sample-size curves. Both use verified source inputs and complete
English masters with bilingual isolated panels. See FigS7/README.txt and
FigS8/README.txt. Figure S3 uses corrected matched-HC3 residual permutations
(99,999, base seed 1707); hepatic epithelial-like p=0.00031, BH q=0.02418.
