APOTOME PAPER — REPRODUCTION RUN-BOOK
=====================================

README.txt and README2.txt are documents, not programs. Read this file with:

  less README.txt

Quick start
-----------
From the directory containing this file:

  ./reproduce_all_figures.sh --check-only
  ./reproduce_all_figures.sh --plan
  ./reproduce_all_figures.sh --output "$PWD" --force

The launcher finds the paper_apotome_repro Python environment, validates the
active repository and renders every established figure at 600 dpi into fresh
temporary stages. A tqdm bar advances only after each staged render passes.
Only after all eleven figures validate are the current analyses and publication
artifacts installed into their canonical analyses/ and FigN/ locations. The
Paper folder itself remains the complete delivered package: active scripts,
raw/raw_data, HIL inputs, source data, figures, panels, legends and provenance.
No top-level apotome_rebuild_* directory is produced. Figures 2 and 3 retain
their mandatory additional isolated computation stages under /tmp.
Canonical installation atomically rewrites every generated master, isolated
panel, legend, source-data file and provenance artifact, even when its bytes
match the previous version. Before the temporary stage can be removed, SHA-256
verification must prove that every installed artifact is byte-for-byte identical
to the output produced by that run.

The optional --output argument is deliberately fail-closed and accepts only the
absolute path of this Paper directory:


If the environment is missing:

  ./Fig1/00_create_conda_envs.sh
  ./reproduce_all_figures.sh

For a nonstandard interpreter location:

  PAPER_PYTHON=/absolute/path/to/python ./reproduce_all_figures.sh

Current numbering
-----------------
  Fig1    behavior, experiment 1
  Fig2    design/anatomy/microscopy/c-FOS
  Fig3    c-FOS/NPY; frozen-reference protected
  Fig4    POMC/c-FOS ARC/ME
  Fig5    GFAP/Iba1 microglia
  FigS1   single-bottle males
  FigS2   single-bottle females
  FigS3   multi-organ H&E WSI kidney/liver/spleen; completed HIL
  FigS4   behavior preference
  FigS5   ACTH/CLIP after 16 h fasting and first-ever exposure
  FigS6   microglial classifier benchmark: morphometric proposals versus DINOv2

There is no active Figure 6 or Figure 7. The multi-organ WSI package
was built as Figure 7, renumbered to Figure 6 on 2026-08-28 and moved into the
supplementary series as Figure S3 on 2026-08-29. Inserting it at S3 shifted
behavior preference from S3 to S4 and ACTH/CLIP from S4 to S5. GFAP/Iba1 is
Figure 5; it was called Figure 6 only while the July ACTH/CLIP experiment
temporarily occupied the main Figure 5 slot. Historical frozen receipts may
retain their original numbers and hashes. Figure S3 is complete; its finalized
HIL receipts are validated before its downstream artifacts are rebuilt.

Language outputs
----------------
Every active multipanel master is English-only. Each renderer exports matching English and
Spanish isolated subpanels, per-panel legends where present, and full-figure
legends in both languages. Retired multipanel masters are not part of the active
figure tree.

Canonical installation
----------------------
Every normal run stages and validates all eleven figures before installing current
outputs into this Paper tree. To include the frozen/certified Figure 3:

  ./reproduce_all_figures.sh --output "$PWD" --force

Without --force, Figure 3 is rebuilt and verified against its frozen reference
but its promoted bytes remain unchanged. The compatibility option --promote is
accepted but no longer needed.

The canonical installer never writes raw_data, raw, scripts, hil_review, legacy
or README.txt. It audits all raw file metadata and all active script hashes both
before and after installation. Superseded generated analyses are moved outside
Paper to a recoverable /tmp path; no new in-tree legacy figure copy is created.
Completed Figure S3 is validated, rebuilt and installed with the other figures.

Supplementary Figure S5
-----------------------
S5 contains nine animals: Water, Sucrose and Allulose n=3 each. Eight are WT;
NPY-M is a male NPY-transgenic Water animal included by author confirmation on
8 September 2026. Its confirmed 546/633 channels identify ACTH/CLIP and c-FOS. ACTH/CLIP objects below the acquisition-specific DAPI
nuclear-area fifth percentile are excluded as likely non-cell fragments. Raw
masks are unchanged.

The master uses six microscopy/cartoon panels in a 2-column × 3-row block, G/H
stacked in a wide right column, and I/J below. G contains ME and ARC; H is
ARC-only with the graph centered. I/J include occurrence maps with a thin light
outer envelope from accepted DAPI-only ARC/ME/VMN anatomy masks and a bold white
HIL 3V contour selected as the observed medoid for each condition. Native
per-acquisition 3V masks exclude nuclear centroids before shell construction and
spatial inference, and density is forced to zero inside the displayed 3V. In
total, 433 DAPI, 0 c-FOS-positive and 0 double-positive objects are excluded.
The masks never define the outer tissue crop or normalization box. I/J also
include six equal-DAPI-density radial-shell profiles, exact exhaustive
condition-label PERMANOVA, BH q values, and separate colorbars.
The multipanel master is English-only. English and Spanish isolated A-J panels
are 600 dpi.

Reproduce only S5:

  export PAPER_PYTHON=/home/server/anaconda3/envs/paper_apotome_repro/bin/python
  ./FigS5/04_run_figure_s5.sh

Check its accepted HIL review without writing:

  ./FigS5/04_run_figure_s5.sh status
  "$PAPER_PYTHON" scripts/utilities/06_import_figure_s5_hil.py --source-only

Check the immutable dedicated HIL DAPI audit:

  ./FigS5/04_run_figure_s5.sh spatial-status

The renderer hash-validates that HIL review, uses each raster only as a negative 3V
exclusion mask, and draws one observed condition-medoid point contour. The masks
never become positive tissue/crop geometry. See
FigS5/provenance/SPATIAL_BOUNDARY_SEMANTICS_CORRECTION.json.

See FigS5/README.txt for endpoints, source tables and provenance.

Supplementary Figure S6
-----------------------
S6 compares existing morphometric proposals with frozen DINOv2 image features
using the same 203 human-reviewed cells from 15 animals. The master reports
balanced accuracy (73.9% versus 57.0%), class recall and confusion matrices.
Intervals resample whole animals and condition on fixed predictions. The DINOv2
classifier was evaluated with nested whole-animal validation. The baseline is
a fixed unsupervised comparator; this is internal validation on a selected
review sample. Figure 5 publication outputs remain unchanged.

The S6 package has frozen, hash-bound prediction/reference tables, an English
600-dpi master, English/Spanish A-D panels and legends, and reproducible source
code. It is included in the all-figure launcher. Rebuild S6 alone:

  "$PAPER_PYTHON" FigS6/01_make_figure_s6_microglia_classifier.py \
    --output-dir /tmp/figs6_benchmark_fresh --dpi 600

See FigS6/README.txt for the complete caption, data and input provenance.

Figure 5 caution
----------------
Figure 5 has complete 61-section HIL ARC/ME/VMN anatomy, but its current
GFAP/Iba1 master remains classifier-provisional. Publication-final use requires
independently HIL-reviewed microglial-state labels or a separately validated frozen
classifier checkpoint. Cells and sections are not biological replicates. See
Fig5/README.txt.

Repository checks
-----------------
These commands read and validate; they do not run figure analyses:

  "$PAPER_PYTHON" scripts/utilities/03_check_reproducibility.py
  "$PAPER_PYTHON" scripts/utilities/01_validate_bundle.py

Expect OVERALL: PASS from both. The first checks syntax, dependencies, sibling
references and byte-identical FigN/scripts-FigN mirrors. The second checks the
publishable tree shape and self-containment.

Licensing
---------
Software code and software-support files are released under the MIT License.
The authors' original documents, figures and data are released under Creative
Commons Attribution 4.0 International (CC BY 4.0). See the root LICENSE file
and LICENSES/ for the exact scope and notices. A file-specific or third-party
notice takes precedence for that file.

Safety and provenance
---------------------
- Read incident.txt before destructive maintenance. The root launcher uses
  absent /tmp stages, validates all figures before canonical installation and
  never treats raw data as a generated replacement target.
- Raw CZI/TIFF files and Cellpose masks are read-only.
- HIL receipts are versioned and hash-bound. Frozen receipts are not rewritten
  merely to update an historical path.
- Historical manuscript completion/HIL receipts retain their pre-sanitization
  DOCX hashes. The Office-package sanitation transition is recorded without
  private review metadata in
  manuscript_sources/PUBLIC_RELEASE_OFFICE_SANITIZATION_20260902.json. The
  segmentation-only Figure S3 transition is recorded in
  manuscript_sources/PUBLIC_RELEASE_MANUSCRIPT_FIGS3_REVISION_20260904.json.
  The authoritative current thesis DOCX/PDF identities after correcting the
  Figure 1/S1/S2 single-bottle terminology are recorded in
  manuscript_sources/PUBLIC_RELEASE_SINGLE_BOTTLE_CONSUMPTION_REVISION_20260904.json.
- A reviewer process must be restarted after its source code changes.
- Active numbered scripts at FigN/ and scripts/FigN/ must remain byte-identical.
- All publication raster outputs are 600 dpi; PDF companions are also produced.

Package layout
--------------
  FigN/              active scripts, promoted graphics, raw/source data,
                     panels, legends, provenance and figure-specific README
  scripts/FigN/      byte-identical mirrors of active numbered scripts
  scripts/shared/    shared HIL/staging helpers
  scripts/utilities/ validation, HIL import and promotion tools
  analyses/          local generated working state; not required in an upload

For an archive, omit machine-local analyses/, recovery_archive/, recovery_support/,
quarantine/, docs/, .claude/, .mypy_cache/, __pycache__/, README2.txt, and
development-only *.orig/*.rej patch remnants. Downloaded third-party
full_text_open_access.pdf article copies are also omitted because an open-access
URL alone is not treated as redistribution evidence. Normal bibliographic
references remain in the manuscript, but the unrelated
literature_webscrap_30_08_2026/ working scrape and unused exploratory inputs are
omitted too. Microscopy tile-stitching receipts used to reconstruct Figures 2–4
remain because they are required to reconstruct microscope acquisitions.
Figure S3 uses its original within-slide organ segmentation masks directly.
Superseded figure graphics are not active publication files. Figure 2's active
HIL readiness helper is
Fig2/09_rebuild_hil_inputs.py, mirrored at scripts/Fig2/. Source-path entries
in manifests are provenance, not runtime dependencies.

The detailed reconstruction history and later corrections are in incident.txt.

Microglial choices (2026-09-07)
-----------------------------
The normal reproduction workflow reuses the current 203 saved human microglial
choices and candidate without prompting. The optional flag
  ./reproduce_all_figures.sh --do_microglial_choices
opens a separate blank 300-cell mask/class review and continues after retraining.
Original choices remain intact. The new data are included directly in GitHub
under Fig5/microglial_review_data/. The review/candidate outputs are retained in
analyses/Fig5/results/current_microglial_choices/. Published Figure 5 and S6
retain their established comparison inputs. See Fig5/README.txt for details.
