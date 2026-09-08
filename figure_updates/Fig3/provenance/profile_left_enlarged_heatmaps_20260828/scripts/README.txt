FIGURE 3 ACCEPTED-FINAL REPRODUCTION
====================================

This package contains the exact English/Spanish analyzers and renders a
receipt-bound Figure 3 candidate from the 75 manifest-bound raw inputs plus the
accepted Water/Allulose ventricle HIL pair. It does not use the
legacy scripts/Fig3/06_run_figure3_staged.sh, whose spatial-analyzer pin belongs
to the superseded pre-human build.

Use Python:
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python

Portable raw-project preparation
--------------------------------
Run only when MinKNOW is stopped or Paper I/O is otherwise approved; this makes
a physical, independent 75-file copy and hashes both sides:

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python \
    Fig3/reproduce_final_20260823/prepare_raw_project.py \
    --project-root /tmp/fig3_reproduction_raw_project_v1/Paper \
    --receipt /tmp/fig3_reproduction_raw_project_v1/RAW_PROJECT_PREPARATION.json \
    --copy-from Fig3/raw_data/legacy_experiment_2025_07_28

Exact fresh reproduction
------------------------
First accept both ventricle polygons in the live reviewer (default port 34247):

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python \
    Fig3/02a_review_cfos_npy_ventricles.py

The raw project, absent output, and accepted-HIL directory are three distinct
arguments. The runner copies and validates only the two masks and two receipts
inside its fresh output stage:

  cd /media/server/STORAGE/Apotome/Paper
  bash Fig3/reproduce_final_20260823/RUN_REPRODUCE_FINAL.sh \
    /tmp/fig3_reproduction_raw_project_v1/Paper \
    /tmp/fig3_reproduction_output_v1 \
    /media/server/STORAGE/Apotome/Paper/Fig3/hil_review/ventricle_cartoon_20260827_v1

For the already certified 2026-08-23 physical raw project, the tested first
argument is /tmp/fig3_human_accepted_raw_project_20260823_v1/Paper. The runner
rehashes all 75 files before analysis and refuses links, missing/extra files,
existing output, or a non-/tmp output.

Expected PASS
-------------
- 75/75 raw files hash-bound and nlink=1.
- Nine animals, exactly three Water, three Sucrose, and three Allulose.
- Accepted spatial analyzer SHA-256:
  67f93ed9cd9c115d43f57729f2de00bd052e86108afa6afdf60c0f981b8d4647.
- One English master PDF/PNG plus A-G isolated PDF/PNG in English and Spanish:
  exactly 30 graphics.
- English/Spanish A-G page and pixel geometry identical; raster resolution
  600 dpi; single-page PDFs, embedded fonts, no Type 3.
- Water and Allulose ventricle masks each validate against an accepted
  reviewer/session/timestamp receipt and the exact registered raw-DAPI source.
- D/E animal-level statistics and F/G exact 1,680-label spatial receipts remain
  bound to the accepted reference receipts.
- Accepted English compositor SHA-256:
  d855867632cdeb02440a7a83935531486522806760e8f8ffc09ddad00bed57b3.
- Panel A uses one 1.56-pt black (#000000) ROI stroke with accepted membership
  and counts, clipped to the microscopy-image extent before black letterboxing.
- F/G use 6.40 x 3.12-inch source canvases: the radial-profile axis is on the
  left at 3.152 inches, and the enlarged condition-means heatmap is on the right
  at 1.541 inches. Each native-aspect panel occupies 7.721 inches in the master;
  the D/E tick-clearance guard stays active.
- final/provenance/README_REPRODUCTION_PASS.json reports status PASS.

Never use --force, reuse an output, write into live Fig3, or replace mandatory
human-accepted Water3 masks. WATER_NPY3 is the authorized stable release
identifier; its unresolved true biological identity remains a provenance
caveat and is not silently changed here.
