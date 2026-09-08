SUPPLEMENTARY FIGURE S5 — FIRST ORAL EXPOSURE AFTER A 16-HOUR FAST

The experiment asks whether oral water, sucrose or allulose preferentially
recruits satiety-associated POMC-related cells in food-deprived mice. Animals
were one year old and had no prior familiarization with the sugars. ACTH/CLIP
F-3 (Santa Cruz Biotechnology, sc-373878; 1:100) identifies the POMC-derived
peptide signal. This cohort was postfixed without vascular perfusion.

The author-requested inclusion on 8 September 2026 adds male NPY-M, an NPY-GFP
transgenic animal in Water. The full cohort is Water/Sucrose/Allulose n=3/3/3:
eight WT animals plus NPY-M. The author confirmed NPY-M channels as DAPI,
488 NPY-GFP, 546 ACTH/CLIP and 633 c-FOS. NPY-GFP is not used as the c-FOS or
ACTH/CLIP channel. The earlier exclusion was a WT-only preparation filter,
not a sex-based or documented image-quality exclusion.

Original microscopy, channel exports, native Cellpose masks and the eight
previously accepted anatomical annotations are preserved. The author reviewed
NPY-M ARC/ME/VMN on DAPI. Each animal also has a separately accepted DAPI-only
third-ventricle contour. Its negative mask excludes nuclear centroids and
intraventricular displayed density; it never defines the outer tissue crop.
The outer envelope comes from the accepted ARC/ME/VMN anatomy masks.

Analysis and display
--------------------
A–F: representative microscopy, segmentation and cellular ACTH/CLIP regions.
G: c-FOS/DAPI fractions in ME and ARC. H: c-FOS-positive fractions among
ACTH/CLIP-positive cells in ARC. I/J: animal-level spatial c-FOS and double-
positive occurrence, using six equal-DAPI-count radial shells and exact
PERMANOVA over 1,680 condition assignments. BH correction covers the two
spatial endpoints. Regional global and pairwise tests are supplied separately.
The animal is the independent unit; cells and shells are measurements within
animals. The original eight-WT results are retained as a sensitivity reference.
ACTH/CLIP objects below the acquisition-specific fifth percentile of DAPI
nuclear area are excluded as likely fragments, with all raw masks retained.

Reproduction from the Paper root
--------------------------------
  export PAPER_PYTHON=/home/server/anaconda3/envs/paper_apotome_repro/bin/python
  ./FigS5/04_run_figure_s5.sh reproduce

This command creates a fresh analysis directory, applies the confirmed manifest
overlay, imports accepted anatomy, quantifies and renders at 600 dpi. It creates
one complete English master and English/Spanish individual panels and legends.
It does not overwrite previously quantified results or promote automatically.

To render an already quantified result set into an absent directory:
  FIGS5_OUTPUT_DIR=/tmp/figs5_candidate_new ./FigS5/04_run_figure_s5.sh render

To check or open the third-ventricle review:
  ./FigS5/04_run_figure_s5.sh spatial-status
  ./FigS5/04_run_figure_s5.sh spatial-review 34250

Accepted inputs and current analysis
-----------------------------------
  provenance/inclusion_NPY_M_20260908/  confirmed overlays and WT-only reference
  hil_review/human_final_including_NPY_M_20260908/  nine anatomy decisions
  hil_review/spatial_3v_including_NPY_M_20260908/    nine ventricular decisions
  analyses/FigS5/results/including_NPY_M_20260908/  quantified current cohort

The original human_final_20260828_v1 and spatial_tissue_20260828_v1 review sets
remain immutable. The historical OUTER_TISSUE label in the latter already
represented a negative ventricular mask; the new metadata migration preserves
its exact points and mask bytes and records their provenance.

Scripts
-------
  01_export_czi_zstacks.py                 original channel export
  02_july_cfos_acth_clip_human_in_loop.py   frozen historical analyzer
  02a_july_cfos_acth_clip_human_in_loop_s5.py  portable quantification and HIL
  03_make_figure_s5_acth_clip_cfos.py        full figure and bilingual panels
  04_run_figure_s5.sh                       prepare/review/build/render/reproduce
  05_review_spatial_tissue_hil.py           DAPI-only ventricular review
  06_include_npy_m.py                      author-confirmed manifest overlay

Figure_S5.pdf/png are the complete English publication figure. Individual
English/Spanish A–J panels are in panels/, captions in legends/, numerical
measurements and statistics in source_data/, and hash-bound receipts in
provenance/. The raw data and original masks are not changed by rendering.
