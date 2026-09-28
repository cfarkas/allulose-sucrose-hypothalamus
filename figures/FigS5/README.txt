FIGURE S5 — ACTH/CLIP-ASSOCIATED c-FOS AFTER FIRST ORAL EXPOSURE

Figure S5 examines oral water, sucrose or allulose after a 16-hour fast in
one-year-old mice without sugar familiarization. ACTH/CLIP F-3 staining
(sc-373878, 1:100) identifies a POMC-derived peptide signal. Tissue was
postfixed without vascular perfusion.

The cohort has three animals per condition: eight wild-type animals and
NPY-M, a male NPY-GFP animal in Water. NPY-M channels are DAPI, 488 NPY-GFP,
546 ACTH/CLIP and 633 c-FOS. NPY-GFP is separate from the two endpoint markers.

Panels
  A–F  Microscopy, segmentation and ACTH/CLIP-associated cellular regions.
  G    c-FOS/DAPI fractions in ME and ARC.
  H    c-FOS-positive fractions among ACTH/CLIP-positive cells in ARC.
  I–J  Spatial c-FOS and double-positive occurrence across radial shells.

Every animal has accepted DAPI-based anatomy and a third-ventricle contour.
The ventricular mask excludes nuclei and displayed density within the lumen;
the outer tissue envelope comes from ARC/ME/VMN anatomy. ACTH/CLIP objects
below the acquisition-specific fifth percentile of DAPI nuclear area are
excluded as likely fragments.

Spatial analysis uses six equal-DAPI-count radial shells and exact PERMANOVA
over 1,680 condition allocations. Benjamini–Hochberg correction covers the two
spatial endpoints. These shells describe spatial sampling, not matched
anatomical subdivisions. Animals are independent units; cells and shells are
within-animal measurements. The wild-type subset provides a sensitivity analysis.

Files
  Figure_S5.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/            Individual A–J panels in English and Spanish, PNG and PDF.
  legends/           Figure and panel captions in both languages.
  source_data/       Animal measurements, shell profiles and statistical tests.
  hil_review/        Accepted anatomy and ventricular contours.
  provenance/        Channel assignments, cohort inputs and validation records.

Reproduce from Paper/ with the paper environment:
  ./FigS5/04_run_figure_s5.sh reproduce
This creates a fresh analysis directory, imports the accepted anatomy,
quantifies the nine-animal cohort, and renders the figure and bilingual panels.
Inspect the generated outputs before installing them into the figure directory.

Render a quantified result set into an absent directory:
  FIGS5_OUTPUT_DIR=/tmp/figs5_render ./FigS5/04_run_figure_s5.sh render
Inspect the ventricular review status with:
  ./FigS5/04_run_figure_s5.sh spatial-status

The complete cohort analysis is in analyses/FigS5/results/including_NPY_M_20260908/.
Raw microscopy and masks are read-only during rendering.
