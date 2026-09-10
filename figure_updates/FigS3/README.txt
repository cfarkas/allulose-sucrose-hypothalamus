FIGURE S3 — KIDNEY/LIVER/SPLEEN WHOLE-SLIDE ANALYSIS
=====================================================

Purpose
-------
This directory is a self-contained, receipt-driven Figure S3 workflow for the
24 Motic H&E whole-slide images in:

  /media/server/STORAGE/Motic_AnatomiaPatologica_2025/
  HE_Liver_Spleen_Kidney_Nancy_2026

Each WSI contains kidney, liver, and spleen. The source directory is read-only;
raw files are referenced and fingerprinted, never renamed or edited. The
workflow exports scanner-native L0 and L2, segments the three organs by computer
vision, selects an equal number of spatially balanced L0 tiles within each
segmented organ, runs HistoPLUS inference, computes model-independent
H&E/nuclear features, performs animal-level statistics, and renders
Figure_S3.png/PDF plus panels and source data.

Control definition
------------------
AGUA and AGUA_Sin_cerebro are both standardized to Water. Specifically,
m26-016 is one of the seven Water controls. "Sin cerebro" is retained only as
sample_note=brain_not_collected; it is not a fourth treatment, exclusion, or
stratum. The cohort is Water n=7, Sucrose n=9, Allulose n=8.

Run or resume
-------------
From the Paper root:

  FigS3/06_run_figure_s3.sh

The run is resumable. Existing TIFF or inference results are reused only from
completed receipts; failed slides are not converted into zeros. Environment
variables allow only scoped operational changes:

  FIGS3_EXPORT_WORKERS=2 FIGS3_GPU_BATCH_SIZE=6 \
    FigS3/06_run_figure_s3.sh

The tested interpreters are:

  /home/server/anaconda3/bin/python
  /home/server/anaconda3/envs/lazyslide311/bin/python

The gated, locally authorized HistoPLUS weight is expected at:

  /home/server/.cache/histoplus/histoplus_cellvit_segmentor_20x.pt

Scripts
-------
00_check_environment.sh
  Verifies 24 MDS sources, Python dependencies, CUDA, GPU identity, and the
  local model weight checksum.

01_prepare_metadata.py
  Reads Hoja1 from EXCEL_ORGANOS.xlsx, forward-fills only the grouped slide
  fields, validates exactly three organs per animal, preserves raw Spanish
  labels, creates safe IDs (m26-001 through m26-024), and writes the 24-row
  slide and 72-row animal-organ manifests.

02_export_mds_levels.py
  Opens the OLE compound document read-only, decodes only DSI0 JPEG pixel
  streams, and writes tiled lossless-deflate BigTIFF for L0 and L2. Each receipt
  distinguishes the scanner's true dimensions from tile-padded dimensions and
  records source/output SHA-256, source and output MPP, tile geometry, and
  scanner metadata. L0 is 0.261780 micrometres/pixel; L2 is 1.047120.

03_segment_organs.py
  Builds a stain-aware tissue mask from each scanner overview, estimates the
  principal specimen axis, partitions tissue into three spatial groups, anchors
  the terminal hematoxylin-rich spleen, and applies the histology-audited
  liver-kidney-spleen placement order (version figs3-three-organ-cv-1.1). Raw
  H&E/stain/edge features and identity confidence remain auditable. It performs
  within-slide segmentation only; no cross-slide image alignment is performed.

04_run_histoplus_gpu.py
  Converts each CV organ mask to WSI coordinates and uses LazySlide's physical-
  resolution tile model at 0.5 micrometres/pixel. It proposes deterministic
  farthest-point tiles plus reserves, checks tissue on the actual L0 pixels,
  and retains exactly 24 tiles per organ. Tile coordinates and analyzed tissue
  area are exported; sampled counts are never extrapolated to whole organs.
  Classical hematoxylin/eosin, focus, watershed nuclear density, area,
  eccentricity, and solidity are computed. HistoPLUS runs on CUDA in float32;
  AMP is intentionally disabled because the installed LazySlide/OpenCV path
  cannot postprocess float16 HoVer maps. Cells, polygons, model classes,
  probabilities, morphology, tiles, and organ assignments are retained.

05_analyze_make_figure.py
  Uses the animal as biological n. Per organ, components must be detected in
  at least 18/24 animals and contain at least 24 total objects to enter CLR,
  PCA, or Water-reference deviation analyses. The tumor-specific raw model
  class is excluded from classification/composition displays and all downstream
  composition analyses, while unaltered raw outputs remain archived. Filtered
  contrasts are marked not tested. Eligible HistoPLUS counts receive a 0.5
  pseudocount and centered-log-ratio transform. Treatment effects use a
  prespecified source-cohort adjustment, HC3 confidence intervals, and 4,999
  Freedman-Lane residual permutations. Benjamini-Hochberg q values are reported
  globally within endpoint family and within organ. Water-reference deviations
  use shrinkage Mahalanobis distance; Water controls are evaluated leave-one-out.
  Classical CV features form a model-independent sensitivity analysis. Paired
  L0 cell insets are generated only for the strongest globally BH-significant
  phenotype contrast, using true model-output polygons and a receipt-driven
  median-proximity selection rule. Independently, the representative overview
  is the Water animal closest to the median whole-slide tissue fraction, with
  sample ID as an exact-tie breaker; this selects m26-015. The compact master
  omits its global title. Panels are lettered in reading order on a 15.0 x 12.3
  inch canvas: workflow/organ Panel A beside Water-deviation Panel B, a
  full-width horizontal histology Panel C, ordination Panel D beside classical
  forest Panel E, and the bottom full row for composition Panel F. Panel F
  replaces the earlier heatmap with organ-specific cell-type CLR effect bars,
  HC3 95% confidence intervals, and high-contrast global-BH q-value stars.
  Panel E marks unadjusted Freedman-Lane permutation p values with stars while
  explicitly stating that no classical feature survives global BH correction.
  Isolated panels are exported in English and Spanish from the same live axes;
  the multipanel master stays English-only.

06_run_figure_s3.sh
  Runs the complete staged workflow and writes a timestamped log.

08_relabel_inference_organs.py
  One-time, idempotent migration used for this completed run after the full-
  resolution histology audit corrected Liver/Kidney identities. It verifies
  every old/new tile against corrected masks, preserves pixels and neural
  predictions, relabels organ/tile keys, refreshes hashes, and writes an audit
  table. Fresh runs use 03_segment_organs.py and do not require migration.

Key outputs
-----------
Figure_S3.png and Figure_S3.pdf
  The 600-dpi raster master and vector PDF.

exports/
  48 BigTIFF files (L0 and L2 for 24 slides), per-slide conversion receipts,
  and the aggregate export manifest.

organ_segmentation/
  Tissue masks, three-organ label masks, colored overlays, computer-vision
  receipts, and the cohort segmentation summary.

histoplus/
  Per-slide candidate/selected tile manifests, direct-pixel QC, classical CV,
  GeoParquet nuclei, phenotype summaries, representative overlays, inference
  receipts, and cohort aggregates.

source_data/
  Metadata, nomenclature audit, animal values, CLR values, PCA scores, all
  statistical contrasts, significant-cell inset manifest, displayed-cell
  barplot values, and model-independent sensitivity values.

panels/, legends/, provenance/
  Isolated panels A-F in English and Spanish (PNG and PDF at 600 dpi), the
  English and Spanish figure legends, an English and a Spanish legend for every
  panel, method sources, hashes, and run receipts. The Spanish panels are
  exported from the same live axes as the English ones through
  scripts/shared/spanish_panel_text.py, so geometry and plotted values are
  identical and only the words differ; the English to Spanish mapping is
  recorded in provenance/Figure_S3_spanish_translation_receipt.json.

Numbering and panel lettering
-----------------------------
This package was built as Figure 7 and renumbered to Figure S3 on 2026-08-28,
taking the slot left empty when the GFAP/Iba1 package returned to Figure 5. The
rename covers directory, file, script, receipt and manifest names as well as the
version tokens inside the receipts; only frozen historical receipts elsewhere in
the paper keep their original Figure 7 wording, by design. The superseded
Figure 7 master, panels and legend are preserved under
legacy/superseded_figure7_numbering_<timestamp>/.

Panels are lettered strictly in reading order, A through F. The earlier master
still carried the letters of two panels that had been merged away, so it ran
A, G, D, F, H, E down the page. The current letters are:

  A  workflow and within-slide organ segmentation
  B  cross-validated Mahalanobis deviation from Water, with ANOVA p values
  C  three-organ histology and the significant segmented-cell examples
  D  organ-specific CLR-PCA
  E  classical H&E forest plots
  F  model-predicted cell-type composition effects

Interpretation boundary
-----------------------
HistoPLUS was trained on human tumor H&E, not normal mouse kidney, liver, or
spleen. The 13 retained display labels are therefore reported as exploratory
model-predicted nuclear phenotypes (displayed as "-like"), not validated mouse
cell identities. The tumor-specific raw output is excluded from the figure and
composition inference and is not evidence of cancer in these tissues. The
manuscript should say
"sugar-associated deviations from Water," not "evoked," causal, or diagnostic
effects. HistoPLUS weights are CC-BY-NC-ND-4.0 and non-commercial; this workflow
is research-use only.

Method sources checked online
-----------------------------
TumorQuantAI: https://github.com/cfarkas/tumorquantai
LazySlide workflow: https://lazyslide.readthedocs.io/en/stable/reference/workflow-cheatsheet.html
LazySlide first analysis: https://lazyslide.readthedocs.io/en/latest/getting-started/first-analysis.html
LazySlide paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC13076205/
HistoPLUS: https://github.com/owkin/histoplus/
LazySlide HistoPLUS implementation/classes:
  https://github.com/rendeirolab/lazyslide-models/blob/main/src/lazyslide_models/segmentation/cellvit_family/histoplus.py
HoVer-Net: https://github.com/vqdang/hover_net
HoVer-Net paper: https://www.sciencedirect.com/science/article/pii/S1361841519301045

Statistical correction, 9 September 2026 (analysis 1.8)
The coefficient-specific Freedman–Lane null retains the other treatment
contrast and source cohort. Both the observed and permuted t statistics use
HC3 standard errors. The previous renderer compared an HC3 observed statistic
with ordinary-OLS permutation statistics and removed both treatment terms
under the null. Counts and effects are unchanged; permutation p and BH q
values were recomputed with 99,999 permutations (base seed 1707).
Hepatic epithelial-like Allulose/Water CLR ratio: 4.05987, 95% CI 2.22130–7.42021,
p=0.00031, global q=0.02418. This is relative phenotype composition, not an
absolute hepatocyte count. Six classical contrasts are nominal; none passes BH.
Numerical test: python FigS3/tests/test_hc3_permutation.py
Reference: Winkler et al., NeuroImage 92, 381–397 (2014),
https://doi.org/10.1016/j.neuroimage.2014.01.060
