FIGURE 2 — ORAL-CHALLENGE DESIGN AND REGIONAL HYPOTHALAMIC c-FOS LABELING

Panels
  A  Experimental timeline for the NPY-GFP/C57BL/6J design.
  B  Investigator-reviewed ARC, ME and VMN anatomy over Cellpose DAPI masks,
     with DAPI-associated c-FOS-positive nuclei highlighted.
  C  Representative native CZI microscopy on a common 600 × 1200 µm field.
  D  Animal-level c-FOS/DAPI fractions in ME, ARC, VMN and non-ROI tissue.

Panel B uses accepted anatomical polygons and saved nuclear segmentation.
Region membership and c-FOS assignments are checked against the source counts,
mask areas and label images. Water E10_FR1-1, Sucrose E7_FR7-1H and Allulose
E8_FR6-4H provide the representative fields in B and C.

Panel C uses calibrated 16-bit native scene-0 mosaics. Adjacent tiles are
registered by integer translations, with common coordinates for DAPI and
c-FOS. Microscopy keeps its native aspect ratio and uses 200 µm scale bars.
Display normalization, registration and channel identity are recorded in the
source tables. Display framing and contrast do not define the quantitative ROIs.

Panel D gives equal weight to animals. Bars show the mean and sample SD;
points identify animals. IQR flags are descriptive and do not exclude animals.
Statistics use ordinary one-way ANOVA and exact two-sided rank-sum comparisons.
The renderer independently checks the statistics against the source tables.

Files
  Figure_2.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/           Individual A–D panels in English and Spanish, PNG and PDF.
  legends/          English and Spanish figure captions.
  source_data/      Animal values, tests, field geometry and displayed counts.
  provenance/       Output hashes, display assets and validation records.
  raw_data/         CZI images, segmentation companions and anatomy inputs.

GitHub supplies the figure graphics, scripts and source tables. The raw
microscopy and accepted masks are supplied by the Zenodo reconstruction.
Panel A has English and Spanish timeline artwork. Figure rendering uses saved
segmentation masks and accepted anatomy; it does not run neural-network inference.

Reproduce from the reconstructed Paper/ directory:
  Fig2/07_run_all.sh /tmp/fig2_stage /tmp/fig2_render

Both arguments must be absolute, absent paths under /tmp. The runner validates
the raw bundle, stages animal statistics and microscopy, checks the anatomy,
then renders the English figure and bilingual panels. Input hashes, field
calibration, channel identity and statistical checks must pass. Outputs remain
in the selected render directory for inspection. Raw microscopy is read-only.
