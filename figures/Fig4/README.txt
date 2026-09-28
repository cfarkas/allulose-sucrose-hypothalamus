FIGURE 4 — POMC/c-FOS LABELING IN ARC AND ME AND SPATIAL CLUSTERING

Figure 4 examines POMC-associated c-FOS labeling after oral exposure in the
familiarized cohort following a four-hour fast.

Panels
  A–C  Representative DAPI, c-FOS and POMC signals with reviewed anatomy.
  D    Cellpose segmentation and cellular assignments.
  E–F  Microscopy and marker-associated signal displays, including NPY-GFP.
  G    Animal-level c-FOS/DAPI fractions in ME and ARC.
  H    Animal-level c-FOS-positive fractions among POMC-positive cells in ME
       and ARC.
  I    Scale-resolved clustering of c-FOS-positive cells in ARC.
  J    Scale-resolved clustering of POMC-associated c-FOS-positive cells in ARC.

Panels G/H show animal means and SDs, ordinary one-way ANOVA and exact
pairwise comparisons. IQR flags do not remove animals. POMC objects smaller
than the mean local DAPI nuclear area are excluded by a uniform size rule;
threshold sensitivity results accompany the source tables. Microscopy displays
measured channel intensities inside accepted marker regions.

For I/J, positive-cell pairs separated by less than r are counted within each
acquisition and tissue side for r from 20 to 150 µm. Conditional random labeling
fixes cell positions, strata and positive counts, with exact null means and
variances. Animal curves express clustering in standard deviations from those
nulls. This comparison conditions on observed field geometry and uses no edge
correction; it does not correct acquisition or detection differences.

The independent unit for I/J is the biological cage. Animal curves receive
equal weight within cages. Conditions are compared through exhaustive cage-label
allocations using a studentized maximum-deviation global envelope over the
whole distance range. Holm adjustment covers the two endpoints. Total c-FOS
uses 14 animals in 2/3/3 Water/Sucrose/Allulose cages; the POMC-associated
endpoint uses 13 animals in 2/2/3 cages. Sparse ME double-positive labeling
does not support the corresponding spatial comparison. G/H include both regions.

Files
  Figure_4.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/           Individual A–J panels in English and Spanish, PNG and PDF.
  legends/          Figure and panel captions in both languages.
  source_data/      Regional values, spatial inputs, cage curves and tests.
  provenance/       Size-filter, geometry and numerical validation records.

Recompute neuronal statistics from the GitHub repository without microscopy:
  python figures/scripts/shared/recheck_current_statistics.py

Stage I/J and assemble the figure from the reconstructed Paper/ directory:
  ./scripts/run_fig4_scale_resolved.sh
Install the validated output into Fig4/ with:
  ./scripts/run_fig4_scale_resolved.sh --publish

The runner reads the packaged POMC cell table and acquisition manifest,
recomputes cage curves and envelope tests, renders bilingual I/J panels, and
assembles the English master with A–H. It needs no GPU or raw microscopy.
Staged outputs are in analyses/fig4_scale_resolved_20260925/assembled/Fig4/.
