FIGURE 3 — NPY-ASSOCIATED c-FOS LABELING AND SPATIAL CLUSTERING

Figure 3 compares 11 animals: 5 Water, 3 Sucrose and 3 Allulose.

Panels
  A    Representative NPY-GFP and c-FOS microscopy, merges and enlarged fields.
  B–D  DAPI nuclear maps with NPY/c-FOS assignments for the three conditions.
  E    c-FOS-positive nuclei as a fraction of DAPI nuclei.
  F    c-FOS/NPY double-positive nuclei as a fraction of NPY-positive nuclei.
  G    Positive-cell pair-counting illustration in Allulose field E8_FR6-4.
  H    Spatial clustering of total c-FOS-positive nuclei.
  I    Spatial clustering of NPY-associated c-FOS-positive nuclei.

Abundance fractions pool counts within each animal and weight animals equally.
The analysis includes an exactly enumerated Welch omnibus statistic, a
Brown–Forsythe variance test, exact studentized contrasts, Hodges–Lehmann shifts
and Cliff's delta. Panel E/F values and comparison definitions are provided in
source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv. IQR flags do not
exclude animals from inference.

Panels H/I count positive-cell pairs within acquisitions and tissue sides at
distances of 20–150 µm. Conditional random labeling fixes cell positions,
strata and positive counts. Each curve expresses observed clustering in
standard deviations from its own null. A studentized maximum-deviation global
envelope compares conditions across the distance range; Holm adjustment covers
the two endpoints. Distances measure proximity, without establishing anatomical
homology or connectivity. The independent unit for these comparisons is the
animal; cage membership is undocumented.

Total c-FOS analysis uses all 11 animals. NPY-associated curves require at least
two double-positive nuclei and use 3/3/3 Water/Sucrose/Allulose animals.
WATER_NPY3 has zero and WATER_NPY4 one such nucleus; both contribute abundance
measurements. WATER_NPY5 has two and contributes to the spatial comparison.

WATER_NPY4–5 and the sugar animals share a Zeiss LSM 780 and 25×/0.8 objective.
The former images use a 2×2 mean reduction to match the sugar images' analysis
pitch of 0.664213 µm. The water group also includes Leica acquisitions.
Physical calibration does not remove staining, detection or sampling differences.

Files
  Figure_3_cFos_NPY.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/       Individual A–I panels in English and Spanish, PNG and PDF.
  legends/      Figure and panel captions in both languages.
  source_data/  Plotted values, tests and portable cell-coordinate inputs.
  provenance/   Cohort, acquisition and figure validation records.

Recompute neuronal statistics from the GitHub repository without microscopy:
  python figures/scripts/shared/recheck_current_statistics.py

Rebuild and install the figure from a complete local Paper/ analysis tree:
  ./scripts/run_fig3_new_water_cohort.sh
The runner requires analyses/fig3_water_extended_20260924 and the accepted
microscopy, mask and panel inputs. --from-raw also extracts CZI channels and
runs GPU Cellpose segmentation. The downloadable Zenodo data do not include
the raw WATER_NPY4–5 acquisitions; GitHub supplies their reviewed cell-coordinate
inputs for numerical reproduction. Fig3/07_rebuild_complete_figure3.sh routes
to the same cohort workflow.
