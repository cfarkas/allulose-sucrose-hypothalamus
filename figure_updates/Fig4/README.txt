FIGURE 4 — POMC/c-FOS IN ARC AND ME, WITH SCALE-RESOLVED SPATIAL PANELS (25 SEPTEMBER 2026)

Figure_4.pdf and its 600-dpi PNG are the complete English master.
Individual A–J panels and their legends are available in English and Spanish.

Panels A–H retain the original microscopy, masks, quantitative summaries and
standalone panel files, and still report ARC and ME separately.

Panels I and J are scale resolved. For every animal the number of pairs of
activated cells closer than r is counted within each acquisition and tissue
side, for r from 20 to 150 micrometres, and expressed in standard deviations
from that animal's own conditional random-labelling null. The null holds every
cell position, stratum and activated-cell count fixed and permutes only which
cells are activated, so its mean and variance are exact and the field geometry
cancels between the observed and the expected count; no edge correction enters.
The independent unit is the biological cage, so animal curves are averaged with
equal weight within a cage. Conditions are compared with a studentized
maximum-deviation global envelope test (Myllymäki et al. 2017, J. R. Stat. Soc.
B 79:381) over an exhaustive enumeration of cage-label allocations. Taking the
maximum over distances controls the error across the whole curve, so the single
p-value printed on each panel needs no further correction across distances;
Holm then covers the two endpoints. This is the analysis Figure 3 panels H and
I carry; only the population, the region and the experimental unit differ.

ME is not analysed in I and J. It retains 74 marker-associated cells with 28
activated across 11 animals, and only 6 of those animals have two or more
activated cells, which is too sparse to define a curve. Panels G and H still
report both regions.

The earlier scalar spatial summaries, core enrichment inside the 80% POMC
density contour and the clustered fraction on the six-neighbour graph, are no
longer drawn. Their tables are unchanged in source_data/, and the panels they
produced are archived with their comparison receipts under
../analyses/fig4_scale_resolved_20260925/archive/Fig4/.

Rebuild and publish Figure 4 from the Paper root:

  ./scripts/run_fig4_scale_resolved.sh

That script re-derives the cage curves and the envelope test from the archived
POMC cell table, re-renders the bilingual I/J panels, re-assembles the master
and installs it here. Without --publish on the final step, outputs stay in
analyses/fig4_scale_resolved_20260925/assembled/Fig4 and nothing here changes.

Scale-resolved methods, per-cage curves, envelope bands and test results:
  ../analyses/fig4_scale_resolved_20260925/README.md
Contour and cluster methods, parameters, condition summaries and tests:
  ../analyses/spatial_shape_20260920/README.md
Historical shell source tables and earlier README:
  ../analyses/spatial_shape_20260920/archive/Fig4/
The original spatial_reanalysis_20260919 package remains the input reference.
