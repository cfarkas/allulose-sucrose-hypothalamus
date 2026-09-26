FIGURE 3 — NPY EXPERIMENT WITH EXTENDED WATER CONTROLS, DENSITY CONTOURS AND CLUSTERS

Figure_3_cFos_NPY.pdf and its 600-dpi PNG are the complete English master.
A–I panels and legends are available individually in English and Spanish.

11 animals are included: 5 Water, 3 Sucrose, 3 Allulose.
The added Water animals (WATER_NPY4, WATER_NPY5; specimens N132-5 and N132-3,
24 September 2026) were acquired on the same Zeiss LSM 780 and the same
LD LCI Plan-Apochromat 25x/0.8 objective as the sugar animals. Their native
0.332106 um sampling was reduced by an exact 2x2 block mean to the sugar
cohort's own 0.664213 um pitch before segmentation, so every animal is
segmented at one physical scale. The original nine animals' counts and spatial
values are byte-identical to the previous revision.

The Water group is far more dispersed than either sugar group, and its spread
coincides with acquisition batch. The abundance endpoints therefore report
Welch's robust F with an exactly enumerated permutation null as the primary
omnibus, a Brown–Forsythe test of equal variance, and exact studentized
permutation tests per contrast, with Hodges–Lehmann shifts and Cliff's delta as
effect sizes. The Welch F approximation, the classical equal-variance ANOVA and
the exact Mann–Whitney tests are retained beside them. Each exact test also
records the smallest p it can return at these group sizes. The E and F panels
print the enumerated permutation p and the exact Mann–Whitney p-values; the
rest of that family is in the source table.
See source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv.

Existing reviewed DAPI/NPY/c-FOS ROIs are fixed. No animal is flagged by the
within-condition 1.5-IQR rule. Water animals with fewer than two NPY/c-FOS
double-positive nuclei stay in abundance while their conditional NPY spatial
metrics are undefined; that endpoint has 3 of 5 estimable Water animals.
A–D retain their reviewed files. G uses Allulose E8_FR6-4, the accepted tissue
and third-ventricle contours, and the median eminence at the bottom.
H and I are scale resolved. Each shows, for every animal, the clustering of
activated cells at distances from 20 to 150 um, expressed in standard
deviations from that animal's own random-labelling null, together with a
studentized maximum-deviation global envelope test comparing conditions across
the whole range. One p-value covers every distance, so no radius is singled
out, and the plot shows at which distances the conditions separate.
See scale_resolved_plan.json in the analysis bundle.

Core enrichment and clustering excess are retained and tested separately in
source_data/Figure_3_cFos_NPY_spatial_permanova.csv. There is no joint
two-measure test, and no multivariate dispersion diagnostic: the single-measure
analogue is the Brown-Forsythe equal-variance test in
source_data/Figure_3_cFos_NPY_abundance_robust_statistics.csv.

E and F print the exact permutation p and the exact Mann-Whitney p-values.
Every displayed value is a number, not a significance symbol. NE means not
estimable.

Water acquisition remains confounded with condition and batch; no test removes
that. The within-batch Allulose–Sucrose contrast is the only unconfounded one.

Rebuild and publish from the Paper root:
  ./scripts/run_fig3_new_water_cohort.sh
Add --from-raw to re-derive the 2026 Water animals from their CZI files.
The usual Fig3/07_rebuild_complete_figure3.sh entry routes to this analysis.

Methods, cohort inputs, IQR ledger and results:
  ../analyses/fig3_water_extended_20260924/README.md
Current receipts: provenance/extended_water_cohort_*.json
The previous nine-animal revision remains in
  ../analyses/fig3_original_cohort_20260920/
Previous canonical files are archived in this analysis bundle's archive/Fig3.
Dated older provenance files describe historical revisions.
