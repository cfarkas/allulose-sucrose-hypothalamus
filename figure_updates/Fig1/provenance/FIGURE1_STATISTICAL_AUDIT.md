# Figure 1 scientific and statistical audit

## Decision

Figure 1 is an exploratory January-cohort analysis. It does not reproduce the legacy pooled-cohort MANOVA, row-level OLS, PCA, LDA, or HC3 mouse-level inference. Repeated rows are used for trajectories; each inferential endpoint contributes exactly one value per cage.

## Cohorts and source integrity

- January cohort: E1-E12, dated January 13-22, 2025. Primary window is Day 1, 3, and 6 (January 13, 15, and 17).
- August add-on: E13/E14, dated August 18, 21, and 24, 2025. It contains 6 consumption-table rows, all Bottle_Weight, and 18 mouse-weight rows. It contains zero measured Bottle_Volume rows and is not pooled with January.
- The raw consumption CSV contains 0 measured zero values. The apparent August zero-volume/intake values in the legacy combined output were join/fill artifacts: absent Bottle_Volume was converted to zero. This pipeline never imputes absent bottle volume.
- E9 is retained: 3 primary rows for one mouse, with genotype explicitly represented as Unknown. The consumption file separately labels its volume record NPY, so the cross-file genotype conflict is preserved rather than guessed away.
- Cages containing multiple mouse-level treatment labels are flagged in cage_design_and_size_qc.csv. For the body-weight endpoint, only valid rows matching the cage bottle treatment enter the cage-balanced primary summary.
- There are 8 nonpositive/missing weight rows. They are flagged as invalid and never interpreted as 0-g mice. None is needed to complete the January Day 1/3/6 primary series.

## Consumption (panels A-B)

- Experimental unit: cage.
- Primary quantity: measured Day-1 bottle volume minus the current measured bottle volume, in mL per cage. It is reported as consumption from Day 1; without leakage/evaporation controls it must not be interpreted as individually ingested volume.
- Primary group sizes: Water 3, Sucrose 5, Allulose 4 cages.
- Day-6 omnibus inference: ordinary equal-variance one-way ANOVA on one value per cage, p=0.0101709. Exact enumeration of all 27,720 cage-label allocations is retained as a sensitivity analysis, p=0.00916306; three exact pairwise tests use Holm correction.
- One primary-window reading (E5, Day 6) increases by more than 2 mL from the prior reading (1 flagged row). It may reflect measurement variation or an undocumented intervention and is retained transparently; the endpoint is still baseline minus Day-6 volume.
- The optional per-mouse sensitivity divides each cage value once by historical nominal occupancy. It never duplicates cage values to mouse rows and is not the primary figure outcome.

## Body weight (panels C-D)

- Panel C displays 23 mice descriptively, uniquely keyed by E|AnimalID. Thick summaries are cage-balanced means, not mouse-level inferential estimates.
- Day-6 inference first averages mouse percent changes within each cage, then compares Water 3, Sucrose 5, and Allulose 4 cages.
- Ordinary one-way ANOVA omnibus p=0.119031; exact cage-label sensitivity p=0.0743867. Three pairwise exact tests use Holm correction. No HC3 mouse-level model is used.
- Descriptive mouse counts are Water 4, Sucrose 9, Allulose 10. They do not replace cage N for inference.

## Design limitations and interpretation

1. Historical nominal cage sizes disagree with observed January weight-record subject counts for E4, E9, E10, E11, E12. Consequently, cage-total consumption is primary; nominal per-mouse results are sensitivity only.
2. Treatment allocation/randomization documentation was not found. Exact label-permutation p values therefore rely on exchangeability and are exploratory, not proof from a documented randomized design.
3. There are only 3-5 cages per treatment. Bootstrap intervals are descriptive, and exact p values are discrete.
4. Sex and genotype are sparse, imbalanced, and partly confounded with cage/treatment. No treatment-by-sex/genotype claim is supported here.
5. Repeated Day 1/3/6 observations are not entered into ordinary independent-row tests. Inference uses one prespecified Day-6 endpoint per cage.
6. Pairwise multiplicity is controlled by Holm across the three treatment contrasts separately for each outcome.
7. No IQR deletion is used. All valid prespecified-window observations are retained.
