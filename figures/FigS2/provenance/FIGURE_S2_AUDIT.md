# Figure S2 scientific and statistical audit

## Decision

Figure S2 is an exploratory sex-stratified companion to Figure 1. It repeats the January single-bottle outcomes within female cages. It is not a treatment-by-sex interaction analysis and must not be interpreted as evidence that a treatment effect differs between sexes.

## Eligibility of cage-level consumption for sex stratification

- All 5 included female cages contain one recorded sex in the bottle-volume file and one recorded sex among primary weight mice.
- Sex and treatment agree across the consumption and primary body-weight sources for every January cage. Therefore cage consumption can be labeled by cage sex without assigning a shared bottle measurement to individual mice.
- Consumption is calculated from cage-bottle volume change and can include unmeasured leakage or evaporation; it is not individual intake.

## Experimental unit and sample sizes

- The experimental unit is the cage for consumption, cage-balanced body-weight summaries, confidence intervals, and all p values.
- Female cage N: Water 1, Sucrose 2, Allulose 2.
- Mouse rows (13 mice total) appear only as pale descriptive points in body-weight trajectories; mice do not replace cages for inference.

## Statistics

- The prespecified endpoint is Day 6. The displayed omnibus p value is an ordinary equal-variance one-way ANOVA on one value per cage. Exact cage-label enumeration of the same F statistic is retained as a sensitivity analysis.
- Pairwise tests exhaustively enumerate treatment labels for the two groups being compared. Holm adjustment is applied to the three pairwise p values separately within each sex and outcome.
- Pointwise trajectory intervals resample complete cage trajectories so the same sampled cage identities are used at Days 1, 3, and 6.
- Female Water has one cage. Its trajectory and endpoint 95% confidence interval are explicitly not estimable; no degenerate one-cage bootstrap interval is reported.

## Interpretation limits

1. Treatment randomization documentation was not found, so label exchangeability is an assumption and all p values are exploratory.
2. Within-sex cage counts are extremely small. Exact p values have low resolution and bootstrap intervals with two or three cages are descriptive and unstable.
3. No treatment-by-sex interaction is estimated or tested. Comparing significance in females with significance in males would be invalid.
4. Sex is assigned at cage level only because every January cage passed the single-sex cross-source audit recorded in cage_sex_assignment_audit.csv.
5. Repeated observations are used for trajectories, not treated as independent rows in the inferential tests.
6. All valid prespecified-window observations are retained; no IQR deletion is used.

## Recorded omnibus results

- Female, Day-6 cumulative cage consumption (mL/cage): one-way ANOVA p=0.031525851198; exact cage-label sensitivity p=0.2 (30 allocations).
- Female, Day-6 cage-mean body-weight change (%): one-way ANOVA p=0.303851826515; exact cage-label sensitivity p=0.4 (30 allocations).
