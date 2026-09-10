# Figure S3: corrected residual-permutation statistics

Analysis version 1.8 uses the same HC3-studentized treatment coefficient in the observed and permuted fits. The previous code used an HC3 observed statistic but ordinary-OLS permutation standard errors. The reduced model now retains the other treatment contrast and cohort as nuisance terms when testing each coefficient.

A fixed 99,999 permutations per contrast (base seed 1707 with the documented endpoint offsets) improves Monte Carlo precision. The `(1 + extreme)/(1 + permutations)` estimator is used; the analysis still relies on residual exchangeability and the nuisance model. HC3 does not by itself establish exchangeability.

Counts, sampling, model coefficients and HC3 confidence intervals are unchanged. The hepatic epithelial-like Allulose/Water CLR ratio remains 4.05987 (95% CI 2.22130–7.42021). Corrected p=0.00031 and global BH q=0.02418, compared with previous p=0.0004 and q=0.0312. Six classical H&E contrasts are nominally below 0.05; none passes global BH. This includes a lower hepatic nuclear-density contrast (p=0.03920; q=0.23520).

`permutation_correction_comparison_20260909.csv` reports every previous and corrected contrast. `tests/test_hc3_permutation.py` compares the vectorized method with independently fitted per-permutation statsmodels HC3 regressions, retaining nuisance variables. The corrected renderer was also rerun at the same fixed seed; all statistical results were unchanged.

Reference: Winkler et al. (2014), NeuroImage 92, 381–397. https://doi.org/10.1016/j.neuroimage.2014.01.060
