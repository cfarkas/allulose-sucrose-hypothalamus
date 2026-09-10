"""Numerical regression test: independent per-draw HC3 OLS reference."""
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('figs3_analysis', ROOT / '05_analyze_make_figure.py')
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)

class HC3PermutationTest(unittest.TestCase):
    def test_matches_independent_scalar_fits_with_nuisance_treatment(self):
        rng = np.random.default_rng(9403)
        frame = pd.DataFrame({'treatment': ['Water']*7+['Sucrose']*9+['Allulose']*8,
                              'cohort': ['NPY','POMC']*12})
        x = np.column_stack([np.ones(24), frame.treatment.eq('Sucrose'),
                             frame.treatment.eq('Allulose'), frame.cohort.eq('POMC')]).astype(float)
        frame['outcome'] = x @ np.array([2, .4, 1.5, .8]) + rng.normal(size=24)*np.linspace(.4, 2.2, 24)
        y = frame.outcome.to_numpy()
        permutations, seed = 199, 723
        got = analysis.regression_one(frame, 'outcome', permutations, seed)
        rng = np.random.default_rng(seed)
        indices = [rng.permutation(24) for _ in range(permutations)]
        observed = sm.OLS(y, x).fit(cov_type='HC3')
        for row, k in zip(got, (1,2)):
            # Explicitly retain the other treatment dummy and cohort under H0.
            nuisance = np.delete(x, k, axis=1)
            fitted = sm.OLS(y, nuisance).fit().fittedvalues
            residual = y - fitted
            statistics = [sm.OLS(fitted + residual[ix], x).fit(cov_type='HC3').tvalues[k]
                          for ix in indices]
            expected = (1 + sum(abs(t) >= abs(observed.tvalues[k]) for t in statistics))/(permutations+1)
            self.assertEqual(row['freedman_lane_permutation_p_value'], expected)
            self.assertAlmostEqual(row['effect'], observed.params[k], places=12)
            self.assertAlmostEqual(row['robust_se'], observed.bse[k], places=12)
            np.testing.assert_allclose([row['ci95_low'], row['ci95_high']], observed.conf_int()[k])

if __name__ == '__main__':
    unittest.main()
