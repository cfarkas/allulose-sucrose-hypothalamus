"""Bounded numerical/oracle checks for the replacement spatial summaries."""
import importlib.util
from itertools import combinations
from pathlib import Path
import unittest
import tempfile
import numpy as np
import pandas as pd

module_path = Path(__file__).with_name('shape_spatial_analysis.py')
if not module_path.exists():
    module_path = Path(__file__).with_name('shape_analysis.py')
spec = importlib.util.spec_from_file_location('shape', module_path)
shape = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shape)


class SpatialTests(unittest.TestCase):
    def test_expected_clustered_matches_exhaustive_labels(self):
        rng = np.random.default_rng(132)
        graphs = [np.array([[0, 1]], int), np.array([[0, 1], [2, 3]], int),
                  np.array([[0, 1], [1, 2], [2, 3], [3, 0]], int)]
        for n in range(2, 8):
            xy = rng.normal(size=(n, 2)) * 20
            graphs_n = [shape.graph_edges(xy, k=2, cap=25)]
            if n == 4:
                graphs_n += graphs[1:]
            if n == 2:
                graphs_n += graphs[:1]
            for edges in graphs_n:
                for m in range(n + 1):
                    observed = []
                    for subset in combinations(range(n), m):
                        labels = np.zeros(n, dtype=bool)
                        labels[list(subset)] = True
                        _, sizes, _, _ = shape.observed_clusters(n, labels, edges)
                        observed.append(np.count_nonzero(sizes))
                    self.assertAlmostEqual(shape.expected_clustered_cells(n, m, edges), np.mean(observed), places=10)

    def test_sparse_pair_and_tissue_gap(self):
        xy = np.array([[5., 5.], [15., 5.]])
        edges = shape.graph_edges(xy)
        self.assertEqual(edges.tolist(), [[0, 1]])
        support = np.ones((25, 25), bool)
        self.assertTrue(shape.segment_valid(xy, edges, support, (1, 1))[0])
        support[:, 10] = False
        self.assertFalse(shape.segment_valid(xy, edges, support, (1, 1))[0])
        self.assertAlmostEqual(shape.expected_clustered_cells(2, 2, edges), 2)

    def test_graph_rigid_motion_invariance(self):
        rng = np.random.default_rng(5)
        xy = rng.uniform(size=(20, 2)) * 100
        expected = shape.graph_edges(xy, 6, 30)
        rotated = xy @ np.array([[0., -1.], [1., 0.]]) + [100, 5]
        np.testing.assert_array_equal(expected, shape.graph_edges(rotated, 6, 30))

    def test_footprint_mass_padding_and_support(self):
        xy = np.array([[20., 20.], [25., 20.], [30., 30.], [35., 25.], [40., 40.], [45., 30.]])
        support = np.ones((75, 75), bool)
        support[:, 50:55] = False
        result = shape.marker_footprint(xy, support, (1., 1.), bandwidth=10., grid=5.)
        self.assertTrue(result['available'])
        self.assertGreaterEqual(result['achieved_mass'], .8 - 1e-12)
        self.assertFalse((result['core_mask'] & ~result['support_mask']).any())
        self.assertFalse(result['core_mask'][[0, -1]].any())
        self.assertFalse(result['core_mask'][:, [0, -1]].any())
        self.assertFalse(shape.marker_footprint(xy[:4], support, (1., 1.))['available'])
        inside = shape.points_in_core(np.array([[-50., -50.], [500., 500.]]), result, 5.)
        np.testing.assert_array_equal(inside, [False, False])

    def test_stratified_aggregation_and_missing(self):
        records = []
        # Strongly different prevalence means pooled expected core count is wrong:
        # stratum expectation 9*.2 + 1*.8 = 2.6, vs pooled 10*.5 = 5.
        for i, (m, core_n, observed) in enumerate([(9, 2, 2), (1, 8, 1)]):
            records.append(dict(cohort='NPY', endpoint='cfos', animal_id='x', region='FIELD', condition='Water', cage='unknown',
                                acquisition_id=str(i), N=10, m=m, shape_N=10, shape_m=m, shape_core_N=core_n,
                                core_observed=observed, core_expected=m*core_n/10, graph_N=10, graph_m=m, graph_edges=3,
                                cluster_count=1, clustered_cells=2, cluster_expected_cells=1, shape_available=True,
                                largest_cluster=2, marker_channel_present=True))
        result = shape.aggregate_animals(pd.DataFrame(records)).iloc[0]
        self.assertAlmostEqual(result.core_expected, 2.6)
        self.assertAlmostEqual(result.core_enrichment, 3/2.6)
        self.assertAlmostEqual(result.cluster_fraction, .4)
        self.assertAlmostEqual(result.cluster_expected_fraction, .2)
        self.assertAlmostEqual(result.cluster_excess, .2)
        for record in records:
            record.update(m=0,shape_m=0,core_observed=0,core_expected=0,graph_m=0,cluster_count=0,clustered_cells=0,cluster_expected_cells=0)
        zero = shape.aggregate_animals(pd.DataFrame(records)).iloc[0]
        self.assertTrue(np.isnan(zero.core_enrichment))
        self.assertTrue(np.isnan(zero.cluster_fraction))
        self.assertEqual(zero.cluster_count, 0)

    def test_exact_and_planned_holm(self):
        result = shape.exact_comparison(np.array([10, 11, 12]), np.array([0, 1, 2]))
        self.assertEqual(result['assignments'], 20)
        self.assertAlmostEqual(result['p_raw'], .1)
        self.assertAlmostEqual(result['minimum_attainable_p'], .1)
        adjusted = shape.holm_planned(np.array([.001, .02] + [np.nan]*10))
        self.assertAlmostEqual(adjusted[0], .012)
        self.assertAlmostEqual(adjusted[1], .22)
        self.assertTrue(np.isnan(adjusted[2:]).all())

    def test_acquisition_bandwidth_sensitivity_reuses_graph(self):
        with tempfile.TemporaryDirectory(prefix='shape_synthetic_', dir='/tmp') as directory:
            root = Path(directory)
            (root / 'contours').mkdir()
            support = np.ones((101, 101), bool)
            geometry_path = root / 'geometry.npz'
            np.savez_compressed(geometry_path, tissue=support, side=np.ones_like(support, int), roi=np.ones_like(support, int))
            xy = np.array([[20.,20.],[25.,20.],[35.,25.],[40.,35.],[45.,45.],[80.,80.]])
            cells = pd.DataFrame(dict(x_um=xy[:,0], y_um=xy[:,1], x_px=xy[:,0], y_px=xy[:,1],
                                      cell_id=[f'c{i}' for i in range(6)], region='FIELD', hemifield='SIDE_A',
                                      in_tissue=True, accepted_legacy=[True]*5+[False], marker_channel_present=True,
                                      is_cfos=[True,True,False,True,False,False]))
            manifest = dict(cohort='NPY', animal_id='SYNTHETIC', condition='Water', cage='unknown',
                            acquisition_id='SYNTHETIC_S01', geometry_path=str(geometry_path), pixel_x_um=1., pixel_y_um=1.)
            config = dict(k=6,edge_cap_um=75.,bandwidth_um=40.,grid_um=5.,contour_mass=.8,
                          minimum_marker_anchors=5,bandwidth_sensitivity=True)
            result = shape.process_acquisition(dict(cells=cells,manifest=manifest,config=config,output=str(root)))
            self.assertEqual(len(result['strata']), 2)
            self.assertEqual(len(result['sensitivity']), 6)
            self.assertEqual(sorted({r['bandwidth_um'] for r in result['sensitivity']}), [20.,40.,80.])
            for endpoint in shape.ENDPOINTS:
                primary = next(r for r in result['strata'] if r['endpoint']==endpoint)
                for row in result['sensitivity']:
                    if row['endpoint']==endpoint:
                        for key in ['graph_edges','cluster_count','clustered_cells','cluster_expected_cells']:
                            self.assertEqual(primary[key], row[key])

    def test_cage_equal_animal_mean(self):
        rows = []
        for i, (n, ratio) in enumerate([(10, 0.), (1000, 2.)]):
            rows.append(dict(cohort='POMC', endpoint='cfos', animal_id=f'a{i}', region='ARC', condition='Water', cage='c1',
                             N=n,m=n,shape_N=n,shape_m=n,cluster_count=1,clustered_cells=n/2,cluster_expected_cells=n/4,
                             core_enrichment=ratio,cluster_fraction=ratio,cluster_expected_fraction=ratio,cluster_excess=ratio))
        units = shape.unit_metrics(pd.DataFrame(rows))
        self.assertEqual(len(units), 1)
        self.assertEqual(units.iloc[0].core_enrichment, 1.)


if __name__ == '__main__':
    unittest.main(verbosity=2)
