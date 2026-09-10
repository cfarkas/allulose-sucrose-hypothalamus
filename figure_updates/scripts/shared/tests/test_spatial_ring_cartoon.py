"""The explanatory geometry must agree with both production shell algorithms."""
import importlib.util
from pathlib import Path
import sys
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'scripts/setup').is_dir())
sys.path.insert(0, str(ROOT/'scripts/shared'))
from spatial_ring_cartoon import illustrative_geometry, ventricle_half_width


def module(figure):
    path=ROOT/figure/'05_analyze_spatial_distributions.py'
    spec=importlib.util.spec_from_file_location('ring_test_'+figure,path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    return mod


def test_illustrated_rings_match_both_production_algorithms():
    points,centre,normal,radius,edges,shells,*_ = illustrative_geometry()
    for figure in ('Fig3','Fig4'):
        result=module(figure).covariance_normalized_radial_shells(points)
        np.testing.assert_array_equal(result[0],shells)
        np.testing.assert_allclose(result[2] if figure=='Fig3' else result[1],edges)
    np.testing.assert_array_equal(np.bincount(shells),[60]*6)
    assert np.count_nonzero(shells<2)==120
    np.testing.assert_allclose(normal.mean(axis=0),0,atol=1e-14)
    np.testing.assert_allclose(np.cov(normal,rowvar=False),np.eye(2),atol=1e-14)


def test_shell_membership_is_invariant_to_translation_rotation_and_stretch():
    points=illustrative_geometry()[0]
    affine=np.array([[1.8,-.5],[.7,.65]])
    for figure in ('Fig3','Fig4'):
        function=module(figure).covariance_normalized_radial_shells
        np.testing.assert_array_equal(function(points)[0],function(points@affine+[81,-32])[0])


def test_illustrative_nuclei_are_outside_the_ventricular_lumen():
    points=illustrative_geometry()[0]
    x,y=points.T
    assert not np.any((y > -1.5) & (np.abs(x) < ventricle_half_width(y)))
