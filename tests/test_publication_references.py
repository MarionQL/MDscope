"""Independent, small numerical reference cases for publication validation.

Oracles use direct sums, combinatorial enumeration and hand-computable values;
they do not call the corresponding production estimator to establish truth.
"""
from itertools import combinations, product
import numpy as np
import pandas as pd
import pytest
from mdscope.statistics.common import (
    autocorrelation, _temporal, _permutation, _fdr,
    _block_bootstrap, _hierarchical, _hierarchical_difference,
)
from mdscope.observables.multi_origin_msd import multi_origin_msd

def test_autocorrelation_direct_biased_oracle():
    x = np.array([1., 2., 4., -1., 3., 8., 0.])
    centered = x - np.mean(x)
    expected = np.array([
        np.dot(centered[:len(x)-k], centered[k:]) / (len(x)-k)
        for k in range(5)
    ])
    expected /= expected[0]
    np.testing.assert_allclose(autocorrelation(x, max_lag=4), expected, atol=1e-12)

def test_auto_block_length_analytic_alternation():
    # Lag-one ACF < 0: the positive sequence ends immediately.
    x = np.array([1., -1.] * 8)
    cfg = {"statistics": {"block": {"length": "auto", "tau_multiplier": 2},
                          "hierarchy": {"block_within_replicate": True}}}
    result = _temporal(x, np.arange(len(x), dtype=float) * 5., cfg)
    assert result["tau_int_frames"] == 0.5
    assert result["block_length_frames"] == 1
    assert result["effective_n"] == len(x)

def test_exact_permutation_independent_vs_enumerated_oracle():
    a = {"r1": 0., "r2": 1.}
    b = {"r1": 3., "r2": 4.}
    cfg = {"statistics": {"hypothesis_tests": {"design": "independent",
          "exact_max_assignments": 1000}}}
    all_values = np.array([0., 1., 3., 4.])
    observed = abs(np.mean([3., 4.]) - np.mean([0., 1.]))
    p = sum(abs(np.mean(all_values[list(c)]) -
                np.mean(np.delete(all_values, list(c)))) >= observed - 1e-15
            for c in combinations(range(4), 2)) / 6
    assert _permutation(a, b, cfg, np.random.default_rng(2)) == p
    assert p == 1/3

def test_exact_paired_sign_flips_independent_oracle():
    a = {"a": 0., "b": 2., "c": 4.}
    b = {"a": 1., "b": 4., "c": 7.}
    cfg = {"statistics": {"hypothesis_tests": {"design": "paired",
          "exact_max_assignments": 1000}}}
    diffs = np.array([1., 2., 3.])
    expected = sum(abs(np.mean(diffs * signs)) >= np.mean(diffs) - 1e-15
                   for signs in product((-1., 1.), repeat=3)) / 8
    assert _permutation(a, b, cfg, np.random.default_rng(3)) == expected
    assert expected == 0.25

def test_bh_reference_with_missing_p_values_and_families():
    data = pd.DataFrame({"family": ["a"]*5 + ["b"],
                         "p": [.01, .04, .03, np.nan, .5, .2]})
    result = _fdr(data, "p", "q", ["family"])
    np.testing.assert_allclose(result["q"].iloc[[0, 1, 2, 4, 5]],
                               [.04, .05333333333333334, .05333333333333334, .5, .2])
    assert np.isnan(result["q"].iloc[3])

def test_block_bootstrap_constant_reference_and_rng_repeatability():
    data = np.ones((12, 2)) * np.array([3., 7.])
    out_a = _block_bootstrap(data, 3, 64, np.random.default_rng(42))
    out_b = _block_bootstrap(data, 3, 64, np.random.default_rng(42))
    np.testing.assert_array_equal(out_a, out_b)
    np.testing.assert_array_equal(out_a, np.tile([3., 7.], (64, 1)))

def test_paired_hierarchical_bootstrap_constant_reference():
    # Each paired difference is exactly +5 regardless of drawn temporal blocks.
    a = {"rep1": np.ones(7)*2, "rep2": np.ones(7)*8}
    b = {"rep1": np.ones(7)*7, "rep2": np.ones(7)*13}
    out = _hierarchical_difference(a, b, 101, np.random.default_rng(5))
    np.testing.assert_allclose(out, 5.)


def test_nonconstant_hierarchy_exact_two_stage_distribution():
    # Each selected replicate contributes ONE temporal draw, regardless of
    # how many original trajectory frames that replicate represented.
    choices = list(product(range(4), repeat=2))
    class Enumeration:
        def __init__(self):
            self.calls = 0
        def integers(self, low, high, size):
            assert (low, high, size) == (0, 2, (16, 2))
            self.calls += 1
            return np.array(choices) // 2 if self.calls == 1 else np.array(choices) % 2
    actual = _hierarchical([[0., 4.], [10., 14.]], 16, Enumeration())
    # Independent enumeration of four equiprobable replicate/draw values.
    expected = [(a+b)/2 for a, b in product([0., 4., 10., 14.], repeat=2)]
    np.testing.assert_array_equal(actual, expected)


def test_nonconstant_paired_hierarchy_preserves_labels_and_draws():
    class Scripted:
        def __init__(self):
            self.script = iter([
                [[0,1], [1,0], [1,1], [0,0]],  # selected pairs
                [[0,1], [1,0], [0,1], [1,1]],  # temporal draws in A
                [[1,0], [1,0], [1,0], [0,1]],  # temporal draws in B
            ])
        def integers(self, low, high, size):
            assert (low, high, size) == (0, 2, (4, 2))
            return np.array(next(self.script))
    a = {"a": np.array([0.,2.]), "b": np.array([10.,14.])}
    b = {"b": np.array([20.,26.]), "a": np.array([3.,7.])}
    actual = _hierarchical_difference(a, b, 4, Scripted())
    np.testing.assert_array_equal(actual, [6.5, 7.5, 11., 3.])

def test_multiorigin_msd_direct_origin_particle_oracle():
    xyz = np.zeros((9, 2, 3))
    xyz[:, 0, 0] = np.arange(9)*2
    xyz[:, 1, 1] = np.arange(9)*3
    lags, msd = multi_origin_msd(xyz, max_lag=6, particle_batch=1)
    expected = [
        np.mean([np.sum((xyz[t+k, p, :2] - xyz[t, p, :2])**2)
                 for p in range(2) for t in range(9-k)])
        for k in range(7)
    ]
    np.testing.assert_array_equal(lags, np.arange(7))
    np.testing.assert_allclose(msd, expected, atol=1e-10)
    np.testing.assert_allclose(msd, 6.5*lags.astype(float)**2, atol=1e-10)


def test_positive_autocorrelation_time_and_effective_sample_reference():
    # Centered ramp [-2,-1,0,1,2] has unbiased ACF [1, .5, -1/6].
    # The first nonpositive lag stops the sum: tau=.5+.5=1, Neff=5/2.
    cfg = {"statistics": {"block": {"max_acf_lag": 2, "tau_multiplier": 2}}}
    meta = _temporal([0.,1.,2.,3.,4.], np.arange(5)*7., cfg)
    np.testing.assert_allclose(meta["acf"], [1., .5, -1/6], atol=1e-14)
    assert meta["tau_int_frames"] == 1.
    assert meta["tau_int_ps"] == 7.
    assert meta["effective_n"] == 2.5 and meta["block_length_frames"] == 2


@pytest.mark.parametrize("design", ["paired", "independent"])
def test_permutation_missing_labels_nan_and_unequal_population_reference(design):
    cfg = {"statistics": {"hypothesis_tests": {"design": design}}}
    # Paired complete cases are x:+2 and y:+4; only two extreme sign flips
    # of four -> p=.5. Unmatched labels are excluded only in paired mode.
    a = {"x": 0., "y": 1., "missing": np.nan}
    b = {"y": 5., "x": 2., "other": 9.}
    if design == "paired":
        expected = .5
    else:
        values = np.array([0.,1.,2.,5.,9.])
        observed = abs(np.mean([2.,5.,9.]) - .5)
        expected = sum(abs(np.mean(np.delete(values, c))-np.mean(values[list(c)])) >= observed-1e-15
                       for c in combinations(range(5), 2))/10
    assert _permutation(a, b, cfg, np.random.default_rng(42)) == pytest.approx(expected)
    assert np.isnan(_permutation({"x": np.nan}, {"x": 1.}, cfg, np.random.default_rng(42)))


def test_bh_ties_nan_families_and_nonconsecutive_indices():
    table = pd.DataFrame({"family": ["a"]*4+["b"]*2,
                          "p": [.01,.01,np.nan,.2,.04,.08]}, index=[9,2,40,5,7,3])
    out = _fdr(table, "p", "q", ["family"])
    np.testing.assert_allclose(out.q.iloc[[0,1,3,4,5]], [.015,.015,.2,.08,.08])
    assert np.isnan(out.q.iloc[2])


def test_moving_block_boundary_and_median_independent_enumeration():
    class Enumeration:
        def integers(self, low, high, size, dtype):
            assert (low,high,size) == (0,2,(4,2))
            return np.array([[0,0],[0,1],[1,0],[1,1]], dtype=dtype)
    np.testing.assert_array_equal(_block_bootstrap([0.,2.,8.], 2, 4, Enumeration(), "median"), [0.,2.,2.,2.])
    # A block longer than the trace is clipped and cannot randomly wrap.
    np.testing.assert_array_equal(_block_bootstrap([0.,2.,8.], 99, 4, np.random.default_rng(9)), [10/3]*4)
