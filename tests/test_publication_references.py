"""Independent, small numerical reference cases for publication validation.

Oracles use direct sums, combinatorial enumeration and hand-computable values;
they do not call the corresponding production estimator to establish truth.
"""
from itertools import combinations, product
import numpy as np
import pandas as pd
from mdscope.analysis.statistics.common import (
    autocorrelation, _temporal, _permutation, _fdr,
    _block_bootstrap, _hierarchical_difference,
)
from mdscope.analysis.observables.multi_origin_msd import multi_origin_msd

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
