"""Prove absent-interaction temporal reuse preserves estimates and RNG."""
import numpy as np
import pandas as pd

from mdscope.statistics.interactions import _interaction_bootstrap


def test_absent_temporal_cache_matches_uncached_science_and_rng():
    cfg = {
        "statistics": {
            "block": {"length": "auto", "tau_multiplier": 2.0},
            "hierarchy": {"block_within_replicate": True},
        }
    }
    times = np.arange(16, dtype=float) * 5
    cases = [
        ("rep1", [], pd.DataFrame()),
        ("rep2", [], pd.DataFrame()),
        ("rep1", [], pd.DataFrame()),
        ("rep1", [10.0, 20.0], pd.DataFrame()),
        ("rep2", [], pd.DataFrame()),
        ("rep1", [35.0], pd.DataFrame()),
    ]

    def execute(reuse):
        rng = np.random.default_rng(2026)
        weights = {}
        zero_meta = {} if reuse else None
        records = []
        for rep_key, presence, lifetimes in cases:
            result = _interaction_bootstrap(
                times, np.asarray(presence), lifetimes, cfg,
                32, rng, weights, rep_key, zero_meta,
            )
            records.append(result)
        return records, rng.integers(0, 1000000, size=20)

    baseline, baseline_rng = execute(False)
    cached, cached_rng = execute(True)
    np.testing.assert_array_equal(baseline_rng, cached_rng)
    for (expected_boot, expected_meta), (actual_boot, actual_meta) in zip(baseline, cached):
        assert expected_meta.keys() == actual_meta.keys()
        for key in expected_meta:
            np.testing.assert_array_equal(
                np.asarray(expected_meta[key]), np.asarray(actual_meta[key]),
                err_msg=f"meta {key}",
            )
        for key in expected_boot:
            np.testing.assert_array_equal(
                expected_boot[key], actual_boot[key], err_msg=f"bootstrap {key}"
            )


def test_absent_temporal_cache_saves_fft_calls(monkeypatch):
    import mdscope.statistics.interactions as interactions

    original = interactions._temporal
    calls = []
    def timed(values, times, cfg):
        calls.append(len(values))
        return original(values, times, cfg)

    monkeypatch.setattr(interactions, "_temporal", timed)
    rng = np.random.default_rng(123)
    weights, metadata = {}, {}
    cfg = {"statistics": {"block": {"length": "auto"},
                          "hierarchy": {"block_within_replicate": True}}}
    times = np.arange(32, dtype=float)
    for _ in range(8):
        interactions._interaction_bootstrap(
            times, [], pd.DataFrame(), cfg, 16, rng,
            weights, "same_rep", metadata,
        )
    assert len(calls) == 1
    assert len(metadata) == 1


def test_disk_backed_weight_cache_preserves_all_estimates_and_rng(monkeypatch):
    import mdscope.statistics.interactions as interactions

    times = np.arange(128, dtype=float) * 5
    lifetimes = pd.DataFrame({'start_ps': [10., 30., 65.],
                              'lifetime_ps': [10., 25., 15.],
                              'left_censored': [False, False, False],
                              'completed': [True, True, True]})
    cfg = {'statistics': {'block': {'length': 4},
                          'hierarchy': {'block_within_replicate': True}}}

    def run(budget):
        monkeypatch.setattr(interactions, '_WEIGHT_CACHE_BYTES', budget)
        rng = np.random.default_rng(456)
        cache, records = {}, []
        for rep in ('one', 'two', 'one', 'three'):
            records.append(_interaction_bootstrap(
                times, times[::3], lifetimes, cfg, 65, rng, cache, rep,
            ))
        resident = sum(v.nbytes for v in cache.values() if not isinstance(v, np.memmap))
        assert resident <= budget
        return records, rng.integers(1000000, size=20)

    expected, expected_rng = run(1024 * 1024)
    actual, actual_rng = run(0)
    np.testing.assert_array_equal(actual_rng, expected_rng)
    for (a, am), (b, bm) in zip(actual, expected):
        for key in a:
            np.testing.assert_array_equal(a[key], b[key])
        for key in am:
            np.testing.assert_array_equal(am[key], bm[key])
