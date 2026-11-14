"""Parallel block-bootstrap results must match the serial estimator and RNG state."""
import numpy as np
import pytest

from mdscope.statistics import common


@pytest.mark.parametrize("reducer", ["mean", "median"])
@pytest.mark.parametrize("draws", [17, 97])
@pytest.mark.parametrize("missing", [False, True])
def test_parallel_bootstrap_preserves_draws_and_rng(monkeypatch, reducer, draws, missing):
    values = np.arange(120, dtype=float).reshape(30, 4)
    if missing:
        values[3, 2] = np.nan
    reference_rng = np.random.default_rng(2026)
    reduce = np.nanmedian if reducer == "median" else np.nanmean
    reference = np.concatenate([
        reduce(values[common._block_indices(len(values), 5, min(32, draws - start), reference_rng)], axis=1)
        for start in range(0, draws, 32)
    ])
    serial_rng = np.random.default_rng(2026)
    parallel_rng = np.random.default_rng(2026)
    monkeypatch.setattr(common, "_bootstrap_workers", lambda: 1)
    serial = common._block_bootstrap(values, 5, draws, serial_rng, reducer)
    monkeypatch.setattr(common, "_bootstrap_workers", lambda: 4)
    parallel = common._block_bootstrap(values, 5, draws, parallel_rng, reducer)
    np.testing.assert_array_equal(serial, reference)
    np.testing.assert_array_equal(parallel, reference)
    assert serial_rng.bit_generator.state == reference_rng.bit_generator.state
    np.testing.assert_array_equal(
        parallel_rng.integers(0, 1000, size=20),
        serial_rng.integers(0, 1000, size=20),
    )


def test_worker_count_respects_slurm_and_affinity(monkeypatch):
    monkeypatch.setattr(common.os, "sched_getaffinity", lambda pid: set(range(6)), raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.delenv("MDSCOPE_STATS_WORKERS", raising=False)
    assert common._bootstrap_workers() == 6
    monkeypatch.setenv("MDSCOPE_STATS_WORKERS", "3")
    assert common._bootstrap_workers() == 3


@pytest.mark.parametrize("n, block", [(0, 2), (89, 1), (89, 4), (89, 89), (89, 100), (65535, 1)])
def test_block_weights_match_expanded_indices_and_rng(n, block):
    expected_rng = np.random.default_rng(103)
    actual_rng = np.random.default_rng(103)
    indices = common._block_indices(n, block, 97, expected_rng)
    expected = np.stack([np.bincount(row, minlength=n) for row in indices])

    weights = common._block_weights(n, block, 97, actual_rng)
    np.testing.assert_array_equal(weights, expected)
    assert weights.dtype == (np.uint32 if n == 65535 else np.uint16)
    assert actual_rng.bit_generator.state == expected_rng.bit_generator.state
