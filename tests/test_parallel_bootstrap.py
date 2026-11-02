"""Parallel block-bootstrap results must match the serial estimator and RNG state."""
import numpy as np
import pytest

from mdscope.statistics import common


@pytest.mark.parametrize("reducer", ["mean", "median"])
@pytest.mark.parametrize("draws", [17, 97])
def test_parallel_bootstrap_preserves_draws_and_rng(monkeypatch, reducer, draws):
    values = np.arange(120, dtype=float).reshape(30, 4)
    values[3, 2] = np.nan
    serial_rng = np.random.default_rng(2026)
    parallel_rng = np.random.default_rng(2026)
    monkeypatch.setattr(common, "_bootstrap_workers", lambda: 1)
    serial = common._block_bootstrap(values, 5, draws, serial_rng, reducer)
    monkeypatch.setattr(common, "_bootstrap_workers", lambda: 4)
    parallel = common._block_bootstrap(values, 5, draws, parallel_rng, reducer)
    np.testing.assert_array_equal(parallel, serial)
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
