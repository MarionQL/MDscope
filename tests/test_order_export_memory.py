"""Profile-memory budget tests for large membrane-order exports."""
from types import SimpleNamespace

import numpy as np
import pytest

from mdscope.pipeline.fused_membrane import _profile_batches


def _profiles(*positions):
    return [
        SimpleNamespace(positions=np.arange(n), profile_id=str(i))
        for i, n in enumerate(positions)
    ]


def test_profile_batches_limit_memory_without_losing_output_order():
    profiles = _profiles(3, 5, 2)
    assert _profile_batches(profiles, 10, 320) == [[0, 1], [2]]
    assert _profile_batches(profiles, 10, 120) == [[0], [1], [2]]
    assert _profile_batches(profiles, 0, 1) == [[0, 1, 2]]
    assert _profile_batches([], 10, 1024) == []
    with pytest.raises(ValueError):
        _profile_batches(profiles, 10, 0)


def test_single_large_order_profile_spills_without_changing_saved_values(tmp_path, monkeypatch):
    from collections import defaultdict

    import pandas as pd

    from mdscope.core.persistence import atomic_npz, atomic_parquet
    from mdscope.pipeline.fused_membrane import FusedMembrane

    observer = FusedMembrane.__new__(FusedMembrane)
    observer.names = ('membrane_order',)
    observer.lipid_contacts = None
    observer.repout = tmp_path / 'result'
    observer.cache = tmp_path / 'cache'
    observer.timings = defaultdict(float)
    observer.anchor_states = {}
    observer.frame_count = 2000
    profile = SimpleNamespace(profile_id='order', positions=np.arange(200),
                              labels=np.array([str(i) for i in range(200)]))
    observer.order = {'profiles': [profile], 'sites': [], 'pairs': []}
    times = np.arange(2000, dtype=float)
    expected = np.arange(400000, dtype=np.float32).reshape(2000, 200)
    for i in range(2):
        selected = slice(i * 1000, (i + 1) * 1000)
        atomic_npz(observer.cache / 'order_profiles' / f'part-{i:06d}.npz',
                   time_ps=times[selected], profile_0000=expected[selected])
        atomic_parquet(observer.cache / 'membrane_order' / f'part-{i:06d}.parquet',
                       pd.DataFrame({'time_ps': times[selected], 'order': [1.] * 1000}))
    monkeypatch.setenv('MDSCOPE_ORDER_FINALIZE_MEMORY_MB', '1')
    observer.finalize()
    with np.load(observer.repout / 'profiles/membrane_order/profile-0000.npz') as saved:
        np.testing.assert_array_equal(saved['time_ps'], times)
        np.testing.assert_array_equal(saved['values'], expected)
        np.testing.assert_array_equal(saved['positions'], profile.positions)
