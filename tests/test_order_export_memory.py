"""Profile-memory budget tests for large membrane-order exports."""
from types import SimpleNamespace

import numpy as np
import pytest

from mdscope.analysis.fused_membrane import _profile_batches


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
