"""Synchronization contract for paired primary/MSD coordinate sources."""
from dataclasses import dataclass
import pytest
from mdscope.core.trajectory_sources import paired_frames


@dataclass
class Frame:
    time: float


def test_paired_frames_preserve_primary_identity_and_order():
    primary = [Frame(t) for t in (0., 10., 20.)]
    secondary = [Frame(t + .05) for t in (0., 10., 20.)]
    observed = list(paired_frames(primary, secondary))
    assert all(a is b for a, b in zip(observed, primary))


def test_paired_frames_reject_time_mismatch_and_length_mismatch():
    with pytest.raises(ValueError, match="time synchronized"):
        list(paired_frames([Frame(0.), Frame(10.)], [Frame(0.), Frame(11.)]))
    with pytest.raises(ValueError, match="sampled frame counts"):
        list(paired_frames([Frame(0.), Frame(10.)], [Frame(0.)]))


def test_single_source_does_not_require_secondary_reader():
    primary = [Frame(0.), Frame(10.)]
    assert list(paired_frames(primary, None)) == primary
