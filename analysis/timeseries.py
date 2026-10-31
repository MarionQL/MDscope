"""Backward-compatible public imports; implementation lives in one-purpose modules.

New code should import directly from alignment/rmsd/rmsf/membrane/density_1d/streaming.
"""
from .core.alignment import _fit_transform, _apply_transform, _valid_box, resolved_selection_counts
from .observables.rmsd import compute_rmsd
from .observables.rmsf import Welford, temporal_rmsf_profiles
from .observables.membrane import MembraneState,initialize_membrane,membrane_metrics
from .observables.density_1d import DensityTarget,DensityProfileAccumulator,temporal_density_profiles
from .streaming import run_streaming_pass
