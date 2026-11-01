"""Replicate-aware MDscope statistics with explicit functional ownership."""
from .common import finite, autocorrelation, _temporal, _block_indices, _block_bootstrap, _block_weights, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta
from .time_series import collect_time_series, analyze_time_series
from .profiles import collect_temporal_profiles, analyze_profiles
from .interactions import analyze_interactions
from .solvation import analyze_solvation_rdfs


def analyze_project(cfg, root, system_repouts):
    analyze_time_series(cfg, root, system_repouts)
    analyze_profiles(cfg, root, system_repouts)
    analyze_interactions(cfg, root, system_repouts)
    analyze_solvation_rdfs(cfg, root, system_repouts)
