"""Replicate-aware MDscope statistics with explicit functional ownership."""
from .common import finite, autocorrelation, _temporal, _block_indices, _block_bootstrap, _block_weights, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta
from .time_series import collect_time_series, analyze_time_series
from .profiles import collect_temporal_profiles, analyze_profiles
from .interactions import analyze_interactions
from .solvation import analyze_solvation_rdfs


def _input_snapshot(root, replicate_paths):
    """Detect concurrent edits during statistics without gating future reads."""
    from pathlib import Path
    root = Path(root)
    records = []
    for relative in replicate_paths:
        rep = root / relative
        for folder in ('summary', 'extensions', 'interactions', 'cache'):
            for path in sorted((rep / folder).rglob('*')):
                if path.suffix not in {'.parquet', '.npz'} or not path.is_file():
                    continue
                if path.name in {'time_series_statistics.parquet', 'autocorrelation.parquet'}:
                    continue  # statistics outputs, not inputs
                stat = path.stat()
                records.append((str(path.relative_to(root)) if path.is_relative_to(root) else str(path.resolve()), stat.st_size, stat.st_mtime_ns))
    return records


def analyze_project(cfg, root, system_repouts):
    from ..config import atomic_json, statistics_family_config
    from ..core.persistence import archive_output
    from pathlib import Path
    import os
    families = ("time_series", "profiles", "interactions", "solvation_rdfs")
    dest = Path(root) / "project_summary/statistics"
    if (dest / "manifest.json").exists():
        archive_output(dest / "manifest.json", Path(root) / ".archives/statistics/manifests",
                       "previous.json")
    manifest = {
        "families": {name: "enabled" if statistics_family_config(cfg, name) is not None else "disabled"
                     for name in families},
        "methods": cfg.get("statistics", {}),
        "status": "running",
        "replicate_paths": [os.path.relpath(rep, root)
                            for reps in system_repouts.values() for rep in reps],
    }
    inputs_before = _input_snapshot(root, manifest['replicate_paths'])
    atomic_json(dest / "manifest.json", manifest)
    for family, folder, analyze in (
        ("time_series", "time_series", analyze_time_series),
        ("profiles", "profiles", analyze_profiles),
        ("interactions", "interactions", analyze_interactions),
        ("solvation_rdfs", "solvation", analyze_solvation_rdfs),
    ):
        if (dest / folder).exists():
            archive_output(dest / folder, Path(root) / ".archives/statistics" / folder,
                           "previous")
        analyze(cfg, root, system_repouts)
    if inputs_before != _input_snapshot(root, manifest['replicate_paths']):
        raise RuntimeError('Persisted statistics inputs changed during finalization; retry on stable outputs')
    manifest["status"] = "complete"
    atomic_json(dest / "manifest.json", manifest)
