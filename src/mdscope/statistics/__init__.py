"""Replicate-aware MDscope statistics with explicit functional ownership."""
from ..pipeline.progress import stage_progress
from .common import finite, autocorrelation, _temporal, _block_indices, _block_bootstrap, _block_weights, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta
from .time_series import collect_time_series, analyze_time_series
from .profiles import collect_temporal_profiles, analyze_profiles
from .interactions import analyze_interactions
from .solvation import analyze_solvation_rdfs
from .diffusion import analyze_multi_origin_diffusion
from .coverage import write_contrast_coverage
from .coupling import analyze_coupling
from .windows import matched_grid


def _input_snapshot(root, replicate_paths):
    """Detect concurrent edits during statistics without gating future reads."""
    from pathlib import Path
    root = Path(root)
    records = []
    for relative in replicate_paths:
        rep = root / relative
        for folder in ('summary', 'extensions', 'interactions', 'cache', 'profiles'):
            for path in sorted((rep / folder).rglob('*')):
                if path.suffix not in {'.parquet', '.npz', '.npy'} or not path.is_file():
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
    families = ("time_series", "profiles", "interactions", "solvation_rdfs", "diffusion", "coupling")
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
    verbose = cfg.get("progress", {}).get("enabled", True)
    with stage_progress("Statistics input inventory", enabled=verbose):
        matched = cfg.get("statistics", {}).get("matched_window_ps") is not None
        bounds = None
        if matched:
            manifest["families"]["matched_time_series"] = manifest["families"]["time_series"]
            if manifest["families"]["matched_time_series"] == "enabled":
                grid = matched_grid(cfg, system_repouts)
                bounds = [float(grid[0]), float(grid[-1])]
                manifest["matched_scalar_window_ps"] = bounds
        inputs_before = _input_snapshot(root, manifest['replicate_paths'])
        atomic_json(dest / "manifest.json", manifest)
    for family, folder, analyze in (
        ("time_series", "time_series", analyze_time_series),
        ("profiles", "profiles", analyze_profiles),
        ("interactions", "interactions", analyze_interactions),
        ("solvation_rdfs", "solvation", analyze_solvation_rdfs),
        ("diffusion", "diffusion", analyze_multi_origin_diffusion),
        ("coupling", "coupling", analyze_coupling),
    ):
        if (dest / folder).exists():
            archive_output(dest / folder, Path(root) / ".archives/statistics" / folder,
                           "previous")
        if manifest["families"][family] == "disabled":
            if verbose:
                print(f"Statistics {family}: skipped (disabled)", flush=True)
            continue
        with stage_progress(f"Statistics {family}", enabled=verbose):
            analyze(cfg, root, system_repouts)
    matched_folder = dest / "matched_time_series"
    if matched_folder.exists():
        archive_output(matched_folder, Path(root) / ".archives/statistics/matched_time_series", "previous")
    if bounds is not None:
        matched_cfg = {**cfg, "statistics": {**cfg["statistics"], "scalar_window_ps": bounds}}
        with stage_progress("Statistics matched_time_series", enabled=verbose):
            analyze_time_series(matched_cfg, root, system_repouts, output_name="matched_time_series")
    coverage_folder = dest / "coverage"
    if coverage_folder.exists():
        archive_output(coverage_folder, Path(root) / ".archives/statistics/coverage", "previous")
    with stage_progress("Statistics contrast coverage", enabled=verbose):
        write_contrast_coverage(dest, system_repouts, manifest["families"])
    with stage_progress("Statistics input verification", enabled=verbose):
        if inputs_before != _input_snapshot(root, manifest['replicate_paths']):
            raise RuntimeError('Persisted statistics inputs changed during finalization; retry on stable outputs')
    manifest["status"] = "complete"
    atomic_json(dest / "manifest.json", manifest)
