"""Replicate-aware MDscope statistics with explicit functional ownership."""
from importlib import import_module

from ..contracts.statistics import STATISTICS_FAMILIES
from ..core.progress import stage_progress

_ANALYZERS = {
    "time_series": ("time_series", "analyze_time_series"),
    "profiles": ("profiles", "analyze_profiles"),
    "interactions": ("interactions", "analyze_interactions"),
    "solvation_rdfs": ("solvation", "analyze_solvation_rdfs"),
    "diffusion": ("diffusion", "analyze_multi_origin_diffusion"),
    "coupling": ("coupling", "analyze_coupling"),
    "lipid_contacts": ("lipid_contacts", "analyze_lipid_contacts"),
}


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
    import os
    from pathlib import Path
    from time import perf_counter

    from ..config import statistics_family_config
    from ..core.persistence import archive_output, atomic_json
    from ..core.profiling import PhaseProfiler
    profiler = PhaseProfiler.from_environment()
    started_total = perf_counter()
    families = STATISTICS_FAMILIES
    dest = Path(root) / "project_summary/statistics"
    if (dest / "manifest.json").exists():
        archive_output(dest / "manifest.json", Path(root) / ".archives/statistics/manifests",
                       "previous.json")
    manifest = {
        "families": {
            name: ("enabled" if statistics_family_config(cfg, name) is not None
                   and (name != "lipid_contacts" or
                        cfg.get("analyses", {}).get("lipid_contacts", {}).get("enabled", False))
                   else "disabled")
            for name in families
        },
        "methods": cfg.get("statistics", {}),
        "input_analyses": cfg.get("analyses", {}),
        "status": "running",
        "replicate_paths": [Path(os.path.relpath(rep, root)).as_posix()
                            for reps in system_repouts.values() for rep in reps],
    }
    verbose = cfg.get("progress", {}).get("enabled", True)
    with stage_progress("Statistics input inventory", enabled=verbose):
        matched = cfg.get("statistics", {}).get("matched_window_ps") is not None
        bounds = None
        if matched:
            manifest["families"]["matched_time_series"] = manifest["families"]["time_series"]
            if manifest["families"]["matched_time_series"] == "enabled":
                from .windows import matched_grid
                grid = matched_grid(cfg, system_repouts)
                bounds = [float(grid[0]), float(grid[-1])]
                manifest["matched_scalar_window_ps"] = bounds
        inputs_before = profiler.call(
            "input_inventory", _input_snapshot, root, manifest["replicate_paths"]
        )
        atomic_json(dest / "manifest.json", manifest)
    for family in STATISTICS_FAMILIES:
        folder = STATISTICS_FAMILIES[family]
        if (dest / folder).exists():
            archive_output(dest / folder, Path(root) / ".archives/statistics" / folder,
                           "previous")
        if manifest["families"][family] == "disabled":
            if verbose:
                print(f"Statistics {family}: skipped (disabled)", flush=True)
            continue
        with stage_progress(f"Statistics {family}", enabled=verbose):
            module, function = _ANALYZERS[family]
            analyze = getattr(import_module(f".{module}", __package__), function)
            profiler.call(family, analyze, cfg, root, system_repouts)
    matched_folder = dest / "matched_time_series"
    if matched_folder.exists():
        archive_output(matched_folder, Path(root) / ".archives/statistics/matched_time_series", "previous")
    if bounds is not None:
        from .time_series import analyze_time_series
        matched_cfg = {**cfg, "statistics": {**cfg["statistics"], "scalar_window_ps": bounds}}
        with stage_progress("Statistics matched_time_series", enabled=verbose):
            profiler.call(
                "matched_time_series", analyze_time_series, matched_cfg, root,
                system_repouts, output_name="matched_time_series",
            )
    coverage_folder = dest / "coverage"
    if coverage_folder.exists():
        archive_output(coverage_folder, Path(root) / ".archives/statistics/coverage", "previous")
    from .coverage import write_contrast_coverage
    with stage_progress("Statistics contrast coverage", enabled=verbose):
        profiler.call(
            "contrast_coverage", write_contrast_coverage,
            dest, system_repouts, manifest["families"],
        )
    with stage_progress("Statistics input verification", enabled=verbose):
        if inputs_before != profiler.call(
            "input_verification", _input_snapshot, root, manifest["replicate_paths"]
        ):
            raise RuntimeError('Persisted statistics inputs changed during finalization; retry on stable outputs')
    manifest["status"] = "complete"
    atomic_json(dest / "manifest.json", manifest)
    if profiler.enabled:
        measured = sum(profiler.seconds.values())
        total = perf_counter() - started_total
        atomic_json(dest / "statistics_phase_profile.json", {
            "elapsed_seconds": total,
            "measured_seconds": measured,
            "other_seconds": max(0.0, total - measured),
            "replicate_count": sum(map(len, system_repouts.values())),
            "family_timings": [
                {"stage": name, "seconds": float(seconds), "calls": profiler.calls[name]}
                for name, seconds in sorted(profiler.seconds.items(),
                                            key=lambda item: -item[1])
            ],
            "note": "Single-run wall timings include NumPy, I/O and processing; "
                    "do not sum nested timers from other project phases.",
        })
