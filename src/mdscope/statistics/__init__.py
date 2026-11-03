"""Replicate-aware MDscope statistics with explicit functional ownership."""
from .common import finite, autocorrelation, _temporal, _block_indices, _block_bootstrap, _block_weights, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta
from .time_series import collect_time_series, analyze_time_series
from .profiles import collect_temporal_profiles, analyze_profiles
from .interactions import analyze_interactions
from .solvation import analyze_solvation_rdfs


def analyze_project(cfg, root, system_repouts):
    from ..config import atomic_json, config_hash, statistics_family_config
    from ..core.persistence import archive_output
    from pathlib import Path
    import json
    families = ("time_series", "profiles", "interactions", "solvation_rdfs")
    dest = Path(root) / "project_summary/statistics"
    previous_hash = "legacy"
    if (dest / "manifest.json").exists():
        try:
            previous_hash = json.loads((dest / "manifest.json").read_text()).get("statistics_hash", "legacy")
        except (ValueError, OSError):
            previous_hash = "unreadable"
        archive_output(dest / "manifest.json", Path(root) / ".archives/statistics/manifests",
                       str(previous_hash) + ".json")
    manifest = {
        "statistics_hash": config_hash(cfg.get("statistics", {})),
        "families": {name: "enabled" if statistics_family_config(cfg, name) is not None else "disabled"
                     for name in families},
        "methods": cfg.get("statistics", {}),
        "status": "running",
    }
    atomic_json(dest / "manifest.json", manifest)
    for family, folder, analyze in (
        ("time_series", "time_series", analyze_time_series),
        ("profiles", "profiles", analyze_profiles),
        ("interactions", "interactions", analyze_interactions),
        ("solvation_rdfs", "solvation", analyze_solvation_rdfs),
    ):
        if (dest / folder).exists():
            archive_output(dest / folder, Path(root) / ".archives/statistics" / folder,
                           str(previous_hash))
        analyze(cfg, root, system_repouts)
    manifest["status"] = "complete"
    atomic_json(dest / "manifest.json", manifest)
