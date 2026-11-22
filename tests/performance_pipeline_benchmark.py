"""Reproducible *performance-only* benchmark on release TPR/XTC test data.

This is NOT a convergence study: inputs are two distinct conditions, each
containing one short trajectory, not independent biological replicates.
The statistics benchmark uses separate output roots with links to the
persisted per-replica data; no production replica is modified.

Run after downloading the test-data-v1 GitHub release:
  MDSCOPE_REAL_DATA=real_test_data MDSCOPE_REAL_OUTPUT=real_benchmark_results \\
    python tests/performance_pipeline_benchmark.py --frames 40 160 --draws 100 500
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import platform
import resource
from time import perf_counter

import yaml

from mdscope.config import load_config
from mdscope.pipeline.runner import run_one
from mdscope.statistics import analyze_project
from real_trajectory_smoke import DATA, ROOT, SYSTEMS, build_config, inspect


def _read(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else {}


def _profile(root):
    result = {}
    for system in SYSTEMS:
        prefix = root / system / "rep1"
        result[system] = {
            "stage_timings": _read(prefix / "stage_timings.json"),
            "streaming": _read(prefix / "streaming_phase_profile.json"),
            "membrane": _read(prefix / "fused_membrane_timings.json"),
            "extensions": _read(prefix / "fused_extensions_timings.json"),
            "distance": _read(prefix / "fused_distance_timings.json"),
        }
    return result


def _trajectory_experiment(root, frames, species):
    if root.exists():
        raise FileExistsError(f"Benchmark output already exists: {root}")
    root.mkdir(parents=True)
    measurements = {}
    for system, (top, xtc) in SYSTEMS.items():
        config = build_config(system, DATA / top, DATA / xtc, species[system])
        data = yaml.safe_load(config.read_text())
        data["project"]["output"] = str(root)
        data["trajectory"]["stop"] = int(frames)
        # Explicitly measure the full cold run. Do not re-use previous output.
        data["statistics"]["enabled"] = False
        config = ROOT / f"performance_{frames}_{system}.yaml"
        config.write_text(yaml.safe_dump(data, sort_keys=False))
        started = perf_counter()
        run_one(config, 0)
        measurements[system] = {"elapsed_seconds": perf_counter() - started}
    return {"wall": measurements, "profiles": _profile(root)}


def _stats_experiment(source, root, frames, draws, occupancy, template):
    if root.exists():
        raise FileExistsError(f"Statistics output already exists: {root}")
    root.mkdir(parents=True)
    repouts = {}
    for system in SYSTEMS:
        destination = root / system
        destination.mkdir()
        link = destination / "rep1"
        link.symlink_to((source / system / "rep1").resolve(), target_is_directory=True)
        repouts[system] = [link]
    cfg = deepcopy(template)
    cfg["progress"]["enabled"] = False
    cfg["statistics"]["enabled"] = True
    cfg["statistics"]["bootstrap_draws"] = int(draws)
    cfg["statistics"]["interaction_min_occupancy"] = float(occupancy)
    cfg["statistics"]["hypothesis_tests"]["permutation_draws"] = 200
    started = perf_counter()
    analyze_project(cfg, root, repouts)
    return {
        "frames": frames,
        "bootstrap_draws": draws,
        "interaction_min_occupancy": occupancy,
        "elapsed_seconds": perf_counter() - started,
        "profile": _read(root / "project_summary/statistics/statistics_phase_profile.json"),
        "interaction_detail": _read(
            root / "project_summary/statistics/interactions/interaction_phase_profile.json"
        ),
        "note": ("One replica per condition: this measures CPU/I/O only, "
                 "not inferential power or scientific uncertainty."),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", nargs="+", type=int, default=[40, 160])
    parser.add_argument("--draws", nargs="+", type=int, default=[100, 500])
    parser.add_argument("--occupancies", nargs="+", type=float, default=[0.0, 0.20],
                        help="Per-identity occupancy filter; 0.2 matches SDC4 production")
    parser.add_argument("--skip-trajectory", action="store_true",
                        help="Reuse existing completed benchmark inputs")
    args = parser.parse_args()
    if any(n < 2 for n in args.frames) or any(n < 2 for n in args.draws):
        parser.error("frames and draws must be at least 2")
    if any(not 0 <= value <= 1 for value in args.occupancies):
        parser.error("occupancy must be between zero and one")
    ROOT.mkdir(parents=True, exist_ok=True)
    os.environ["MDSCOPE_PROFILE"] = "1"
    os.environ["MDSCOPE_STATS_WORKERS"] = "1"
    species = {system: inspect(system, DATA / top, DATA / xtc)
               for system, (top, xtc) in SYSTEMS.items()}
    results = {
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "stats_workers": 1,
            "source": "test-data-v1 release trajectories",
            "per_condition_replicates": 1,
        },
        "trajectory": {},
        "statistics": {},
    }
    for frames in args.frames:
        root = ROOT / f"benchmark_frames_{frames}"
        if not args.skip_trajectory:
            results["trajectory"][str(frames)] = _trajectory_experiment(
                root, frames, species)
        elif not all((root / system / "rep1/summary/timeseries.parquet").is_file()
                     for system in SYSTEMS):
            raise FileNotFoundError(f"Missing precomputed benchmark results: {root}")

    longest = max(args.frames)
    source = ROOT / f"benchmark_frames_{longest}"
    template = load_config(ROOT / f"performance_{longest}_{next(iter(SYSTEMS))}.yaml")
    for occupancy in args.occupancies:
        for draws in args.draws:
            label = f"{draws}draws_{int(occupancy * 100):02d}pct"
            root = ROOT / f"benchmark_statistics_{longest}f_{label}"
            results["statistics"][label] = _stats_experiment(
                source, root, longest, draws, occupancy, template)

    results["peak_rss_mb"] = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
        if platform.system() == "Linux" else None
    )
    output = ROOT / "performance_benchmark.json"
    output.write_text(json.dumps(results, indent=2))
    print(f"Performance report: {output}", flush=True)
    for frames, record in results["trajectory"].items():
        for system, entry in record["wall"].items():
            print(f"BENCH trajectory frames={frames} system={system} "
                  f"wall={entry['elapsed_seconds']:.3f}s", flush=True)
    for label, entry in results["statistics"].items():
        detail = entry["interaction_detail"]
        identities = detail.get("total_unique_identities_across_kinds")
        print(f"BENCH statistics frames={longest} scenario={label} "
              f"identities={identities} wall={entry['elapsed_seconds']:.3f}s", flush=True)


if __name__ == "__main__":
    main()
