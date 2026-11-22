"""Benchmark real finalization (statistics, plots, PDF tables) on the two test trajectories.

    python tests/real_finalization_benchmark.py

Uses MDSCOPE_REAL_DATA (default examples/data) and MDSCOPE_REAL_OUTPUT (default
real_test_results). Per-replicate outputs are generated first and are NOT part
of the timed finalization. Five thousand bootstrap draws, one replica per
system: a runtime/resource benchmark, not an inferential result.
"""
from __future__ import annotations

import cProfile
import json
import os
import pstats
import subprocess
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("MDSCOPE_REAL_DATA", "examples/data")).resolve()
OUT = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve() / "finalization_benchmark"
# Expensive production mechanisms/solvation/shared-state kernels are excluded:
# this benchmarks finalization on representative real outputs.
ENABLED = {"alignment", "rmsd", "radius_of_gyration", "density_1d", "membrane",
           "membrane_order", "lipid_enrichment"}


def build_config() -> Path:
    cfg = yaml.safe_load((ROOT / "synd4_full_mechanisms.yaml").read_text(encoding="utf-8"))
    cfg["project"] = {"name": "release_finalize_benchmark", "output": str(OUT / "results")}
    cfg["systems"] = {
        "synd4_glycolipids": {
            "label": "Synd4 + Glycolipids",
            "topology": str(DATA / "carb-sdc4-step7_production.tpr"),
            "replicates": [{"label": "rep1", "trajectory": str(DATA / "carb-sdc4-rep1_test.xtc")}],
        },
        "synd4_no_glycolipids": {
            "label": "Synd4 - Glycolipids",
            "topology": str(DATA / "nocarb-sdc4-step7_production.tpr"),
            "replicates": [{"label": "rep1", "trajectory": str(DATA / "nocarb-sdc4-rep1_test.xtc")}],
        },
    }
    # ~200-frame samples, not the 1000+ production frames the YAML's start skips.
    cfg["trajectory"] = {"start": 0, "stop": None, "stride": 1, "chunk_frames": 50}
    for name, spec in list(cfg["analyses"].items()):
        if name in ENABLED:
            cfg["analyses"][name] = {**spec, "enabled": True} if isinstance(spec, dict) else {"enabled": True}
        else:
            # Drop disabled analyses' options: some (e.g. SDC4 hs_sugar_water_cutoff_A)
            # are only valid while their analysis is enabled.
            cfg["analyses"][name] = {"enabled": False}
    cfg["cache"]["aligned_coordinates"] = False
    cfg["statistics"]["enabled"] = True
    cfg["statistics"]["bootstrap_draws"] = 5000
    cfg["plots"]["enabled"] = True
    cfg["plots"]["table_reports_pdf"] = True  # replaces the retired report section
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "benchmark.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return path


def mdscope(*args: str) -> None:
    subprocess.run([sys.executable, "-u", "-m", "mdscope.cli", *args], check=True)


def main() -> None:
    if (OUT / "results").exists():
        raise SystemExit(f"{OUT / 'results'} exists; remove it for a cold benchmark")
    config = build_config()
    mdscope("validate", str(config))
    for index in ("0", "1"):
        mdscope("replicate", str(config), "--index", index)
    from mdscope.cli import main as cli_main
    profiler = cProfile.Profile()
    sys.argv = ["mdscope", "finalize", str(config)]
    started = time.perf_counter()
    profiler.runcall(cli_main)
    elapsed = time.perf_counter() - started
    profiler.dump_stats(str(OUT / "finalize.prof"))
    with (OUT / "profile_summary.txt").open("w", encoding="utf-8") as handle:
        pstats.Stats(profiler, stream=handle).sort_stats("cumulative").print_stats(75)
    timing = OUT / "results/project_summary/finalization_timings.json"
    stages = json.loads(timing.read_text()) if timing.is_file() else {}
    for row in stages.get("stages", []):
        print(f"FINALIZE {row['stage']}: {row['seconds']:.2f} s")
    print(f"TOTAL FINALIZE: {elapsed:.2f} s wall (profiled)")
    print(f"FINALIZATION BENCHMARK: PASS ({OUT})")


if __name__ == "__main__":
    main()
