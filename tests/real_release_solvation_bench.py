"""Benchmark solvation on BOTH published MDscope release trajectories.

Execute in the GitHub Actions release-benchmark workflow. Each variant runs
in a fresh process so VmHWM / ru_maxrss is meaningful across comparisons.
Results are machine-readable, and the script fails on missing input or a
failed benchmark instead of silently reporting synthetic timings.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import MDAnalysis as mda
import yaml

SYSTEMS = {
    "carb_sdc4": ("carb-sdc4-step7_production.tpr",
                  "carb-sdc4-rep1_test.xtc"),
    "nocarb_sdc4": ("nocarb-sdc4-step7_production.tpr",
                    "nocarb-sdc4-rep1_test.xtc"),
}
BENCH = Path(__file__).resolve().parents[1] / "scripts" / "benchmark_solvation_real.py"


def _configuration(topology, trajectory, system, data, *, n_frames):
    u = mda.Universe(str(topology), str(trajectory))
    # Check the actual released topology rather than assuming ion species.
    ions = {name: selection for name, selection in {
        "K": "resname POT", "Cl": "resname CLA", "Na": "resname SOD",
    }.items() if len(u.select_atoms(selection))}
    targets = {
        "protein": {"acidic_oxygen": {
            "selection": "protein and ((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))",
            "group_by": "residue"}},
        "hs": {"glycan_oxygen": {
            "selection": "resname AGLCNA AIDOA BGLCA BXYL and name O*",
            "group_by": "residue"}},
        "lipid": {"phosphate_oxygen": {
            "selection": "resname POPC DPPC POPE POPI POPS PSM DOPC and name O11 O12 O13 O14",
            "group_by": "residue"}},
    }
    sites = {}
    for cls, options in targets.items():
        query = next(iter(options.values()))["selection"]
        count = len(u.select_atoms(query))
        if not count:
            raise AssertionError(f"{system}: required {cls} site selection is empty")
        sites[cls] = count
    water = u.select_atoms("resname TIP3 and name OH2")
    if not len(water):
        raise AssertionError(f"{system}: released TPR has no TIP3 OH2 water oxygens")
    summary = {"atoms": len(u.atoms), "water_oxygen_atoms": len(water),
               "target_atom_counts": sites, "ion_species": list(ions),
               "trajectory_frames_available": len(u.trajectory)}
    del u
    cfg = {
        "project": {"name": "release_solvation_benchmark",
                    "output": str(data)},
        "trajectory": {"start": 0, "stop": None, "stride": 1,
                       "chunk_frames": min(10, n_frames)},
        "systems": {system: {"topology": str(topology),
                             "replicates": [{
                                 "trajectory": str(trajectory), "label": "rep1"}]}},
        "analyses": {
            "alignment": {"enabled": False},
            "solvation": {
                "enabled": True,
                "preset": "standard",
                "water_selection": "resname TIP3 and name OH2",
                "ion_selections": ions,
                "targets": targets,
                "search_radius_A": 8.0,
                "rdf_bin_width_A": 0.2,
                "shell_cutoffs_A": {"water": 3.5, "ions": 4.0},
                "options": {"track_target_classes": ["protein"]},
            },
        },
        "statistics": {"enabled": False},
        "plots": {"enabled": False},
    }
    return cfg, summary


def run(data_dir, output_dir, *, standard_frames, kinetics_frames, sustained_frames=120):
    data_dir = Path(data_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for system, (tpr, xtc) in SYSTEMS.items():
        top, traj = data_dir / tpr, data_dir / xtc
        for path in (top, traj):
            if not path.is_file():
                raise FileNotFoundError(path)
        cfg, topology = _configuration(
            top, traj, system, output_dir,
            n_frames=max(standard_frames, kinetics_frames))
        config_path = output_dir / f"{system}.yaml"
        config_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        print(f"RELEASE_TOPOLOGY {system} {json.dumps(topology, sort_keys=True)}",
              flush=True)
        for preset, n, methods in (
            ("standard", standard_frames, ("auto", "nsgrid")),
            ("kinetics", kinetics_frames, ("auto",)),
            ("standard", sustained_frames, ("auto",)),
        ):
            for method in methods:
                report_path = output_dir / f"{system}-{preset}-{method}-{n}frames.json"
                command = [
                    sys.executable, str(BENCH), str(config_path),
                    "--system", system, "--rep", "rep1",
                    "--frames", str(n), "--preset", preset,
                    "--neighbor-method", method,
                    "--output", str(report_path),
                ]
                print(f"BENCHMARK_START {system} {preset} {method} frames={n}",
                      flush=True)
                proc = subprocess.run(
                    command, text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, check=False)
                if proc.returncode:
                    print(proc.stdout[-14000:], flush=True)
                    raise RuntimeError(
                        f"benchmark failed: {system} {preset} {method} "
                        f"exit={proc.returncode}")
                report = json.loads(report_path.read_text())
                # Save compact summaries to the CI output; detailed per-target
                # selection audits remain in the attached JSON reports.
                compact = {key: report.get(key) for key in (
                    "system", "preset", "neighbor_method",
                    "analyzed_frames", "time_s_total", "time_s_updates",
                    "time_s_finalize", "frames_per_second_updates",
                    "peak_rss_MiB", "rss_MiB_after_universe",
                    "rss_MiB_after_observer_init", "rss_MiB_after_updates",
                    "rss_MiB_after_finalize", "cache_bytes",
                    "finalized_bytes", "total_output_bytes")}
                print("BENCHMARK_RESULT " + json.dumps(compact, sort_keys=True),
                      flush=True)
                reports.append(compact)
    (output_dir / "summary.json").write_text(
        json.dumps({"released_assets": "test-data-v1",
                    "method": "fresh Python process per run",
                    "runs": reports}, indent=2))
    return reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="real_test_data")
    parser.add_argument("--output-dir", default="release_solvation_benchmark")
    parser.add_argument("--standard-frames", type=int, default=24)
    parser.add_argument("--kinetics-frames", type=int, default=12)
    parser.add_argument("--sustained-frames", type=int, default=120)
    args = parser.parse_args()
    if min(args.standard_frames, args.kinetics_frames, args.sustained_frames) < 1:
        parser.error("frame counts must be positive")
    run(args.data_dir, args.output_dir,
        standard_frames=args.standard_frames,
        kinetics_frames=args.kinetics_frames,
        sustained_frames=args.sustained_frames)


if __name__ == "__main__":
    main()
