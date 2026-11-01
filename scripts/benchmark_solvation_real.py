"""Run a *real SolvationObserver* preflight/benchmark on HPC trajectory frames.

Does not edit existing MDscope outputs. Measures fused observer update/flush/
finalize costs and produces atom-selection audit and output-volume statistics.
It is opt-in: no trajectory paths or CHARMM atom identities are invented.
Usage:
  python scripts/benchmark_solvation_real.py synd4_solvation.yaml --system synd4_glycolipids --rep rep1 --frames 100 --preset standard
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import tempfile
import time

import MDAnalysis as mda
import numpy as np

from mdscope.config import iter_replicates, load_config, time_ps
from mdscope.observables.solvation import SolvationObserver


def audit(observer):
    """Validate selections against the loaded, actual topology."""
    rows = []
    for target, kind, label, group in observer.targets:
        rows.append({"target_id": target, "target_class": kind,
                     "site_class": label, "n_atoms": len(group),
                     "residue_indices": [int(x) for x in group.resindices],
                     "resnames": sorted(map(str, set(group.resnames))),
                     "atom_names": sorted(map(str, set(group.names)))})
    return {"targets": rows,
            "solvents": {
                name: {"atoms": len(group),
                       "resnames": sorted(map(str, set(group.resnames))),
                       "atom_names": sorted(map(str, set(group.names)))}
                for name, group in observer.solvents.items()}}


def benchmark(cfg, config_path, *, system, replicate, frames, preset, output=None):
    choices = [rep for rep in iter_replicates(cfg, config_path)
               if rep.system == system and rep.label == replicate]
    if len(choices) != 1:
        raise ValueError(f"No unique system/replicate: {system}/{replicate}")
    rep = choices[0]
    if not cfg["analyses"]["solvation"].get("enabled", False):
        raise ValueError("Enable analyses.solvation in the configuration")
    cfg["analyses"]["solvation"]["preset"] = preset
    # Select just the specified frames; preserve original stride and start.
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    stop = tr.get("stop")
    stride = int(tr["stride"])
    if stride <= 0 or frames <= 0:
        raise ValueError("Positive stride and frame count required")
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    if start >= len(u.trajectory):
        raise ValueError("trajectory.start is beyond available frames")
    u.trajectory[start]
    with tempfile.TemporaryDirectory(prefix="mdscope-solvation-bench-") as scratch:
        obs = SolvationObserver(u, cfg, scratch)
        topology_audit = audit(obs)
        init_t = time.perf_counter()
        update_time = 0.
        completed = 0
        for ts in u.trajectory[start:stop:stride]:
            if completed >= frames:
                break
            stamp = time_ps(ts, completed, cfg)
            begin = time.perf_counter()
            obs.update(ts, completed, stamp)
            update_time += time.perf_counter()-begin
            completed += 1
            if completed % int(tr.get("chunk_frames", 1000)) == 0:
                obs.flush((completed-1)//int(tr["chunk_frames"]))
        if completed % int(tr.get("chunk_frames", 1000)):
            obs.flush(completed//int(tr["chunk_frames"]))
        before_finalize = time.perf_counter()
        obs.finalize()
        elapsed = time.perf_counter()-init_t
        finalize_time = time.perf_counter()-before_finalize
        sizes = {str(p.relative_to(scratch)): p.stat().st_size
                 for p in Path(scratch).rglob("*") if p.is_file()}
        outcome = {
            "topology": str(rep.topology),
            "trajectory": str(rep.trajectory),
            "system": rep.system, "replicate": rep.label, "preset": preset,
            "analyzed_frames": completed,
            "time_s_total": round(elapsed, 4),
            "time_s_updates": round(update_time, 4),
            "time_s_finalize": round(finalize_time, 4),
            "frames_per_second_updates": (round(completed/update_time, 3)
                                          if update_time else None),
            "peak_rss_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
                            if os.uname().sysname == "Linux" else None,
            "output_bytes": sizes, "total_output_bytes": sum(sizes.values()),
            "audit": topology_audit,
            "limitations": [
                "Python process peak RSS includes Universe and module imports.",
                "SASA and other MDscope observers are not included in this benchmark.",
                "Short-frame results may underestimate long-run events or finalization cost.",
                "Existing source trajectory sampling limits residence-time interpretation."
            ],
        }
    if output:
        Path(output).write_text(json.dumps(outcome, indent=2))
    return outcome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--system", required=True)
    parser.add_argument("--rep", required=True)
    parser.add_argument("--frames", type=int, default=100)
    parser.add_argument("--preset", choices=("standard", "detailed", "kinetics"),
                        default="standard")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cfg = load_config(args.config)
    print(json.dumps(benchmark(cfg, args.config, system=args.system,
          replicate=args.rep, frames=args.frames, preset=args.preset,
          output=args.output), indent=2))


if __name__ == "__main__":
    main()
