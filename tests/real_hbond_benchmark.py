"""Release-data parity and speed benchmark for opt-in protein H-bond selection.

Run: python tests/real_hbond_benchmark.py
Requires both paired TPR/XTC files from GitHub release test-data-v1.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from mdscope.pipeline.fused_extensions import FusedExtensions

DATA = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data"))
OUT = Path(os.environ.get("MDSCOPE_HBOND_BENCH_OUTPUT", "real_test_results/hbond_benchmark"))
SYSTEMS = {
    "carb": ("carb-sdc4-step7_production.tpr", "carb-sdc4-rep1_test.xtc"),
    "nocarb": ("nocarb-sdc4-step7_production.tpr", "nocarb-sdc4-rep1_test.xtc"),
}


def detector(universe, optimized, path):
    cfg = {
        "selections": {"protein": "protein"},
        "analyses": {
            "hydrogen_bonds": {
                "enabled": True, "mode": "protein_protein",
                "update_selections": False,
                "restrict_search_to_protein": optimized,
                "d_h_cutoff_A": 1.2, "d_a_cutoff_A": 3.0,
                "angle_cutoff_deg": 150.0,
            }
        },
    }
    observer = FusedExtensions(universe, cfg, path)
    return observer.hb["hydrogen_bonds"]["hba"]


def run_frames(hba, u, indices):
    start = time.perf_counter()
    for frame in indices:
        ts = u.trajectory[int(frame)]
        hba._ts = ts
        hba._single_frame()
    elapsed = time.perf_counter() - start
    hba._conclude()
    events = np.asarray(hba.results.hbonds, dtype=float).reshape(-1, 6)
    if len(events):
        events = events[np.lexsort((events[:, 3], events[:, 2], events[:, 1], events[:, 0]))]
    return elapsed, events


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    n = int(os.environ.get("MDSCOPE_HBOND_BENCH_FRAMES", "8"))
    assert n >= 2
    reports = {}
    for name, (top, xtc) in SYSTEMS.items():
        u = mda.Universe(str(DATA / top), str(DATA / xtc))
        indices = np.linspace(0, len(u.trajectory) - 1, min(n, len(u.trajectory)), dtype=int)
        timings = {}
        events = {}
        for optimized in (False, True):
            key = "protein_only" if optimized else "baseline"
            hba = detector(u, optimized, OUT / name / key)
            timings[key], events[key] = run_frames(hba, u, indices)
        baseline, candidate = events["baseline"], events["protein_only"]
        # Exact atom/frame identity and numerical geometry parity are mandatory.
        assert baseline.shape == candidate.shape, (name, baseline.shape, candidate.shape)
        np.testing.assert_array_equal(baseline[:, :4], candidate[:, :4],
                                      err_msg=f"{name}: event identities differ")
        np.testing.assert_allclose(baseline[:, 4:], candidate[:, 4:],
                                   rtol=1e-6, atol=1e-5,
                                   err_msg=f"{name}: distances or angles differ")
        reports[name] = {
            "frames": indices.tolist(),
            "events": int(len(baseline)),
            "seconds": timings,
            "speedup": timings["baseline"] / timings["protein_only"],
            "parity": "passed",
        }
        print(json.dumps({name: reports[name]}, indent=2), flush=True)
    (OUT / "results.json").write_text(json.dumps(reports, indent=2) + "\n")


if __name__ == "__main__":
    main()
