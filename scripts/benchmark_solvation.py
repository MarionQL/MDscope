"""Opt-in bounded solvation neighbor-search microbenchmark (no input trajectories).

This measures periodic pair search + histogram/unique-pair grouping, not full
MDscope pipeline I/O or physically representative contacts. Run on HPC with
realistic atom counts and compare with actual replicate walltime/memory.
"""
import argparse
import json
import time
import tracemalloc

import numpy as np
from MDAnalysis.lib.distances import capped_distance


def run(*, sites, waters, ions, frames, box_A, cutoff_A, seed=2026):
    rng = np.random.default_rng(seed)
    box = np.array([box_A]*3 + [90.]*3)
    site_xyz = rng.uniform(0, box_A, (sites, 3)).astype(np.float32)
    solvents = {
        "water": rng.uniform(0, box_A, (waters, 3)).astype(np.float32),
        "ions": rng.uniform(0, box_A, (ions, 3)).astype(np.float32),
    }
    tracemalloc.start()
    started = time.perf_counter()
    pairs_seen = 0
    for _ in range(frames):
        # Move positions slightly to avoid timing an identical array.
        for xyz in solvents.values():
            xyz[:] = (xyz + .1) % box_A
            pairs, distances = capped_distance(
                site_xyz, xyz, max_cutoff=cutoff_A, box=box,
                return_distances=True)
            pairs_seen += len(pairs)
            if len(pairs):
                # Representative class-bin aggregation, without dataframe I/O.
                histogram = np.histogram(distances, bins=40, range=(0, cutoff_A))
                assert histogram[0].sum() == len(distances)
    seconds = time.perf_counter()-started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return {
        "frames": frames, "sites": sites, "waters": waters, "ions": ions,
        "cutoff_A": cutoff_A, "pairs": pairs_seen,
        "seconds": round(seconds, 4),
        "frames_per_second": round(frames/seconds, 3),
        "python_tracemalloc_peak_MB": round(peak / 1e6, 2),
        "caution": "tracemalloc excludes some native allocations and input trajectory I/O"
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sites", type=int, default=2000)
    parser.add_argument("--waters", type=int, default=100000)
    parser.add_argument("--ions", type=int, default=500)
    parser.add_argument("--frames", type=int, default=25)
    parser.add_argument("--box-A", type=float, default=160.)
    parser.add_argument("--cutoff-A", type=float, default=8.)
    args = parser.parse_args()
    if min(args.sites, args.waters, args.frames, args.box_A, args.cutoff_A) <= 0 or args.ions < 0:
        parser.error("all positive dimensions required (ions may be zero)")
    print(json.dumps(run(
        sites=args.sites, waters=args.waters, ions=args.ions,
        frames=args.frames, box_A=args.box_A, cutoff_A=args.cutoff_A),
        indent=2))
