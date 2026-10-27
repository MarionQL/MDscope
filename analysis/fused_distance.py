"""Chunked residue-to-membrane minimum distances inside the shared frame loop.

Exact distance search uses the existing hybrid capped-distance/blockwise
algorithm; only trajectory iteration and checkpoint storage are changed.
"""
from __future__ import annotations

from pathlib import Path
import re
import time

import MDAnalysis as mda
import numpy as np

from .config import analysis_runs, time_ps
from .membrane_order import _atomic_npz
from .progress import track
from .residue_membrane_distance import (
    DEFAULT_LIPID_RESNAMES, _minimum_by_residue_hybrid,
    _valid_box, _without_hydrogens, run_residue_membrane_distance,
)


def enabled_names(cfg):
    return tuple(run.run_id for run in analysis_runs(cfg, "residue_membrane_distance"))


class FusedDistances:
    def __init__(self, universe, cfg, repout):
        self.u, self.cfg, self.repout = universe, cfg, Path(repout)
        self.cache = self.repout / "cache" / "fused_distance"
        self.cache.mkdir(parents=True, exist_ok=True)
        self.runs = {}
        self.frames = 0
        self.compute_seconds = 0.
        for run in analysis_runs(cfg, "residue_membrane_distance"):
            spec = run.spec
            selection = str(spec.get("selection") or cfg["selections"]["protein"])
            selected = universe.select_atoms(selection)
            if not len(selected):
                raise ValueError(
                    "residue_membrane_distance.selection matched zero atoms: " + selection
                )
            query = selected.residues.atoms if spec.get("expand_residue_atoms", True) else selected
            query = _without_hydrogens(query, bool(spec.get("exclude_hydrogens", True)))
            membrane_selection = str(spec.get("membrane_selection") or (
                "resname " + " ".join(map(str, spec.get("lipid_resnames") or DEFAULT_LIPID_RESNAMES))
            ))
            membrane = _without_hydrogens(
                universe.select_atoms(membrane_selection),
                bool(spec.get("exclude_hydrogens", True))
            )
            if not len(query) or not len(membrane):
                raise ValueError("residue_membrane_distance query and membrane selections must both contain atoms")
            residues = query.residues
            lookup = {int(res.ix): i for i, res in enumerate(residues)}
            mapping = np.fromiter(
                (lookup[int(ix)] for ix in query.resindices), dtype=int, count=len(query)
            )
            self.runs[run.run_id] = {
                "spec": spec, "query": query, "membrane": membrane,
                "mapping": mapping, "nres": len(residues),
                "values": [], "times": [],
            }

    def update(self, ts, frame_index):
        start = time.perf_counter()
        t = time_ps(ts, frame_index, self.cfg)
        box = _valid_box(self.u)
        for entry in self.runs.values():
            s = entry["spec"]
            values = _minimum_by_residue_hybrid(
                entry["query"], entry["membrane"], entry["mapping"],
                entry["nres"], box,
                s.get("initial_search_cutoff_A", 6.0),
                s.get("search_growth_factor", 2.0),
                s.get("query_block_atoms", 512),
                s.get("fallback_block_atoms", 4096),
            )
            entry["values"].append(values.astype(np.float32))
            entry["times"].append(t)
        self.frames += 1
        self.compute_seconds += time.perf_counter() - start

    def flush(self, part):
        for run_id, entry in self.runs.items():
            directory = self.cache / run_id
            directory.mkdir(parents=True, exist_ok=True)
            _atomic_npz(
                directory / f"part-{part:06d}.npz",
                time_ps=np.asarray(entry["times"], dtype=float),
                values=np.asarray(entry["values"], dtype=np.float32).reshape(
                    (-1, entry["nres"])
                ),
            )
            entry["times"].clear()
            entry["values"].clear()

    def prune_from(self, next_part):
        for path in self.cache.rglob("part-*.npz"):
            match = re.search(r"part-(\d+)", path.name)
            if match and int(match.group(1)) >= next_part:
                path.unlink()

    def finalize(self, rep):
        result = {}
        for run_id, entry in self.runs.items():
            result[run_id] = run_residue_membrane_distance(
                rep, self.cfg, self.repout,
                spec=entry["spec"], run_id=run_id,
                precomputed_parts=self.cache / run_id,
            )
        from .config import atomic_json
        atomic_json(self.repout / "fused_distance_timings.json", {
            "enabled_instances": list(self.runs),
            "frames_processed_this_invocation": self.frames,
            "compute_seconds": self.compute_seconds,
        })
        return result


def run_fused_distances(rep, cfg, repout):
    """Fallback for systems with no regular streaming stage."""
    if not enabled_names(cfg):
        return {}
    tr = cfg["trajectory"]
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observer = FusedDistances(u, cfg, repout)
    observer.prune_from(0)
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    for i, ts in enumerate(track(
        frames, total=len(frames), desc=f"{rep.system}/{rep.label}: fused distances",
        enabled=cfg.get("progress", {}).get("enabled", True),
    )):
        observer.update(ts, i)
        if (i + 1) % chunk == 0:
            observer.flush(i // chunk)
    if observer.frames % chunk:
        observer.flush(observer.frames // chunk)
    return observer.finalize(rep)
