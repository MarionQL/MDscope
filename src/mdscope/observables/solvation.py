"""Shared periodic solvation neighbor search, compact RDFs and optional identity events.

Bulk-normalized radial distributions are NOT corrected for inaccessible volume.
All distances are Å and times ps; the actual water/ion atom selections are
explicit to avoid guessing CHARMM lipid and glycan atom chemistry.
"""
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors

from ..config import atomic_json, atomic_parquet, ensure_dir

PRESETS = {"standard": (False, False), "detailed": (True, True),
           "kinetics": (True, True)}
COUNT_COLUMNS = ["frame", "time_ps", "solvent", "target_id", "target_class",
                 "site_class", "n_solvent"]
EVENT_COLUMNS = ["frame", "time_ps", "solvent", "target_id", "solvent_id"]
BRIDGE_COLUMNS = ["frame", "time_ps", "solvent", "solvent_id", "target_a", "target_b"]


def rdf_from_counts(counts, edges, density_exposure):
    """Number of observed pairs / expected pairs at BULK solvent density."""
    shell_volume = (4 * np.pi / 3) * np.diff(np.asarray(edges, float) ** 3)
    if density_exposure <= 0:
        return np.full(len(counts), np.nan)
    return np.asarray(counts, float) / (shell_volume * density_exposure)


def contact_episodes(events, observed_frames):
    """Continuous *sampled-frame* contact runs, preserving boundary censoring.

    observed_frames is the ordered (frame, time_ps) sequence. Missing contact
    records represent absent contacts, not unknown observations.
    """
    columns = ["solvent", "target_id", "solvent_id", "start_frame", "end_frame",
               "duration_ps", "left_censored", "right_censored"]
    if not observed_frames:
        return pd.DataFrame(columns=columns)
    frames = [int(f) for f, _ in observed_frames]
    times = {int(f): float(t) for f, t in observed_frames}
    order = {f: i for i, f in enumerate(frames)}
    pairs = defaultdict(set)
    for f, solvent, target, molecule in events:
        pairs[(solvent, target, molecule)].add(int(f))
    rows = []
    for (solvent, target, molecule), active in pairs.items():
        previous = start = None
        for frame in sorted(active, key=order.__getitem__):
            if start is None:
                start = previous = frame
            elif order[frame] == order[previous] + 1:
                previous = frame
            else:
                rows.append((solvent, target, molecule, start, previous,
                             times[previous] - times[start], start == frames[0],
                             previous == frames[-1]))
                start = previous = frame
        rows.append((solvent, target, molecule, start, previous,
                     times[previous] - times[start], start == frames[0],
                     previous == frames[-1]))
    return pd.DataFrame(rows, columns=columns)


def _targets(u, spec):
    result = []
    seen = set()
    for kind, groups in (spec.get("targets") or {}).items():
        for label, definition in groups.items():
            selection = definition if isinstance(definition, str) else definition["selection"]
            grouping = "residue" if isinstance(definition, str) else definition.get("group_by", "residue")
            atoms = u.select_atoms(selection)
            if not len(atoms):
                raise ValueError(f"solvation target {kind}/{label} has no atoms")
            if grouping == "residue":
                chunks = [(str(res.ix), atoms[atoms.resindices == res.ix])
                          for res in atoms.residues]
            elif grouping == "selection":
                chunks = [("all", atoms)]
            elif grouping == "fragment":
                try:
                    chunks = [(str(min(f.indices)), f) for f in atoms.fragments]
                except Exception as exc:
                    raise ValueError("solvation fragment grouping requires topology bonds") from exc
            else:
                raise ValueError("solvation group_by: use residue, fragment or selection")
            for ident, group in chunks:
                overlap = seen.intersection(map(int, group.indices))
                if overlap:
                    raise ValueError("solvation target selections overlap")
                seen.update(map(int, group.indices))
                result.append((f"{kind}/{label}/{ident}", kind, label, group))
    if not result:
        raise ValueError("solvation.targets must contain nonempty target groups")
    return result


class SolvationObserver:
    """Fused observer. Histograms and compact tables are checkpointed per chunk."""

    def __init__(self, universe, cfg, repout):
        self.u = universe
        self.spec = cfg["analyses"]["solvation"]
        self.repout = Path(repout)
        preset = self.spec.get("preset", "standard")
        if preset not in PRESETS:
            raise ValueError("solvation.preset: expected standard, detailed or kinetics")
        events, bridges = PRESETS[preset]
        opts = self.spec.get("options") or {}
        self.events_on = bool(opts.get("events", events))
        self.bridges_on = bool(opts.get("bridges", bridges))
        if self.bridges_on:
            self.events_on = True
        self.targets = _targets(universe, self.spec)
        indices = np.concatenate([g.indices for _, _, _, g in self.targets])
        self.site_atoms = universe.atoms[indices]
        self.target_of_site = np.concatenate(
            [np.full(len(g), i, dtype=np.int32)
             for i, (_, _, _, g) in enumerate(self.targets)])
        self.solvents = {}
        for name, setting in (("water", "water_selection"), ("ions", "ion_selection")):
            if self.spec.get(setting):
                group = universe.select_atoms(self.spec[setting])
                if not len(group) or len(np.unique(group.resindices)) != len(group):
                    raise ValueError(f"solvation.{setting}: choose exactly one atom per water/ion")
                self.solvents[name] = group
        if not self.solvents:
            raise ValueError("solvation requires water_selection and/or ion_selection")
        radius = float(self.spec.get("search_radius_A", 8.0))
        width = float(self.spec.get("rdf_bin_width_A", 0.2))
        if not np.isfinite(radius) or not np.isfinite(width) or min(radius, width) <= 0:
            raise ValueError("solvation RDF radius and bin width must be positive")
        self.edges = np.linspace(0, radius, int(np.ceil(radius / width)) + 1)
        self.shells = {"water": 3.5, "ions": 4.0, **(self.spec.get("shell_cutoffs_A") or {})}
        if any(not np.isfinite(float(x)) or not 0 < float(x) <= radius for x in self.shells.values()):
            raise ValueError("solvation shell cutoffs must be positive and <= search radius")
        self.hist = defaultdict(lambda: np.zeros(len(self.edges) - 1, np.int64))
        self.exposure = defaultdict(float)
        self.rows, self.events, self.bridges, self.frame_rows = [], [], [], []
        self.frames = 0

    def update(self, ts, frame_index, time_ps):
        box = np.asarray(self.u.dimensions, float)
        if box.shape != (6,) or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0):
            raise ValueError("solvation needs a valid periodic box")
        volume = abs(float(np.linalg.det(triclinic_vectors(box))))
        if volume <= 0:
            raise ValueError("solvation box volume must be positive")
        self.frame_rows.append((int(frame_index), float(time_ps)))
        for solvent, atoms in self.solvents.items():
            pair, distances = capped_distance(
                self.site_atoms.positions, atoms.positions,
                max_cutoff=float(self.edges[-1]), box=box, return_distances=True)
            indices = self.target_of_site[pair[:, 0]] if len(pair) else np.empty(0, int)
            shell = distances <= float(self.shells[solvent])
            contacts = defaultdict(set)
            for ti, si in zip(indices[shell], pair[shell, 1]):
                contacts[int(ti)].add(int(atoms.indices[si]))
            for i, (target, kind, label, group) in enumerate(self.targets):
                key = (solvent, kind, label)
                self.hist[key] += np.histogram(distances[indices == i], bins=self.edges)[0]
                self.exposure[key] += len(group) * len(atoms) / volume
                self.rows.append((int(frame_index), float(time_ps), solvent,
                                  target, kind, label, len(contacts[i])))
                if self.events_on:
                    self.events.extend((int(frame_index), float(time_ps), solvent,
                                        target, solvent_id)
                                       for solvent_id in sorted(contacts[i]))
            if self.bridges_on:
                inverse = defaultdict(set)
                for target_index, solvent_ids in contacts.items():
                    for sid in solvent_ids:
                        inverse[sid].add(target_index)
                for sid, target_indices in inverse.items():
                    for left, right in combinations(sorted(target_indices), 2):
                        self.bridges.append((int(frame_index), float(time_ps), solvent, sid,
                                             self.targets[left][0], self.targets[right][0]))
        self.frames += 1

    def flush(self, part):
        root = ensure_dir(self.repout / "cache" / "solvation")
        for name, values, columns in (
            ("counts", self.rows, COUNT_COLUMNS),
            ("events", self.events, EVENT_COLUMNS),
            ("bridges", self.bridges, BRIDGE_COLUMNS),
            ("frames", self.frame_rows, ["frame", "time_ps"]),
        ):
            if values or name in {"counts", "frames"}:
                atomic_parquet(root / f"{name}-{part:06d}.parquet",
                               pd.DataFrame(values, columns=columns))
            values.clear()
        history = [(solvent, kind, label, *map(int, bins), float(self.exposure[(solvent, kind, label)]))
                   for (solvent, kind, label), bins in sorted(self.hist.items())]
        atomic_parquet(root / f"hist-{part:06d}.parquet",
            pd.DataFrame(history, columns=["solvent", "target_class", "site_class"]
                         + [f"bin_{i}" for i in range(len(self.edges)-1)] + ["exposure"]))
        self.hist.clear()
        self.exposure.clear()

    def prune_from(self, part):
        for path in (self.repout / "cache" / "solvation").glob("*.parquet"):
            try:
                number = int(path.stem.rsplit("-", 1)[1])
            except (IndexError, ValueError):
                continue
            if number >= part:
                path.unlink()

    def finalize(self):
        cache = self.repout / "cache" / "solvation"
        dest = ensure_dir(self.repout / "extensions" / "solvation")
        def read(name):
            files = sorted(cache.glob(f"{name}-*.parquet"))
            return pd.concat((pd.read_parquet(p) for p in files), ignore_index=True) if files else pd.DataFrame()
        for name in ("counts", "events", "bridges"):
            df = read(name)
            if not df.empty or name == "counts":
                atomic_parquet(dest / f"{name}.parquet", df)
        hist = read("hist")
        rows = []
        if not hist.empty:
            for (solvent, kind, label), data in hist.groupby(["solvent", "target_class", "site_class"]):
                counts = data[[f"bin_{i}" for i in range(len(self.edges)-1)]].sum(axis=0).to_numpy(dtype=float)
                g = rdf_from_counts(counts, self.edges, float(data.exposure.sum()))
                rows.extend(zip([solvent]*len(g), [kind]*len(g), [label]*len(g),
                                self.edges[:-1], self.edges[1:], counts, g))
        atomic_parquet(dest / "rdf.parquet", pd.DataFrame(rows, columns=[
            "solvent", "target_class", "site_class", "r_inner_A", "r_outer_A",
            "pair_count", "g_bulk"]))
        if self.events_on:
            frames = read("frames").drop_duplicates("frame").sort_values("frame")
            events = read("events")
            episodes = contact_episodes(
                events[["frame", "solvent", "target_id", "solvent_id"]].itertuples(index=False, name=None)
                if not events.empty else [],
                list(frames[["frame", "time_ps"]].itertuples(index=False, name=None)))
            atomic_parquet(dest / "episodes.parquet", episodes)
        atomic_json(dest / "metadata.json", {
            "preset": self.spec.get("preset", "standard"), "edges_A": self.edges.tolist(),
            "shell_cutoffs_A": self.shells, "events": self.events_on,
            "bridges": self.bridges_on,
            "rdf_limitation": "Bulk-density normalized; no excluded-volume correction. Evaluate alongside SASA.",
            "residence_limitation": "Episodes represent consecutive sampled contacts; temporal resolution is trajectory frame interval.",
        })
        return {"frames": len(read("frames")), "fused": True}
