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



def water_donor_hbonds(site_xyz, water_xyz, water_hydrogens, pairs, box,
                       max_da_A=3.5, min_angle_deg=150.):
    """Return water-index/site-index/H-index triples satisfying D-H...A.

    Explicit water O-H bonds and target-acceptor atom selections are required.
    This deliberately does not classify target-donor/water-acceptor contacts.
    """
    from MDAnalysis.lib.distances import minimize_vectors
    pairs = np.asarray(pairs, dtype=int).reshape(-1, 2)
    if not len(pairs):
        return np.empty((0, 3), dtype=int)
    donor = np.asarray(water_xyz, float)[pairs[:, 1]]
    acceptor = np.asarray(site_xyz, float)[pairs[:, 0]]
    da = np.linalg.norm(minimize_vectors(acceptor-donor, box), axis=1)
    output = []
    for hi in (0, 1):
        hyd = np.asarray(water_hydrogens, float)[pairs[:, 1], hi]
        dh = minimize_vectors(donor-hyd, box)
        ah = minimize_vectors(acceptor-hyd, box)
        norms = np.linalg.norm(dh, axis=1)*np.linalg.norm(ah, axis=1)
        cosine = np.divide(np.einsum("ij,ij->i", dh, ah), norms,
                           out=np.zeros(len(da)), where=norms > 1e-12)
        angle = np.degrees(np.arccos(np.clip(cosine, -1., 1.)))
        output.extend((int(pairs[i,1]), int(pairs[i,0]), hi)
                      for i in np.flatnonzero((da <= max_da_A) &
                                               (angle >= min_angle_deg)))
    return np.asarray(output, dtype=int).reshape(-1, 3)


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
        self.water_hbond_bridges = bool(opts.get("water_hbond_bridges", False))
        if self.water_hbond_bridges:
            self.bridges_on = True
        if self.bridges_on:
            self.events_on = True
        self.targets = _targets(universe, self.spec)
        indices = np.concatenate([g.indices for _, _, _, g in self.targets])
        self.site_atoms = universe.atoms[indices]
        self.target_of_site = np.concatenate(
            [np.full(len(g), i, dtype=np.int32)
             for i, (_, _, _, g) in enumerate(self.targets)])
        self.solvents = {}
        selectors = {}
        if self.spec.get("water_selection"):
            selectors["water"] = self.spec["water_selection"]
        ion_specs = self.spec.get("ion_selections") or {}
        if ion_specs and self.spec.get("ion_selection"):
            raise ValueError("solvation: specify ion_selection OR ion_selections")
        if ion_specs:
            for label, selection in ion_specs.items():
                if not label or "/" in str(label):
                    raise ValueError("solvation ion species labels must be nonempty and exclude '/'")
                selectors["ions/" + str(label)] = selection
        elif self.spec.get("ion_selection"):
            selectors["ions"] = self.spec["ion_selection"]
        assigned = set()
        for name, query in selectors.items():
            group = universe.select_atoms(query)
            if not len(group) or len(np.unique(group.resindices)) != len(group):
                raise ValueError(f"solvation.{name}: select exactly one atom per solvent entity")
            overlap = set(map(int, group.indices)) & assigned
            if overlap:
                raise ValueError(f"solvation solvent species selections overlap: {name}")
            assigned.update(map(int, group.indices))
            self.solvents[name] = group
        if not self.solvents:
            raise ValueError("solvation requires water_selection and/or ion_selection")
        self.hbond_rows = []
        self.water_hydrogens = None
        self.hbond_da_A = float(opts.get("hbond_da_cutoff_A", 3.5))
        self.hbond_angle = float(opts.get("hbond_angle_deg", 150.))
        if self.water_hbond_bridges:
            if "water" not in self.solvents:
                raise ValueError("water H-bond bridges require water_selection")
            if not (0 < self.hbond_da_A <= float(self.spec.get("search_radius_A", 8.))
                    and 0 < self.hbond_angle <= 180):
                raise ValueError("invalid water H-bond bridge thresholds")
            ids = []
            for oxygen in self.solvents["water"]:
                try:
                    bonded = oxygen.bonded_atoms
                except Exception as exc:
                    raise ValueError("water H-bond bridging requires O-H topology bonds") from exc
                hydrogen = [int(atom.index) for atom in bonded
                            if str(atom.name).upper().lstrip("0123456789").startswith("H")]
                if len(hydrogen) != 2:
                    raise ValueError("water oxygen must have precisely two bonded H atoms")
                ids.extend(hydrogen)
            self.water_hydrogens = self.u.atoms[ids]

        radius = float(self.spec.get("search_radius_A", 8.0))
        width = float(self.spec.get("rdf_bin_width_A", 0.2))
        if not np.isfinite(radius) or not np.isfinite(width) or min(radius, width) <= 0:
            raise ValueError("solvation RDF radius and bin width must be positive")
        self.edges = np.linspace(0, radius, int(np.ceil(radius / width)) + 1)
        self.shells = {"water": 3.5, "ions": 4.0, **(self.spec.get("shell_cutoffs_A") or {})}
        self.shells.update({key: self.shells.get(key, self.shells["ions"])
                            for key in self.solvents if key.startswith("ions/")})
        if any(not np.isfinite(float(x)) or not 0 < float(x) <= radius for x in self.shells.values()):
            raise ValueError("solvation shell cutoffs must be positive and <= search radius")
        # Optional local 3D WATER density maps. The 3 anchor atoms define a
        # co-moving orthonormal frame for a named target, not a global grid.
        from MDAnalysis.lib.distances import minimize_vectors
        self._minimize_vectors = minimize_vectors
        self.grid_frames = {}
        self.grid_width = float(self.spec.get("grid_width_A", 12.0))
        self.grid_spacing = float(self.spec.get("grid_spacing_A", 0.75))
        if self.spec.get("grid_frames"):
            if "water" not in self.solvents:
                raise ValueError("solvation 3D hydration maps require water_selection")
            if not 0 < self.grid_spacing <= self.grid_width <= 24.0:
                raise ValueError("solvation grid width must be <=24 Å and spacing positive")
            n = int(np.ceil(self.grid_width / self.grid_spacing))
            if n > 64:
                raise ValueError("solvation grid exceeds 64 bins per axis")
            self.grid_edges = np.linspace(-self.grid_width/2, self.grid_width/2, n+1)
            known = {target for target, _, _, _ in self.targets}
            for target, selection in self.spec["grid_frames"].items():
                if target not in known:
                    raise ValueError(f"unknown 3D hydration target {target}")
                frame = universe.select_atoms(selection)
                if len(frame) != 3:
                    raise ValueError(f"3D hydration frame {target} needs exactly 3 ordered atoms")
                target_atoms = next(g for name, _, _, g in self.targets if name == target)
                if int(frame.indices[0]) not in set(map(int, target_atoms.indices)):
                    raise ValueError("3D grid origin must belong to the target site selection")
                self.grid_frames[target] = frame
        self.grid_counts = {key: np.zeros((len(self.grid_edges)-1,)*3, dtype=np.int64)
                            for key in self.grid_frames} if self.grid_frames else {}
        self.grid_exposure = 0.0
        self.hist = defaultdict(lambda: np.zeros(len(self.edges) - 1, np.int64))
        self.exposure = defaultdict(float)
        self.rows, self.events, self.bridges, self.frame_rows = [], [], [], []
        self.frames = 0
        self.class_names = list(dict.fromkeys(
            (kind, label) for _, kind, label, _ in self.targets))
        class_idx = {name: i for i, name in enumerate(self.class_names)}
        self.class_of_target = np.array(
            [class_idx[(kind, label)] for _, kind, label, _ in self.targets],
            dtype=np.int32)
        self.class_atoms = np.bincount(
            self.class_of_target,
            weights=[len(group) for _, _, _, group in self.targets],
            minlength=len(self.class_names))

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
            # Bin all site–solvent pairs in a single vectorized pass.
            # Per-target boolean masking here would scale as O(targets*pairs)
            # and is particularly expensive for hundreds of lipid molecules.
            radial_bin = np.searchsorted(self.edges, distances, side="right")-1
            valid = (radial_bin >= 0) & (radial_bin < len(self.edges)-1)
            n_bins = len(self.edges)-1
            flat = self.class_of_target[indices[valid]] * n_bins + radial_bin[valid]
            histogram = np.bincount(flat, minlength=len(self.class_names)*n_bins).reshape(
                len(self.class_names), n_bins)
            for class_index, (kind, label) in enumerate(self.class_names):
                key = (solvent, kind, label)
                self.hist[key] += histogram[class_index]
                self.exposure[key] += self.class_atoms[class_index] * len(atoms)/volume
            for i, (target, kind, label, group) in enumerate(self.targets):
                self.rows.append((int(frame_index), float(time_ps), solvent,
                                  target, kind, label, len(contacts[i])))
                if self.events_on:
                    self.events.extend((int(frame_index), float(time_ps), solvent,
                                        target, solvent_id)
                                       for solvent_id in sorted(contacts[i]))
            if solvent == "water" and self.grid_frames:
                # Reuse candidate pairs from the same periodic neighbor search.
                # For each target, map water positions to its local anchor frame.
                water_xyz = atoms.positions
                for target_i, (target_id, _, _, _) in enumerate(self.targets):
                    anchors = self.grid_frames.get(target_id)
                    if anchors is None:
                        continue
                    origin, second, third = anchors.positions
                    v1, v2 = self._minimize_vectors(
                        np.asarray([second-origin, third-origin], dtype=float), box)
                    norm1 = np.linalg.norm(v1)
                    if norm1 < 1e-7:
                        raise ValueError(f"degenerate 3D grid frame for {target_id}")
                    e1 = v1/norm1
                    v2 = v2-np.dot(v2,e1)*e1
                    norm2 = np.linalg.norm(v2)
                    if norm2 < 1e-7:
                        raise ValueError(f"collinear 3D grid frame for {target_id}")
                    e2 = v2/norm2
                    rotation = np.column_stack((e1, e2, np.cross(e1,e2)))
                    js = np.unique(pair[indices == target_i, 1])
                    if len(js):
                        delta = self._minimize_vectors(
                            np.asarray(water_xyz[js]-origin,dtype=float),box)
                        local = delta @ rotation
                        # Inscribed sphere is fully covered by the original
                        # site-to-water search because origin belongs to target.
                        # Cube corners outside that sphere have no sampling
                        # guarantee and must never be interpreted as zero density.
                        local = local[np.linalg.norm(local, axis=1) <= self.grid_width/2]
                        self.grid_counts[target_id] += np.histogramdd(
                            local, bins=(self.grid_edges,)*3)[0].astype(np.int64)
                self.grid_exposure += len(atoms)/volume
            if solvent == "water" and self.water_hbond_bridges:
                candidate = pair[distances <= self.hbond_da_A]
                qualified = water_donor_hbonds(
                    self.site_atoms.positions, atoms.positions,
                    self.water_hydrogens.positions.reshape(-1, 2, 3),
                    candidate, box, max_da_A=self.hbond_da_A,
                    min_angle_deg=self.hbond_angle)
                by_water = defaultdict(lambda: defaultdict(set))
                for water_idx, site_idx, hydrogen_idx in qualified:
                    by_water[int(water_idx)][int(self.target_of_site[site_idx])].add(
                        int(hydrogen_idx))
                for water_idx, target_hydrogens in by_water.items():
                    for a, b in combinations(sorted(target_hydrogens), 2):
                        # Two distinct donor H atoms must participate.
                        if any(h1 != h2 for h1 in target_hydrogens[a]
                               for h2 in target_hydrogens[b]):
                            self.hbond_rows.append(
                                (int(frame_index), float(time_ps),
                                 int(atoms.indices[water_idx]),
                                 self.targets[a][0], self.targets[b][0]))
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
            ("water_hbond_bridges", self.hbond_rows,
             ["frame", "time_ps", "water_id", "target_a", "target_b"]),
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
        if self.grid_frames:
            records = []
            for target, grid in self.grid_counts.items():
                i,j,k = np.nonzero(grid)
                records.extend((target,int(x),int(y),int(z),int(grid[x,y,z]))
                               for x,y,z in zip(i,j,k))
                grid.fill(0)
            atomic_parquet(root / f"grid-{part:06d}.parquet", pd.DataFrame(
                records, columns=["target_id","ix","iy","iz","count"]))
            atomic_parquet(root / f"grid_exposure-{part:06d}.parquet",
                           pd.DataFrame({"exposure":[self.grid_exposure]}))
            self.grid_exposure = 0.
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
        for name in ("counts", "events", "bridges", "water_hbond_bridges"):
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
        if self.grid_frames:
            grids = read("grid")
            exposure = read("grid_exposure").exposure.sum()
            if not grids.empty:
                grids = grids.groupby(["target_id","ix","iy","iz"],as_index=False).agg(
                    count=("count","sum"))
                width = float(self.grid_edges[1]-self.grid_edges[0])
                for dim in "xyz":
                    grids[dim+"_A"] = self.grid_edges[grids["i"+dim].to_numpy(int)] + width/2
                grids["g_bulk"] = grids["count"]/(width**3*exposure)
                grids["valid_spherical_domain"] = True
            atomic_parquet(dest / "water_density_3d.parquet", grids)
        # Persist target-to-topology-residue membership for SASA reuse.
        membership = [(target, int(res.ix)) for target, _, _, group in self.targets
                      for res in group.residues]
        atomic_parquet(dest / "target_residues.parquet", pd.DataFrame(
            membership, columns=["target_id", "topology_residue_index"]))
        if self.events_on:
            frames = read("frames").drop_duplicates("frame").sort_values("frame")
            events = read("events")
            episodes = contact_episodes(
                events[["frame", "solvent", "target_id", "solvent_id"]].itertuples(index=False, name=None)
                if not events.empty else [],
                list(frames[["frame", "time_ps"]].itertuples(index=False, name=None)))
            atomic_parquet(dest / "episodes.parquet", episodes)
        from .solvation_postprocess import write_solvation_derivatives
        write_solvation_derivatives(dest, self.repout,
            sasa_run=self.spec.get("sasa_run", "default"))
        atomic_json(dest / "metadata.json", {
            "preset": self.spec.get("preset", "standard"), "edges_A": self.edges.tolist(),
            "shell_cutoffs_A": self.shells, "events": self.events_on,
            "3d_grid_frames": list(self.grid_frames),
            "3d_grid_domain": "Only observed voxels within an inscribed sphere around a target origin; other voxels are not sampled.",
            "3d_grid_edges_A": self.grid_edges.tolist() if self.grid_frames else None,
            "bridges": self.bridges_on,
            "water_hbond_bridges": self.water_hbond_bridges,
            "rdf_limitation": "Bulk-density normalized; no excluded-volume correction. Evaluate alongside SASA.",
            "residence_limitation": "Episodes represent consecutive sampled contacts; temporal resolution is trajectory frame interval.",
        })
        return {"frames": len(read("frames")), "fused": True}
