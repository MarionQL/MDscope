"""Position- and chemistry-resolved water coordination of attached HS."""
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from ..config import atomic_parquet, time_ps
from ..observables.solvation import unique_target_solvent_pairs
from .hs_topology import discover_covalent_hs

OUTPUT = Path("extensions/sdc4_mechanisms/hs_sugar_hydration.parquet")


def sugar_map(universe, spec):
    """Number each linear glycan from its protein bond, independent of resids."""
    chains, attachments, _ = discover_covalent_hs(universe, spec["hs_auto"])
    bonds = universe.bonds.indices
    resindices = universe.atoms.resindices
    left, right = bonds.T
    edges = resindices[bonds]
    edges = edges[edges[:, 0] != edges[:, 1]]
    rows = []
    atoms = []
    owners = []
    for chain, group in sorted(chains.items()):
        residues = set(group.resindices)
        attachment = attachments[chain].indices
        roots = set(resindices[right[np.isin(left, attachment)]]) | set(
            resindices[left[np.isin(right, attachment)]])
        roots &= residues
        chain_edges = edges[np.isin(edges, list(residues)).all(axis=1)]
        neighbors = {i: set() for i in residues}
        for a, b in chain_edges:
            neighbors[a].add(b)
            neighbors[b].add(a)
        if len(roots) != 1:
            raise ValueError(f"{chain}: expected one attachment sugar")
        current = roots.pop()
        visited = set()
        while current not in visited:
            visited.add(current)
            residue = universe.residues[current]
            sugar = group[group.resindices == current]
            oxygen = sugar[np.char.startswith(sugar.names.astype(str), "O")]
            if not len(oxygen):
                raise ValueError(f"{chain}: no oxygen atoms in {residue.resname} {residue.resid}")
            owners.extend([len(rows)] * len(oxygen))
            atoms.extend(oxygen.indices)
            rows.append((chain, len(visited), residue.resname, int(residue.resid), len(oxygen)))
            onward = neighbors[current] - visited
            if len(neighbors[current]) > (1 if len(visited) == 1 else 2):
                raise ValueError(f"{chain}: expected a linear HS chain")
            if not onward:
                break
            current = onward.pop()
        if visited != residues:
            raise ValueError(f"{chain}: incomplete glycan traversal")
    columns = ["hs_chain", "sugar_position", "sugar_type", "sugar_resid", "oxygen_count"]
    return pd.DataFrame(rows, columns=columns), universe.atoms[atoms], np.asarray(owners)


def run_hs_hydration(rep, cfg, repout, expected_times):
    spec = cfg["analyses"]["sdc4_mechanisms"]
    cutoff = float(spec["hs_sugar_water_cutoff_A"])
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    sugars, oxygens, owner = sugar_map(u, spec)
    water = u.select_atoms(spec["water_oxygen_selection"])
    if not len(water) or len(np.unique(water.resindices)) != len(water):
        raise ValueError("HS hydration requires one oxygen per water molecule")
    tr = cfg["trajectory"]
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    if not len(frames):
        raise ValueError("No sampled frames for HS hydration")
    times = np.empty(len(frames))
    counts = np.empty((len(frames), len(sugars)), dtype=np.int32)
    for i, ts in enumerate(frames):
        box = ts.dimensions
        if box is None or not np.isfinite(box).all() or np.any(box <= 0):
            raise ValueError(f"Missing or invalid periodic cell at frame {ts.frame}")
        pairs = capped_distance(oxygens.positions, water.positions, cutoff,
                                box=box, method="nsgrid", return_distances=False)
        counts[i], _ = unique_target_solvent_pairs(owner[pairs[:, 0]], pairs[:, 1], len(sugars))
        times[i] = time_ps(ts, i, cfg)
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("HS hydration requires finite increasing frame times")
    if times.shape != expected_times.shape or not np.allclose(times, expected_times, rtol=0, atol=1e-6):
        raise ValueError("HS hydration frames differ from completed replica times")
    table = pd.DataFrame({col: np.tile(sugars[col], len(frames)) for col in sugars})
    table["time_ps"] = np.repeat(times, len(sugars))
    table["water_cutoff_A"] = cutoff
    table["water_count"] = counts.ravel()
    atomic_parquet(repout / OUTPUT, table)


def temporal_hs_hydration_profiles(repout, cfg, expected_times=None):
    cutoff = cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("hs_sugar_water_cutoff_A")
    if cutoff is None:
        return {}
    path = repout / OUTPUT
    if not path.is_file():
        raise ValueError(f"Missing HS sugar hydration at {path}; run mdscope replicate first")
    rows = pd.read_parquet(path)
    if rows.empty or not rows.water_cutoff_A.eq(float(cutoff)).all():
        raise ValueError(f"HS sugar hydration cutoff changed or data empty at {path}")
    wide = rows.pivot(index="time_ps", columns=["hs_chain", "sugar_position", "sugar_type"],
                      values="water_count")
    times = wide.index.to_numpy(float)
    values = wide.to_numpy(float)
    if (not np.isfinite(times).all() or not np.isfinite(values).all()
            or np.any(values < 0) or np.any(np.diff(times) <= 0)
            or (expected_times is not None and (times.shape != expected_times.shape or not np.allclose(times, expected_times, rtol=0, atol=1e-6)))):
        raise ValueError(f"Invalid HS sugar hydration frame grid at {path}")
    out = {}
    for chain in wide.columns.get_level_values("hs_chain").unique():
        matrix = wide[chain]
        out[f"hs_water_by_position__{chain}"] = {
            "times": times, "values": matrix.to_numpy(float), "reducer": "mean",
            "positions": np.asarray([pos for pos, _ in matrix.columns]),
            "labels": np.asarray([f"{pos}:{kind}" for pos, kind in matrix.columns]),
        }
        types = matrix.columns.get_level_values("sugar_type")
        labels = sorted(set(types))
        out[f"hs_water_by_sugar_type__{chain}"] = {
            "times": times, "reducer": "mean", "positions": np.arange(len(labels)),
            "values": np.column_stack([matrix.loc[:, types == kind].mean(axis=1) for kind in labels]),
            "labels": np.asarray([
                f"{kind} ({','.join(str(pos) for pos, name in matrix.columns if name == kind)})"
                for kind in labels
            ]),
        }
    return out


class HydrationObserver:
    """Fused, restartable hydration counts for HS, protein, and membrane sites.

    A target is an individual HS sugar, protein residue or selected membrane
    residue. Aggregate counts are unions of water identities, never sums of
    target counts. All distances are in Å and frame times in ps.
    """

    TARGET_COLUMNS = [
        "compartment", "target_id", "chain", "sugar_position", "resname",
        "resid", "topology_residue_index", "n_sites", "water_cutoff_A",
    ]
    COUNT_COLUMNS = [
        "frame", "time_ps", "compartment", "target_id", "water_count",
        "water_cutoff_A",
    ]
    AGGREGATE_COLUMNS = [
        "frame", "time_ps", "compartment", "level", "group_id",
        "n_targets", "water_count", "water_cutoff_A",
    ]
    SUGAR_COLUMNS = [
        "hs_chain", "sugar_position", "sugar_type", "sugar_resid",
        "oxygen_count", "time_ps", "water_cutoff_A", "water_count",
    ]

    def __init__(self, universe, cfg, repout):
        self.u, self.cfg, self.repout = universe, cfg, Path(repout)
        analyses = cfg.get("analyses", {})
        spec = analyses.get("hydration", {})
        hs_spec = analyses.get("sdc4_mechanisms", {})
        hs_cutoff = hs_spec.get("hs_sugar_water_cutoff_A")
        self.hs_enabled = hs_cutoff is not None
        general = bool(spec.get("enabled", False))
        if not self.hs_enabled and not general:
            raise ValueError("HydrationObserver requires an enabled hydration analysis")
        water_sel = (spec.get("water_oxygen_selection")
                     or hs_spec.get("water_oxygen_selection")
                     or analyses.get("solvation", {}).get("water_selection"))
        if not water_sel:
            raise ValueError("Hydration requires water_oxygen_selection")
        self.water = universe.select_atoms(water_sel)
        if not len(self.water) or len(np.unique(self.water.resindices)) != len(self.water):
            raise ValueError("Hydration requires exactly one oxygen per water molecule")

        self.targets = []
        self.site_indices = []
        self.site_owner = []
        self.groups = {}
        self.sugars = None
        self.hs_indices = []
        self.hs_cutoff = float(hs_cutoff) if self.hs_enabled else None
        cutoff = float(spec.get("cutoff_A", 3.5))
        if general and (not np.isfinite(cutoff) or cutoff <= 0):
            raise ValueError("hydration.cutoff_A must be finite and positive")
        if self.hs_enabled and (not np.isfinite(self.hs_cutoff) or self.hs_cutoff <= 0):
            raise ValueError("hs_sugar_water_cutoff_A must be finite and positive")

        if self.hs_enabled:
            self.sugars, sites, owners = sugar_map(universe, hs_spec)
            for i, row in self.sugars.iterrows():
                indices = sites.indices[owners == i]
                target = self._add("hs", f"{row.hs_chain}:{int(row.sugar_position)}",
                                   row.hs_chain, int(row.sugar_position),
                                   row.sugar_type, int(row.sugar_resid),
                                   int(universe.atoms[indices[0]].resindex),
                                   indices, self.hs_cutoff)
                self.hs_indices.append(target)
                self._group("hs", "chain", row.hs_chain, target)
                self._group("hs", "sugar_type", f"{row.hs_chain}:{row.sugar_type}", target)
                self._group("hs", "total", "all", target)

        if general:
            protein_sel = spec.get("protein_selection") or cfg.get("selections", {}).get("protein", "protein")
            protein = universe.select_atoms(protein_sel).select_atoms("not name H*")
            if not len(protein):
                raise ValueError("Hydration protein selection has no heavy atoms")
            for res in protein.residues:
                atoms = protein[protein.resindices == res.ix]
                target = self._add("protein", str(res.ix), "", -1, res.resname,
                                   int(res.resid), int(res.ix), atoms.indices, cutoff)
                self._group("protein", "total", "all", target)
            membrane_sel = spec.get("membrane_selection")
            if membrane_sel:
                membrane = universe.select_atoms(membrane_sel).select_atoms("not name H*")
                if not len(membrane):
                    raise ValueError("Hydration membrane selection has no heavy atoms")
                for res in membrane.residues:
                    atoms = membrane[membrane.resindices == res.ix]
                    target = self._add("membrane", str(res.ix), "", -1, res.resname,
                                       int(res.resid), int(res.ix), atoms.indices, cutoff)
                    self._group("membrane", "lipid_type", res.resname, target)
                    self._group("membrane", "total", "all", target)
        if not self.targets:
            raise ValueError("Hydration selected no targets")
        if set(self.water.indices).intersection(self.site_indices):
            raise ValueError("Hydration targets overlap water selection")
        self.sites = universe.atoms[self.site_indices]
        self.owner = np.asarray(self.site_owner, dtype=np.intp)
        self.cutoffs = np.asarray([r["water_cutoff_A"] for r in self.targets], float)
        self.rows, self.aggregate_rows, self.hs_rows = [], [], []
        self.frames = 0
        self.last_time = None

    def _add(self, compartment, identity, chain, position, resname, resid,
             resindex, indices, cutoff):
        if not len(indices):
            raise ValueError(f"Hydration target {compartment}/{identity} is empty")
        target = len(self.targets)
        self.targets.append(dict(compartment=compartment, target_id=str(identity),
            chain=str(chain), sugar_position=int(position), resname=str(resname),
            resid=int(resid), topology_residue_index=int(resindex),
            n_sites=len(indices), water_cutoff_A=float(cutoff)))
        self.site_indices.extend(map(int, indices))
        self.site_owner.extend([target] * len(indices))
        return target

    def _group(self, compartment, level, group_id, target):
        self.groups.setdefault((compartment, level, str(group_id)), []).append(target)

    def update(self, ts, frame_index, time_ps):
        box = ts.dimensions
        if (box is None or np.asarray(box).shape != (6,)
                or not np.isfinite(box).all() or np.any(np.asarray(box) <= 0)):
            raise ValueError(f"Invalid hydration periodic box at frame {ts.frame}")
        if not np.isfinite(time_ps) or (
                self.last_time is not None and time_ps <= self.last_time):
            raise ValueError("Hydration requires finite, increasing sampled frame times")
        self.last_time = float(time_ps)
        pairs, distances = capped_distance(
            self.sites.positions, self.water.positions, float(self.cutoffs.max()),
            box=box, method="nsgrid", return_distances=True)
        if len(pairs):
            owners = self.owner[pairs[:, 0]]
            valid = distances <= self.cutoffs[owners]
            target_ids = owners[valid]
            water_ids = pairs[valid, 1]
        else:
            target_ids = np.empty(0, dtype=np.intp)
            water_ids = np.empty(0, dtype=np.intp)
        counts, unique_pairs = unique_target_solvent_pairs(
            target_ids, water_ids, len(self.targets))
        frame = int(frame_index)
        t = float(time_ps)
        self.rows.extend((frame, t, target["compartment"], target["target_id"],
                          int(counts[i]), target["water_cutoff_A"])
                         for i, target in enumerate(self.targets))
        for (compartment, level, label), indices in self.groups.items():
            # A water shared by two targets belongs to the union only once.
            mask = np.isin(unique_pairs[:, 0], indices)
            n_unique = len(np.unique(unique_pairs[mask, 1]))
            self.aggregate_rows.append(
                (frame, t, compartment, level, label, len(indices),
                 n_unique, float(self.cutoffs[indices[0]])))
        if self.hs_enabled:
            for index, sugar in zip(self.hs_indices, self.sugars.itertuples(index=False)):
                self.hs_rows.append((
                    sugar.hs_chain, int(sugar.sugar_position), sugar.sugar_type,
                    int(sugar.sugar_resid), int(sugar.oxygen_count),
                    t, self.hs_cutoff, int(counts[index])))
        self.frames += 1

    def flush(self, part):
        from ..config import ensure_dir
        cache = ensure_dir(self.repout / "cache" / "hydration")
        for name, rows, columns in (
            ("counts", self.rows, self.COUNT_COLUMNS),
            ("aggregates", self.aggregate_rows, self.AGGREGATE_COLUMNS),
            ("hs", self.hs_rows, self.SUGAR_COLUMNS),
        ):
            if rows:
                atomic_parquet(cache / f"{name}-{part:06d}.parquet",
                               pd.DataFrame(rows, columns=columns))
                rows.clear()

    def prune_from(self, part):
        cache = self.repout / "cache" / "hydration"
        for path in cache.glob("*.parquet"):
            try:
                number = int(path.stem.rsplit("-", 1)[1])
            except (IndexError, ValueError):
                continue
            if number >= part:
                path.unlink()

    @staticmethod
    def _consolidate(parts, destination):
        """Stream record batches to one atomic public Parquet product."""
        import os
        import tempfile
        import pyarrow as pa
        import pyarrow.parquet as pq
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not parts:
            raise ValueError(f"Missing hydration cache chunks for {destination}")
        with tempfile.NamedTemporaryFile(
                dir=destination.parent, suffix=".parquet", delete=False) as tmp:
            temporary = Path(tmp.name)
        try:
            schema = pq.read_schema(parts[0])
            with pq.ParquetWriter(temporary, schema) as writer:
                for part in parts:
                    for batch in pq.ParquetFile(part).iter_batches(batch_size=65536):
                        writer.write_table(pa.Table.from_batches([batch], schema=schema))
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)

    def finalize(self):
        cache = self.repout / "cache" / "hydration"
        dest = self.repout / "extensions" / "hydration"
        for name, output in (("counts", "counts.parquet"),
                             ("aggregates", "aggregate_counts.parquet")):
            self._consolidate(sorted(cache.glob(f"{name}-*.parquet")), dest / output)
        atomic_parquet(dest / "target_metadata.parquet",
                       pd.DataFrame(self.targets, columns=self.TARGET_COLUMNS))
        if self.hs_enabled:
            self._consolidate(
                sorted(cache.glob("hs-*.parquet")), self.repout / OUTPUT)
        return {"frames": self.frames, "fused": True,
                "targets": len(self.targets)}
