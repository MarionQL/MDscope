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
    cutoff = cfg["analyses"].get("sdc4_mechanisms", {}).get("hs_sugar_water_cutoff_A")
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
