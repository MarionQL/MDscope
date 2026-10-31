"""Exact residue-to-membrane minimum distance, memory-mapped frame profiles."""
from __future__ import annotations
import hashlib
from pathlib import Path
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import USED_OPENMP,capped_distance,distance_array
from ..profile_common import _residue_label
from ..pipeline.progress import track
from ..config import DEFAULT_INSTANCE_ID,analysis_groups,analysis_runs,atomic_parquet,ensure_dir,profile_dir,time_ps
DEFAULT_LIPID_RESNAMES = [
    "POPC",
    "POPE",
    "POPI",
    "DOPC",
    "DPPC",
    "PSM",
    "CHL1",
    "BGAL",
    "GM1a",
    "CER1"
]


def _valid_box(universe):
    dims = np.asarray(
        universe.dimensions,
        float
    ) if universe.dimensions is not None else np.array([])
    return universe.dimensions if len(
        dims
    ) >= 3 and np.all(np.isfinite(dims[:3])) and np.all(dims[:3] > 0) else None


def _without_hydrogens(atoms, exclude: bool):
    if not exclude:
        return atoms
    names = np.char.upper(np.asarray(atoms.names, dtype=str))
    is_h = np.char.startswith(names, "H")
    for prefix in ("1H", "2H", "3H"):
        is_h |= np.char.startswith(names, prefix)
    return atoms[~is_h]


def _index_hash(atoms) -> str:
    return hashlib.sha256(
        np.asarray(atoms.indices, dtype=np.int64).tobytes()
    ).hexdigest()


def _minimum_by_residue_hybrid(
    query,
    membrane,
    atom_to_residue,
    n_residues,
    box,
    initial_cutoff_A,
    search_growth_factor,
    query_block_atoms,
    membrane_block_atoms,
):
    # First resolve residues that are already close to the membrane with
    # one inexpensive capped neighbor search. Any residue with a neighbor
    # inside this cutoff already has its exact minimum because all omitted
    # distances are necessarily larger than the cutoff.
    minima = np.full(n_residues, np.inf)
    cutoff = max(float(initial_cutoff_A), 0.1)
    query_block = max(1, int(query_block_atoms))
    membrane_block = max(1, int(membrane_block_atoms))
    query_indices = np.arange(len(query), dtype=int)

    # Bound the initial capped search too. This prevents a large temporary
    # pair list from being created for every query atom at once.
    for qstart in range(0, len(query), query_block):
        qindices = query_indices[qstart:qstart + query_block]
        pairs, distances = capped_distance(
            query.positions[qindices],
            membrane.positions,
            max_cutoff=cutoff,
            box=box,
            return_distances=True
        )
        if len(pairs):
            np.minimum.at(
                minima,
                atom_to_residue[qindices[pairs[:, 0].astype(int)]],
                np.asarray(distances, float)
            )

    # Give unresolved residues one expanded capped search before falling
    # back to the exact blockwise distance matrix. This makes the configured
    # search_growth_factor functional without allowing an unbounded pair list.
    unresolved = ~np.isfinite(minima)
    growth = max(float(search_growth_factor), 1.0)
    if np.any(unresolved) and growth > 1.0:
        expanded_cutoff = cutoff * growth
        active = query_indices[unresolved[atom_to_residue]]
        for qstart in range(0, len(active), query_block):
            qindices = active[qstart:qstart + query_block]
            pairs, distances = capped_distance(
                query.positions[qindices],
                membrane.positions,
                max_cutoff=expanded_cutoff,
                box=box,
                return_distances=True,
            )
            if len(pairs):
                np.minimum.at(
                    minima,
                    atom_to_residue[qindices[pairs[:, 0].astype(int)]],
                    np.asarray(distances, float),
                )

    # For residues still farther from the membrane, calculate the exact
    # minimum with MDAnalysis' optimized distance_array in fixed-size blocks.
    # Reuse the full-block result array to avoid repeated large allocations.
    unresolved = ~np.isfinite(minima)
    active = query_indices[unresolved[atom_to_residue]]
    backend = "OpenMP" if USED_OPENMP else "serial"
    for qstart in range(0, len(active), query_block):
        qindices = active[qstart:qstart + query_block]
        atom_min = np.full(len(qindices), np.inf)
        full_result = np.empty(
            (len(qindices), membrane_block),
            dtype=np.float64,
        )
        for mstart in range(0, len(membrane), membrane_block):
            mstop = min(mstart + membrane_block, len(membrane))
            membrane_count = mstop - mstart
            if membrane_count == membrane_block:
                result = full_result
            else:
                result = np.empty(
                    (len(qindices), membrane_count),
                    dtype=np.float64,
                )
            d = distance_array(
                query.positions[qindices],
                membrane.positions[mstart:mstop],
                box=box,
                result=result,
                backend=backend,
            )
            atom_min = np.minimum(atom_min, np.min(d, axis=1))
        np.minimum.at(minima, atom_to_residue[qindices], atom_min)
    return minima


def run_residue_membrane_distance(
    rep,
    cfg: dict,
    repout: Path,
    spec: dict | None = None,
    run_id: str = DEFAULT_INSTANCE_ID,
    precomputed_parts: Path | None = None,
) -> dict:
    spec = spec if spec is not None else cfg["analyses"]["residue_membrane_distance"]
    if not spec.get("enabled", False):
        return {"enabled": False}
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    selected = universe.select_atoms(selection)
    if not len(selected):
        raise ValueError(
            "residue_membrane_distance.selection matched zero atoms: "
            f"{selection}"
        )
    query = (
        selected.residues.atoms
        if spec.get("expand_residue_atoms", True)
        else selected
    )
    query = _without_hydrogens(query, bool(spec.get("exclude_hydrogens", True)))
    membrane_selection = str(
        spec.get("membrane_selection")
        or (
            "resname "
            + " ".join(
                map(
                    str,
                    spec.get("lipid_resnames") or DEFAULT_LIPID_RESNAMES,
                )
            )
        )
    )
    membrane = _without_hydrogens(
        universe.select_atoms(membrane_selection),
        bool(spec.get("exclude_hydrogens", True))
    )
    if not len(query) or not len(membrane):
        raise ValueError(
            "residue_membrane_distance query and membrane selections "
            "must both contain atoms"
        )
    residues = query.residues
    residue_lookup = {int(res.ix): i for i, res in enumerate(residues)}
    atom_to_residue = np.fromiter(
        (residue_lookup[int(ix)] for ix in query.resindices),
        dtype=int,
        count=len(query)
    )
    out = ensure_dir(profile_dir(repout, "residue_membrane_distance", run_id))
    atomic_parquet(out / "selection_audit.parquet", pd.DataFrame([{
        "analysis_instance": str(run_id),
        "analysis_label": str(spec.get("label", run_id)),
        "protein_selection": selection, "membrane_selection": membrane_selection,
        "protein_atom_count": len(query), "protein_residue_count": len(residues),
        "membrane_atom_count": len(membrane),
        "membrane_residue_count": len(membrane.residues),
        "protein_atom_index_hash": _index_hash(query),
        "membrane_atom_index_hash": _index_hash(membrane),
        "membrane_resnames": ", ".join(sorted(set(map(str, membrane.resnames)))),
    }]))
    np.save(
        out / "membrane_atom_indices.npy",
        np.asarray(membrane.indices, dtype=np.int64),
    )
    tr = cfg["trajectory"]
    frame_slice = slice(tr.get("start"), tr.get("stop"), int(tr.get("stride", 1)))
    if precomputed_parts is None:
        nframes = len(universe.trajectory[frame_slice])
        part_paths = []
    else:
        part_paths = sorted(Path(precomputed_parts).glob("part-*.npz"))
        nframes = 0
        for path in part_paths:
            with np.load(path, allow_pickle=False) as data:
                nframes += len(data["time_ps"])
    if not nframes:
        raise ValueError(
            "residue_membrane_distance trajectory slice contains zero frames"
        )
    tmp = out / "distance_matrix_A.tmp.npy"
    matrix_path = out / "distance_matrix_A.npy"
    # The framexresidue matrix can be large, so a float32 memmap writes
    # rows directly to disk instead of holding the full matrix in RAM.
    matrix = np.lib.format.open_memmap(
        tmp,
        mode="w+",
        dtype=np.float32,
        shape=(nframes, len(residues))
    )
    times = np.empty(nframes, float)
    if precomputed_parts is None:
        for i, ts in enumerate(track(universe.trajectory[frame_slice],total=nframes,
                    desc=f"{rep.system}/{rep.label}: residue distances",
                    enabled=cfg.get("progress",{}).get("enabled",True))):
            times[i] = time_ps(ts, i, cfg)
            matrix[i] = _minimum_by_residue_hybrid(
                query, membrane, atom_to_residue, len(residues), _valid_box(universe),
                spec.get("initial_search_cutoff_A", 6.0),
                spec.get("search_growth_factor", 2.0),
                spec.get("query_block_atoms", 512),
                spec.get("fallback_block_atoms", 4096),
            ).astype(np.float32)
    else:
        # Rebuild the exact original memory-mapped output from checkpointed
        # per-frame distance chunks; never iterate the XTC a second time.
        offset = 0
        for path in part_paths:
            with np.load(path, allow_pickle=False) as data:
                values = data["values"]
                t = data["time_ps"]
                if values.ndim != 2 or values.shape != (len(t), len(residues)):
                    raise ValueError(f"Invalid fused residue distance chunk: {path}")
                matrix[offset:offset + len(t)] = values
                times[offset:offset + len(t)] = t
                offset += len(t)
        if offset != nframes:
            raise RuntimeError("Fused distance chunk count mismatch")
    matrix.flush()
    del matrix
    tmp.replace(matrix_path)
    np.save(out / "time_ps.npy", times)
    selected_counts = np.bincount(atom_to_residue, minlength=len(residues))
    metadata = pd.DataFrame({
        "analysis_instance": str(run_id),
        "analysis_label": str(spec.get("label", run_id)),
        "residue_order": np.arange(len(residues), dtype=int),
        "topology_residue_index": [int(res.ix) for res in residues],
        "segid": [str(getattr(res, "segid", "") or "?") for res in residues],
        "resid": [int(res.resid) for res in residues],
        "resname": [str(res.resname) for res in residues],
        "residue": [_residue_label(res) for res in residues],
        "n_selected_atoms": selected_counts,
    })
    atomic_parquet(out / "residues.parquet", metadata)
    distances = np.load(matrix_path, mmap_mode="r")
    q1, median, q3 = np.nanquantile(distances, [0.25, 0.5, 0.75], axis=0)
    profile = metadata.copy()
    profile["n_frames"] = np.sum(np.isfinite(distances), axis=0)
    profile["mean_min_distance_A"] = np.nanmean(distances, axis=0)
    profile["sd_min_distance_A"] = (
        np.nanstd(distances, axis=0, ddof=1)
        if nframes > 1
        else np.nan
    )
    profile["median_min_distance_A"] = median
    profile["q1_min_distance_A"], profile["q3_min_distance_A"] = q1, q3
    profile["iqr_min_distance_A"] = q3 - q1
    profile["min_distance_A"], profile["max_distance_A"] = np.nanmin(
        distances,
        axis=0
    ), np.nanmax(distances, axis=0)
    if spec.get("contact_cutoff_A") is not None:
        cutoff = float(spec["contact_cutoff_A"])
        profile["contact_cutoff_A"] = cutoff
        profile["fraction_frames_within_cutoff"] = np.nanmean(
            distances <= cutoff,
            axis=0,
        )
    atomic_parquet(out / "median_profile.parquet", profile)
    if spec.get("groups", False):
        group_rows = []
        for group_id, group in analysis_groups(cfg).items():
            group_residues = set(
                map(
                    int,
                    universe.select_atoms(group["selection"]).residues.ix,
                )
            )
            subset = profile[
                profile["topology_residue_index"].isin(group_residues)
            ].copy()
            if not subset.empty:
                subset.insert(0, "group_residue_index", np.arange(1, len(subset) + 1))
                subset.insert(0, "group", group_id)
                group_rows.append(subset)
        if group_rows:
            atomic_parquet(
                out / "group_median_profiles.parquet",
                pd.concat(group_rows, ignore_index=True)
            )
    return {
        "enabled": True,
        "run_id": run_id,
        "frames": nframes,
        "residues": len(residues),
        "selected_atoms": len(query),
        "membrane_atoms": len(membrane)
    }


def temporal_residue_membrane_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Expose frame x residue distance matrices as standard temporal
    profiles.
    """
    out = {}
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        directory = profile_dir(repout, "residue_membrane_distance", run.run_id)
        matrix_path, time_path, residue_path = (
            directory / "distance_matrix_A.npy", directory / "time_ps.npy",
            directory / "residues.parquet"
        )
        if not (matrix_path.exists() and time_path.exists() and residue_path.exists()):
            continue
        values = np.asarray(np.load(matrix_path, mmap_mode="r"), float)
        times = np.asarray(np.load(time_path), float)
        residues = pd.read_parquet(residue_path).sort_values("residue_order")
        positions = residues["residue_order"].to_numpy(float)
        labels = residues["residue"].astype(str).to_numpy()
        prefix = f"residue_membrane_distance__{run.run_id}"
        out[f"{prefix}__mean_min_distance_A"] = {
            "times": times,
            "values": values,
            "positions": positions,
            "labels": labels,
            "reducer": "mean",
        }
        out[f"{prefix}__median_min_distance_A"] = {
            "times": times,
            "values": values,
            "positions": positions,
            "labels": labels,
            "reducer": "median",
        }
        cutoff = run.spec.get("contact_cutoff_A")
        if cutoff is not None:
            out[f"{prefix}__fraction_frames_within_cutoff"] = {
                "times": times,
                "values": (values <= float(cutoff)).astype(float),
                "positions": positions,
                "labels": labels,
                "reducer": "mean",
            }
        if run.spec.get("groups", False):
            group_path = directory / "group_median_profiles.parquet"
            if group_path.exists():
                grouped = pd.read_parquet(group_path)
                index_by_topology = dict(zip(
                    residues.topology_residue_index.astype(int),
                    range(len(residues))
                ))
                for group_id, group in grouped.groupby("group", sort=False):
                    order = group.sort_values("group_residue_index")
                    columns = [
                        index_by_topology[
                            int(ix)
                        ]
                        for ix in order.topology_residue_index
                        if int(ix) in index_by_topology
                    ]
                    if not columns:
                        continue
                    gvalues = values[:, columns]
                    gpos = order["group_residue_index"].to_numpy(float)[:len(columns)]
                    glabel = order["residue"].astype(str).to_numpy()[:len(columns)]
                    gp = f"{prefix}__group_{group_id}"
                    out[f"{gp}__mean_min_distance_A"] = {
                        "times": times,
                        "values": gvalues,
                        "positions": gpos,
                        "labels": glabel,
                        "reducer": "mean",
                    }
                    out[f"{gp}__median_min_distance_A"] = {
                        "times": times,
                        "values": gvalues,
                        "positions": gpos,
                        "labels": glabel,
                        "reducer": "median",
                    }
                    if cutoff is not None:
                        out[f"{gp}__fraction_frames_within_cutoff"] = {
                            "times": times,
                            "values": (gvalues <= float(cutoff)).astype(float),
                            "positions": gpos,
                            "labels": glabel,
                            "reducer": "mean",
                        }
    return out
