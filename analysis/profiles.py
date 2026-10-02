"""
Frame-resolved structural profiles: DSSP, SASA, and residue-to-membrane
distances.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.dssp import DSSP
from MDAnalysis.lib.distances import capped_distance, distance_array

from .config import (
    DEFAULT_INSTANCE_ID,
    analysis_groups,
    analysis_runs,
    atomic_parquet,
    ensure_dir,
    profile_dir,
    time_ps,
)

def _residue_label(res) -> str:
    chain = getattr(res, "chainID", "") or getattr(res, "segid", "") or "?"
    return f"{chain}:{res.resname}{res.resid}"


def run_structure_profiles(rep, cfg: dict, repout: Path):
    if not cfg["analyses"]["dssp"].get(
        "enabled",
        False
    ) and not cfg["analyses"]["sasa"].get("enabled", False):
        return
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    if cfg["analyses"]["dssp"].get("enabled", False):
        _run_dssp(universe, cfg, repout)
    if cfg["analyses"]["sasa"].get("enabled", False):
        universe.trajectory[0]
        _run_sasa(universe, cfg, repout)


def _run_dssp(universe, cfg: dict, repout: Path):
    spec, tr = cfg["analyses"]["dssp"], cfg["trajectory"]
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    atoms = universe.select_atoms(selection)
    if not len(atoms):
        raise ValueError(f"dssp.selection matched zero atoms: {selection}")
    analysis = DSSP(
        atoms,
        guess_hydrogens=bool(spec.get("guess_hydrogens", True)),
        heavyatom_names=tuple(spec.get("heavyatom_names", ("N", "CA", "C", "O O1 OT1"))),
        hydrogen_name=str(spec.get("hydrogen_name", "H HN HT1 HT2 HT3")),
    )
    kwargs = {"start": tr.get("start"), "stop": tr.get("stop"), "step": int(tr["stride"])}
    backend = str(spec.get("backend", "serial"))
    if backend != "serial":
        kwargs["backend"] = backend
        if spec.get("n_workers") is not None:
            kwargs["n_workers"] = int(spec["n_workers"])
        if spec.get("n_parts") is not None:
            kwargs["n_parts"] = int(spec["n_parts"])
    analysis.run(**kwargs)

    residues = list(atoms.select_atoms("protein").residues)
    resids = np.asarray(analysis.results.resids)
    if len(residues) != len(resids):
        mapped, used = [], set()
        for resid in resids:
            candidates = [r for r in residues if r.resid == resid and r.ix not in used]
            if not candidates:
                raise RuntimeError("Could not map DSSP residues back to topology residues")
            mapped.append(candidates[0])
            used.add(candidates[0].ix)
        residues = mapped
    codes = np.asarray(analysis.results.dssp).astype(str)
    codes[codes == "-"] = "C"
    times = np.asarray(analysis.times, float)
    if cfg.get("time", {}).get("dt_ps") is not None:
        # Keep DSSP on the same user-defined analyzed-frame clock as the
        # streamed analyses.
        times = float(cfg["time"].get(
            "origin_ps",
            0.0
        )) + np.arange(len(times)) * float(cfg["time"]["dt_ps"])
    nframes, nres = codes.shape
    out = ensure_dir(repout / "profiles" / "dssp")
    atomic_parquet(out / "assignments.parquet", pd.DataFrame({
        "time_ps": np.repeat(times, nres),
        "topology_residue_index": np.tile([int(r.ix) for r in residues], nframes),
        "residue": np.tile([_residue_label(r) for r in residues], nframes),
        "dssp": codes.ravel(),
    }))


def _freesasa_structure(atoms):
    try:
        import freesasa
    except ImportError as exc:
        raise RuntimeError("SASA requires FreeSASA (`pip install freesasa`).") from exc
    structure = freesasa.Structure()
    for atom in atoms:
        x, y, z = map(float, atom.position)
        structure.addAtom(
            str(atom.name),
            str(atom.resname),
            int(atom.resid),
            str(atom.segid),
            x,
            y,
            z
        )
    return freesasa, structure


def _run_sasa(universe, cfg: dict, repout: Path):
    spec, tr = cfg["analyses"]["sasa"], cfg["trajectory"]
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    atoms = universe.select_atoms(selection)
    if not len(atoms):
        raise ValueError(f"sasa.selection matched zero atoms: {selection}")
    residues = list(atoms.residues)
    residue_lookup = {int(res.ix): i for i, res in enumerate(residues)}
    atom_residue = np.fromiter(
        (residue_lookup[int(ix)] for ix in atoms.resindices),
        dtype=int,
        count=len(atoms)
    )
    rows = []
    for frame_index, ts in enumerate(
        universe.trajectory[tr.get("start"):tr.get("stop"):int(tr["stride"])]
    ):
        # FreeSASA's structure stores coordinates, so it must be rebuilt
        # after each trajectory update.
        freesasa, structure = _freesasa_structure(atoms)
        params = freesasa.Parameters()
        algorithm = str(spec.get("algorithm", "shrake_rupley")).lower()
        if algorithm in {"shrake_rupley", "shrake-rupley", "sr"}:
            params.setAlgorithm(freesasa.ShrakeRupley)
            if spec.get("n_points") is not None:
                params.setNPoints(int(spec["n_points"]))
        elif algorithm in {"lee_richards", "lee-richards", "lr"}:
            params.setAlgorithm(freesasa.LeeRichards)
            if spec.get("n_slices") is not None:
                params.setNSlices(int(spec["n_slices"]))
        else:
            raise ValueError("sasa.algorithm must be shrake_rupley or lee_richards")
        if spec.get("probe_radius_A") is not None:
            params.setProbeRadius(float(spec["probe_radius_A"]))
        if spec.get("n_threads") is not None:
            params.setNThreads(int(spec["n_threads"]))
        result = freesasa.calc(structure, params)
        areas = np.fromiter(
            (result.atomArea(i) for i in range(len(atoms))),
            dtype=float,
            count=len(atoms)
        )
        residue_areas = np.bincount(atom_residue, weights=areas, minlength=len(residues)) / 100.0
        rows.extend(
            (time_ps(ts, frame_index, cfg), int(res.ix), _residue_label(res), float(area)) for res,
            area in zip(residues, residue_areas)
        )
    out = ensure_dir(repout / "profiles" / "sasa")
    atomic_parquet(
        out / "residue_sasa.parquet",
        pd.DataFrame(rows, columns=["time_ps", "topology_residue_index", "residue", "sasa_nm2"])
    )


def postprocess_structure_profiles(
    repout: Path,
    cfg: dict
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """
    Return additional whole-system and group time-series produced by
    profile analyses.
    """
    whole, grouped = [], []
    raw_group = repout / "raw" / "group_residues.parquet"
    group_map = pd.read_parquet(raw_group)[[
        "group",
        "topology_residue_index"
    ]].drop_duplicates() if raw_group.exists() else None

    dssp_file = repout / "profiles" / "dssp" / "assignments.parquet"
    if dssp_file.exists():
        data = pd.read_parquet(dssp_file)
        if not data.empty:
            fractions = data.groupby(["time_ps", "dssp"]).size().unstack(fill_value=0)
            fractions = fractions.div(
                fractions.sum(axis=1),
                axis=0
            ).reset_index().rename(columns={
                "H": "helix_fraction",
                "E": "sheet_fraction",
                "C": "coil_fraction"
            })
            for col in ("helix_fraction", "sheet_fraction", "coil_fraction"):
                if col not in fractions:
                    fractions[col] = 0.0
            fractions = fractions[["time_ps", "helix_fraction", "sheet_fraction", "coil_fraction"]]
            atomic_parquet(repout / "profiles" / "dssp" / "fraction_timeseries.parquet", fractions)
            whole.append(fractions)
            occupancy = data.groupby([
                "topology_residue_index",
                "residue",
                "dssp"
            ]).size().rename("frames").reset_index()
            occupancy["occupancy"] = occupancy["frames"] / occupancy.groupby([
                "topology_residue_index",
                "residue"
            ])["frames"].transform("sum")
            atomic_parquet(repout / "profiles" / "dssp" / "residue_occupancy.parquet", occupancy)
            if group_map is not None:
                gd = data.merge(group_map, on="topology_residue_index", how="inner")
                if not gd.empty:
                    gfrac = gd.groupby(["group", "time_ps", "dssp"]).size().unstack(fill_value=0)
                    gfrac = gfrac.div(
                        gfrac.sum(axis=1),
                        axis=0
                    ).reset_index().rename(columns={
                        "H": "helix_fraction",
                        "E": "sheet_fraction",
                        "C": "coil_fraction"
                    })
                    for col in ("helix_fraction", "sheet_fraction", "coil_fraction"):
                        if col not in gfrac:
                            gfrac[col] = 0.0
                    gfrac = gfrac[[
                        "group",
                        "time_ps",
                        "helix_fraction",
                        "sheet_fraction",
                        "coil_fraction"
                    ]]
                    atomic_parquet(
                        repout / "profiles" / "dssp" / "group_fraction_timeseries.parquet",
                        gfrac
                    )
                    grouped.append(gfrac)

    sasa_file = repout / "profiles" / "sasa" / "residue_sasa.parquet"
    if sasa_file.exists():
        data = pd.read_parquet(sasa_file)
        if not data.empty:
            total = data.groupby(
                "time_ps",
                as_index=False
            )["sasa_nm2"].sum().rename(columns={"sasa_nm2": "total_sasa_nm2"})
            atomic_parquet(repout / "profiles" / "sasa" / "total_timeseries.parquet", total)
            whole.append(total)
            summary = data.groupby(["topology_residue_index", "residue"])["sasa_nm2"].agg(
                mean="mean",
                sd="std",
                median="median",
                min="min",
                max="max"
            ).reset_index()
            atomic_parquet(repout / "profiles" / "sasa" / "residue_summary.parquet", summary)
            if group_map is not None:
                gd = data.merge(group_map, on="topology_residue_index", how="inner")
                if not gd.empty:
                    gtotal = gd.groupby(
                        ["group", "time_ps"],
                        as_index=False
                    )["sasa_nm2"].sum().rename(columns={"sasa_nm2": "total_sasa_nm2"})
                    atomic_parquet(
                        repout / "profiles" / "sasa" / "group_total_timeseries.parquet",
                        gtotal
                    )
                    grouped.append(gtotal)
    return whole, grouped


def temporal_structure_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Return per-frame residue matrices for SASA and DSSP occupancy
    profiles.
    """
    out = {}
    sasa_file = Path(repout) / "profiles" / "sasa" / "residue_sasa.parquet"
    if sasa_file.exists():
        data = pd.read_parquet(sasa_file)
        if not data.empty:
            residues = data[[
                "topology_residue_index",
                "residue"
            ]].drop_duplicates().sort_values("topology_residue_index")
            pivot = data.pivot(
                index="time_ps",
                columns="topology_residue_index",
                values="sasa_nm2"
            ).reindex(columns=residues.topology_residue_index)
            out["sasa_mean"] = {
                "times": pivot.index.to_numpy(float), "values": pivot.to_numpy(float),
                "positions": np.arange(
                    len(residues),
                    dtype=float
                ), "labels": residues.residue.to_numpy(str), "reducer": "mean",
            }

    dssp_file = Path(repout) / "profiles" / "dssp" / "assignments.parquet"
    if dssp_file.exists():
        data = pd.read_parquet(dssp_file)
        if not data.empty:
            residues = data[[
                "topology_residue_index",
                "residue"
            ]].drop_duplicates().sort_values("topology_residue_index")
            for code, label in (("H", "helix"), ("E", "sheet"), ("C", "coil")):
                indicator = data.assign(value=(data["dssp"] == code).astype(float))
                pivot = indicator.pivot(
                    index="time_ps",
                    columns="topology_residue_index",
                    values="value"
                ).reindex(columns=residues.topology_residue_index).fillna(0.0)
                out[f"dssp_{label}_occupancy"] = {
                    "times": pivot.index.to_numpy(float), "values": pivot.to_numpy(float),
                    "positions": np.arange(
                        len(residues),
                        dtype=float
                    ), "labels": residues.residue.to_numpy(str), "reducer": "mean",
                }
    return out


DEFAULT_LIPID_RESNAMES =[
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
    return hashlib.sha256(np.asarray(atoms.indices, dtype=np.int64).tobytes()).hexdigest()


def _minimum_by_residue_adaptive(
    query,
    membrane,
    atom_to_residue,
    n_residues,
    box,
    initial_cutoff_A,
    growth_factor,
    fallback_block_atoms
):
    # Start with capped neighbor searches because they avoid an
    # O(N_query*N_membrane) distance matrix in typical frames.
    # The cutoff grows only for residues still missing a neighbor; a
    # blocked full-distance fallback guarantees an exact minimum.
    minima = np.full(n_residues, np.inf)
    unresolved = np.ones(n_residues, bool)
    cutoff, growth = max(float(initial_cutoff_A), 0.1), max(float(growth_factor), 1.1)
    if box is not None:
        from MDAnalysis.lib.mdamath import triclinic_vectors
        tri = triclinic_vectors(np.asarray(box, float))
        alpha, beta, gamma = np.deg2rad(np.asarray(box, float)[3:6])
        ceiling = min(
            tri[1, 1] * np.sin(alpha),
            tri[2, 2] * np.sin(alpha),
            tri[0, 0] * np.sin(beta),
            tri[2, 2] * np.sin(beta),
            tri[0, 0] * np.sin(gamma),
            tri[1, 1] * np.sin(gamma)
        ) / 2.0
        ceiling = max(0.1, ceiling * 0.999999)
    else:
        xyz = np.vstack((query.positions, membrane.positions))
        ceiling = max(cutoff, float(np.linalg.norm(np.ptp(xyz, axis=0))) + cutoff)
    query_indices = np.arange(len(query), dtype=int)
    while np.any(unresolved):
        active = query_indices[unresolved[atom_to_residue]]
        if not len(active):
            break
        current = min(cutoff, ceiling)
        pairs, distances = capped_distance(
            query.positions[active],
            membrane.positions,
            max_cutoff=current,
            box=box,
            return_distances=True
        )
        if len(pairs):
            np.minimum.at(
                minima,
                atom_to_residue[active[pairs[:, 0].astype(int)]],
                np.asarray(distances, float)
            )
            unresolved[np.isfinite(minima)] = False
        if current >= ceiling:
            break
        cutoff = min(cutoff * growth, ceiling)
    if np.any(unresolved):
        active = query_indices[unresolved[atom_to_residue]]
        atom_min = np.full(len(active), np.inf)
        block = max(1, int(fallback_block_atoms))
        for start in range(0, len(membrane), block):
            d = distance_array(
                query.positions[active],
                membrane.positions[start:start + block],
                box=box
            )
            atom_min = np.minimum(atom_min, np.min(d, axis=1))
        np.minimum.at(minima, atom_to_residue[active], atom_min)
    return minima


def run_residue_membrane_distance(
    rep,
    cfg: dict,
    repout: Path,
    spec: dict | None = None,
    run_id: str = DEFAULT_INSTANCE_ID
) -> dict:
    spec = spec if spec is not None else cfg["analyses"]["residue_membrane_distance"]
    if not spec.get("enabled", False):
        return {"enabled": False}
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    selected = universe.select_atoms(selection)
    if not len(selected):
        raise ValueError(f"residue_membrane_distance.selection matched zero atoms: {selection}")
    query = selected.residues.atoms if spec.get("expand_residue_atoms", True) else selected
    query = _without_hydrogens(query, bool(spec.get("exclude_hydrogens", True)))
    membrane_selection = str(spec.get("membrane_selection") or ("resname " + " ".join(map(
        str,
        spec.get("lipid_resnames") or DEFAULT_LIPID_RESNAMES
    ))))
    membrane = _without_hydrogens(
        universe.select_atoms(membrane_selection),
        bool(spec.get("exclude_hydrogens", True))
    )
    if not len(query) or not len(membrane):
        raise ValueError(
            "residue_membrane_distance query and membrane selections must both contain atoms"
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
        "analysis_instance": str(run_id), "analysis_label": str(spec.get("label", run_id)),
        "protein_selection": selection, "membrane_selection": membrane_selection,
        "protein_atom_count": len(query), "protein_residue_count": len(residues),
        "membrane_atom_count": len(membrane), "membrane_residue_count": len(membrane.residues),
        "protein_atom_index_hash": _index_hash(query),
        "membrane_atom_index_hash": _index_hash(membrane),
        "membrane_resnames": ", ".join(sorted(set(map(str, membrane.resnames)))),
    }]))
    np.save(out / "membrane_atom_indices.npy", np.asarray(membrane.indices, dtype=np.int64))

    tr = cfg["trajectory"]
    frame_slice = slice(tr.get("start"), tr.get("stop"), int(tr.get("stride", 1)))
    nframes = len(universe.trajectory[frame_slice])
    if not nframes:
        raise ValueError("residue_membrane_distance trajectory slice contains zero frames")
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
    for i, ts in enumerate(universe.trajectory[frame_slice]):
        times[i] = time_ps(ts, i, cfg)
        matrix[i] = _minimum_by_residue_adaptive(
            query, membrane, atom_to_residue, len(residues), _valid_box(universe),
            spec.get(
                "initial_search_cutoff_A",
                6.0
            ), spec.get("search_growth_factor", 2.0), spec.get("fallback_block_atoms", 4096)
        ).astype(np.float32)
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
    profile["sd_min_distance_A"] = np.nanstd(distances, axis=0, ddof=1) if nframes > 1 else np.nan
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
        profile["fraction_frames_within_cutoff"] = np.nanmean(distances <= cutoff, axis=0)
    atomic_parquet(out / "median_profile.parquet", profile)

    if spec.get("groups", False):
        group_rows = []
        for group_id, group in analysis_groups(cfg).items():
            group_residues = set(map(int, universe.select_atoms(group["selection"]).residues.ix))
            subset = profile[profile["topology_residue_index"].isin(group_residues)].copy()
            if not subset.empty:
                subset.insert(0, "group_residue_index", np.arange(1, len(subset) + 1))
                subset.insert(0, "group", group_id)
                group_rows.append(subset)
        if group_rows:
            atomic_parquet(
                out / "group_median_profiles.parquet",
                pd.concat(group_rows, ignore_index=True)
            )
    return{
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
                    columns =[
                        index_by_topology[
                            int(ix)
                        ]for ix in order.topology_residue_index if int(ix) in index_by_topology
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
