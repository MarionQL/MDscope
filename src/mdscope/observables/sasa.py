"""FreeSASA residue solvent-accessible surface area."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from ..core.profiling import _residue_label
from ..pipeline.progress import track
from ..config import ensure_dir,atomic_parquet,time_ps
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


# Explicit van der Waals radii (Angstrom) for heterogeneous membrane atoms.
# This is a geometric approximation, not a force-field SASA parameterization.
ELEMENT_RADII_A = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
                   "P": 1.80, "S": 1.80, "F": 1.47, "CL": 1.75,
                   "BR": 1.85, "I": 1.98}


def sasa_instances(universe, cfg):
    """Build independent target/occluder SASA selections from one Universe."""
    from ..config import analysis_runs
    try:
        import freesasa
    except ImportError as exc:
        raise RuntimeError("SASA requires freesasa") from exc
    instances = {}
    for run in analysis_runs(cfg, "sasa"):
        spec = run.spec
        target_sel = str(spec.get("target_selection") or spec.get("selection")
                         or cfg["selections"]["protein"])
        occ_sel = str(spec.get("occluder_selection") or target_sel)
        target = universe.select_atoms(target_sel)
        occluders = universe.select_atoms(occ_sel)
        if not len(target) or not len(occluders):
            raise ValueError(f"SASA run {run.run_id}: target/occluder selection is empty")
        target_idx = np.flatnonzero(np.isin(occluders.indices, target.indices))
        if len(target_idx) != len(target):
            raise ValueError(f"SASA run {run.run_id}: occluders must include every target atom")
        atom_selection = spec.get("atom_sasa_selection")
        requested = (universe.select_atoms(str(atom_selection)).indices
                     if atom_selection else np.empty(0, dtype=int))
        atom_positions = np.flatnonzero(np.isin(target.indices, requested))
        if atom_selection and (not len(atom_positions) or
                               len(atom_positions) != len(requested)):
            raise ValueError("atom_sasa_selection must be a nonempty subset of SASA target atoms")
        motif_sequences = spec.get("atom_sasa_sequences") or []
        if motif_sequences:
            from ..mechanisms.mechanistic_common import sequence_matches
            matched = sequence_matches(
                target.select_atoms("name CA"),
                {f"atom_sasa_{i}": sequence for i, sequence in
                 enumerate(motif_sequences)},
            )
            expected = spec.get("atom_sasa_expected_matches")
            if not matched or (expected is not None and len(matched) != int(expected)):
                raise ValueError(
                    f"atom_sasa_sequences mapped {len(matched)} occurrences; "
                    f"expected {expected}"
                )
            motif_resindices = np.unique(np.concatenate(list(matched.values())))
            atom_positions = atom_positions[np.isin(
                target.resindices[atom_positions], motif_resindices
            )]
            if not len(atom_positions):
                raise ValueError("No selected SASA sidechain atoms in mapped motifs")
        residues = list(target.residues)
        residue_ids = np.asarray([r.ix for r in residues], dtype=int)
        mapping = np.searchsorted(residue_ids, target.resindices)
        # Indices from MDAnalysis residues are sorted by topology index.
        if not np.array_equal(residue_ids[mapping], target.resindices):
            raise ValueError("SASA target residue mapping failed")
        radii_spec = spec.get("element_radii_A", {})
        radii_model = {**ELEMENT_RADII_A, **{str(k).upper(): float(v)
                                         for k, v in radii_spec.items()}}
        elements = getattr(occluders, "elements", None)
        if elements is None:
            # GROMACS/PSF topologies frequently omit an element field. Infer
            # only common biopolymer/lipid elements from CHARMM-style names;
            # atypical names still require explicit elements in the topology.
            def _symbol(name):
                name = str(name).lstrip("0123456789").upper()
                if name.startswith(("CL", "BR")):
                    return name[:2]
                if name and name[0] in "HCNOPSFI":
                    return name[0]
                raise ValueError(
                    f"Cannot infer SASA element for atom name {name!r}; "
                    "provide topology elements"
                )
            symbols = [_symbol(name) for name in occluders.names]
        else:
            symbols = [str(e).upper() for e in elements]
        missing = sorted(set(symbols) - radii_model.keys())
        if missing:
            raise ValueError(f"SASA run {run.run_id}: no radii for elements {missing}")
        radii = [float(radii_model[e]) for e in symbols]
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
        instances[run.run_id] = {
            "occluders": occluders, "radii": radii, "target_idx": target_idx,
            "atom_positions": atom_positions,
            "selected_atom_indices": target.indices[atom_positions],
            "selected_atom_resindices": target.resindices[atom_positions],
            "selected_atom_names": target.names[atom_positions],
            "residues": residues, "mapping": mapping, "params": params,
            "library": freesasa, "periodic": bool(spec.get("periodic", False)),
            "reference_run": spec.get("reference_run"),
            "target_atoms": np.asarray(target.indices, dtype=int),
            "target_radii": np.asarray(radii, dtype=float)[target_idx],
            "probe_radius": float(spec.get("probe_radius_A", 1.4)),
            "settings": (algorithm, spec.get("n_points"), spec.get("n_slices"),
                         spec.get("probe_radius_A", 1.4)),
        }
    for run_id, state in instances.items():
        reference = state["reference_run"]
        if reference is None:
            continue
        if reference == run_id or reference not in instances:
            raise ValueError(f"SASA run {run_id}: unknown/self reference_run {reference!r}")
        baseline = instances[reference]
        if (not np.array_equal(state["target_atoms"], baseline["target_atoms"])
                or not np.array_equal(state["target_radii"], baseline["target_radii"])
                or state["settings"] != baseline["settings"]):
            raise ValueError(
                f"SASA run {run_id}: reference_run must have identical target atoms, "
                "target radii, sampling algorithm, and probe radius"
            )
    return instances


def write_sasa_burial(repout, run_id, reference_run):
    """Join identical residue/time records to quantify environmental occlusion."""
    root = Path(repout) / "profiles" / "sasa"
    environment = pd.read_parquet(root / run_id / "residue_sasa.parquet")
    ref_dir = root if reference_run == "default" else root / reference_run
    baseline = pd.read_parquet(ref_dir / "residue_sasa.parquet")
    keys = ["time_ps", "topology_residue_index"]
    paired = baseline.merge(environment, on=keys, how="inner",
                            suffixes=("_isolated", "_environment"),
                            validate="one_to_one")
    if len(paired) != len(baseline) or len(paired) != len(environment):
        raise ValueError(f"SASA run {run_id}: reference and environment frames/residues differ")
    if not (paired["residue_isolated"] == paired["residue_environment"]).all():
        raise ValueError(f"SASA run {run_id}: residue identity mismatch")
    paired = paired.rename(columns={
        "residue_isolated": "residue",
        "sasa_nm2_isolated": "isolated_sasa_nm2",
        "sasa_nm2_environment": "environment_sasa_nm2",
    }).drop(columns=["residue_environment"])
    paired["buried_sasa_nm2"] = (
        paired["isolated_sasa_nm2"] - paired["environment_sasa_nm2"]
    )
    isolated = paired["isolated_sasa_nm2"].to_numpy(float)
    paired["buried_fraction"] = np.divide(
        paired["buried_sasa_nm2"].to_numpy(float),
        isolated, out=np.full(len(isolated), np.nan), where=isolated > 0,
    )
    target = root / run_id
    atomic_parquet(target / "residue_burial.parquet", paired)
    summary = paired.groupby(["topology_residue_index", "residue"], as_index=False).agg(
        mean_buried_sasa_nm2=("buried_sasa_nm2", "mean"),
        std_buried_sasa_nm2=("buried_sasa_nm2", "std"),
        mean_buried_fraction=("buried_fraction", "mean"),
        mean_environment_sasa_nm2=("environment_sasa_nm2", "mean"),
        mean_isolated_sasa_nm2=("isolated_sasa_nm2", "mean"),
    )
    atomic_parquet(target / "residue_burial_summary.parquet", summary)


def sasa_residue_rows(state, time, atom_sink=None):
    """Return per-residue protein SASA in nm², including occlusion by other atoms."""
    atoms = state["occluders"]
    coords = atoms.positions.astype(np.float64, copy=False)
    if state["periodic"]:
        from MDAnalysis.lib.distances import minimize_vectors
        box = np.asarray(atoms.universe.dimensions, dtype=float)
        if box.shape != (6,) or np.any(~np.isfinite(box)) or np.any(box[:3] <= 0):
            raise ValueError("Periodic SASA requires valid box dimensions")
        target_coords = coords[state["target_idx"]]
        if np.any(np.ptp(target_coords, axis=0) >= box[:3] / 2):
            raise ValueError(
                "Periodic SASA requires a whole compact protein; unwrap the "
                "protein before SASA if it crosses the periodic boundary"
            )
        anchor = target_coords.mean(axis=0)
        coords = anchor + minimize_vectors(
            np.asarray(coords - anchor, dtype=np.float64), box
        )
    target_idx = state["target_idx"]
    # A distant atom cannot occlude a target SASA sphere. Restrict large
    # environments using a conservative bound on both van der Waals radii
    # plus the solvent probe; this preserves target areas exactly.
    if len(coords) > len(target_idx):
        from MDAnalysis.lib.distances import capped_distance
        max_cutoff = (
            float(np.max(state["target_radii"])) +
            float(np.max(state["radii"])) +
            2.0 * state["probe_radius"] + 1e-5
        )
        pairs = capped_distance(coords[target_idx], coords, max_cutoff,
                                method="nsgrid", return_distances=False)
        candidate_idx = np.union1d(target_idx, pairs[:, 1])
    else:
        candidate_idx = np.arange(len(coords))
    candidate_xyz = coords[candidate_idx]
    candidate_radii = np.asarray(state["radii"])[candidate_idx]
    # FreeSASA 2.2.1 Shrake-Rupley can crash in native code when an atom
    # has no overlapping probe-inflated neighbor. Its area is exactly a
    # sphere: remove such atoms from the native input and calculate it
    # analytically. They cannot occlude any remaining atom by definition.
    from MDAnalysis.lib.distances import capped_distance
    inflated = candidate_radii + state["probe_radius"]
    pairs, distances = capped_distance(
        candidate_xyz, candidate_xyz, 2*float(inflated.max()),
        return_distances=True,
    )
    overlap = ((pairs[:, 0] != pairs[:, 1])
               & (distances < inflated[pairs[:, 0]] + inflated[pairs[:, 1]]))
    connected = np.unique(pairs[overlap].ravel())
    candidate_areas = 4*np.pi*inflated**2
    if len(connected):
        result = state["library"].calcCoord(
            candidate_xyz[connected].ravel().tolist(),
            candidate_radii[connected].tolist(), state["params"],
        )
        candidate_areas[connected] = np.fromiter(
            (result.atomArea(int(i)) for i in range(len(connected))),
            dtype=float, count=len(connected),
        )
    candidate_target_idx = np.searchsorted(candidate_idx, target_idx)
    areas = candidate_areas[candidate_target_idx]
    if atom_sink is not None:
        for j, target_position in enumerate(state["atom_positions"]):
            atom_sink.append((
                time, int(state["selected_atom_indices"][j]),
                int(state["selected_atom_resindices"][j]),
                str(state["selected_atom_names"][j]),
                float(areas[target_position] / 100.),
            ))
    values = np.bincount(
        state["mapping"], weights=areas, minlength=len(state["residues"])
    ) / 100.0
    return [(time, int(res.ix), _residue_label(res), float(area))
            for res, area in zip(state["residues"], values)]


def _run_sasa(universe, cfg: dict, repout: Path):
    """Standalone SASA path using the same implementation as the fused observer."""
    states = sasa_instances(universe, cfg)
    tr = cfg["trajectory"]
    rows = {key: [] for key in states}
    atom_rows = {key: [] for key in states}
    frames = universe.trajectory[tr.get("start"):tr.get("stop"):int(tr["stride"])]
    for frame_index, ts in enumerate(track(
            frames, total=len(frames), desc="SASA",
            enabled=cfg.get("progress", {}).get("enabled", True))):
        t = time_ps(ts, frame_index, cfg)
        for key, state in states.items():
            rows[key].extend(sasa_residue_rows(state, t, atom_rows[key]))
    for key, values in rows.items():
        output_dir = (Path(repout) / "profiles" / "sasa" /
                      key if key != "default" else Path(repout) / "profiles" / "sasa")
        atomic_parquet(ensure_dir(output_dir) / "residue_sasa.parquet",
                       pd.DataFrame(values, columns=[
                           "time_ps", "topology_residue_index", "residue", "sasa_nm2"
                       ]))
        if atom_rows[key]:
            atomic_parquet(ensure_dir(output_dir) / "atom_sasa.parquet",
                           pd.DataFrame(atom_rows[key], columns=[
                               "time_ps", "atom_index", "topology_residue_index",
                               "atom_name", "sasa_nm2"
                           ]))
    for key, state in states.items():
        if state["reference_run"] is not None:
            write_sasa_burial(repout, key, state["reference_run"])

