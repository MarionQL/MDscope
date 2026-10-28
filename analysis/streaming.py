"""Single-pass trajectory coordinator: analysis-specific metrics live in separate modules.

Compatibility: raw/summary output schema and checkpoint format retained.
"""
from __future__ import annotations
from pathlib import Path
from time import perf_counter
import MDAnalysis as mda
import numpy as np
import pandas as pd
from .config import analysis_groups, analysis_runs, atomic_parquet,ensure_dir,time_ps
from .progress import track
from .profiling import PhaseProfiler
from .core.trajectory_sources import open_msd_universe, paired_frames
from .msd import initialize_msd
from .fused_membrane import FusedMembrane, enabled_names as fused_enabled_names
from .fused_extensions import FusedExtensions, enabled_names as extension_enabled_names
from .fused_distance import FusedDistances, enabled_names as distance_enabled_names
from .sdc4_mechanisms import SDC4Mechanisms, enabled as sdc4_enabled
from .alignment import _apply_transform,_fit_transform,_valid_box
from .rmsd import compute_rmsd
from .rmsf import Welford
from .radius_of_gyration import compute_rg
from .density_1d import DensityProfileAccumulator
from .membrane import initialize_membrane,membrane_metrics
from .interactions import chain_label,initialize_ligand_contacts,initialize_salt_bridges,ligand_contact_events,residue_label,salt_bridge_events
def _checkpoint_payload(path, cfg_hash, processed, next_part, rmsf, groups, densities, membrane):
    payload = {
        "config_hash": np.asarray(cfg_hash), "processed": np.asarray(processed),
        "next_part": np.asarray(next_part),
        "w_n": np.asarray(rmsf.n), "w_mean": rmsf.mean, "w_M2": rmsf.M2,
    }
    for group_id, group in groups.items():
        safe = group_id.replace("/", "_")
        payload[f"g_{safe}_n"] = np.asarray(group["rmsf"].n)
        payload[f"g_{safe}_mean"] = group["rmsf"].mean
        payload[f"g_{safe}_M2"] = group["rmsf"].M2
    for density in densities:
        payload.update(density.checkpoint_payload())
    if membrane.enabled and membrane.msd:
        payload.update({
            "membrane_msd_last_positions": membrane.msd_last_positions,
            "membrane_msd_displacements": membrane.msd_displacements,
            "membrane_msd_resindices": membrane.msd_resindices,
        })
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **payload)
    tmp.replace(path)


def _restore_checkpoint(checkpoint, cfg_hash, rmsf, groups, densities, membrane):
    cp = np.load(checkpoint, allow_pickle=False)
    if str(cp["config_hash"].item()) != cfg_hash or cp["w_mean"].shape != rmsf.mean.shape:
        raise ValueError("checkpoint does not match current run")
    rmsf.n, rmsf.mean, rmsf.M2 = (
        int(cp["w_n"]), cp["w_mean"].astype(float), cp["w_M2"].astype(float)
    )
    for group_id, group in groups.items():
        safe = group_id.replace("/", "_")
        group["rmsf"].n = int(cp[f"g_{safe}_n"])
        group["rmsf"].mean = cp[f"g_{safe}_mean"].astype(float)
        group["rmsf"].M2 = cp[f"g_{safe}_M2"].astype(float)
    for density in densities:
        density.restore_checkpoint(cp)
    if membrane.enabled and membrane.msd:
        if not np.array_equal(
            np.asarray(cp["membrane_msd_resindices"], int),
            membrane.msd_resindices
        ):
            raise ValueError("membrane MSD checkpoint residue selection changed")
        membrane.msd_last_positions = np.asarray(cp["membrane_msd_last_positions"], float)
        membrane.msd_displacements = np.asarray(cp["membrane_msd_displacements"], float)
    return int(cp["processed"]), int(cp["next_part"])


def run_streaming_pass(rep, cfg: dict, outroot: Path, cfg_hash: str):
    # Protein, membrane, density, salt-bridge, and ligand-contact
    # calculations share one trajectory pass.
    # Reading a large XTC once is substantially cheaper than reopening
    # it independently for each simple metric.
    repout = ensure_dir(outroot / rep.system / rep.label)
    raw, cache = ensure_dir(repout / "raw"), ensure_dir(repout / "cache")
    checkpoint = cache / "stream_checkpoint.npz"
    profiler = PhaseProfiler.from_environment()
    setup_started = perf_counter() if profiler.enabled else None
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    selections = cfg["selections"]
    protein = universe.select_atoms(selections["protein"])
    if not len(protein):
        raise ValueError(f"Protein selection matched zero atoms: {selections['protein']}")
    alignment_spec = cfg["analyses"].get("alignment", {})
    alignment_enabled = bool(alignment_spec.get("enabled", True))
    alignment_selection = (
        str(
            alignment_spec.get("selection") or selections.get("alignment") or "protein and name CA"
        )
    )
    features = universe.select_atoms(selections.get("features") or alignment_selection)
    if not len(features):
        raise ValueError("Feature selection matched zero atoms")

    # Frame 0 is the structural reference for all alignment/RMSD calculations.
    universe.trajectory[0]
    alignment_atoms = (
        universe.select_atoms(alignment_selection)
        if alignment_enabled
        else universe.atoms[[]]
    )
    if alignment_enabled and not len(alignment_atoms):
        raise ValueError(
            f"Alignment selection matched zero atoms: {alignment_selection}"
        )
    reference_alignment = (
        alignment_atoms.positions.astype(np.float64, copy=True)
        if alignment_enabled
        else None
    )
    reference_features = features.positions.astype(np.float64, copy=True)

    ligand = (
        universe.select_atoms(selections["ligand"])
        if selections.get("ligand")
        else None
    )
    reference_ligand = (
        ligand.positions.astype(np.float64, copy=True)
        if ligand is not None and len(ligand)
        else None
    )
    salt_state = initialize_salt_bridges(universe, cfg)
    ligand_state = initialize_ligand_contacts(universe, protein, ligand, cfg)
    atomic_parquet(
        raw / "protein_chains.parquet",
        pd.DataFrame({"chain": sorted({chain_label(atom) for atom in protein})})
    )

    rg_spec = cfg["analyses"].get("radius_of_gyration", {})
    rg_atoms = universe.select_atoms(str(rg_spec.get("selection") or selections["protein"]))
    if rg_spec.get("enabled", True) and not len(rg_atoms):
        raise ValueError("radius_of_gyration.selection matched zero atoms")

    densities = [DensityProfileAccumulator(
        universe,
        cfg,
        run.spec,
        run.run_id
    ) for run in analysis_runs(cfg, "density_1d")]
    tr = cfg["trajectory"]
    start, stop, stride, chunk = (
        tr.get("start"), tr.get("stop"), int(tr["stride"]), int(tr["chunk_frames"])
    )
    rmsf_acc = Welford((len(features), 3))

    groups = {}
    # The universe is still on frame 0, so all copied coordinates below
    # are reference coordinates.
    for group_id, group_spec in analysis_groups(cfg).items():
        atoms = universe.select_atoms(group_spec["selection"])
        measured = features.intersection(atoms)
        if not len(atoms) or not len(measured):
            raise ValueError(f"Analysis group {group_id!r} resolved to zero atoms/features")
        align_selection = group_spec.get("alignment_selection") if alignment_enabled else None
        if align_selection:
            group_alignment = universe.select_atoms(str(align_selection))
            if not len(group_alignment):
                raise ValueError(
                    f"Group {group_id!r} alignment selection matched zero atoms: "
                    f"{align_selection}"
                )
            group_alignment_reference = group_alignment.positions.astype(
                np.float64,
                copy=True,
            )
        else:
            group_alignment = None
            group_alignment_reference = None
        groups[group_id] = {
            "atoms": atoms,
            "measure": measured,
            "reference": measured.positions.astype(np.float64, copy=True),
            "alignment": group_alignment,
            "alignment_reference": group_alignment_reference,
            "rmsf": Welford((len(measured), 3)),
            "rg": rg_atoms.intersection(atoms),
        }

    if start is not None:
        universe.trajectory[int(start)]
    membrane = initialize_membrane(universe, cfg)
    msd_universe = open_msd_universe(
        rep, cfg, membrane, universe, 0 if start is None else int(start)
    )
    if msd_universe is not None:
        print(f"MSD original coordinate source: {rep.msd_trajectory}", flush=True)
    fused = FusedMembrane(universe, cfg, repout) if fused_enabled_names(cfg) else None
    extensions = FusedExtensions(universe, cfg, repout) if extension_enabled_names(cfg) else None
    distances = FusedDistances(universe, cfg, repout) if distance_enabled_names(cfg) else None
    sdc4 = SDC4Mechanisms(universe, cfg, repout) if sdc4_enabled(cfg) else None

    processed = part = 0
    if checkpoint.exists():
        try:
            processed, part = _restore_checkpoint(
                checkpoint,
                cfg_hash,
                rmsf_acc,
                groups,
                densities,
                membrane
            )
        except Exception:
            for path in raw.glob("*-part-*.parquet"):
                path.unlink(missing_ok=True)
            for path in cache.glob("aligned-part-*.npz"):
                path.unlink(missing_ok=True)
            checkpoint.unlink(missing_ok=True)
            rmsf_acc = Welford((len(features), 3))
            for group in groups.values():
                group["rmsf"] = Welford((len(group["measure"]), 3))
            densities = [DensityProfileAccumulator(
                universe,
                cfg,
                run.spec,
                run.run_id
            ) for run in analysis_runs(cfg, "density_1d")]
            for density in densities:
                density.clear_frame_parts(repout)
            if fused is not None:
                fused.prune_from(0)
            if extensions is not None:
                extensions.prune_from(0)
            if distances is not None:
                distances.prune_from(0)
            if sdc4 is not None:
                sdc4.prune_from(0)
            if start is not None:
                universe.trajectory[int(start)]
            membrane = initialize_membrane(universe, cfg)
            if msd_universe is not None:
                msd_universe.trajectory[0 if start is None else int(start)]
                initialize_msd(msd_universe, membrane, cfg["analyses"]["membrane"])
            processed = part = 0

    if extensions is not None:
        extensions.prune_from(part)
    if distances is not None:
        distances.prune_from(part)
    if sdc4 is not None:
        sdc4.prune_from(part)
    if fused is not None:
        # Keep only chunks committed by the streaming checkpoint. An interrupted
        # chunk may have been written just before a crash and must be replayed.
        fused.prune_from(part)
    first = (0 if start is None else int(start)) + processed * stride
    rows, group_rows, salt_rows, ligand_rows, coords, coord_times = [], [], [], [], [], []
    frames=universe.trajectory[slice(first,stop,stride)]
    raw_msd_frames = (
        msd_universe.trajectory[slice(first,stop,stride)]
        if msd_universe is not None else None
    )
    input_frames = paired_frames(frames, raw_msd_frames)
    if setup_started is not None:
        profiler.add('initialization', perf_counter() - setup_started)
    for ts in track(profiler.frames(input_frames) if profiler.enabled else input_frames,
                    total=len(frames),
                    desc=f"{rep.system}/{rep.label}: streaming",
                    enabled=cfg.get("progress",{}).get("enabled",True)):
        if profiler.enabled:
            frame_started = perf_counter()
            measured_before = sum(profiler.seconds.values())
        tps = time_ps(ts, processed, cfg)
        # Membrane observers must consume raw wrapped coordinates, not aligned
        # protein coordinates. All three analyses share this exact frame.
        if fused is not None:
            profiler.call("fused_membrane_update", fused.update, ts, processed)
        if extensions is not None:
            profiler.call("fused_extensions_update", extensions.update, ts, processed)
        if distances is not None:
            profiler.call("fused_distances_update", distances.update, ts, processed)
        if sdc4 is not None:
            profiler.call("sdc4_mechanisms_update", sdc4.update, ts, processed)
        for density in densities:
            profiler.call("density_update", density.update, universe, tps)
        membrane_values = profiler.call(
            "membrane_metrics", membrane_metrics, universe, membrane, msd_universe
        )
        if alignment_enabled:
            global_transform = profiler.call(
                "protein_alignment_fit", _fit_transform,
                alignment_atoms.positions, reference_alignment,
            )
            aligned_features = profiler.call(
                "feature_alignment_apply", _apply_transform,
                features.positions, global_transform,
            )
        else:
            global_transform = None
            aligned_features = features.positions.astype(
                np.float64,
                copy=True,
            )

        rmsd_value = (
            float(
                compute_rmsd(
                    aligned_features,
                    reference_features
                )
            ) if cfg["analyses"]["rmsd"].get("enabled", True) else np.nan
        )
        rg_value = compute_rg(rg_atoms) if rg_spec.get("enabled", True) else np.nan
        if reference_ligand is not None:
            if alignment_enabled:
                aligned_ligand = _apply_transform(
                    ligand.positions,
                    global_transform,
                )
            else:
                aligned_ligand = ligand.positions.astype(
                    np.float64,
                    copy=True,
                )
            ligand_rmsd = float(compute_rmsd(
                aligned_ligand,
                reference_ligand
            ))
        else:
            ligand_rmsd = np.nan
        rmsf_acc.update(aligned_features)

        for group_id, group in groups.items():
            if group["alignment"] is not None:
                group_transform = _fit_transform(
                    group["alignment"].positions,
                    group["alignment_reference"],
                )
                positions = _apply_transform(
                    group["measure"].positions,
                    group_transform,
                )
            elif alignment_enabled:
                positions = _apply_transform(
                    group["measure"].positions,
                    global_transform,
                )
            else:
                positions = group["measure"].positions.astype(
                    np.float64,
                    copy=True,
                )
            group["rmsf"].update(positions)
            group_rows.append({
                "time_ps": tps,
                "group": group_id,
                "rmsd_A": (
                    float(
                        compute_rmsd(
                            positions,
                            group["reference"]
                        )
                    )
                    if cfg["analyses"]["rmsd"].get("enabled", True)
                    else np.nan
                ),
                "radius_of_gyration_A": (
                    compute_rg(group["rg"])
                    if rg_spec.get("enabled", True) and len(group["rg"])
                    else np.nan
                ),
            })

        salt_rows.extend(profiler.call(
            "salt_bridges", salt_bridge_events, salt_state, tps, _valid_box(universe)
        ))
        ligand_min, ligand_count, new_contacts = profiler.call(
            "ligand_contacts", ligand_contact_events, ligand_state, tps,
            _valid_box(universe),
        )
        ligand_rows.extend(new_contacts)
        row = {"time_ps": tps, "rmsd_A": rmsd_value, "radius_of_gyration_A": rg_value}
        if reference_ligand is not None:
            row["ligand_rmsd_A"] = ligand_rmsd
        if cfg["analyses"]["ligand_contacts"].get("enabled", False):
            row["ligand_min_distance_A"] = ligand_min
            row["ligand_contact_residue_count"] = ligand_count
        row.update(membrane_values)
        rows.append(row)
        if cfg["cache"].get("aligned_coordinates", True):
            coords.append(aligned_features.astype(
                cfg["cache"].get("dtype", "float32"),
                copy=True
            ))
            coord_times.append(tps)
        processed += 1
        if len(rows) >= chunk:
            _flush(
                raw,
                cache,
                rows,
                group_rows,
                salt_rows,
                ligand_rows,
                coords,
                coord_times,
                part,
                cfg
            )
            for density in densities:
                density.flush_frames(repout, part)
            if fused is not None:
                fused.flush(part)
            if extensions is not None:
                extensions.flush(part)
            if distances is not None:
                distances.flush(part)
            if sdc4 is not None:
                sdc4.flush(part)
            rows, group_rows, salt_rows, ligand_rows, coords, coord_times = [], [], [], [], [], []
            part += 1
            _checkpoint_payload(
                checkpoint,
                cfg_hash,
                processed,
                part,
                rmsf_acc,
                groups,
                densities,
                membrane
            )

        if profiler.enabled:
            measured_during = sum(profiler.seconds.values()) - measured_before
            profiler.add('other_frame_work', max(0., perf_counter() - frame_started - measured_during))

    final_started = perf_counter() if profiler.enabled else None
    if rows or group_rows or salt_rows or ligand_rows or coords:
        _flush(
            raw,
            cache,
            rows,
            group_rows,
            salt_rows,
            ligand_rows,
            coords,
            coord_times,
            part,
            cfg
        )
        for density in densities:
            density.flush_frames(repout, part)
        if fused is not None:
            fused.flush(part)
        if extensions is not None:
            extensions.flush(part)
        if distances is not None:
            distances.flush(part)
        if sdc4 is not None:
            sdc4.flush(part)
        part += 1
        _checkpoint_payload(
            checkpoint,
            cfg_hash,
            processed,
            part,
            rmsf_acc,
            groups,
            densities,
            membrane
        )

    fused_results = fused.finalize() if fused is not None else {}
    extension_results = extensions.finalize() if extensions is not None else {}
    distance_results = distances.finalize(rep) if distances is not None else {}
    sdc4_results = sdc4.finalize() if sdc4 is not None else {}

    if cfg["cache"].get("aligned_coordinates", True):
        atomic_parquet(
            cache / "feature_atoms.parquet",
            pd.DataFrame({
                "atom_index": features.indices,
                "topology_residue_index": features.resindices,
                "residue": [residue_label(atom) for atom in features],
                "atom_name": features.names,
            })
        )

    rmsf_values = rmsf_acc.rmsf()
    rmsf_enabled = cfg["analyses"]["rmsf"].get("enabled", True)
    if rmsf_enabled:
        atomic_parquet(
            raw / "rmsf.parquet",
            pd.DataFrame({
                "atom_index": features.indices,
                "topology_residue_index": features.resindices,
                "residue": [residue_label(atom) for atom in features],
                "atom_name": features.names,
                "rmsf_A": rmsf_values,
                "equivalent_bfactor_A2": (
                    (8.0 * np.pi ** 2 / 3.0) * rmsf_values ** 2
                ),
            }),
        )
    group_rmsf = []
    group_residues = []
    for group_id, group in groups.items():
        if rmsf_enabled:
            group_rmsf.extend(
                (
                    group_id,
                    int(atom.index),
                    int(atom.resindex),
                    residue_label(atom),
                    atom.name,
                    float(value),
                )
                for atom, value in zip(group["measure"], group["rmsf"].rmsf())
            )
        group_residues.extend(
            (group_id, int(res.ix), str(res.resname), int(res.resid))
            for res in group["atoms"].residues
        )
    if group_rmsf:
        atomic_parquet(
            raw / "group_rmsf.parquet",
            pd.DataFrame(
                group_rmsf,
                columns=[
                    "group",
                    "atom_index",
                    "topology_residue_index",
                    "residue",
                    "atom_name",
                    "rmsf_A",
                ],
            ),
        )
    if group_residues:
        atomic_parquet(
            raw / "group_residues.parquet",
            pd.DataFrame(
                group_residues,
                columns=[
                    "group",
                    "topology_residue_index",
                    "resname",
                    "resid",
                ],
            ).drop_duplicates(),
        )
    for density in densities:
        density.finalize(repout)

    if rmsf_enabled:
        try:
            bfactor = (8.0 * np.pi ** 2 / 3.0) * rmsf_values ** 2
            if not hasattr(universe.atoms, "tempfactors"):
                universe.add_TopologyAttr(
                    "tempfactors",
                    np.zeros(len(universe.atoms)),
                )
            universe.atoms.tempfactors = 0.0
            features.tempfactors = bfactor
            protein.write(str(raw / "rmsf_bfactor.pdb"))
        except Exception as exc:
            print(f"Warning: could not write RMSF B-factor PDB: {exc}")
    if final_started is not None:
        profiler.add('finalization', perf_counter() - final_started)
    profiler.write(repout / "streaming_phase_profile.json", frames=processed)
    return {
        "frames": processed,
        "topology_atoms": len(universe.atoms),
        "protein_atoms": len(protein),
        "feature_atoms": len(features),
        "alignment_atoms": len(alignment_atoms),
        "fused_membrane": fused_results,
        "fused_extensions": extension_results,
        "fused_distances": distance_results,
        "sdc4_mechanisms": sdc4_results,
    }


def _flush(raw, cache, rows, group_rows, salt_rows, ligand_rows, coords, times, part, cfg):
    # Chunked writes cap RAM use and make completed trajectory segments
    # recoverable after wall-time/job failures.
    tag = f"part-{part:06d}"
    if rows:
        atomic_parquet(raw / f"timeseries-{tag}.parquet", pd.DataFrame(rows))
    if group_rows:
        atomic_parquet(raw / f"group_timeseries-{tag}.parquet", pd.DataFrame(group_rows))
    if salt_rows:
        columns = [
            "time_ps",
            "residue_1",
            "residue_2",
            "chain_1",
            "chain_2",
            "scope",
            "interaction_class",
            "min_distance_A",
        ]
        atomic_parquet(
            raw / f"salt_bridges-{tag}.parquet",
            pd.DataFrame(salt_rows, columns=columns),
        )
    if ligand_rows and cfg["analyses"]["ligand_contacts"].get("save_atom_events", True):
        columns = [
            "time_ps",
            "ligand_atom",
            "protein_atom",
            "protein_residue",
            "distance_A",
        ]
        atomic_parquet(
            raw / f"ligand_contacts-{tag}.parquet",
            pd.DataFrame(ligand_rows, columns=columns),
        )
    if coords:
        tmp = cache / f"aligned-{tag}.tmp.npz"
        # float32 aligned coordinates are accurate enough for RMSF while
        # roughly halving this cache's storage footprint.
        np.savez_compressed(tmp, xyz=np.stack(coords), time_ps=np.asarray(times, float))
        tmp.replace(cache / f"aligned-{tag}.npz")
