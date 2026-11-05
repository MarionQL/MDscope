"""Single-pass trajectory coordinator: analysis-specific metrics live in separate modules.

Compatibility: raw/summary output schema retained; atomic NPZ checkpoints
include committed artifact fingerprints and replicate identity.
"""
from __future__ import annotations
from pathlib import Path
from time import perf_counter
import json
from zipfile import BadZipFile
import MDAnalysis as mda
import numpy as np
import pandas as pd
from ..config import analysis_groups, analysis_runs, atomic_parquet,ensure_dir,time_ps
from .progress import track
from ..core.profiling import PhaseProfiler
from ..core.trajectory_sources import open_msd_universe, paired_frames
from ..core.checkpoints import (validate_progress, validate_committed_parts,
                                chunk_artifacts, validate_artifact_inventory,
                                numbered_artifacts, prune_artifact_orphans)
from ..core.coordinate_contracts import CoordinateSpace, require_coordinate_space, msd_source_space
from ..observables.msd import initialize_msd
from .fused_membrane import FusedMembrane, enabled_names as fused_enabled_names
from .fused_extensions import FusedExtensions, enabled_names as extension_enabled_names
from .fused_distance import FusedDistances, enabled_names as distance_enabled_names
from ..mechanisms.sdc4_mechanisms import SDC4Mechanisms, enabled as sdc4_enabled
from ..core.alignment import _apply_transform,_fit_transform,_valid_box
from ..observables.rmsd import compute_rmsd
from ..observables.rmsf import Welford
from ..observables.radius_of_gyration import compute_rg
from ..observables.density_1d import DensityProfileAccumulator
from ..observables.membrane import initialize_membrane,membrane_metrics
from ..interactions import chain_label,initialize_ligand_contacts,initialize_salt_bridges,ligand_contact_events,residue_label,salt_bridge_events
def _checkpoint_payload(path, processed, next_part, rmsf, groups, densities, membrane,
                        rep, artifacts):
    artifacts.update(chunk_artifacts(path.parent.parent, next_part - 1))
    payload = {
        "processed": np.asarray(processed),
        "next_part": np.asarray(next_part),
        "w_n": np.asarray(rmsf.n), "w_mean": rmsf.mean, "w_M2": rmsf.M2,
        "artifact_inventory_version": np.asarray(1),
        "artifact_inventory": np.asarray(json.dumps(artifacts, sort_keys=True)),
        "replicate_identity": np.asarray([rep.system, rep.label]),
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


def _restore_checkpoint(checkpoint, rmsf, groups, densities, membrane):
    with np.load(checkpoint, allow_pickle=False) as cp:
        if cp["w_mean"].shape != rmsf.mean.shape or cp["w_M2"].shape != rmsf.M2.shape:
            raise ValueError("checkpoint accumulator shape does not match the current selection")
        counts = (["processed", "next_part", "w_n"]
                  + [f"g_{key.replace('/', '_')}_n" for key in groups]
                  + [key for key in cp.files if key.startswith("density_") and key.endswith("_frames")])
        for key in counts:
            value = cp[key]
            if value.shape != () or value.dtype.kind not in "iu":
                raise ValueError(f"Checkpoint counter {key} must be an integer scalar")
        processed = int(cp["processed"])
        if any(int(cp[key]) != processed for key in counts[2:]):
            raise ValueError("Checkpoint accumulator counts disagree with processed frames")
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
        return processed, int(cp["next_part"])


def validate_stream_inventory(rep, repout: Path, *, completed=False) -> dict:
    """Validate committed chunks without opening a topology or trajectory."""
    checkpoint = repout / "cache" / "stream_checkpoint.npz"
    try:
        with np.load(checkpoint, allow_pickle=False) as cp:
            counter = cp["next_part"]
            if counter.shape != () or counter.dtype.kind not in "iu" or int(counter) < 1:
                raise ValueError("Invalid checkpoint part counter")
            next_part = int(counter)
            if "replicate_identity" in cp and not np.array_equal(
                    cp["replicate_identity"], [rep.system, rep.label]):
                raise ValueError("Checkpoint replicate identity does not match requested system/label")
            if "artifact_inventory_version" not in cp:
                if completed and "artifact_inventory" not in cp:
                    # Completed legacy runs need not resume from this state.
                    raw = repout / "raw"
                    if next_part > sum(1 for _ in raw.glob("timeseries-part-*.parquet")):
                        raise ValueError("Legacy checkpoint references missing committed timeseries parts")
                    validate_committed_parts(raw, next_part)
                    return {}
                raise ValueError("Legacy checkpoint has no committed artifact inventory; archive this replicate before rerunning")
            if int(cp["artifact_inventory_version"].item()) != 1:
                raise ValueError("Unsupported committed artifact inventory version")
            if "replicate_identity" not in cp:
                raise ValueError("Missing checkpoint replicate identity")
            artifacts = json.loads(str(cp["artifact_inventory"].item()))
        if not isinstance(artifacts, dict) or next_part > len(artifacts):
            raise ValueError("Invalid committed artifact counters")
        validate_committed_parts(repout / "raw", next_part)
        validate_artifact_inventory(repout, artifacts, next_part)
    except (OSError, ValueError, KeyError, EOFError, TypeError, IndexError, BadZipFile) as exc:
        raise RuntimeError(
            f"Invalid streaming checkpoint {checkpoint}: {exc}. Existing output was preserved."
        ) from exc
    return artifacts


def run_streaming_pass(rep, cfg: dict, outroot: Path, cfg_hash=None):
    # Protein, membrane, density, salt-bridge, and ligand-contact
    # calculations share one trajectory pass.
    # Reading a large XTC once is substantially cheaper than reopening
    # it independently for each simple metric.
    repout = ensure_dir(outroot / rep.system / rep.label)
    raw, cache = ensure_dir(repout / "raw"), ensure_dir(repout / "cache")
    checkpoint = cache / "stream_checkpoint.npz"
    artifacts = {}
    # Validate persisted file identity before constructors, output writes or
    # orphan pruning. Older callers may still supply the ignored cfg_hash.
    if checkpoint.exists():
        artifacts = validate_stream_inventory(rep, repout)
    elif any(numbered_artifacts(repout)):
        raise RuntimeError(
            f"Unsafe streaming restart in {repout}: checkpoint is absent but numbered artifacts exist. "
            "Existing output was preserved."
        )
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

    # Resolve metric-specific group filters once; grouping itself is shared.
    metric_groups = {}
    for metric in ("rmsd", "radius_of_gyration"):
        requested = cfg["analyses"][metric].get("groups")
        if requested is not None:
            if not isinstance(requested, list) or not all(isinstance(x, str) for x in requested):
                raise ValueError(f"analyses.{metric}.groups must be a list of group IDs")
            unknown = set(requested) - set(groups)
            if unknown:
                raise ValueError(f"Unknown {metric} group IDs: {sorted(unknown)}")
            metric_groups[metric] = set(requested)
        else:
            metric_groups[metric] = None

    processed = part = 0
    if checkpoint.exists():
        try:
            processed, part = _restore_checkpoint(
                checkpoint,
                rmsf_acc,
                groups,
                densities,
                membrane,
            )
            total_frames = len(universe.trajectory[slice(start,stop,stride)])
            validate_progress(processed, part, chunk, total_frames)
        except (OSError, ValueError, KeyError, EOFError, TypeError, IndexError, BadZipFile) as exc:
            # A corrupted checkpoint should not silently destroy committed data.
            # Operators can archive/remove the replicate to start a clean run.
            raise RuntimeError(
                f"Invalid streaming checkpoint {checkpoint}: {exc}. "
                "Existing output was preserved; move/archive the replicate "
                "directory before restarting from frame zero."
            ) from exc

    # Never prune existing chunks when the checkpoint is absent or references
    # missing committed parts. Both cases require operator recovery.
    try:
        validate_committed_parts(raw, part)
    except ValueError as exc:
        raise RuntimeError(
            f"Unsafe streaming restart in {raw}: {exc}. "
            "Existing output was preserved."
        ) from exc

    # Parts created after the last durable checkpoint are uncommitted.
    # Keep earlier parts intact so a resume never double-counts frames.
    prune_artifact_orphans(repout, part)
    atomic_parquet(
        raw / "protein_chains.parquet",
        pd.DataFrame({"chain": sorted({chain_label(atom) for atom in protein})})
    )
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
        # Validate coordinate provenance before unaligned membrane observers.
        source_space = CoordinateSpace.SOURCE
        require_coordinate_space("membrane observers", source_space, CoordinateSpace.SOURCE)
        # The MSD reader is independent when explicitly configured.
        msd_space = msd_source_space(msd_universe is not None)
        # All raw observers run before fitted feature arrays are constructed.
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

        feature_space = CoordinateSpace.FITTED if alignment_enabled else CoordinateSpace.SOURCE
        if alignment_enabled:
            require_coordinate_space("RMSD", feature_space, CoordinateSpace.FITTED)
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
            # Group IDs are shared by RMSD, RMSF and Rg; each metric may
            # independently select a subset without creating another observer.
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
                    and (metric_groups["rmsd"] is None or group_id in metric_groups["rmsd"])
                    else np.nan
                ),
                "radius_of_gyration_A": (
                    compute_rg(group["rg"])
                    if rg_spec.get("enabled", True) and len(group["rg"])
                    and (metric_groups["radius_of_gyration"] is None or group_id in metric_groups["radius_of_gyration"])
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
                processed,
                part,
                rmsf_acc,
                groups,
                densities,
                membrane,
                rep,
                artifacts,
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
            processed,
            part,
            rmsf_acc,
            groups,
            densities,
            membrane,
            rep,
            artifacts,
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
