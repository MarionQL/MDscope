"""Read-only persisted-input inspection shared by Studio and finalization.

Registry owns artifact paths/schemas. This module owns readiness facts; browsing
requires saved provenance/markers while compatible direct finalization retains
its existing structural validation policy. No topology/trajectory is opened.
"""
import json
from copy import deepcopy
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import numpy as np
import pandas as pd

from ..config import analysis_runs, config_hash
from ..mechanisms.hs_hydration import temporal_hs_hydration_profiles
from ..mechanisms.residue_occlusion import (
    FRACTION_COLUMNS,
    IDENTITY_COLUMNS,
    OCCLUSION_FILE,
    WATER_COLUMNS,
    residue_environment_pending,
)
from ..observables.density_1d import missing_membrane_density_runs
from ..observables.lipid_enrichment import footprint_pending
from ..registry import analysis_output_plan, artifact_schema, core_saved_contracts, required_artifacts


def check_artifact(path: Path) -> int:
    """Check a required product using its metadata only."""
    size = path.stat().st_size
    if not path.is_file() or not size:
        raise ValueError(f"Missing or empty artifact: {path}")
    if path.suffix == ".parquet":
        from pyarrow.parquet import read_metadata
        read_metadata(path)
    elif path.suffix == ".npy":
        np.load(path, allow_pickle=False, mmap_mode="r")
    elif path.suffix == ".npz":
        # ZIP CRC verification streams payload blocks instead of loading a
        # potentially large order-profile array into RAM.
        with ZipFile(path) as archive:
            if not archive.namelist() or archive.testzip() is not None:
                raise ValueError(f"Invalid NPZ payload: {path}")
    elif path.suffix == ".json" and not isinstance(json.loads(path.read_text()), dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return size


def validate_analysis_output(entry, cfg, *, allow_missing_occlusion=False):
    """Check enabled primary/companion products before reuse or finalization."""
    try:
        occlusion = entry.path.parent / OCCLUSION_FILE
        spec = cfg["analyses"].get(entry.name, {})
        may_add_occlusion = (
            entry.name == "sdc4_mechanisms" and spec.get("residue_occlusion")
            and allow_missing_occlusion and residue_environment_pending(occlusion, spec)
        )
        for path in required_artifacts(entry, cfg):
            if may_add_occlusion and path == occlusion:
                continue
            check_artifact(path)
            if path != entry.path and path.suffix == ".parquet":
                # Companions must satisfy the registry's existing column
                # contracts too; metadata readability alone misses damaged tables.
                from pyarrow.parquet import read_schema
                expected_columns = artifact_schema(entry, path)["columns"]
                if expected_columns is not None:
                    missing = set(expected_columns) - set(read_schema(path).names)
                    if missing:
                        raise ValueError(f"{path}: missing columns {sorted(missing)}")
        if entry.path.suffix == ".parquet":
            from pyarrow.parquet import read_schema
            columns = set(read_schema(entry.path).names)
            expected = {
                "solvation": {"frame", "time_ps", "target_id", "n_solvent"},
                "dssp": {"time_ps", "topology_residue_index", "dssp"},
                "sasa": {"time_ps", "topology_residue_index", "sasa_nm2"},
                "hydrogen_bonds": {"time_ps", "donor", "acceptor"},
                "lipid_hydrogen_bonds": {"time_ps", "donor", "acceptor"},
            }.get(entry.name, {"time_ps"})
            if not expected <= columns:
                raise ValueError(f"{entry.path}: missing columns {sorted(expected - columns)}")
        if entry.name == "membrane_order":
            index = pd.read_parquet(entry.path.parent / "profile_index.parquet")
            for filename in index["file"]:
                path = entry.path.parent / str(filename)
                if path.parent != entry.path.parent:
                    raise ValueError(f"Invalid profile filename: {filename}")
                check_artifact(path)
                with np.load(path, allow_pickle=False) as profile:
                    if not {"time_ps", "values", "positions", "labels"} <= set(profile.files):
                        raise ValueError(f"{path}: missing order-profile arrays")
        if entry.name == "sdc4_mechanisms":
            if not may_add_occlusion and cfg["analyses"]["sdc4_mechanisms"].get("residue_occlusion"):
                columns = set(read_schema(occlusion).names)
                required = {*IDENTITY_COLUMNS, *FRACTION_COLUMNS}
                if spec.get("residue_water_cutoff_A") is not None:
                    required.update(WATER_COLUMNS)
                if not required <= columns:
                    raise ValueError(f"{occlusion}: missing residue occlusion columns")
            path = entry.path.parent / "lipid_density_maps.npz"
            with np.load(path, allow_pickle=False) as maps:
                expected = {"counts", "density_per_nm2", "enrichment_fold",
                            "bulk_density_per_nm2", "species", "leaflet",
                            "x_edges_A", "y_edges_A", "nframes"}
                if not expected <= set(maps.files):
                    raise ValueError(f"{path}: missing density-map arrays")
        if entry.name == "residue_membrane_distance":
            directory = entry.path.parent
            matrix = np.load(entry.path, allow_pickle=False, mmap_mode="r")
            times = np.load(directory / "time_ps.npy", allow_pickle=False)
            from pyarrow.parquet import read_metadata
            residues = read_metadata(directory / "residues.parquet").num_rows
            if (matrix.ndim != 2 or times.ndim != 1 or matrix.shape != (len(times), residues)
                    or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)):
                raise ValueError(f"{directory}: distance matrix/time/residue shape or time-axis mismatch")
    except (OSError, ValueError, KeyError, EOFError, TypeError, IndexError, BadZipFile) as exc:
        raise RuntimeError(
            f"invalid required {entry.name} artifact: {exc}. Existing outputs were preserved; "
            "restore the required product or recompute in a separate output directory."
        ) from exc


def replicate_state_prefix(rep) -> str:
    # Keep the historical replica marker filename. It encodes identity only.
    identity_hash = config_hash({
        "system": rep.system,
        "replicate": rep.label,
    })
    return f"rep-{identity_hash}"


def validate_replica_summary(repout, rep):
    """Validate the common frame-grid and available checkpoint identity facts."""
    timeseries = repout / "summary" / "timeseries.parquet"
    if not timeseries.is_file():
        raise ValueError(f"missing {timeseries}; finish the replica before offline statistics or plotting")
    times = pd.read_parquet(timeseries, columns=["time_ps"])["time_ps"]
    if (times.empty or not np.isfinite(times.to_numpy(dtype=float)).all()
            or not times.is_monotonic_increasing or not times.is_unique):
        raise ValueError(f"invalid time_ps in {timeseries} (empty, nonfinite, duplicated or nonmonotonic)")
    checkpoint = repout / "cache" / "stream_checkpoint.npz"
    if checkpoint.is_file():
        with np.load(checkpoint, allow_pickle=False) as cp:
            if "replicate_identity" in cp and not np.array_equal(
                    cp["replicate_identity"], [rep.system, rep.label]):
                raise ValueError("checkpoint replicate identity mismatch")
            if "processed" in cp:
                count = cp["processed"]
                if count.shape != () or count.dtype.kind not in "iu" or len(times) != int(count):
                    raise ValueError("summary frame count disagrees with checkpoint")
    return times.to_numpy(float)


def validate_saved_analysis(entry, cfg, times):
    """Check registered products and their applicable saved-data prerequisites."""
    core = core_saved_contracts(cfg, entry.name)
    if core:
        for item in core:
            path = entry.path.parent.parent / item["path_within_replica"]
            check_artifact(path)
            table = pd.read_parquet(path, columns=item["columns"])
            if "time_ps" in table:
                selected = table.loc[table["group"] == item["group"]] if item.get("group") is not None else table
                if not np.array_equal(selected["time_ps"].to_numpy(float), times):
                    raise ValueError(f"{path}: time grid or group coverage disagrees with summary")
        return
    validate_analysis_output(entry, cfg)
    if entry.name == "lipid_enrichment" and footprint_pending(
            entry.path, cfg["analyses"]["lipid_enrichment"]):
        raise ValueError("lipid footprint settings changed; run MDscope with this config")
    if entry.name == "hydration":
        temporal_hs_hydration_profiles(entry.path.parents[2], cfg, times)


def validate_additional_products(repout, cfg, times):
    """Check derived prerequisites not represented by registered primary products."""
    temporal_hs_hydration_profiles(repout, cfg, times)
    if missing_membrane_density_runs(repout, cfg):
        raise ValueError("missing membrane-relative density; run MDscope with this config")
    msd = cfg.get("analyses", {}).get("membrane", {}).get("msd", {})
    if msd.get("enabled") and msd.get("multi_origin"):
        for name in ("multi_origin_diffusion.parquet", "multi_origin_msd.parquet"):
            check_artifact(repout / "summary" / name)


def validate_replica_products(cfg, root, replicates):
    # Finalization consumes persisted replica products, not the trajectory.
    # A changed reporting configuration or source revision must not force
    # expensive completed trajectories to be rerun. Validate the actual
    # required input using its persisted structure.
    missing = []
    verbose = cfg.get("progress", {}).get("enabled", True)
    for index, rep in enumerate(replicates, 1):
        if verbose:
            print(f"Validating replica {index}/{len(replicates)} {rep.system}/{rep.label}", flush=True)
        repout = root / rep.system / rep.label
        try:
            times = validate_replica_summary(repout, rep)
            for entry in analysis_output_plan(cfg, repout):
                validate_saved_analysis(entry, cfg, times)
            validate_additional_products(repout, cfg, times)
            if requires_temporal_rmsf_statistics(cfg):
                from ..observables.rmsf import temporal_rmsf_readiness
                ready, reason = temporal_rmsf_readiness(repout, times)
                if not ready:
                    raise ValueError(reason)
        except (OSError, ValueError, KeyError, ImportError, EOFError, TypeError, IndexError, BadZipFile, RuntimeError) as exc:
            missing.append(
                f"{rep.system}/{rep.label}: {exc}; invalid persisted inputs under {repout}; "
                "inspect the affected saved artifact before retrying"
            )
    if missing:
        raise RuntimeError(
            "Cannot finalize: required persisted replica data are missing or invalid. "
            "Existing outputs were not replaced. "
            + "; ".join(missing)
        )


def _verified_saved_configuration(current, saved, system, replica, analysis, run_id="default"):
    """Compare source identity and analysis settings, not reporting options."""
    for key in ("trajectory", "time", "selections", "cache"):
        a, b = current.get(key), saved.get(key)
        if key == "trajectory":
            a = {k: v for k, v in (a or {}).items() if k != "chunk_frames"}
            b = {k: v for k, v in (b or {}).items() if k != "chunk_frames"}
        if key == "selections":
            def scientific_selections(value):
                value = deepcopy(value or {})
                value["groups"] = {name: {k: v for k, v in spec.items()
                                          if k not in {"label", "linestyle"}}
                                   for name, spec in value.get("groups", {}).items()}
                return value
            a, b = scientific_selections(a), scientific_selections(b)
        if a != b:
            return False, key + " configuration differs"
    if system not in saved["systems"]:
        return False, "saved condition identity differs"
    if current["systems"][system].get("topology") != saved["systems"][system].get("topology"):
        return False, "topology source differs"
    def by_label(cfg):
        return {str(item.get("label", f"rep{i}")): {k: v for k, v in item.items() if k != "plot_label"}
                for i, item in enumerate(cfg["systems"][system]["replicates"], 1)}
    if replica not in by_label(saved) or by_label(current)[replica] != by_label(saved)[replica]:
        return False, "replicate input identity differs"
    if analysis in {"rmsd", "rmsf"} and (
            current.get("analyses", {}).get("alignment") != saved.get("analyses", {}).get("alignment")):
        return False, "alignment configuration differs"
    active = {r.run_id: r.spec for r in analysis_runs(current, analysis)}
    previous = {r.run_id: r.spec for r in analysis_runs(saved, analysis)}
    if run_id in active:
        # Labels are presentation; compare the actual selected run, not siblings.
        clean = lambda spec: {k: v for k, v in spec.items() if k not in {"label", "enabled"}}
        if run_id not in previous or clean(active[run_id]) != clean(previous[run_id]):
            return False, "analysis settings differ"
    return True, ""


def inspect_results_directory(document, base, results_directory):
    """Read-only inventory of registered per-replica analysis artifacts.

    Registered products require matching saved provenance and a completed stage
    plus shared summary-time/checkpoint and artifact checks. These are readiness
    facts, not evidence of chemical validity or sampling convergence.
    """
    from ..config import iter_replicates, resolve_config
    from ..registry import core_saved_contracts, required_artifacts, saved_output_plan

    root = Path(results_directory).expanduser()
    if not root.is_absolute():
        root = Path(base) / root
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"Results directory does not exist: {root}")
    cfg = resolve_config(document, base=Path(base))
    saved = None
    saved_config = root / "resolved_config.yaml"
    if saved_config.is_file():
        import yaml
        try:
            saved = resolve_config(yaml.safe_load(saved_config.read_text()), base=Path(base))
        except (OSError, ValueError, TypeError, KeyError, yaml.YAMLError):
            saved = None
    inventory = []
    for replica in iter_replicates(cfg, Path(base) / "analysis.yaml"):
        replica_root = root / replica.system / replica.label
        try:
            times = validate_replica_summary(replica_root, replica)
            validate_additional_products(replica_root, cfg, times)
            summary_error = None
        except (OSError, ValueError, KeyError, TypeError, EOFError, BadZipFile, RuntimeError) as exc:
            times, summary_error = None, str(exc)
        current_entries = list(saved_output_plan(cfg, replica_root))
        saved_entries = []
        if saved is not None:
            try:
                saved_entries = list(saved_output_plan(saved, replica_root))
            except (KeyError, ValueError, TypeError):
                pass
        entries = {entry.checkpoint_key: (entry, cfg) for entry in current_entries}
        for entry in saved_entries:
            entries.setdefault(entry.checkpoint_key, (entry, saved))
        for entry, source_cfg in entries.values():
            required = required_artifacts(entry, source_cfg)
            present = []
            for path in required:
                try:
                    if path.is_file() and path.stat().st_size > 0:
                        present.append(path)
                except OSError:
                    pass  # An inaccessible product is unavailable, never verified.

            verified = False
            reason = "Missing or unreadable saved resolved configuration"
            provenance_unknown = saved is None
            if saved is None and len(present) == len(required):
                try:
                    for path in required:
                        check_artifact(path)
                except (RuntimeError, OSError, ValueError, KeyError, TypeError) as exc:
                    reason, provenance_unknown = str(exc), False
            if saved is not None:
                try:
                    compatible, reason = _verified_saved_configuration(
                        cfg, saved, replica.system, replica.label, entry.name, entry.run_id or "default")
                    if compatible and len(present) == len(required):
                        import json
                        key = (replicate_state_prefix(replica) + "-" +
                               ("summary" if core_saved_contracts(source_cfg, entry.name)
                                else entry.checkpoint_key))
                        marker = root / ".state" / f"{key}.json"
                        try:
                            completed = json.loads(marker.read_text()).get("status") == "complete"
                        except (OSError, ValueError, AttributeError):
                            completed = False
                        if not completed:
                            reason = "No completed analysis stage marker"
                        else:
                            if summary_error:
                                raise ValueError(summary_error)
                            validate_saved_analysis(entry, source_cfg, times)
                            verified = True
                            reason = "Validated stage marker, saved configuration, summary and registered artifacts"
                except (RuntimeError, OSError, ValueError, KeyError, TypeError) as exc:
                    reason = str(exc)
            statistics_ready, statistics_reason = verified, reason
            if verified and entry.name == "rmsf":
                from ..observables.rmsf import temporal_rmsf_readiness
                statistics_ready, statistics_reason = temporal_rmsf_readiness(replica_root, times)
            inventory.append({
                "reason": reason,
                "statistics_ready": statistics_ready,
                "statistics_reason": statistics_reason,
                "summary_valid": summary_error is None,
                "source_spec": deepcopy(next(
                    (run.spec for run in analysis_runs(source_cfg, entry.name)
                     if run.run_id == (entry.run_id or "default")),
                    source_cfg["analyses"].get(entry.name, {}))),
                "expected_replicates": [r.label for r in iter_replicates(cfg, Path(base) / "analysis.yaml")
                                        if r.system == replica.system],
                "scheduled": entry.checkpoint_key in {item.checkpoint_key for item in current_entries},
                "analysis": entry.name,
                "run_id": entry.run_id or "default",
                "condition": replica.system,
                "replica": replica.label,
                "present": len(present),
                "required": len(required),
                "status": ("verified_existing" if verified else
                           "unknown" if provenance_unknown and len(present) == len(required) else
                           "incompatible" if len(present) == len(required) else
                           "partial" if present else "missing"),
                "missing": [str(path.relative_to(root)) for path in required
                            if path not in present],
            })
    return inventory




def requires_temporal_rmsf_statistics(cfg):
    """An exact request is stronger than legacy family-wide skip semantics."""
    from ..config import statistics_family_config
    if statistics_family_config(cfg, "profiles") is None:
        return False
    profiles = cfg.get("statistics", {}).get("metric_selection", {}).get("profiles") if isinstance(
        cfg.get("statistics", {}).get("metric_selection"), dict) else None
    return bool(cfg.get("statistics", {}).get("enabled", False) and profiles and
                any(item["profile"] == "rmsf" or item["profile"].startswith("rmsf__group_")
                    for item in profiles))


def unavailable_saved_statistics(cfg, inventory, *, allow_scheduled=True):
    """Exact temporal RMSF needs more persisted inputs than a static raw plot."""
    from ..config import analysis_enabled
    if not requires_temporal_rmsf_statistics(cfg):
        return []
    if allow_scheduled and analysis_enabled(cfg, "rmsf"):
        return []
    from ..config import iter_replicates
    rows = {(row["condition"], row["replica"]): row for row in inventory or []
            if row["analysis"] == "rmsf"}
    for rep in iter_replicates(cfg, Path("analysis.yaml")):
        row = rows.get((rep.system, rep.label))
        if not row or row["status"] != "verified_existing" or not row.get("statistics_ready", False):
            return ["profiles"]
    return []


def configuration_with_saved_sources(cfg, inventory, *, statistics=False):
    """Use fully verified saved-only sources in downstream readers, never scheduling.

    Returns a detached resolved configuration. The user's enable switches and
    saved resolved_config.yaml are preserved; no replicate task runs here.
    """
    from ..config import iter_replicates
    result = deepcopy(cfg)
    replicas = list(iter_replicates(cfg, Path("analysis.yaml")))
    by_key = {(r["analysis"], r["run_id"], r["condition"], r["replica"]): r
              for r in inventory}
    ready = {}
    for analysis, run_id in sorted({(r["analysis"], r["run_id"]) for r in inventory}):
        rows = [by_key.get((analysis, run_id, rep.system, rep.label)) for rep in replicas]
        if rows and all(row and row["status"] == "verified_existing" and "source_spec" in row
                        and (not statistics or row.get("statistics_ready", True))
                        for row in rows):
            ready.setdefault(analysis, {})[run_id] = rows[0]["source_spec"]
    for analysis, runs in ready.items():
        current = result["analyses"].get(analysis, {})
        if not current.get("enabled", False):
            source = deepcopy(next(iter(runs.values())))
            source["enabled"] = True
            if "default" not in runs:
                source["runs"] = deepcopy(runs)
            result["analyses"][analysis] = source
        elif current.get("runs") is not None:
            active_runs = {run.run_id for run in analysis_runs(result, analysis)}
            for name, spec in runs.items():
                if name != "default" and name not in active_runs:
                    # Saved run specs already include their original base values.
                    current["runs"][name] = deepcopy(spec)
    return result
