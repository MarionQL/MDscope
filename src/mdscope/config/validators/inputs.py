"""Independent analysis validators, in historical evaluation order."""
from pathlib import Path
import numpy as np
from ...errors import InputPathError
from ..runs import _validate_path_component, analysis_enabled


def validate_potassium_bridge_audit(cfg):
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        analyses = cfg["analyses"]
        solvation = analyses["solvation"]
        motif = solvation["targets"].get("motif", {}).get("acidic_patch", {})
        hs = analyses["sdc4_mechanisms"]["hs_auto"]
        if not (solvation.get("water_selection") and isinstance(motif, dict)
                and motif.get("sequence") and motif.get("selection")
                and hs.get("glycan_selection") and hs.get("attachment_selection")):
            raise ValueError("potassium_bridge_audit requires solvation water and acidic_patch "
                             "selectors and sdc4_mechanisms.hs_auto selections")


def validate_lipid_contacts(cfg):
    contacts = cfg.get('analyses', {}).get('lipid_contacts', {})
    if contacts.get('enabled', False):
        if not contacts.get('selection') or not contacts.get('species_selections'):
            raise ValueError('lipid_contacts requires selection and species_selections')
        if float(contacts.get('cutoff_A', 5.)) <= 0 or int(contacts.get('batch_atoms', 4096)) < 1:
            raise ValueError('lipid_contacts cutoff_A and batch_atoms must be positive')


def validate_input_systems(cfg, base, *, check_paths):
    systems = cfg.get("systems")
    if not systems:
        raise ValueError("Config must define at least one system under systems:")
    for field in ("stride", "chunk_frames", "start", "stop"):
        value = cfg["trajectory"].get(field)
        if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
            raise ValueError(f"trajectory.{field} must be an integer")
    if cfg["trajectory"].get("stride", 1) < 1:
        raise ValueError("trajectory.stride must be >= 1")
    if cfg["trajectory"].get("chunk_frames", 1000) < 1:
        raise ValueError("trajectory.chunk_frames must be >= 1")
    start, stop = cfg['trajectory'].get('start'), cfg['trajectory'].get('stop')
    # Preserve supported Python negative indexing, which needs trajectory
    # length to resolve. Reject an obviously empty forward slice early.
    if start is not None and stop is not None and int(start) >= 0 and int(stop) >= 0 and int(stop) <= int(start):
        raise ValueError('trajectory.stop must exceed start for a forward frame window')
    time_spec = cfg.get("time", {}) or {}
    if time_spec.get("dt_ps") is not None:
        dt_ps = float(time_spec["dt_ps"])
        if not np.isfinite(dt_ps) or dt_ps <= 0:
            raise ValueError("time.dt_ps must be a finite value > 0")
    origin_ps = float(time_spec.get("origin_ps", 0.0))
    if not np.isfinite(origin_ps):
        raise ValueError("time.origin_ps must be finite")
    for name, spec in systems.items():
        if not isinstance(name, str):
            raise ValueError("System names must be strings")
        _validate_path_component(name, "system name")
        if not isinstance(spec, dict):
            raise ValueError(f"System {name!r} must be a mapping")
        if not spec.get("topology") or not spec.get("replicates"):
            raise ValueError(
                f"System {name!r} must define topology and at least one replicate"
            )
        top = Path(spec["topology"])
        top = top if top.is_absolute() else (base / top).resolve()
        if check_paths and not top.is_file():
            raise InputPathError(f"System {name!r} topology does not exist: {top}")
        labels = []
        for index, rep in enumerate(spec["replicates"], 1):
            if isinstance(rep, dict):
                if not rep.get("trajectory"):
                    raise ValueError(
                        f"System {name!r} replicate {index} is missing trajectory"
                    )
                trajectory_value = rep["trajectory"]
                label = str(rep.get("label", f"rep{index}"))
            else:
                trajectory_value = rep
                label = f"rep{index}"
            _validate_path_component(
                label,
                f"System {name!r} replicate {index} label",
            )
            labels.append(label)
            traj = Path(trajectory_value)
            traj = traj if traj.is_absolute() else (base / traj).resolve()
            if check_paths and not traj.is_file():
                raise InputPathError(
                    f"System {name!r} replicate {label!r} trajectory does not exist: {traj}"
                )
            if isinstance(rep, dict) and rep.get("msd_trajectory"):
                msd_path = Path(rep["msd_trajectory"])
                if not msd_path.is_absolute():
                    msd_path = (base / msd_path).resolve()
                if check_paths and not msd_path.is_file():
                    raise InputPathError(
                        f"Unmodified membrane MSD trajectory does not exist: {msd_path}"
                    )
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        if duplicates:
            raise ValueError(
                f"System {name!r} has duplicate replicate labels: "
                + ", ".join(duplicates)
            )


def validate_inputs(cfg, base, *, check_paths):
    validate_potassium_bridge_audit(cfg)
    validate_lipid_contacts(cfg)
    validate_input_systems(cfg, base, check_paths=check_paths)
