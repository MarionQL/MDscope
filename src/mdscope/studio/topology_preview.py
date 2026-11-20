"""Nonblocking, topology-only metadata for Studio authoring.

Worker threads never access Streamlit state or trajectories.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json

from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError

_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="mdscope-topology")


def topology_signature(systems, base):
    """Fingerprint topology paths and file revisions, not unrelated analysis settings."""
    base = Path(base).resolve()
    entries = []
    for system_id, spec in systems.items():
        path = Path(spec.get("topology") or "")
        path = path if path.is_absolute() else base / path
        path = path.resolve()
        try:
            stat = path.stat()
            revision = [stat.st_size, stat.st_mtime_ns]
        except OSError:
            revision = None
        entries.append([system_id, str(path), revision])
    return hashlib.sha256(json.dumps(entries).encode()).hexdigest()


def inspect_topologies(systems, base):
    """Return bounded summaries per system, without retaining MDAnalysis universes."""
    base = Path(base).resolve()
    summaries = {}
    for system_id, spec in systems.items():
        path = Path(spec.get("topology") or "")
        path = (path if path.is_absolute() else base / path).resolve()
        if not path.is_file():
            summaries[system_id] = {"error": f"Topology file not found: {path}"}
            continue
        try:
            universe = Universe(str(path))
            atoms = universe.atoms
            identifiers = {}
            for key, attribute in (("SegIDs", "segids"), ("Chain IDs", "chainIDs"),
                                   ("Residue names", "resnames")):
                try:
                    values = getattr(atoms, attribute)
                    identifiers[key] = sorted({str(v).strip() for v in values if str(v).strip()})
                except (AttributeError, KeyError):
                    identifiers[key] = []
            summaries[system_id] = {
                "atoms": len(atoms),
                "residues": len(universe.residues),
                "segments": len(universe.segments),
                "identifiers": identifiers,
                # Retain only a topology location; group selection checks run separately.
                "path": str(path),
            }
        except (OSError, ValueError, TypeError, AttributeError) as exc:
            summaries[system_id] = {"error": f"Could not read {path}: {exc}"}
    return summaries


def submit_topologies(systems, base):
    return _POOL.submit(inspect_topologies, systems, str(base))


def inspect_groups(systems, base, groups):
    """Validate group matches independently for every topology, without analysis preflight."""
    base = Path(base).resolve()
    result = {}
    for system_id, spec in systems.items():
        path = Path(spec.get("topology") or "")
        path = (path if path.is_absolute() else base / path).resolve()
        try:
            universe = Universe(str(path))
        except (OSError, ValueError, TypeError) as exc:
            result[system_id] = {"error": str(exc)}
            continue
        matches = {}
        for group_id, spec in groups.items():
            selection = spec if isinstance(spec, str) else spec.get("selection", "")
            if not selection or not selection.strip():
                matches[group_id] = {"error": "Enter an atom selection"}
                continue
            try:
                atoms = universe.select_atoms(selection)
                matches[group_id] = {"atoms": len(atoms), "residues": len(atoms.residues)}
            except (SelectionError, ValueError, AttributeError, TypeError) as exc:
                matches[group_id] = {"error": str(exc)}
        result[system_id] = matches
    return result


def submit_groups(systems, base, groups):
    return _POOL.submit(inspect_groups, systems, str(base), groups)
