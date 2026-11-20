"""Nonblocking, topology-only metadata and group validation for Studio.

Each system topology is parsed once per preview request. Workers never access
Streamlit session state, trajectories, or analysis outputs.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json

from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError

_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="mdscope-topology")


def topology_signature(systems, base):
    """Fingerprint topology paths and revisions, not unrelated analysis settings."""
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


def inspect_preview(systems, base, groups):
    """Read each topology once, extracting metadata and optional group matches."""
    base = Path(base).resolve()
    metadata = {}
    matches_by_system = {}
    for system_id, spec in systems.items():
        path = Path(spec.get("topology") or "")
        path = (path if path.is_absolute() else base / path).resolve()
        if not path.is_file():
            error = f"Topology file not found: {path}"
            metadata[system_id] = {"error": error}
            matches_by_system[system_id] = {"error": error}
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
            metadata[system_id] = {
                "atoms": len(atoms),
                "residues": len(universe.residues),
                "segments": len(universe.segments),
                "identifiers": identifiers,
                "path": str(path),
            }
        except (OSError, ValueError, TypeError, AttributeError) as exc:
            error = f"Could not read {path}: {exc}"
            metadata[system_id] = {"error": error}
            matches_by_system[system_id] = {"error": error}
            continue
        matches = {}
        for group_id, group_spec in groups.items():
            selection = (group_spec if isinstance(group_spec, str)
                         else group_spec.get("selection", ""))
            if not selection or not selection.strip():
                matches[group_id] = {"error": "Enter an atom selection"}
                continue
            try:
                selected_atoms = universe.select_atoms(selection)
                matches[group_id] = {
                    "atoms": len(selected_atoms), "residues": len(selected_atoms.residues),
                }
            except (SelectionError, ValueError, AttributeError, TypeError) as exc:
                matches[group_id] = {"error": str(exc)}
        matches_by_system[system_id] = matches
    return {"metadata": metadata, "groups": matches_by_system}


def inspect_topologies(systems, base):
    """Compatibility helper for metadata-only callers."""
    return inspect_preview(systems, base, {})["metadata"]


def inspect_groups(systems, base, groups):
    """Compatibility helper for group-only callers."""
    return inspect_preview(systems, base, groups)["groups"]


def submit_preview(systems, base, groups):
    """One background job per requested configuration, never two competing reads."""
    return _POOL.submit(inspect_preview, systems, str(base), groups)
