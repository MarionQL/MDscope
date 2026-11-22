"""Nonblocking topology metadata and reusable-group validation for Studio.

Topologies are cached by filesystem revision in a bounded worker-local cache.
Group edits reuse their loaded Universe and never invalidate metadata.
"""
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path
import hashlib
import json

from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError

# One worker serializes access to cached MDAnalysis Universes. Keep only
# two loaded revisions to bound memory on large multi-condition projects.
_POOL = ThreadPoolExecutor(max_workers=1, thread_name_prefix="mdscope-topology")


def _topology_path(spec, base):
    path = Path(spec.get("topology") or "")
    return (path if path.is_absolute() else Path(base).resolve() / path).resolve()


def _topology_revision(path):
    try:
        stat = path.stat()
        return (stat.st_size, stat.st_mtime_ns)
    except OSError:
        return None


def topology_signature(systems, base):
    """Only a topology path or file revision change invalidates metadata."""
    entries = [
        [system_id, str(_topology_path(spec, base)),
         _topology_revision(_topology_path(spec, base))]
        for system_id, spec in systems.items()
    ]
    return hashlib.sha256(json.dumps(entries).encode()).hexdigest()


@lru_cache(maxsize=2)
def _cached_universe(path, _revision):  # _revision only keys the cache: an edited file is re-read
    # Called only by the single worker for asynchronous requests.
    return Universe(path)


def _read_universe(path):
    if not path.is_file():
        raise FileNotFoundError(f"Topology file not found: {path}")
    return _cached_universe(str(path), _topology_revision(path))


def inspect_topologies(systems, base):
    """Get stable metadata once per topology version, without trajectories."""
    metadata = {}
    for system_id, spec in systems.items():
        path = _topology_path(spec, base)
        try:
            universe = _read_universe(path)
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
            metadata[system_id] = {"error": f"Could not read {path}: {exc}"}
    return metadata


def inspect_groups(systems, base, groups):
    """Resolve groups against cached topologies, independently of metadata."""
    matches_by_system = {}
    for system_id, spec in systems.items():
        path = _topology_path(spec, base)
        try:
            universe = _read_universe(path)
        except (OSError, ValueError, TypeError) as exc:
            matches_by_system[system_id] = {"error": str(exc)}
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
    return matches_by_system


def inspect_preview(systems, base, groups):
    """Combined synchronous entry point for testing and scripted consumers."""
    return {
        "metadata": inspect_topologies(systems, base),
        "groups": inspect_groups(systems, base, groups),
    }


def submit_topologies(systems, base):
    return _POOL.submit(inspect_topologies, systems, str(base))


def submit_groups(systems, base, groups):
    return _POOL.submit(inspect_groups, systems, str(base), groups)
