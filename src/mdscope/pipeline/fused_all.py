"""One trajectory pass for observer-only configurations (no protein streaming)."""
from __future__ import annotations

from pathlib import Path
import MDAnalysis as mda

from .fused_membrane import FusedMembrane, enabled_names as membrane_enabled
from .fused_extensions import FusedExtensions, enabled_names as extension_enabled
from .fused_distance import FusedDistances, enabled_names as distance_enabled
from ..mechanisms.sdc4_mechanisms import SDC4Mechanisms, enabled as sdc4_enabled
from .progress import track


def run_fused_only(rep, cfg, repout, *, only=None):
    """Share source frames among requested observer families only.

    A missing new companion must not rewrite unrelated completed analyses.
    """
    requested = set(only) if only is not None else None
    membrane = membrane_enabled(cfg)
    extensions = extension_enabled(cfg)
    distance = distance_enabled(cfg)
    want_membrane = bool(membrane and (requested is None or set(membrane) & requested))
    selected_extensions = tuple(name for name in extensions
                                if requested is None or name in requested)
    want_distance = bool(distance and (requested is None or "residue_membrane_distance" in requested))
    want_sdc4 = sdc4_enabled(cfg) and (requested is None or "sdc4_mechanisms" in requested)
    if not (want_membrane or selected_extensions or want_distance or want_sdc4):
        return {}
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observers = []
    if want_membrane:
        observers.append(FusedMembrane(u, cfg, repout))
    if selected_extensions:
        observers.append(FusedExtensions(u, cfg, repout, names=selected_extensions))
    if want_distance:
        observers.append(FusedDistances(u, cfg, repout))
    if want_sdc4:
        observers.append(SDC4Mechanisms(u, cfg, repout))
    for observer in observers:
        observer.prune_from(0)
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    count = 0
    for count, ts in enumerate(track(
        frames, total=len(frames), desc=f"{rep.system}/{rep.label}: fused observers",
        enabled=cfg.get("progress", {}).get("enabled", True),
    ), start=1):
        for observer in observers:
            observer.update(ts, count - 1)
        if count % chunk == 0:
            for observer in observers:
                observer.flush((count - 1) // chunk)
    if count % chunk:
        for observer in observers:
            observer.flush(count // chunk)
    return {
        type(observer).__name__: (
            observer.finalize(rep) if isinstance(observer, FusedDistances)
            else observer.finalize()
        )
        for observer in observers
    }
