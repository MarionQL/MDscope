"""One trajectory pass for observer-only configurations (no protein streaming)."""
from __future__ import annotations

from pathlib import Path
import MDAnalysis as mda

from .fused_membrane import FusedMembrane, enabled_names as membrane_enabled
from .fused_extensions import FusedExtensions, enabled_names as extension_enabled
from .fused_distance import FusedDistances, enabled_names as distance_enabled
from ..mechanisms.sdc4_mechanisms import SDC4Mechanisms, enabled as sdc4_enabled
from ..progress import track


def run_fused_only(rep, cfg, repout):
    """Share the same Universe/frame for *all* optional observer families."""
    if not (membrane_enabled(cfg) or extension_enabled(cfg) or distance_enabled(cfg) or sdc4_enabled(cfg)):
        return {}
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observers = []
    if membrane_enabled(cfg):
        observers.append(FusedMembrane(u, cfg, repout))
    if extension_enabled(cfg):
        observers.append(FusedExtensions(u, cfg, repout))
    if distance_enabled(cfg):
        observers.append(FusedDistances(u, cfg, repout))
    if sdc4_enabled(cfg):
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
