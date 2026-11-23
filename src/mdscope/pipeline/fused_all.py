"""One source-coordinate trajectory pass for pending fused observer outputs."""
from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory

import MDAnalysis as mda
import pandas as pd

from ..config import time_ps

from ..config import analysis_runs
from ..config.runs import distance_enabled_names as distance_enabled
from ..config.runs import extensions_enabled_names as extension_enabled
from ..config.runs import membrane_enabled_names as membrane_enabled
from ..core.persistence import ensure_dir
from ..core.progress import track
from ..systems.sdc4.mechanisms import SDC4Mechanisms
from ..systems.sdc4.mechanisms import enabled as sdc4_enabled
from .fused_distance import FusedDistances
from .fused_extensions import FusedExtensions
from .fused_membrane import FusedMembrane


def run_fused_only(rep, cfg, repout, *, only=None, validate_summary_grid=False,
                   on_finalized=None):
    """Replay only pending outputs in a shared source-frame pass.

    `only` accepts registry FusedOutput entries (including their run_id)
    or legacy family-name strings. None means all enabled fused outputs.
    Selected families use short-lived chunk staging so existing completed
    observer chunks are not overwritten during an independent replay.
    When recovering a primary-stream run, validate the saved time grid before
    finalizing any newly produced output. The optional callback records each
    successfully finalized observer before a later observer can fail.
    """
    requested = None if only is None else {
        (entry.name, entry.run_id) if hasattr(entry, "run_id") else (str(entry), "*")
        for entry in only
    }

    def wants(name, run_id=None):
        return (requested is None or (name, "*") in requested
                or (name, run_id) in requested)

    membrane = tuple(name for name in membrane_enabled(cfg) if wants(name))
    extensions = extension_enabled(cfg)
    sasa_run_ids = tuple(
        run.run_id for run in analysis_runs(cfg, "sasa")
        if wants("sasa", None if run.run_id == "default" else run.run_id)
    ) if "sasa" in extensions else ()
    selected_extensions = tuple(
        name for name in extensions
        if (name == "sasa" and sasa_run_ids) or (name != "sasa" and wants(name))
    )
    distance = tuple(
        run_id for run_id in distance_enabled(cfg)
        if wants("residue_membrane_distance", run_id)
    )
    want_sdc4 = sdc4_enabled(cfg) and wants("sdc4_mechanisms")
    if not (membrane or selected_extensions or distance or want_sdc4):
        return {}

    repout = ensure_dir(Path(repout))
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    observers = []
    with ExitStack() as stack:
        stack.callback(u.trajectory.close)  # an open reader locks the file on Windows
        tr = cfg["trajectory"]
        start = 0 if tr.get("start") is None else int(tr["start"])
        u.trajectory[start]
        def staging(name):
            return Path(stack.enter_context(
                TemporaryDirectory(dir=repout, prefix=f".fused-only-{name}-")
            ))

        if membrane:
            options = ({"names": membrane, "cache_dir": staging("membrane")}
                       if requested is not None else {})
            observers.append(FusedMembrane(u, cfg, repout, **options))
        if selected_extensions:
            options = ({"names": selected_extensions, "sasa_run_ids": sasa_run_ids,
                        "cache_dir": staging("extensions")}
                       if requested is not None else {})
            observers.append(FusedExtensions(u, cfg, repout, **options))
        if distance:
            options = ({"run_ids": distance, "cache_dir": staging("distance")}
                       if requested is not None else {})
            observers.append(FusedDistances(u, cfg, repout, **options))
        if want_sdc4:
            options = ({"cache_dir": staging("sdc4")}
                       if requested is not None else {})
            observers.append(SDC4Mechanisms(u, cfg, repout, **options))
        for observer in observers:
            observer.prune_from(0)
        frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
        expected = None
        if validate_summary_grid:
            summary = repout / "summary" / "timeseries.parquet"
            if summary.is_file():
                expected = pd.read_parquet(summary, columns=["time_ps"])["time_ps"]
        if not len(frames) or (expected is not None and len(expected) != len(frames)):
            raise ValueError("Trajectory frame count does not match the replicate frame grid")
        chunk = max(1, int(tr.get("chunk_frames", 1000)))
        count = 0
        for count, ts in enumerate(track(
            frames, total=len(frames), desc=f"{rep.system}/{rep.label}: fused observers",
            enabled=cfg.get("progress", {}).get("enabled", True),
        ), start=1):
            if expected is not None and time_ps(ts, count - 1, cfg) != expected.iat[count - 1]:
                raise ValueError("Trajectory times do not match the replicate frame grid")
            for observer in observers:
                observer.update(ts, count - 1)
            if count % chunk == 0:
                for observer in observers:
                    observer.flush((count - 1) // chunk)
        if count % chunk:
            for observer in observers:
                observer.flush(count // chunk)
        results = {}
        for observer in observers:
            result = (observer.finalize(rep) if isinstance(observer, FusedDistances)
                      else observer.finalize())
            results[type(observer).__name__] = result
            if on_finalized is not None:
                on_finalized()
        return results
