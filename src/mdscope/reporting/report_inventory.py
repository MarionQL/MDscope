"""Pure persisted-table inventory and grouping for offline evidence reports.

No trajectory readers, figures or PDF dependencies belong in this module.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd

CATEGORIES = (
    ("SDC4 ectodomain orientation and membrane contacts", ("sdc4_ectodomain_",)),
    ("NxIP accessibility and occlusion", ("sdc4_motif_",)),
    ("Individual HS orientation, contacts and hydration", ("sdc4_hs_",)),
    ("TM dimer crossing and GXXXG interface", ("sdc4_tm_", "sdc4_gxxxg_")),
    ("Local membrane thickness, curvature and order", ("sdc4_local_", "sdc4_bulk_")),
    ("Protein-centered lipid enrichment", ("lipid_enrichment_",)),
    ("Bilayer normal deviations", ("bilayer_normal_",)),
    ("Lipid tail order and other membrane metrics", ("membrane_order_", "lipid_order_")),
)

def _read(path):
    return pd.read_parquet(path) if path.is_file() else pd.DataFrame()


def _sample_frame_table(df, max_frames):
    if len(df) <= max_frames:
        return df
    indices = np.linspace(0, len(df)-1, max_frames).astype(int)
    return df.iloc[np.unique(indices)].copy()


def _summary_inventory(root, cfg, system_repouts):
    records = []
    for system, paths in system_repouts.items():
        for repout in paths:
            path = repout / "summary/timeseries.parquet"
            if not path.is_file():
                path = repout / "extensions/sdc4_mechanisms/timeseries.parquet"
            data = _read(path)
            t = pd.to_numeric(data.get("time_ps", pd.Series(dtype=float)),
                              errors="coerce").dropna()
            mechanism = repout / "extensions/sdc4_mechanisms"
            audit_file = mechanism / "selection_audit.json"
            audit = json.loads(audit_file.read_text()) if audit_file.exists() else {}
            records.append({
                "system": system, "replicate": repout.name,
                "frames": len(t),
                "first_ps": float(t.min()) if len(t) else np.nan,
                "last_ps": float(t.max()) if len(t) else np.nan,
                "mechanisms": (mechanism / "timeseries.parquet").exists(),
                "lipid_maps": (mechanism / "lipid_density_maps.npz").exists(),
                "hs_chains": len(audit.get("hs_atom_counts", {})),
                "nxip_motifs": len(audit.get("motif_resindices", {})),
            })
    return pd.DataFrame(records)


def _metric_groups(table):
    if table.empty or "metric" not in table:
        return []
    names = table["metric"].dropna().astype(str).unique().tolist()
    result = []
    used = set()
    for label, prefixes in CATEGORIES:
        found = sorted(m for m in names if m.startswith(prefixes))
        if found:
            result.append((label, found))
            used.update(found)
    other = sorted(set(names) - used)
    if other:
        result.append(("Other available scalar observables", other))
    return result


