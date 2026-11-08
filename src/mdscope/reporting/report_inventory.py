"""Pure persisted-table inventory and grouping for offline evidence reports.

No trajectory readers, figures or PDF dependencies belong in this module.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from pathlib import Path

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





def inference_status(q_value, p_value=None):
    """Adjusted support is determined by q, never by bootstrap interval signs."""
    try:
        if np.isfinite(float(q_value)):
            return 'Adjusted support' if float(q_value) <= .05 else 'Exploratory'
    except (TypeError, ValueError):
        pass
    return 'Descriptive / unadjusted'


def report_table_status(path, cfg, root):
    """One eligibility policy for PDF tables and current bundle exports."""
    from ..config import presentation_enabled
    from .figure_inventory import figure_family, statistics_status
    path = Path(path)
    family = figure_family(path)
    if not presentation_enabled(cfg, 'report', family):
        return 'disabled'
    if 'statistics' in path.relative_to(root).parts:
        return statistics_status(cfg, root, {'solvation': 'solvation_rdfs'}.get(family, family))
    return 'persisted' if path.is_file() else 'missing_source'


def read_report_table(path, cfg, root):
    if report_table_status(path, cfg, root) not in {'valid', 'legacy_unverified', 'persisted'}:
        return pd.DataFrame()
    return _read(path)


def report_table_inventory(root, cfg):
    """Audit exported project tables, including rejected stale/partial products."""
    import pyarrow.parquet as pq
    root = Path(root)
    rows = []
    for path in sorted((root/'project_summary').rglob('*.parquet')):
        rows.append({'table': str(path.relative_to(root)),
                     'status': report_table_status(path, cfg, root),
                     'rows': pq.ParquetFile(path).metadata.num_rows})
    return pd.DataFrame(rows, columns=['table', 'status', 'rows'])


def scalar_evidence_overview(comparisons):
    """Count reported contrasts by existing category; never rank by CI sign."""
    rows = []
    for category, metrics in _metric_groups(comparisons):
        data = comparisons.loc[comparisons.metric.isin(metrics)]
        estimate = pd.to_numeric(data['mean_difference_system2_minus_system1'], errors='coerce')
        q = pd.to_numeric(data.get('mean_q_value', pd.Series(np.nan, index=data.index)), errors='coerce')
        finite = np.isfinite(estimate)
        supported = finite & np.isfinite(q) & (q <= .05)
        exploratory = finite & np.isfinite(q) & (q > .05)
        descriptive = finite & ~np.isfinite(q)
        rows.append((category, len(data), int(supported.sum()), int(exploratory.sum()),
                     int(descriptive.sum()), int((~finite).sum())))
    return rows
