"""SDC4 saved-data reporting exports, separate from replica scientific products.

These CSVs are presentation summaries regenerated from persisted Parquet; they
are not primary MDscope numerical products or statistical inference inputs.
Historical extension CSVs are intentionally left untouched.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..core.persistence import archive_output, atomic_text


def export_region_replica_summaries(root, system, reps, tables):
    """Write independent per-replica reporting summaries with legacy columns."""
    output_root = (Path(root) / "project_summary" / "reporting_exports"
                   / "sdc4" / str(system))
    for rep in reps:
        repname = Path(rep).name
        frame = tables.get(repname)
        destination = output_root / repname / "lipid_region_replicate_summary.csv"
        if frame is None or "moiety" not in frame:
            _archive_unavailable_export(root, destination)
            continue
        contacts = frame.loc[frame.moiety.eq("any")].copy()
        if contacts.empty:
            _archive_unavailable_export(root, destination)
            continue
        for column in ("matched_local_contact_fraction",
                       "matched_local_composition_enrichment"):
            if column not in contacts:
                contacts[column] = np.nan
        summary = contacts.groupby(
            ["region", "species", "leaflet"], as_index=False
        ).agg(
            mean_contact_fraction=("contact_fraction", "mean"),
            mean_local_enrichment=("local_composition_enrichment", "mean"),
            mean_matched_local_contact_fraction=("matched_local_contact_fraction", "mean"),
            mean_matched_local_enrichment=("matched_local_composition_enrichment", "mean"),
            n_matched_local_frames=("matched_local_contact_fraction", "count"),
            n_matched_enrichment_frames=("matched_local_composition_enrichment", "count"),
            n_frames_with_species=("contact_fraction", "count"),
        )
        atomic_text(destination, summary.to_csv(index=False))


def _archive_unavailable_export(root, destination):
    """Keep an older presentation summary out of current empty/missing results."""
    if destination.is_file():
        archive = Path(root) / '.archives' / destination.parent.relative_to(root)
        archive_output(destination, archive, destination.name)
