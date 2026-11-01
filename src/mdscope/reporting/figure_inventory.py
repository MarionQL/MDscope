"""Inventory generated scientific figures and audit coverage of reportable analyses.

Only persisted artifacts are inspected; never reopen trajectories. Primary
figures are selected deterministically, with all figures preserved in the ZIP.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


PRIMARY_SUFFIXES = (
    "__timeseries.png", "__profile.png", "__system_difference.png",
    "occupancy_heatmap.png", "lifetime_distribution.png",
    "distance_heatmap.png", "lipid_enrichment_heatmap.png",
    "state_populations.png", "state_sampling_diagnostics.png",
    "lipid_density_maps.png", "lipid_enrichment_maps.png",
)
SUPPLEMENT_SUFFIXES = (
    "__distribution.png", "__autocorrelation.png",
    "__bootstrap_table.png", "__comparison_table.png",
)


def figure_inventory(root, system_repouts):
    """List every generated PNG with its origin and preferred presentation."""
    root=Path(root)
    paths=list((root/"project_summary"/"figures").rglob("*.png"))
    paths+=list((root/"project_summary"/"shared_states").glob("*.png"))
    for system,reps in system_repouts.items():
        for rep in reps:
            paths+=list((Path(rep)/"extensions").rglob("*.png"))
            paths+=list((Path(rep)/"figures").rglob("*.png"))
    rows=[]
    for path in sorted(set(paths)):
        name=path.name
        kind=("primary" if name.endswith(PRIMARY_SUFFIXES) else
              "supplement" if name.endswith(SUPPLEMENT_SUFFIXES) else "other")
        rows.append({"figure":str(path.relative_to(root)),
                     "section":str(path.parent.relative_to(root)),
                     "kind":kind,"bytes":path.stat().st_size})
    return pd.DataFrame(rows,columns=["figure","section","kind","bytes"])


def report_coverage(root, system_repouts):
    """Audit which statistical analysis families have corresponding figures.

    This checks family-level presence, not scientific correctness or the
    completeness of every plotted metric.
    """
    root=Path(root)
    fig=figure_inventory(root,system_repouts)
    paths=fig.figure.tolist()
    rows=[]
    specs=(
        ("Time series","statistics/time_series/system_statistics.parquet",
         "project_summary/figures/time_series/"),
        ("Profiles","statistics/profiles/system_profile_statistics.parquet",
         "project_summary/figures/profiles/"),
        ("Interactions","statistics/interactions/",
         "project_summary/figures/interactions/"),
    )
    for family,source,visual in specs:
        path=root/"project_summary"/source
        has_data=(path.is_dir() and any(path.rglob("*.parquet"))) if source.endswith("/") else path.exists()
        has_figure=any(p.startswith(visual) for p in paths)
        rows.append({"analysis":family,"data_available":bool(has_data),
                     "figure_available":bool(has_figure),
                     "status":"covered" if has_data and has_figure else
                              "missing_figure" if has_data else "not_generated"})
    for family,filename in (
        ("Solvation","species_timeseries.parquet"),
        ("SDC4 mechanisms","timeseries.parquet"),
    ):
        folder="solvation" if family=="Solvation" else "sdc4_mechanisms"
        found=any((Path(rep)/"extensions"/folder/filename).exists()
                  for reps in system_repouts.values() for rep in reps)
        has_figure=any("solvation" in p for p in paths) if family=="Solvation" else any(
            "sdc4" in p or "lipid_" in p or "shared_states" in p for p in paths)
        rows.append({"analysis":family,"data_available":found,
                     "figure_available":has_figure,
                     "status":"covered" if found and has_figure else
                              "missing_figure" if found else "not_generated"})
    return pd.DataFrame(rows)
