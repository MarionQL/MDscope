"""Bounded review archive export; no PDF rendering or trajectory dependencies."""
from __future__ import annotations
from pathlib import Path
import tempfile
import zipfile
import pandas as pd
from .report_inventory import _sample_frame_table, report_table_status
from .figure_inventory import figure_inventory, report_coverage, metric_coverage

def _make_bundle(root, pdf, system_repouts, options, cfg=None):
    """Pack sufficient machine-readable evidence without massive trajectories."""
    dest = pdf.parent / "MDscope_review_bundle.zip"
    max_frames = int(options.get("sample_frames_per_replicate", 250))
    with tempfile.NamedTemporaryFile(
        dir=dest.parent, prefix=".review-", suffix=".zip", delete=False
    ) as f:
        tmp = Path(f.name)
    try:
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=6) as archive:
            archive.write(pdf, arcname=pdf.name)
            # Complete machine-readable figure manifest and coverage audit.
            figures=figure_inventory(root,system_repouts,cfg)
            coverage=report_coverage(root,system_repouts,cfg)
            if not (pdf.parent/'figure_manifest.csv').exists():
                archive.writestr("project_summary/report/figure_manifest.csv",
                                 figure_inventory(root,system_repouts,cfg,include_rejected=True).to_csv(index=False))
            archive.writestr("project_summary/report/figure_coverage.csv",
                             coverage.to_csv(index=False))
            archive.writestr("project_summary/report/figure_metric_coverage.csv",
                             metric_coverage(root,cfg,system_repouts).to_csv(index=False))
            for item in figures.itertuples():
                path_=root/item.figure
                archive.write(path_,arcname=item.figure)
            # Use this export's staged manifests, never previous report metadata.
            for sidecar in sorted(pdf.parent.iterdir()):
                if sidecar.suffix in {'.json', '.csv'}:
                    archive.write(sidecar, arcname='project_summary/report/'+sidecar.name)
            for name in ("resolved_config.yaml", "run_metadata.json"):
                path = root / name
                if path.is_file():
                    archive.write(path, arcname=name)
            project = root / "project_summary"
            if project.is_dir():
                for path in sorted(project.rglob("*")):
                    if 'report' in path.relative_to(project).parts:
                        continue
                    if path.is_file() and path.suffix.lower() in {".parquet", ".json"}:
                        if path.suffix == ".parquet":
                            if cfg is not None and report_table_status(path,cfg,root) not in {'valid','legacy_unverified','persisted'}:
                                continue
                            archive.write(path, arcname=str(path.relative_to(root)))
                            if path.stat().st_size > 50_000_000:
                                continue  # preserve full Parquet without materializing huge CSVs
                            df = pd.read_parquet(path)
                            archive.writestr(
                                str(path.relative_to(root).with_suffix(".csv")),
                                df.to_csv(index=False),
                            )
                        else:
                            archive.write(path, arcname=str(path.relative_to(root)))
            for system, paths in system_repouts.items():
                for repout in paths:
                    relative = Path(system) / repout.name
                    for name in (
                        "extensions/sdc4_mechanisms/selection_audit.json",
                        "extensions/sdc4_mechanisms/lipid_density_maps.npz",
                        "extensions/sdc4_mechanisms/hs_conditioned_solvation.parquet",
                    ):
                        path = repout / name
                        if not path.is_file():
                            continue
                        if path.suffix == ".parquet":
                            archive.writestr(
                                str(relative / Path(name).with_suffix(".csv")),
                                pd.read_parquet(path).to_csv(index=False),
                            )
                        else:
                            archive.write(path, arcname=str(relative / name))
                    for name in ("lipid_contacts", "potassium_bridge_audit"):
                        for path in sorted((repout / "extensions" / name).glob("*")):
                            if path.suffix in {".parquet", ".json"}:
                                archive.write(path, arcname=str(relative / path.relative_to(repout)))
                    # Preserve mechanistic event/residence tables as
                    # compressed Parquet. Avoid expanding large sparse
                    # molecule-contact maps into enormous CSVs in memory.
                    mechanism = repout / "extensions" / "sdc4_mechanisms"
                    for name in (
                        "lipid_region_specificity.parquet",
                        "lipid_residue_contacts.parquet",
                        "motif_lipid_residence.parquet",
                        "motif_lipid_residence_intermittent.parquet",
                        "motif_lipid_contact_events.parquet",
                    ):
                        path = mechanism / name
                        if path.is_file():
                            archive.write(
                                path, arcname=str(relative / "extensions" /
                                                  "sdc4_mechanisms" / name)
                            )
                    # Include full compact solvation summaries and bounded
                    # samples of large event/voxel/frame-level tables.
                    sol_dir=repout/"extensions"/"solvation"
                    if sol_dir.is_dir():
                        for sol_file in sorted(sol_dir.glob("*.parquet")):
                            # Preserve every result exactly, without loading
                            # full million-row event tables into pandas.
                            archive.write(sol_file,arcname=str(
                                relative/"extensions"/"solvation"/sol_file.name))
                            if sol_file.stat().st_size > 50_000_000:
                                continue  # full Parquet is already archived
                            sol_data=pd.read_parquet(sol_file)
                            large={"counts","events","episodes","bridges",
                                   "water_hbond_bridges","hydration_sasa",
                                   "water_density_3d","species_timeseries"}
                            if sol_file.stem in large:
                                sol_data=_sample_frame_table(sol_data,max_frames)
                                suffix="_sample.csv"
                            else:
                                suffix=".csv"
                            archive.writestr(
                                str(relative/"extensions"/"solvation"/
                                    (sol_file.stem+suffix)),
                                sol_data.to_csv(index=False))
                        metadata=sol_dir/"metadata.json"
                        if metadata.exists():
                            archive.write(metadata,arcname=str(
                                relative/"extensions"/"solvation"/"metadata.json"))
                    path = repout / "summary/timeseries.parquet"
                    if path.is_file():
                        sampled = _sample_frame_table(pd.read_parquet(path), max_frames)
                        archive.writestr(
                            str(relative / "summary/timeseries_sample.csv"),
                            sampled.to_csv(index=False),
                        )
            archive.writestr(
                "README.txt",
                "MDscope scientific review bundle\n"
                "PDF: MDscope_evidence_report.pdf\n"
                "project_summary/statistics: full-precision exported statistical tables\n"
                "project_summary/report/figure_manifest.csv: all generated PNGs\n"
                "project_summary/report/figure_coverage.csv: family-level visual audit\n"
                "project_summary/shared_states: model and state tables\n"
                "system/replicate/extensions: selection audits, lipid maps, conditioned HS summaries\n"
                "system/replicate/summary/timeseries_sample.csv: uniformly sampled frames "
                f"(at most {max_frames} per replicate)\n"
                "Raw trajectories, TPR files, large caches and full frame-level "
                "time-series tables are intentionally excluded.\n",
            )
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)
    return dest
