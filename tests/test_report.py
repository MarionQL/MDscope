"""End-to-end PDF/ZIP export tests with small persisted MDscope-style tables."""
from __future__ import annotations

import json
import zipfile

import numpy as np
import pandas as pd

from mdscope.reporting.report import generate_report, _sample_frame_table


def test_report_pdf_and_review_bundle(tmp_path):
    root = tmp_path / "results"
    summary = root / "project_summary"
    stats = summary / "statistics/time_series"
    stats.mkdir(parents=True)
    cfg = {
        "project": {"name": "synthetic_sdc4", "output": str(root)},
        "trajectory": {"start": 1000, "stop": None, "stride": 10},
        "systems": {
            "carb": {"replicates": [{"label": "rep1"}]},
            "nocarb": {"replicates": [{"label": "rep1"}]},
        },
        "plots": {"enabled": False},
        "report": {
            "enabled": True, "bundle_enabled": True,
            "max_metrics_per_section": 3,
            "sample_frames_per_replicate": 5,
            "include_lipid_maps": False,
        },
    }
    (root / "run_metadata.json").write_text(json.dumps({
        "mdscope_version": "test", "config_hash": "abc123"
    }))
    (root / "resolved_config.yaml").write_text("project:\n  name: synthetic_sdc4\n")
    system_paths = {}
    for system, shift in (("carb", 0.), ("nocarb", 3.)):
        folder = root / system / "rep1"
        (folder / "summary").mkdir(parents=True)
        mechanism = folder / "extensions/sdc4_mechanisms"
        mechanism.mkdir(parents=True)
        pd.DataFrame({
            "time_ps": np.arange(12, dtype=float) * 100,
            "sdc4_ectodomain_A_height_A": np.arange(12, dtype=float) + shift,
        }).to_parquet(folder / "summary/timeseries.parquet", index=False)
        (mechanism / "selection_audit.json").write_text(json.dumps({
            "hs_atom_counts": {"HS_A": 12},
            "motif_resindices": {"NxIP_A": [1, 2, 3, 4]},
            "gxxxg_resindices": {"A": [5, 6, 7, 8, 9]},
            "lipid_counts": {"POPC": 12, "CHL1": 7},
        }))
        system_paths[system] = [folder]
    for system, shift in (("carb",0.),("nocarb",1.)):
        sol=system_paths[system][0]/"extensions"/"solvation"
        sol.mkdir(parents=True)
        pd.DataFrame({
            "time_ps":np.arange(12,dtype=float)*100,
            "solvent":"water","target_class":"protein","site_class":"acidic",
            "mean_n_solvent":np.arange(12)*.1+shift,
            "hydrated_fraction":.5
        }).to_parquet(sol/"species_timeseries.parquet")
        pd.DataFrame({
            "solvent":["water"],"target_id":["acidic_1"],
            "target_class":["protein"],"site_class":["acidic"],
            "n_frames":[12],"mean_n_solvent":[2.+shift],
            "sd_n_solvent":[.2],"hydrated_fraction":[.5]
        }).to_parquet(sol/"target_summary.parquet")
        pd.DataFrame({"frame":np.arange(12),"time_ps":np.arange(12)*100.,
            "solvent":"water","target_id":"acidic_1","solvent_id":np.arange(12)}
        ).to_parquet(sol/"events.parquet")
    rdfdir=summary/"statistics"/"solvation"
    rdfdir.mkdir(parents=True)
    pd.DataFrame([
        {"system":system,"solvent":"water","target_class":"protein",
         "site_class":"acidic","r_inner_A":float(i),
         "r_outer_A":float(i+1),"mean_g_bulk":1.+shift+.1*i,
         "sd_g_bulk":.2,"n_replicates":3}
        for system,shift in (("carb",0.),("nocarb",1.))
        for i in range(5)
    ]).to_parquet(rdfdir/"system_rdf.parquet")
    metric = "sdc4_ectodomain_A_height_A"
    pd.DataFrame([
        {"system": system, "group": None, "metric": metric,
         "n_replicates": 1, "mean": shift + 5.5,
         "mean_ci_low": shift + 4., "mean_ci_high": shift + 7.}
        for system, shift in (("carb", 0.), ("nocarb", 3.))
    ]).to_parquet(stats / "system_statistics.parquet", index=False)
    pd.DataFrame([{
        "system_1": "carb", "system_2": "nocarb",
        "group": None, "metric": metric,
        "mean_difference_system2_minus_system1": 3.,
        "mean_difference_ci_low": 2., "mean_difference_ci_high": 4.,
        "mean_p_value": .05, "mean_q_value": .07,
    }]).to_parquet(stats / "system_comparisons.parquet", index=False)
    pd.DataFrame([
        {"system": system, "replicate": "rep1", "group": None,
         "metric": metric, "mean": shift+5.5,
         "mean_ci_low": shift+4., "mean_ci_high": shift+7.}
        for system, shift in (("carb", 0.), ("nocarb", 3.))
    ]).to_parquet(stats / "replicate_statistics.parquet", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir=summary/"figures"/"time_series"
    figdir.mkdir(parents=True)
    fig,ax=plt.subplots(figsize=(4,2))
    ax.plot([0,1,2],[1,2,1])
    ax.set_xlabel("Time (ns)")
    fig.savefig(figdir/"rmsd_A__timeseries.png")
    plt.close(fig)
    output = generate_report(cfg, root, system_paths)
    pdf = summary / "report/MDscope_evidence_report.pdf"
    bundle = summary / "report/MDscope_review_bundle.zip"
    assert output["pdf"] == str(pdf)
    assert output["bundle"] == str(bundle)
    assert pdf.read_bytes().startswith(b"%PDF-")
    assert pdf.stat().st_size > 3000
    with zipfile.ZipFile(bundle) as z:
        names = set(z.namelist())
        assert "MDscope_evidence_report.pdf" in names
        assert "project_summary/statistics/time_series/system_comparisons.csv" in names
        assert "carb/rep1/summary/timeseries_sample.csv" in names
        assert "carb/rep1/extensions/sdc4_mechanisms/selection_audit.json" in names
        assert "resolved_config.yaml" in names
        assert "carb/rep1/extensions/solvation/events.parquet" in names
        assert "carb/rep1/extensions/solvation/target_summary.csv" in names
        assert "project_summary/statistics/solvation/system_rdf.csv" in names
        assert "project_summary/figures/time_series/rmsd_A__timeseries.png" in names
        assert "project_summary/report/figure_manifest.csv" in names
        assert "project_summary/report/figure_coverage.csv" in names
        manifest=pd.read_csv(z.open("project_summary/report/figure_manifest.csv"))
        assert "project_summary/figures/time_series/rmsd_A__timeseries.png" in manifest.figure.tolist()

        sampled = pd.read_csv(z.open("carb/rep1/summary/timeseries_sample.csv"))
        assert len(sampled) == 5
        assert "sdc4_ectodomain_A_height_A" in sampled
        assert not any(name.endswith(".xtc") for name in names)
    # The command must be safely repeatable without archiving or rerunning.
    again = generate_report(cfg, root, system_paths, force=True)
    assert again["pdf"] == output["pdf"]


def test_bounded_sample_preserves_endpoints():
    df = pd.DataFrame({"time_ps": np.arange(100), "value": np.arange(100)})
    sampled = _sample_frame_table(df, 7)
    assert len(sampled) == 7
    assert sampled.iloc[0].time_ps == 0
    assert sampled.iloc[-1].time_ps == 99
