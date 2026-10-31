"""Offline report-data layer remains usable without PDF/trajectory dependencies."""
import json
import numpy as np
import pandas as pd
from mdscope.analysis.reporting.report_inventory import (
    _read, _sample_frame_table, _summary_inventory, _metric_groups,
)
from mdscope.analysis.reporting.report import _metric_groups as report_metric_groups

def test_inventory_prefers_finalized_summary_without_reading_trajectories(tmp_path):
    rep = tmp_path / "condition" / "rep1"
    (rep / "summary").mkdir(parents=True)
    (rep / "extensions/sdc4_mechanisms").mkdir(parents=True)
    pd.DataFrame({"time_ps": [100.0, 200.0, 300.0]}).to_parquet(
        rep / "summary/timeseries.parquet"
    )
    (rep / "extensions/sdc4_mechanisms/selection_audit.json").write_text(
        json.dumps({"hs_atom_counts": {"a": 5}, "motif_resindices": {"a": [1]}})
    )
    inventory = _summary_inventory(tmp_path, {}, {"condition": [rep]})
    row = inventory.iloc[0]
    assert row["frames"] == 3
    assert row["first_ps"] == 100.0 and row["last_ps"] == 300.0
    assert row["hs_chains"] == 1 and row["nxip_motifs"] == 1
    assert row["mechanisms"] == False

def test_metric_groups_consistent_with_report_public_facade():
    data = pd.DataFrame({"metric": [
        "sdc4_tm_crossing_angle_deg", "some_unclassified_metric",
        "sdc4_tm_crossing_angle_deg", np.nan
    ]})
    assert _metric_groups(data) == report_metric_groups(data)
    assert _metric_groups(data)[-1] == (
        "Other available scalar observables", ["some_unclassified_metric"]
    )
    assert _read(__import__("pathlib").Path("/nonexistent/mdscope.parquet")).empty
    assert len(_sample_frame_table(pd.DataFrame({"x": range(101)}), 5)) == 5
