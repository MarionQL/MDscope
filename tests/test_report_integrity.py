"""Preserve the completed-statistics boundary when inputs change."""
import pandas as pd
import pytest

from mdscope.reporting.figure_inventory import statistics_status
from mdscope.statistics import analyze_project
from test_presentation_controls import _fixture


def test_inputs_changed_during_statistics_leave_an_incomplete_manifest(tmp_path, monkeypatch):
    import mdscope.statistics as statistics
    cfg, paths, _ = _fixture(tmp_path)
    cfg['statistics']['bootstrap_draws']=20
    source=paths['condition'][0]/'summary/timeseries.parquet'
    def change_input(*args):
        data=pd.read_parquet(source);data['rmsd_A']+=20.;data.to_parquet(source)
    monkeypatch.setattr(statistics,'analyze_profiles',change_input)
    with pytest.raises(RuntimeError,match='inputs changed'):
        analyze_project(cfg,tmp_path,paths)
    assert statistics_status(cfg,tmp_path,'time_series')=='incomplete'
