"""Statistics profiling is opt-in and leaves inference results unchanged."""
import json
import pytest

from mdscope.contracts.statistics import STATISTICS_FAMILIES
from mdscope.statistics.dispatch import analyze_project


def test_statistics_family_profile_records_inventory_and_stage_timings(
        tmp_path, monkeypatch):
    monkeypatch.setenv("MDSCOPE_PROFILE", "1")
    cfg = {
        "progress": {"enabled": False},
        "statistics": {"enabled": True, "families": {
            name: False for name in STATISTICS_FAMILIES
        }},
    }
    analyze_project(cfg, tmp_path, {})
    dest = tmp_path / "project_summary/statistics"
    assert json.loads((dest / "manifest.json").read_text())["status"] == "complete"
    prof = json.loads((dest / "statistics_phase_profile.json").read_text())
    stages = {item["stage"]: item for item in prof["family_timings"]}
    assert {"input_inventory", "contrast_coverage", "input_verification"} <= stages.keys()
    assert all(item["seconds"] >= 0 for item in prof["family_timings"])
    assert prof["replicate_count"] == 0
    assert prof["elapsed_seconds"] >= 0


def test_statistics_profile_disabled_has_no_extra_output(tmp_path, monkeypatch):
    monkeypatch.delenv("MDSCOPE_PROFILE", raising=False)
    cfg = {
        "progress": {"enabled": False},
        "statistics": {"enabled": True, "families": {
            name: False for name in STATISTICS_FAMILIES
        }},
    }
    analyze_project(cfg, tmp_path, {})
    assert not (tmp_path / "project_summary/statistics/statistics_phase_profile.json").exists()


@pytest.mark.parametrize("change", ["modify", "create", "delete"])
def test_statistics_refuses_publication_when_grouped_rmsf_input_changes(tmp_path, monkeypatch, change):
    import pandas as pd

    from mdscope.statistics import coverage

    rep = tmp_path / "condition" / "rep1"
    raw = rep / "raw" / "group_rmsf.parquet"
    raw.parent.mkdir(parents=True)
    if change != "create":
        pd.DataFrame({"group": ["A"], "atom_index": [0]}).to_parquet(raw)
    cfg = {"progress": {"enabled": False}, "statistics": {"families": {
        name: False for name in STATISTICS_FAMILIES
    }}}
    original = coverage.write_contrast_coverage

    def concurrent_edit(*args, **kwargs):
        result = original(*args, **kwargs)
        if change == "delete":
            raw.unlink()
        else:
            pd.DataFrame({"group": ["A"], "atom_index": [1]}).to_parquet(raw)
        return result

    monkeypatch.setattr(coverage, "write_contrast_coverage", concurrent_edit)
    with pytest.raises(RuntimeError, match="inputs changed"):
        analyze_project(cfg, tmp_path, {"condition": [rep]})
    manifest = tmp_path / "project_summary/statistics/manifest.json"
    assert json.loads(manifest.read_text())["status"] != "complete"
