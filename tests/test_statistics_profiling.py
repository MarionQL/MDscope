"""Statistics profiling is opt-in and leaves inference results unchanged."""
import json

from mdscope.statistics import analyze_project


def test_statistics_family_profile_records_inventory_and_stage_timings(
        tmp_path, monkeypatch):
    import mdscope.statistics as statistics
    monkeypatch.setenv("MDSCOPE_PROFILE", "1")
    monkeypatch.setattr(statistics, "write_contrast_coverage",
                        lambda *args: None)
    cfg = {
        "progress": {"enabled": False},
        "statistics": {"enabled": True, "families": {
            name: False for name in statistics.STATISTICS_FAMILIES
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
    import mdscope.statistics as statistics
    monkeypatch.delenv("MDSCOPE_PROFILE", raising=False)
    monkeypatch.setattr(statistics, "write_contrast_coverage",
                        lambda *args: None)
    cfg = {
        "progress": {"enabled": False},
        "statistics": {"enabled": True, "families": {
            name: False for name in statistics.STATISTICS_FAMILIES
        }},
    }
    analyze_project(cfg, tmp_path, {})
    assert not (tmp_path / "project_summary/statistics/statistics_phase_profile.json").exists()
