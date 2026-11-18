"""Individual-figure allowlists must not change metric semantics or load excluded families."""
from copy import deepcopy
from pathlib import Path
from runpy import run_path

from PIL import Image
import pytest
import yaml

from mdscope.config import DEFAULTS, figure_enabled, load_config, plot_enabled
from mdscope.reporting.plot_ids import id_from_figure_path, validate_plot_ids
from mdscope.reporting.plotting import plot_all, plot_interactions
from mdscope.reporting.figure_inventory import figure_status


def _fake_data(tmp_path):
    root = Path(__file__).resolve().parents[1]
    functions = run_path(str(root / "scripts" / "plot_preview.py"))
    return functions["create_fixture"](tmp_path)


def test_exact_figure_allowlist_preserves_legacy_defaults():
    cfg = deepcopy(DEFAULTS)
    cfg["plots"]["enabled"] = True
    assert plot_enabled(cfg, "time_series")
    assert figure_enabled(cfg, "time_series:rmsd_A:timeseries")
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    assert plot_enabled(cfg, "time_series")
    assert not plot_enabled(cfg, "profiles")
    assert figure_enabled(cfg, "time_series:rmsd_A:timeseries")
    assert not figure_enabled(cfg, "time_series:rmsd_A:distribution")
    cfg["plots"]["families"] = {"time_series": False}
    assert not figure_enabled(cfg, "time_series:rmsd_A:timeseries")


@pytest.mark.parametrize("item", [
    "time_series:rmsd_A:bogus", "profiles:foo:bad",
    "interactions:kind:unregistered", "sdc4:anything",
    "solvation:no_such_plot", "diffusion:nope",
    "unknown:foo:timeseries", "time_series::timeseries",
])
def test_unsupported_plot_ids_are_rejected(item):
    with pytest.raises(ValueError, match="unsupported plot ID"):
        validate_plot_ids([item])


def test_named_plot_allowlist_is_validated_before_trajectory_work(tmp_path, project):
    path, root = project("selected_fig")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"] = {"enabled": True, "include": [
        "time_series:rmsd_A:timeseries", "solvation:shell_rdf",
        "profiles:residue_membrane_distance__membrane:tiny:distance_heatmap",
    ]}
    path.write_text(yaml.safe_dump(cfg))
    result = load_config(path)
    assert result["plots"]["include"] == cfg["plots"]["include"]
    assert not root.exists()


def test_only_selected_time_series_is_saved(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, tmp_path, repouts)
    produced = sorted((tmp_path / "project_summary/figures").rglob("*.png"))
    assert [path.name for path in produced] == ["rmsd_A__timeseries.png"]
    assert id_from_figure_path(produced[0]) == "time_series:rmsd_A:timeseries"
    with Image.open(produced[0]) as figure:
        figure.verify()


def test_selected_distribution_does_not_write_time_series(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:distribution"]
    cfg["plots"]["distributions"] = True
    plot_all(cfg, tmp_path, repouts)
    produced = sorted((tmp_path / "project_summary/figures").rglob("*.png"))
    assert [path.name for path in produced] == ["rmsd_A__distribution.png"]


def test_interaction_filter_is_pushed_before_loading_tables(tmp_path, monkeypatch):
    import mdscope.reporting.plotting as module

    seen = []
    monkeypatch.setattr(module, "interaction_tables", lambda *args, **kwargs:
                        (seen.append(kwargs["include_kinds"]) or {}))
    cfg = {"plots": {"enabled": True, "include": [
        "interactions:hydrogen_bonds:occupancy_heatmap",
    ]}, "statistics": {"enabled": False}, "systems": {
        "condition": {"replicates": [{"label": "rep1"}]},
    }}
    plot_interactions(cfg, tmp_path, {"condition": [tmp_path / "rep1"]})
    assert seen == [{"hydrogen_bonds"}]


def test_old_generated_figure_excluded_when_no_longer_selected(tmp_path):
    rep = tmp_path / "condition" / "rep1"
    rep.mkdir(parents=True)
    dest = tmp_path / "project_summary/figures/time_series/rmsd_A__timeseries.png"
    dest.parent.mkdir(parents=True)
    Image.new("RGB", (80, 50), "white").save(dest)
    cfg = {"plots": {"enabled": True, "include": [
        "time_series:radius_of_gyration_A:timeseries",
    ]}, "statistics": {"enabled": False}}
    assert figure_status(dest, cfg, tmp_path, {"condition": [rep]}) == "disabled"


def test_solvation_not_selected_never_opens_rdf(tmp_path, monkeypatch):
    import mdscope.reporting.solvation_figures as module

    def forbidden(*args, **kwargs):
        raise AssertionError("solvation RDF should not be loaded")
    monkeypatch.setattr(module, "solvation_rdf_figure", forbidden)
    monkeypatch.setattr(module, "solvation_timeseries_figure", lambda *args: None)
    cfg = {"plots": {"enabled": True, "include": ["solvation:water_and_ions"]}}
    module.plot_solvation(cfg, tmp_path, {"condition": []})
