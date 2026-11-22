"""User-facing YAML and saved-statistics controls for exploratory regions."""
import json
from copy import deepcopy
from pathlib import Path

import pandas as pd
import pytest
import yaml
from test_residue_response import _inputs

from mdscope.config import load_config, resolve_config
from mdscope.statistics.residue_regions import scan_saved_residue_regions
from mdscope.studio.model import new_document, review


def test_discovery_yaml_survives_studio_round_trip_and_plan(tmp_path):
    doc = new_document()
    settings = {"enabled": True, "profiles": ["sasa_mean__isolated"],
                "chains": ["PROA"], "min_length": 4, "minimum_cosine": 0.8}
    doc["statistics"]["residue_response"] = settings
    text, plan = review(doc, tmp_path)
    assert yaml.safe_load(text) == doc
    assert plan["residue_response"]["settings"] == settings
    assert plan["residue_response"]["status"] == "disabled_statistics"
    assert resolve_config(doc)["statistics"]["residue_response"] == settings


@pytest.mark.parametrize("settings", [
    {"enabled": "yes"}, {"profiles": ["density_1d"]}, {"profiles": "rmsf"}, {"profiles": ["rmsf", "rmsf"]},
    {"chains": [""]}, {"min_length": True}, {"min_length": 2.5},
    {"min_length": 0}, {"minimum_cosine": float("nan")},
    {"minimum_cosine": 1.1}, {"unknown_option": True},
])
def test_invalid_discovery_yaml_rejected(settings):
    doc = new_document()
    doc["statistics"]["residue_response"] = settings
    with pytest.raises(ValueError, match="residue_response"):
        resolve_config(doc)


def _saved(tmp_path, settings=None):
    reps, contrasts = _inputs()
    folder = tmp_path / "profiles"
    folder.mkdir()
    reps.to_parquet(folder / "replicate_profiles.parquet", index=False)
    contrasts.to_parquet(folder / "system_comparisons.parquet", index=False)
    methods = {} if settings is None else {"residue_response": settings}
    (tmp_path / "manifest.json").write_text(json.dumps({
        "status": "complete", "families": {"profiles": "enabled"},
        "methods": methods,
        "replicate_paths": [f"{s}/rep{i}" for s in ("control", "treated")
                            for i in range(1, 4)],
    }))
    return folder


def test_saved_discovery_uses_yaml_settings_and_archives_disabled_products(tmp_path):
    folder = _saved(tmp_path, {"enabled": True, "min_length": 4})
    scan_saved_residue_regions(tmp_path)
    regions = pd.read_parquet(folder / "residue_response_regions.parquet")
    assert [(r.start_residue, r.end_residue) for r in regions.itertuples()] == [(10, 13)]
    first = (folder / "residue_response_matrix.parquet").read_bytes()
    before = (folder / "system_comparisons.parquet").read_bytes()
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["methods"]["residue_response"] = {"enabled": False}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    scan_saved_residue_regions(tmp_path)
    assert not (folder / "residue_response_matrix.parquet").exists()
    archived = list((folder / ".archives").rglob("residue_response_matrix.parquet"))
    assert len(archived) == 1 and archived[0].read_bytes() == first
    assert (folder / "system_comparisons.parquet").read_bytes() == before
    meta = json.loads((folder / "residue_response_manifest.json").read_text())
    assert meta["status"] == "disabled"
    assert meta["settings"]["enabled"] is False


@pytest.mark.parametrize("selection", [{"profiles": []}, {"chains": []}])
def test_explicit_empty_discovery_selection_does_not_mean_all(tmp_path, selection):
    folder = _saved(tmp_path, selection)
    scan_saved_residue_regions(tmp_path)
    assert not (folder / "residue_response_matrix.parquet").exists()
    assert json.loads((folder / "residue_response_manifest.json").read_text())["status"] == "empty_selection"


def test_omitted_settings_preserve_screen_and_explicit_subset_keeps_units(tmp_path):
    folder = _saved(tmp_path)
    scan_saved_residue_regions(tmp_path)
    assert len(pd.read_parquet(folder / "residue_response_regions.parquet")) == 2
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    manifest["methods"]["residue_response"] = {"profiles": ["sasa_mean__isolated"]}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    scan_saved_residue_regions(tmp_path)
    matrix = pd.read_parquet(folder / "residue_response_matrix.parquet")
    assert set(matrix.profile) == {"sasa_mean__isolated"}
    assert set(matrix.difference_system2_minus_system1) == {-2., 2.}
    assert pd.read_parquet(folder / "residue_response_regions.parquet").empty


@pytest.mark.parametrize("name", ["synd4_baseline.yaml", "synd4_full_mechanisms.yaml",
                                  "synd4_full_mechanisms_finalize_fast.yaml"])
def test_maintained_sdc4_configurations_retain_semantics(name):
    path = Path(__file__).parents[1] / name
    original = load_config(path, check_paths=False)
    raw = yaml.safe_load(path.read_text())
    raw.setdefault("statistics", {})["residue_response"] = {"enabled": True}
    explicit = resolve_config(raw, base=path.parent)
    assert original == explicit


def test_plot_inventory_rejects_disabled_or_changed_discovery_settings(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    # No PNG decoding is necessary for disabled/stale inputs.
    path = tmp_path / "project_summary/figures/profiles/residue_response/a__heatmap.png"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"old figure")
    folder = tmp_path / "project_summary/statistics/profiles"
    folder.mkdir(parents=True)
    settings = {"enabled": True, "profiles": None, "chains": None,
                "min_length": 3, "minimum_cosine": 0.7}
    (folder / "residue_response_manifest.json").write_text(json.dumps({"status": "complete", "settings": settings}))
    cfg = {"statistics": {"residue_response": {"enabled": False}}, "plots": {"enabled": True}}
    assert figure_status(path, cfg, tmp_path, {}) == "disabled"
    cfg["statistics"]["residue_response"] = {"min_length": 5}
    assert figure_status(path, cfg, tmp_path, {}) == "stale_configuration"


def test_profile_entry_point_honors_discovery_without_changing_primary_statistics(tmp_path, monkeypatch):
    import mdscope.statistics.profiles as module
    from mdscope.config import DEFAULTS
    cfg = deepcopy(DEFAULTS)
    cfg["progress"]["enabled"] = False
    cfg["statistics"]["bootstrap_draws"] = 10
    cfg["statistics"]["hypothesis_tests"]["enabled"] = False
    paths = {s: [tmp_path / s / f"rep{i}" for i in range(3)] for s in ("left", "right")}
    def profiles(rep, _cfg):
        import numpy as np
        times = np.arange(6.)
        value = np.tile(np.array([1., 2., 3., 4.]), (6, 1)) + int(rep.name[-1])*.1
        if rep.parent.name == "right":value += 1.
        return {name: {"times": times, "values": value, "positions": np.arange(4),
                       "labels": np.array([f"OTHER:ALA{i}" for i in range(4)])}
                for name in ("sasa_mean__isolated", "residue_occlusion__any")}
    monkeypatch.setattr(module, "collect_temporal_profiles", profiles)
    cfg["statistics"]["residue_response"] = {"enabled": True, "chains": ["OTHER"]}
    module.analyze_profiles(cfg, tmp_path, paths)
    folder = tmp_path / "project_summary/statistics/profiles"
    baseline = pd.read_parquet(folder / "system_comparisons.parquet")
    assert not pd.read_parquet(folder / "residue_response_regions.parquet").empty
    cfg["statistics"]["residue_response"] = {"enabled": False}
    module.analyze_profiles(cfg, tmp_path, paths)
    pd.testing.assert_frame_equal(pd.read_parquet(folder / "system_comparisons.parquet"), baseline)
    assert not (folder / "residue_response_matrix.parquet").exists()
