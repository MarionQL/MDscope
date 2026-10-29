"""Catalog contracts: enabled planning, stable paths and checkpoint aliases."""
from copy import deepcopy
from pathlib import Path
import pytest
from mdscope.analysis.config import DEFAULTS
from mdscope.analysis.registry import FUSED_OUTPUTS, fused_output, fused_output_plan


@pytest.mark.parametrize("name,relative,marker", [
    ("membrane_order", "profiles/membrane_order/timeseries.parquet", "membrane-order"),
    ("bilayer_normal", "extensions/bilayer_normal/timeseries.parquet", "bilayer_normal"),
    ("lipid_enrichment", "extensions/lipid_enrichment/timeseries.parquet", "lipid_enrichment"),
    ("hydrogen_bonds", "interactions/hydrogen_bonds/events.parquet", "hydrogen-bonds"),
    ("lipid_hydrogen_bonds", "interactions/lipid_hydrogen_bonds/events.parquet", "lipid_hydrogen_bonds"),
    ("dssp", "profiles/dssp/assignments.parquet", "profiles"),
    ("sasa", "profiles/sasa/residue_sasa.parquet", "profiles"),
    ("local_environment", "extensions/local_environment/timeseries.parquet", "local_environment"),
    ("sdc4_mechanisms", "extensions/sdc4_mechanisms/timeseries.parquet", "sdc4-mechanisms"),
])
def test_fused_contracts_preserve_legacy_paths(name, relative, marker, tmp_path):
    out = fused_output(tmp_path, name)
    assert out.path == tmp_path / relative
    assert out.checkpoint_key == marker


def test_distance_runs_are_identified_independently(tmp_path):
    a = fused_output(tmp_path, "residue_membrane_distance", run_id="r1")
    b = fused_output(tmp_path, "residue_membrane_distance", run_id="r2")
    assert a.path != b.path
    assert a.checkpoint_key == "residue-membrane-r1"
    with pytest.raises(ValueError, match="requires a run_id"):
        fused_output(tmp_path, "residue_membrane_distance")


def test_disabled_catalog_is_empty(tmp_path):
    cfg = deepcopy(DEFAULTS)
    for settings in cfg["analyses"].values():
        if isinstance(settings, dict):
            settings["enabled"] = False
    assert fused_output_plan(cfg, tmp_path) == ()


@pytest.mark.parametrize("name", [
    "membrane_order", "bilayer_normal", "lipid_enrichment",
    "hydrogen_bonds", "lipid_hydrogen_bonds", "dssp",
    "sasa", "local_environment", "sdc4_mechanisms",
])
def test_single_enabled_analysis_has_one_output(tmp_path, name):
    cfg = deepcopy(DEFAULTS)
    for settings in cfg["analyses"].values():
        if isinstance(settings, dict):
            settings["enabled"] = False
    cfg["analyses"][name]["enabled"] = True
    names = {x.name for x in fused_output_plan(cfg, tmp_path)}
    assert names == {name}


def test_parallel_dssp_stays_specialized(tmp_path):
    cfg = deepcopy(DEFAULTS)
    for settings in cfg["analyses"].values():
        if isinstance(settings, dict):
            settings["enabled"] = False
    cfg["analyses"]["dssp"].update({"enabled": True, "backend": "parallel"})
    assert not any(x.name == "dssp" for x in fused_output_plan(cfg, tmp_path))
