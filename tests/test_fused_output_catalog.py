"""Catalog contracts: enabled planning, persisted paths and instance IDs."""
from copy import deepcopy
import pytest
from mdscope.config import DEFAULTS
from mdscope.registry import fused_output, fused_output_plan


@pytest.mark.parametrize("name,relative,marker", [
    ("membrane_order", "profiles/membrane_order/timeseries.parquet", "membrane_order"),
    ("bilayer_normal", "extensions/bilayer_normal/timeseries.parquet", "bilayer_normal"),
    ("lipid_enrichment", "extensions/lipid_footprint/timeseries.parquet", "lipid_enrichment"),
    ("hydrogen_bonds", "interactions/hydrogen_bonds/events.parquet", "hydrogen_bonds"),
    ("lipid_hydrogen_bonds", "interactions/lipid_hydrogen_bonds/events.parquet", "lipid_hydrogen_bonds"),
    ("dssp", "profiles/dssp/assignments.parquet", "dssp"),
    ("sasa", "profiles/sasa/residue_sasa.parquet", "sasa"),
    ("local_environment", "extensions/local_environment/timeseries.parquet", "local_environment"),
    ("sdc4_mechanisms", "extensions/sdc4_mechanisms/timeseries.parquet", "sdc4_mechanisms"),
])
def test_fused_contracts(name, relative, marker, tmp_path):
    out = fused_output(tmp_path, name)
    assert out.path == tmp_path / relative
    assert out.checkpoint_key == marker


def test_distance_runs_are_identified_independently(tmp_path):
    a = fused_output(tmp_path, "residue_membrane_distance", run_id="r1")
    b = fused_output(tmp_path, "residue_membrane_distance", run_id="r2")
    assert a.path != b.path
    assert a.checkpoint_key == "residue_membrane_distance-r1"
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
    cfg["analyses"].setdefault(name, {})["enabled"] = True
    names = {x.name for x in fused_output_plan(cfg, tmp_path)}
    assert names == {name}


def test_parallel_dssp_stays_specialized(tmp_path):
    cfg = deepcopy(DEFAULTS)
    for settings in cfg["analyses"].values():
        if isinstance(settings, dict):
            settings["enabled"] = False
    cfg["analyses"]["dssp"].update({"enabled": True, "backend": "parallel"})
    assert not any(x.name == "dssp" for x in fused_output_plan(cfg, tmp_path))


def test_fused_distance_selective_run_keeps_other_checkpoint_parts(tmp_path):
    """A missing named distance output cannot prune a completed sibling."""
    from mdscope.pipeline.fused_distance import FusedDistances
    from mdscope_test_support import membrane_universe

    cfg = {"selections": {"protein": "name CA"}, "analyses": {
        "residue_membrane_distance": {"enabled": True, "runs": {
            "near": {"membrane_selection": "resname POPC"},
            "far": {"membrane_selection": "resname POPC"},
        }},
    }}
    cache = tmp_path / "cache/fused_distance"
    preserved = cache / "far/part-000000.npz"
    preserved.parent.mkdir(parents=True)
    preserved.write_bytes(b"do not replace")
    observer = FusedDistances(membrane_universe(), cfg, tmp_path, run_ids=("near",))
    assert list(observer.runs) == ["near"]
    near = cache / "near/part-000000.npz"
    near.parent.mkdir(parents=True)
    near.write_bytes(b"stale")
    observer.prune_from(0)
    assert not near.exists()
    assert preserved.read_bytes() == b"do not replace"


def test_fused_distance_rejects_unknown_requested_run_id(tmp_path):
    from mdscope.pipeline.fused_distance import FusedDistances
    from mdscope_test_support import membrane_universe

    cfg = {"selections": {"protein": "name CA"}, "analyses": {
        "residue_membrane_distance": {"enabled": True, "runs": {
            "near": {"membrane_selection": "resname POPC"},
        }},
    }}
    with pytest.raises(ValueError, match="Unknown residue-distance run IDs"):
        FusedDistances(membrane_universe(), cfg, tmp_path, run_ids=("wrong",))
