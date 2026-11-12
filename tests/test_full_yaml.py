"""Protect original simulation paths while expanding the full analysis YAML."""
from pathlib import Path
import yaml


def test_full_yaml_preserves_all_original_simulations():
    repo = Path(__file__).resolve().parents[1]
    original = yaml.safe_load((repo / "synd4_baseline.yaml").read_text())
    full = yaml.safe_load((repo / "synd4_full_mechanisms.yaml").read_text())
    assert full["systems"] == original["systems"]
    assert len(full["systems"]) == 2
    assert sum(len(s["replicates"]) for s in full["systems"].values()) == 6
    for field in ("start", "stop", "stride"):
        assert full["trajectory"][field] == original["trajectory"][field]
    assert full["project"]["output"] != original["project"]["output"]
    for name in (
        "membrane_order", "bilayer_normal", "lipid_enrichment",
        "lipid_hydrogen_bonds", "sdc4_mechanisms", "shared_states",
        "sasa", "lipid_contacts", "local_environment", "coupling_screen",
        "potassium_bridge_audit", "solvation",
    ):
        assert full["analyses"][name]["enabled"], name
    hs = full["analyses"]["sdc4_mechanisms"]["hs_auto"]
    assert hs["enabled"] and hs["expected_chains"] == 6
    assert "report" not in full and full["plots"]["enabled"] is False

    analyses = full["analyses"]
    assert analyses["membrane"]["msd"]["multi_origin"]
    assert analyses["solvation"]["options"]["bridge_spatial_ion"] == "ions/K"
    assert len(analyses["solvation"]["options"]["bridge_target_pairs"]) == 12
    assert analyses["sasa"]["runs"]["environment"]["reference_run"] == "isolated"
    assert not analyses["ligand_contacts"]["enabled"] and full["selections"]["ligand"] is None
    fast = yaml.safe_load((repo / "synd4_full_mechanisms_finalize_fast.yaml").read_text())
    assert fast["project"] == full["project"]
    assert fast["systems"] == full["systems"]
    assert fast["trajectory"] == full["trajectory"]
    assert fast["statistics"]["interaction_exclude_kinds"] == ["lipid_hydrogen_bonds"]
    assert fast["analyses"]["shared_states"]["min_complete_fraction"] == .80
    assert "common_phosphate" not in fast["analyses"]["residue_membrane_distance"]["runs"]
    assert fast["analyses"]["sasa"]["runs"]["isolated"]["atom_sasa_sequences"] == ["VVHPLVPL"]
    # The smaller finalization config targets older persisted replica products.
    variant = {"shared_states", "sasa", "sdc4_mechanisms", "residue_membrane_distance"}
    assert {key: value for key, value in fast["analyses"].items() if key not in variant} == {
        key: value for key, value in analyses.items() if key not in variant
    }
