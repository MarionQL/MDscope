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
    assert full["trajectory"] == original["trajectory"]
    assert full["project"]["output"] != original["project"]["output"]
    for name in (
        "membrane_order", "bilayer_normal", "lipid_enrichment",
        "lipid_hydrogen_bonds", "sdc4_mechanisms", "shared_states",
        "sasa",
    ):
        assert full["analyses"][name]["enabled"], name
    hs = full["analyses"]["sdc4_mechanisms"]["hs_auto"]
    assert hs["enabled"] and hs["expected_chains"] == 6
    assert full["report"]["enabled"] and full["report"]["bundle_enabled"]
