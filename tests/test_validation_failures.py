"""Actionable config/topology failures precede output creation or replacement."""
import re
from copy import deepcopy
from pathlib import Path

import MDAnalysis as mda
import pytest
import yaml
from mdscope_test_support import freesasa_backend
from mdscope_test_support import invoke_cli as invoke

from mdscope.config import load_config
from mdscope.pipeline.preflight import validate_selections
from mdscope.pipeline.runner import run_one


def _set(cfg, keys, value):
    target = cfg
    for key in keys[:-1]: target = target.setdefault(key, {})
    target[keys[-1]] = value


@pytest.mark.parametrize("keys,value,diagnostic", [
    (("analyses", "lipid_enrichment"), {"enabled": True, "species_selections": {"POPC": "all"}}, "footprint_selection"),
    (("analyses", "lipid_enrichment"), {"enabled": True, "species_selections": {"POPC": "all"}, "footprint_selection": "protein", "shell_edges_A": [0, 5, 5]}, "shell_edges_A"),
    (("analyses", "lipid_enrichment"), {"enabled": True, "species_selections": {"POPC": "all"}, "footprint_selection": "protein", "area_grid_size": 16.5}, "area_grid_size"),
    (("trajectory","stride"), 1.5, "stride"),
    (("trajectory","chunk_frames"), True, "chunk_frames"),
    (("trajectory","start"), .5, "start"),
    (("trajectory","stride"), 0, "stride"),
    (("trajectory","chunk_frames"), 0, "chunk_frames"),
    (("time","dt_ps"), -1., "dt_ps"),
    (("time","origin_ps"), float("inf"), "origin_ps"),
    (("analyses","rmsd","runs"), {"extra": {}}, "runs"),
    (("analyses","rmsd","groups"), ["unknown"], "group IDs"),
    (("analyses","density_1d"), {"enabled":True,"runs":{"../bad":{"selections":{"x":"all"}}}}, "run ID"),
    (("analyses","density_1d"), {"enabled":True,"runs":{"default":{"selections":{"x":"all"}}}}, "reserved"),
    (("analyses","density_1d"), {"enabled":True,"selections":{"x":{"selection":"all","mode":"bogus"}}}, "mode"),
    (("analyses","sasa"), {"enabled":True,"n_points":0}, "n_points"),
    (("analyses","sasa"), {"enabled":True,"probe_radius_A":float("nan")}, "radii"),
    (("analyses","solvation"), {"enabled":True,"preset":"bogus"}, "preset"),
    (("analyses","solvation"), {"enabled":True,"targets":{"protein":{"site":"all"}},"water_selection":"all","ion_selection":"all","ion_selections":{"K":"all"}}, "OR"),
    (("selections","groups"), {"bad":{}}, "missing selection"),
    (("statistics","families"), {"profiles":{"bootstrap":"no"}}, "bootstrap"),
    (("statistics","families"), {"profiles":{"hierarchy":{"enabled":"no"}}}, "hierarchy"),
    (("statistics","hierarchy"), {"resample_replicates":"no"}, "hierarchy"),
    (("statistics","hypothesis_tests","design"), "unrecognized", "design"),
    (("statistics","hypothesis_families"), {"first":["x"],"second":["x"]}, "multiple hypothesis"),
    (("statistics","hypothesis_families"), {1:["x"]}, "named nonempty"),
])
def test_invalid_config_categories_are_rejected_without_mutating_outputs(tmp_path, project, keys, value, diagnostic):
    path, root = project("invalid")
    cfg = yaml.safe_load(path.read_text()); _set(cfg, keys, value)
    path.write_text(yaml.safe_dump(cfg))
    root.mkdir(); sentinel = root/"keep.txt"; sentinel.write_text("previous scientific results")
    before = sentinel.read_bytes()
    with pytest.raises((ValueError,FileNotFoundError), match=re.escape(diagnostic)):
        load_config(path)
    assert sentinel.read_bytes() == before and list(root.iterdir()) == [sentinel]


@pytest.mark.parametrize("content", ["[1, 2]", "analyses: []", "statistics: true", "analyses: {rmsd: invalid}"])
def test_yaml_shape_errors_are_explicit(tmp_path, content):
    path = tmp_path/"bad.yaml"; path.write_text(content)
    with pytest.raises(ValueError,match="mapping|boolean"):
        load_config(path)


@pytest.mark.parametrize("field", ["topology", "trajectory"])
def test_missing_files_and_directories_are_not_accepted_as_inputs(tmp_path, project, field):
    path, _ = project("missing")
    cfg = yaml.safe_load(path.read_text())
    target = cfg["systems"]["tiny"] if field == "topology" else cfg["systems"]["tiny"]["replicates"][0]
    for value in ("missing-file", str(tmp_path)):
        target[field] = value; path.write_text(yaml.safe_dump(cfg))
        with pytest.raises(FileNotFoundError): load_config(path)


def _bad_selection(cfg, case):
    if case == "protein": cfg["selections"]["protein"] = "name NONEXISTENT"
    elif case == "alignment": cfg["selections"]["alignment"] = "index 0"
    elif case == "group": cfg["selections"]["groups"] = {"bad":{"selection":"name NONEXISTENT"}}
    elif case == "species": cfg["analyses"]["bilayer_normal"] = {"enabled":True,"species_selections":{"absent":"resname NONE"}}
    elif case == "overlapping_anchors": cfg["analyses"]["bilayer_normal"] = {"enabled":True,"species_selections":{"one":"protein","two":"protein"}}
    elif case == "sasa_reference": cfg["analyses"]["sasa"] = {"enabled":True,"runs":{
        "isolated":{"target_selection":"index 0:1"},
        "environment":{"target_selection":"index 1:2","occluder_selection":"all","reference_run":"isolated"}}}
    elif case == "sasa_occluder": cfg["analyses"]["sasa"] = {"enabled":True,"target_selection":"all","occluder_selection":"index 0"}


@pytest.mark.parametrize("case,diagnostic", [("protein","zero-atom"),("alignment","alignment"),("group","zero atoms"),
        ("species","matched no atoms"),("overlapping_anchors","overlap"),("sasa_reference","reference_run"),("sasa_occluder","every target")])
def test_selection_preflight_and_run_share_errors_and_preserve_outputs(tmp_path, project, monkeypatch, case, diagnostic):
    freesasa_backend(monkeypatch, tmp_path, available=True)
    path, root = project("selection")
    raw = yaml.safe_load(path.read_text()); _bad_selection(raw,case); path.write_text(yaml.safe_dump(raw))
    cfg = load_config(path)
    root.mkdir(); sentinel = root/"keep.txt"; sentinel.write_bytes(b"unchanged results")
    original = mda.Universe
    def topology_only(*a, **k):
        assert len(a) == 1, "preflight opened trajectory before rejecting invalid selections"
        return original(*a, **k)
    monkeypatch.setattr(mda,"Universe",topology_only)
    with pytest.raises(ValueError,match=diagnostic): validate_selections(cfg,tmp_path)
    with pytest.raises(ValueError,match=diagnostic): run_one(path,0)
    assert list(root.iterdir()) == [sentinel] and sentinel.read_bytes() == b"unchanged results"


@pytest.mark.parametrize("case,diagnostic", [("protein","zero-atom"),("alignment","alignment"),
        ("sasa_reference","reference_run"),("overlapping_anchors","overlap")])
def test_cli_failures_return_nonzero_without_replacing_existing_outputs(tmp_path, project, monkeypatch, case, diagnostic):
    freesasa_backend(monkeypatch, tmp_path, available=True)
    path, root = project("cli_bad")
    cfg = yaml.safe_load(path.read_text()); _bad_selection(cfg,case); path.write_text(yaml.safe_dump(cfg))
    root.mkdir(); sentinel = root/"keep.txt"; sentinel.write_bytes(b"saved report")
    result = invoke(tmp_path,"replicate",path,"--index","0",ok=False)
    assert result.returncode == 2 and result.stderr.startswith("mdscope: error: ")
    assert diagnostic in result.stderr and "Traceback" not in result.stderr
    assert list(root.iterdir()) == [sentinel] and sentinel.read_bytes() == b"saved report"


def test_duplicate_and_unmatched_paired_replica_labels_fail_early(tmp_path, project):
    path, _ = project("pairing"); cfg = yaml.safe_load(path.read_text())
    cfg["systems"]["tiny"]["replicates"] *= 2
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError,match="duplicate replicate"): load_config(path)
    cfg["systems"]["tiny"]["replicates"] = cfg["systems"]["tiny"]["replicates"][:1]
    cfg["systems"]["second"] = deepcopy(cfg["systems"]["tiny"])
    cfg["systems"]["second"]["replicates"][0]["label"] = "different"
    cfg["statistics"]["hypothesis_tests"] = {"enabled":True,"design":"paired"}
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError,match="identical replicate labels"): load_config(path)
    cfg["statistics"]["hypothesis_tests"]["enabled"] = False
    cfg["statistics"]["families"] = {"time_series":{"hypothesis_tests":True}}
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError,match="identical replicate labels"): load_config(path)


def test_overlapping_analysis_groups_are_valid_distinct_outputs(tmp_path, project):
    path, _ = project("groups"); raw = yaml.safe_load(path.read_text())
    raw["selections"]["groups"] = {"whole":{"selection":"protein"},"subset":{"selection":"index 0:2"}}
    path.write_text(yaml.safe_dump(raw)); cfg = load_config(path)
    validate_selections(cfg,tmp_path)


@pytest.mark.parametrize("keys,expected", [
    (("projeckt",), "project"),
    (("analyses", "rmsff"), "rmsf"),
    (("analyses", "rmsd", "enabld"), "enabled"),
    (("analyses", "density_1d", "runs", "new_density", "bin_wdth_A"), "bin_width_A"),
    (("analyses", "sasa", "runs", "environment", "probe_radus_A"), "probe_radius_A"),
    (("analyses", "membrane", "msd", "multi_origin_max_lag_frame"), "multi_origin_max_lag_frames"),
    (("analyses", "membrane_order", "tails", "sn1", "carbon_pattrn"), "carbon_pattern"),
    (("analyses", "solvation", "targets", "protein", "site", "selecton"), "selection"),
    (("analyses", "sdc4_mechanisms", "hs_auto", "glycn_selection"), "glycan_selection"),
    (("statistics", "block", "tau_multiplir"), "tau_multiplier"),
    (("statistics", "adaptive_hierarchy", "min_daws"), "min_draws"),
    (("statistics", "families", "profiles", "bootstrp"), "bootstrap"),
    (("statistics", "families", "time_series", "adaptive_hierarchy", "min_daws"), "min_draws"),
    (("selections", "groups", "A", "alignment_selecton"), "alignment_selection"),
    (("plots", "font_sizes", "bsae"), "base"),
])
def test_unknown_yaml_keys_fail_without_modifying_outputs(tmp_path, project, keys, expected):
    path, root = project("unknown_key")
    cfg = yaml.safe_load(path.read_text())
    _set(cfg, keys, True)
    path.write_text(yaml.safe_dump(cfg))
    root.mkdir()
    sentinel = root / "keep.txt"
    sentinel.write_bytes(b"existing scientific output")
    with pytest.raises(ValueError, match=r"Unknown YAML option .*" + re.escape(expected)):
        load_config(path)
    assert list(root.iterdir()) == [sentinel]
    assert sentinel.read_bytes() == b"existing scientific output"


def test_unknown_replica_file_key_reports_exact_index(tmp_path, project):
    path, _ = project("unknown_rep")
    cfg = yaml.safe_load(path.read_text())
    cfg["systems"]["tiny"]["replicates"][0]["trajectry"] = "mistyped"
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match=r"systems.tiny.replicates\[0\].trajectry.*trajectory"):
        load_config(path)


def test_dynamic_yaml_identifiers_remain_supported():
    from mdscope.config import DEFAULTS
    from mdscope.config_schema import validate_known_keys

    spec = {
        "systems": {"custom_condition": {"topology": "system.tpr",
                    "replicates": [{"trajectory": "run.xtc", "plot_label": "R1"}]}},
        "selections": {"groups": {"chain_C": {"selection": "segid C",
                      "alignment_selection": "segid C and name CA"}}},
        "analyses": {
            "density_1d": {"runs": {"membrane_normal": {
                "selections": {"custom_species": {"selection": "resname X",
                                                  "mode": "number"}}}}},
            "sasa": {"runs": {"environment": {"target_selection": "protein",
                "element_radii_A": {"C": 1.7}}}},
            "membrane_order": {"tails": {"new_tail": {"carbon_pattern": "C[0-9]+"}},
                               "lipid_overrides": {"NEWLIPID": {"new_tail": {
                                   "carbon_names": ["C1", "C2"], "positions": [1, 2]}}}},
            "solvation": {"targets": {"custom_molecule": {"custom_site": {
                "selection": "resname ABC", "group_by": "residue"}}},
                "ion_selections": {"Mg": "resname MG"}},
            "sdc4_mechanisms": {"motifs": {"custom_regex": "G...G"},
                "lipid_anchors": {"new_lipid": "resname NEWLIPID"}},
            "coupling_screen": {"enabled": False, "pairs": [["a", "b"]]},
        },
        "statistics": {"hypothesis_families": {"planned": ["metric"]},
            "families": {"time_series": {"adaptive_hierarchy": {
                "min_draws": 100, "max_draws": 500,
                "batch_size": 100, "stable_checks": 2,
                "relative_tolerance": 0.05}}}},
        "plots": {"families": {"profiles": True}},
    }
    validate_known_keys(spec, DEFAULTS)


@pytest.mark.parametrize("filename", [
    "synd4_baseline.yaml", "synd4_full_mechanisms.yaml",
    "synd4_full_mechanisms_finalize_fast.yaml",
])
def test_maintained_yaml_uses_only_supported_option_names(filename):
    from mdscope.config import DEFAULTS
    from mdscope.config_schema import validate_known_keys

    path = Path(__file__).resolve().parents[1] / filename
    validate_known_keys(yaml.safe_load(path.read_text()), DEFAULTS)
