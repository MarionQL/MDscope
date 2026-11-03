"""Actionable config/topology failures precede output creation or replacement."""
from copy import deepcopy
import re
import MDAnalysis as mda
import pytest
import yaml
from mdscope.config import load_config
from mdscope.core.validation import validate_selections
from mdscope.pipeline.runner import run_one
from test_cli_end_to_end import project, invoke


def _set(cfg, keys, value):
    target = cfg
    for key in keys[:-1]: target = target.setdefault(key, {})
    target[keys[-1]] = value


@pytest.mark.parametrize("keys,value,diagnostic", [
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


@pytest.mark.parametrize("case", ["protein", "alignment", "sasa_reference", "overlapping_anchors"])
def test_cli_failures_return_nonzero_without_replacing_existing_outputs(tmp_path, project, case):
    path, root = project("cli_bad")
    cfg = yaml.safe_load(path.read_text()); _bad_selection(cfg,case); path.write_text(yaml.safe_dump(cfg))
    root.mkdir(); sentinel = root/"keep.txt"; sentinel.write_bytes(b"saved report")
    result = invoke(tmp_path,"replicate",path,"--index","0",ok=False)
    assert "ValueError" in result.stderr
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
