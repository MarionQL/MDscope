"""Regression checks for #225, #226 and #227 without trajectory reads."""
from pathlib import Path
import json
import subprocess

import pytest
import yaml


def test_source_identity_does_not_use_launch_directory(tmp_path, monkeypatch):
    from mdscope.core.provenance import _source_identity
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    subprocess.run(["git", "init", "-q", str(foreign)], check=True)
    monkeypatch.chdir(foreign)
    identity = _source_identity()
    assert identity["git_commit"] is None or len(identity["git_commit"]) == 40
    # A foreign repo has no commit, and must never be mistaken for MDscope.
    assert identity["git_commit"] != "foreign"


def test_source_identity_without_checkout(tmp_path, monkeypatch):
    from mdscope.core import provenance
    original = provenance.__file__
    module = tmp_path / "src" / "mdscope" / "core" / "provenance.py"
    module.parent.mkdir(parents=True)
    module.write_text("placeholder")
    monkeypatch.setattr(provenance, "__file__", str(module))
    identity = provenance._source_identity()
    assert identity["git_commit"] is None
    assert identity["git_dirty"] is None
    assert "checkout" in identity["git_unavailable_reason"]
    monkeypatch.setattr(provenance, "__file__", original)


def test_initialization_rejects_changed_config_without_modifying_saved_files(tmp_path, monkeypatch):
    from mdscope.pipeline.runner import initialize
    from mdscope.core import provenance
    cfg = {"project": {"output": str(tmp_path / "out")}, "systems": {"a": {}}}
    config_file = tmp_path / "analysis.yaml"
    config_file.write_text("project: {}")
    monkeypatch.setattr("mdscope.pipeline.runner.environment_metadata", lambda: {"git_commit": None})
    root, _ = initialize(config_file, cfg)
    before = (root / "resolved_config.yaml").read_bytes()
    metadata = (root / "run_metadata.json").read_bytes()
    initialize(config_file, cfg)
    assert (root / "resolved_config.yaml").read_bytes() == before
    assert (root / "run_metadata.json").read_bytes() == metadata
    changed = {**cfg, "systems": {"b": {}}}
    with pytest.raises(ValueError, match="config hash"):
        initialize(config_file, changed)
    assert (root / "resolved_config.yaml").read_bytes() == before


@pytest.mark.parametrize(("counts", "expected"), [
    ({"a": 3, "b": 3}, .1),
    ({"a": 4, "b": 4}, 2 / 70),
])
def test_independent_exact_resolution(counts, expected):
    from mdscope.reporting.report import _design_resolution
    manifest = {"methods": {"hypothesis_tests": {"design": "independent"}},
                "replicate_paths": [f"{system}/{index}" for system, n in counts.items()
                                    for index in range(n)]}
    _, floor = _design_resolution(manifest)
    actual, kind = floor("a", "b")
    assert kind == "exact"
    assert actual == pytest.approx(expected)


def test_paired_resolution_uses_intersection_not_smallest_count():
    from mdscope.reporting.report import _design_resolution
    manifest = {"methods": {"hypothesis_tests": {"design": "paired"}},
                "replicate_paths": ["a/r1", "a/r2", "a/r3", "b/r1", "b/r2", "b/r3", "b/other"]}
    _, floor = _design_resolution(manifest)
    assert floor("a", "b") == (.25, "exact")


def test_monte_carlo_resolution_is_explicitly_not_exact():
    from mdscope.reporting.report import _design_resolution
    manifest = {"methods": {"hypothesis_tests": {
        "design": "independent", "exact_max_assignments": 10, "permutation_draws": 99}},
        "replicate_paths": [f"{system}/{index}" for system in ("a", "b") for index in range(3)]}
    _, floor = _design_resolution(manifest)
    assert floor("a", "b") == (.01, "monte_carlo_resolution")
