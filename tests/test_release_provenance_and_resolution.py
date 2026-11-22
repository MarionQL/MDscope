"""Regression checks for #225, #226 and #227 without trajectory reads."""
import subprocess

import pytest


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
    module = tmp_path / "src" / "mdscope" / "core" / "provenance.py"
    module.parent.mkdir(parents=True)
    module.write_text("placeholder")
    monkeypatch.setattr(provenance, "__file__", str(module))
    identity = provenance._source_identity()
    assert identity["git_commit"] is None
    assert identity["git_dirty"] is None
    assert "checkout" in identity["git_unavailable_reason"]


def test_initialization_rejects_changed_config_without_modifying_saved_files(tmp_path, monkeypatch):
    from mdscope.pipeline.runner import initialize
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
    from mdscope.errors import ConfigError
    with pytest.raises(ConfigError, match="config hash"):
        initialize(config_file, changed)
    assert (root / "resolved_config.yaml").read_bytes() == before


def test_cli_reports_changed_project_config_as_user_error(tmp_path, project):
    import yaml
    from mdscope_test_support import invoke_cli
    path, root = project("hash_guard")
    invoke_cli(tmp_path, "replicate", path, "--index", "0")
    saved = (root / "resolved_config.yaml").read_bytes()
    raw = yaml.safe_load(path.read_text())
    raw["trajectory"]["stride"] = 2
    path.write_text(yaml.safe_dump(raw))
    result = invoke_cli(tmp_path, "replicate", path, "--index", "0", ok=False)
    assert result.returncode == 2 and result.stderr.startswith("mdscope: error: ")
    assert "config hash" in result.stderr and "Traceback" not in result.stderr
    assert (root / "resolved_config.yaml").read_bytes() == saved


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


@pytest.mark.parametrize(("design", "n_1", "n_2", "expected"), [
    ("independent", 3, 3, .1),
    ("independent", 4, 4, 2 / 70),
    ("independent", 3, 4, 1 / 35),  # unequal groups: no mirrored assignment
    ("independent", 2, 1, 1 / 3),
    ("paired", 3, 3, 2 / 8),
])
def test_reported_floor_equals_production_permutation_minimum(design, n_1, n_2, expected):
    """Oracle: maximally separated replicas reach the estimator's smallest p."""
    import numpy as np

    from mdscope.statistics.permutation_resolution import permutation_floor
    from mdscope.statistics.common import _permutation
    tests = {"enabled": True, "design": design}
    a = {f"r{i}": float(i) for i in range(n_1)}
    b = ({key: value + 100. + 10 * i for i, (key, value) in enumerate(a.items())}
         if design == "paired" else {f"r{i}": 100. + i for i in range(n_2)})
    observed = _permutation(a, b, {"statistics": {"hypothesis_tests": tests}},
                            np.random.default_rng(0))
    assert observed == pytest.approx(expected)
    matched = len(a.keys() & b.keys())
    assert permutation_floor(n_1, n_2, n_matched=matched, tests=tests) == (
        pytest.approx(expected), "exact")


def test_packaged_logo_gracefully_handles_missing_or_unreadable_assets(tmp_path):
    from mdscope.studio.app import _studio_logo_markup
    logo = tmp_path / "logo.png"
    assert _studio_logo_markup(logo) == ""
    logo.mkdir()
    assert _studio_logo_markup(logo) == ""
    logo.rmdir()
    logo.write_bytes(b"test-image-bytes")
    assert 'data:image/png;base64,' in _studio_logo_markup(logo)


def test_source_identity_dirty_checkout_and_foreign_cwd(monkeypatch, tmp_path):
    from mdscope.core import provenance
    checkout = tmp_path / "checkout"
    module = checkout / "src/mdscope/core/provenance.py"
    module.parent.mkdir(parents=True)
    module.write_text("test")
    subprocess.run(["git", "init", "-q", str(checkout)], check=True)
    subprocess.run(["git", "-C", str(checkout), "config", "user.email", "test@example.org"], check=True)
    subprocess.run(["git", "-C", str(checkout), "config", "user.name", "Test"], check=True)
    subprocess.run(["git", "-C", str(checkout), "add", "."], check=True)
    subprocess.run(["git", "-C", str(checkout), "commit", "-qm", "initial"], check=True)
    expected = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    foreign = tmp_path / "other"
    foreign.mkdir()
    subprocess.run(["git", "init", "-q", str(foreign)], check=True)
    monkeypatch.chdir(foreign)
    monkeypatch.setattr(provenance, "__file__", str(module))
    assert provenance._source_identity()["git_commit"] == expected
    assert provenance._source_identity()["git_dirty"] is False
    module.write_text("changed")
    assert provenance._source_identity()["git_dirty"] is True


def test_findings_support_counts_only_replicas_with_metric_rows(tmp_path):
    import pandas as pd

    from mdscope.reporting.report import _replicate_support, _support_floor
    from mdscope.statistics.coverage import _identity
    path = tmp_path / "time_series"
    path.mkdir()
    table = pd.DataFrame([
        {"system": "a", "replicate": "r1", "group": "", "metric": "rmsd"},
        {"system": "a", "replicate": "r2", "group": "", "metric": "rmsd"},
        {"system": "b", "replicate": "r1", "group": "", "metric": "rmsd"},
        {"system": "a", "replicate": "r3", "group": "", "metric": "rg"},
    ])
    table.to_parquet(path / "replicate_statistics.parquet")
    identity = _identity(table.iloc[:1], ("group", "metric"))[0]
    assert _replicate_support(tmp_path, lambda family: family == "time_series") == {}
    support = _replicate_support(tmp_path, lambda family: False)
    assert support["matched_time_series", identity, "a"] == {"r1", "r2"}
    assert ("matched_time_series", identity, "a") not in _replicate_support(
        tmp_path, lambda family: family == "matched_time_series")
    a = support["time_series", identity, "a"]
    b = support["time_series", identity, "b"]
    assert a == {"r1", "r2"}
    assert b == {"r1"}
    minimum, kind = _support_floor(a, b, {"hypothesis_tests": {"design": "independent"}})
    assert minimum == pytest.approx(1 / 3)  # 2 vs 1 unequal groups; see oracle test
    assert kind == "exact"
