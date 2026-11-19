"""Failure-injection checks for actionable, non-destructive stage boundaries (#183)."""

from types import SimpleNamespace

import MDAnalysis as mda
import pytest

from mdscope.config import load_config
from mdscope.core.validation import validate_selections
from mdscope.pipeline.runner import _validate_replica_products
from mdscope.pipeline.streaming import validate_stream_inventory


def test_topology_load_failure_reports_system_path_and_original_cause(
    tmp_path, project, monkeypatch
):
    config, root = project("boundary_topology")
    cfg = load_config(config)
    def invalid_topology(*args, **kwargs):
        raise ValueError("unsupported topology encoding")

    monkeypatch.setattr(mda, "Universe", invalid_topology)
    with pytest.raises(ValueError, match="System 'tiny'.*topology") as failure:
        validate_selections(cfg, config.parent)
    assert "unsupported topology encoding" in str(failure.value)
    assert str(config.parent) in str(failure.value)
    assert isinstance(failure.value.__cause__, ValueError)
    assert not root.exists()


def test_missing_replica_summary_error_names_file_and_preserves_outputs(tmp_path):
    from mdscope.config import DEFAULTS
    root = tmp_path / "saved"
    root.mkdir()
    sentinel = root / "keep.txt"
    sentinel.write_bytes(b"saved findings")
    rep = SimpleNamespace(system="example", label="rep1")
    with pytest.raises(RuntimeError, match="example/rep1") as failure:
        _validate_replica_products(DEFAULTS, root, [rep])
    assert "summary/timeseries.parquet" in str(failure.value)
    assert "finish the replica" in str(failure.value)
    assert sentinel.read_bytes() == b"saved findings"
    assert sorted(root.iterdir()) == [sentinel]


def test_corrupt_checkpoint_reports_replica_and_preserves_raw_bytes(tmp_path):
    root = tmp_path / "saved"
    repout = root / "example" / "rep1"
    checkpoint = repout / "cache" / "stream_checkpoint.npz"
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_bytes(b"corrupted checkpoint")
    sentinel = repout / "summary" / "timeseries.parquet"
    sentinel.parent.mkdir()
    sentinel.write_bytes(b"previous results")
    rep = SimpleNamespace(system="example", label="rep1")
    with pytest.raises(RuntimeError, match="example/rep1") as failure:
        validate_stream_inventory(rep, repout)
    assert str(checkpoint) in str(failure.value)
    assert isinstance(failure.value.__cause__, (OSError, ValueError, EOFError))
    assert checkpoint.read_bytes() == b"corrupted checkpoint"
    assert sentinel.read_bytes() == b"previous results"


@pytest.mark.parametrize("manifest_content,expected", [
    (None, "missing"),
    ("{invalid", "invalid"),
    ('{"status":"running"}', "running"),
])
def test_statistical_plot_manifest_failure_preserves_existing_figures(
    tmp_path, project, manifest_content, expected
):
    import yaml
    from mdscope.pipeline.runner import finalize_project

    config, root = project("plot_manifest_diagnostic")
    raw = yaml.safe_load(config.read_text())
    raw["plots"] = {
        "enabled": True,
        "include": ["time_series:rmsd_A:bootstrap_table"],
    }
    config.write_text(yaml.safe_dump(raw))
    root.mkdir()
    existing = root / "project_summary" / "figures" / "old.png"
    existing.parent.mkdir(parents=True)
    existing.write_bytes(b"prior figure")
    manifest = root / "project_summary" / "statistics" / "manifest.json"
    if manifest_content is not None:
        manifest.parent.mkdir(parents=True)
        manifest.write_text(manifest_content)
    with pytest.raises(RuntimeError, match="completed statistics") as failure:
        finalize_project(
            config, cfg=load_config(config, check_paths=False),
            root=root, phase="plotting"
        )
    assert expected in str(failure.value)
    assert str(manifest) in str(failure.value)
    assert existing.read_bytes() == b"prior figure"
