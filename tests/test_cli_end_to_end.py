"""Exercise the installed CLI across a real (tiny) PDB/DCD project.

These tests deliberately spawn the console script rather than calling its Python
implementations: argparse, packaging, configuration paths and exit codes matter.
"""
from __future__ import annotations


import numpy as np
import pandas as pd
import pytest
import yaml

from mdscope_test_support import invoke_cli as invoke


def assert_summary(root):
    table = root / "tiny" / "rep1" / "summary" / "timeseries.parquet"
    assert table.is_file(), f"Missing replicate timeseries: {table}"
    df = pd.read_parquet(table)
    assert len(df) == 6
    assert "rmsd_A" in df
    assert np.isfinite(df["rmsd_A"]).all()
    assert (root / "project_summary").is_dir()


def test_init_and_validate_cli(tmp_path, project):
    initialized = tmp_path / "initialized.yaml"
    invoke(
        tmp_path, "init", "--output", initialized,
        stdin="example\ntiny.pdb\ntiny.dcd\n\n",
    )
    generated = yaml.safe_load(initialized.read_text())
    assert generated["project"]["name"] == "example"
    assert generated["systems"]["system1"]["replicates"][0]["trajectory"] == "tiny.dcd"
    cfg, _ = project("validation")
    result = invoke(tmp_path, "validate", cfg)
    assert "Valid config" in result.stdout
    assert "protein=4" in result.stdout


def test_replicate_finalize_and_report_cli(tmp_path, project):
    cfg, root = project("separate", statistics=True)
    invoke(tmp_path, "replicate", cfg, "--index", "0")
    assert (root / "tiny" / "rep1" / "summary" / "timeseries.parquet").is_file()
    invoke(tmp_path, "finalize", cfg)
    assert_summary(root)
    statistics_dir = root / "project_summary" / "statistics"
    findings = statistics_dir / "findings.md"
    content = findings.read_bytes()
    assert b"MDscope findings" in content
    assert not list(root.rglob("*.png"))
    invoke(tmp_path, "report", statistics_dir)
    assert findings.read_bytes() == content
    resolved = (root / "resolved_config.yaml").read_bytes()
    options = yaml.safe_load(cfg.read_text())
    options["plots"] = {"enabled": True, "families": {"profiles": False, "interactions": False}}
    cfg.write_text(yaml.safe_dump(options))
    invoke(tmp_path, "plotting", cfg)
    assert list((root / "project_summary" / "figures").rglob("*__timeseries.png"))
    assert (root / "resolved_config.yaml").read_bytes() == resolved
    assert findings.read_bytes() == content


def test_run_cli(tmp_path, project):
    cfg, root = project("all_at_once")
    invoke(tmp_path, "run", cfg)
    assert_summary(root)


@pytest.mark.parametrize("args", [
    ("replicate", "missing.yaml", "--index", "0"),
    ("validate", "missing.yaml"),
    ("finalize", "missing.yaml"),
    ("report", "missing.yaml"),
])
def test_missing_config_exits_nonzero(tmp_path, args):
    invoke(tmp_path, *args, ok=False)


def test_invalid_replicate_index_and_incomplete_finalization(tmp_path, project):
    cfg, root = project("failures")
    invoke(tmp_path, "replicate", cfg, "--index", "99", ok=False)
    invoke(tmp_path, "finalize", cfg, ok=False)
    assert not (root / "project_summary").exists()
