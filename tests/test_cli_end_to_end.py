"""Exercise the installed CLI across a real (tiny) PDB/DCD project.

These tests deliberately spawn the console script rather than calling its Python
implementations: argparse, packaging, configuration paths and exit codes matter.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
import yaml


def invoke(tmp_path: Path, *args: str, stdin: str | None = None, ok: bool = True):
    command = ["mdscope", *map(str, args)]
    result = subprocess.run(
        command, input=stdin, text=True, capture_output=True,
        cwd=tmp_path, timeout=150, check=False,
    )
    if ok:
        assert result.returncode == 0, (
            f"{command!r} failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    else:
        assert result.returncode != 0, f"{command!r} unexpectedly succeeded"
    return result


@pytest.fixture
def project(tmp_path):
    pdb = tmp_path / "tiny.pdb"
    xyz = np.array([
        [10., 10., 10.], [13., 10., 10.],
        [10., 13., 10.], [10., 10., 13.],
    ])
    lines = ["CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1"]
    for i, (x, y, z) in enumerate(xyz, 1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}           C"
        )
    pdb.write_text("\n".join([*lines, "TER", "END", ""]))
    universe = mda.Universe(str(pdb))
    trajectory = tmp_path / "tiny.dcd"
    with mda.Writer(str(trajectory), n_atoms=len(universe.atoms)) as writer:
        for frame in range(6):
            positions = xyz.copy()
            positions[-1, 2] += frame * 0.12
            universe.atoms.positions = positions
            universe.dimensions = np.array([50., 50., 50., 90., 90., 90.])
            writer.write(universe.atoms)

    def config(label="results", statistics=False):
        analyses = {
            name: {"enabled": False}
            for name in (
                "membrane", "membrane_order", "dssp", "sasa",
                "salt_bridges", "hydrogen_bonds", "ligand_contacts",
                "density_1d", "residue_membrane_distance",
                "bilayer_normal", "lipid_enrichment",
                "lipid_hydrogen_bonds", "sdc4_mechanisms",
                "lifetimes",
            )
        }
        analyses.update({
            "alignment": {"enabled": True},
            "rmsd": {"enabled": True},
            "rmsf": {"enabled": True},
            "radius_of_gyration": {"enabled": True},
        })
        specification = {
            "project": {"name": "cli_smoke", "output": label},
            "trajectory": {"stride": 1, "chunk_frames": 3},
            "selections": {
                "protein": "protein", "alignment": "protein and name CA",
                "features": "protein and name CA", "ligand": None,
            },
            "systems": {
                "tiny": {
                    "topology": pdb.name,
                    "replicates": [{"trajectory": trajectory.name, "label": "rep1"}],
                },
            },
            "analyses": analyses,
            "statistics": {"enabled": statistics, "bootstrap_draws": 100},
            "plots": {"enabled": False},
        }
        path = tmp_path / f"{label}.yaml"
        path.write_text(yaml.safe_dump(specification, sort_keys=False))
        return path, tmp_path / label
    return config


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
