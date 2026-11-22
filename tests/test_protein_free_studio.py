"""Independent protein-free topology/production entry-point acceptance cases."""
import MDAnalysis as mda
import pandas as pd
import pytest
import yaml

from mdscope.config import load_config
from mdscope.pipeline.preflight import validate_selections
from mdscope.pipeline.runner import run_one
from mdscope.studio.model import new_document
from mdscope_test_support import freesasa_backend, membrane_universe


def lipid_project(tmp_path):
    u = membrane_universe()
    lipids = u.select_atoms("resname POPC")
    top = tmp_path / "lipids.pdb"
    traj = tmp_path / "lipids.xtc"
    lipids.write(str(top))
    with mda.Writer(str(traj), n_atoms=len(lipids)) as writer:
        for _ in u.trajectory:
            writer.write(lipids)
    doc = new_document()
    doc["systems"] = {"bilayer": {"topology": str(top), "replicates": [{"label": "rep1", "trajectory": str(traj)}]}}
    doc["project"]["output"] = str(tmp_path / "results")
    doc["time"] = {"dt_ps": 1.0}
    # Unused invalid defaults must never even reach the selection parser.
    doc["selections"].update(protein="invalid syntax", features="invalid syntax", alignment="invalid syntax")
    doc["selections"]["groups"] = {"lipids": {"selection": "resname POPC"}}
    doc["analyses"]["membrane"].update(enabled=True, headgroup_selection="name P")
    doc["analyses"]["density_1d"] = {"enabled": True, "selections": {"heads": "name P"}}
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(doc))
    return path, doc


def test_protein_free_preflight_and_production_share_requirements(tmp_path):
    path, _ = lipid_project(tmp_path)
    cfg = load_config(path)
    counts = validate_selections(cfg, tmp_path)["bilayer"]
    assert counts["required_by"]["protein"] == []
    assert counts["protein"] is None
    assert counts["details"]["protein"] == {"status": "not_required"}
    assert counts["groups"]["lipids"]["group_atoms"] == 8
    run_one(path, 0)
    output = tmp_path / "results/bilayer/rep1"
    series = pd.read_parquet(output / "summary/timeseries.parquet")
    # Independent geometry: 100x100 area, four anchors in each leaflet.
    assert series.area_per_lipid_A2.tolist() == pytest.approx([2500.0] * 5)
    assert series.membrane_thickness_A.tolist() == pytest.approx([40.0] * 5)
    assert len(pd.read_parquet(output / "profiles/density_1d/default/heads.parquet")) > 0


@pytest.mark.parametrize("analysis", ["rmsd", "rmsf", "alignment"])
def test_required_structural_selection_still_rejects_bad_defaults(tmp_path, analysis):
    path, doc = lipid_project(tmp_path)
    doc["analyses"][analysis]["enabled"] = True
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match="selection"):
        validate_selections(load_config(path), tmp_path)
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("analysis", ["sasa", "radius_of_gyration", "rmsd"])
def test_required_default_protein_cannot_be_empty(tmp_path, monkeypatch, analysis):
    freesasa_backend(monkeypatch, tmp_path, available=True)
    path, doc = lipid_project(tmp_path)
    doc["selections"].update(protein="protein", features="all")
    doc["analyses"][analysis]["enabled"] = True
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match="zero-atom protein"):
        validate_selections(load_config(path), tmp_path)
    assert not (tmp_path / "results").exists()


def test_required_alignment_cannot_use_insufficient_nonprotein_atoms(tmp_path):
    path, doc = lipid_project(tmp_path)
    doc["selections"].update(features="all", alignment="index 0 1")
    doc["analyses"]["alignment"]["enabled"] = True
    path.write_text(yaml.safe_dump(doc))
    with pytest.raises(ValueError, match="at least three atoms"):
        validate_selections(load_config(path), tmp_path)


def test_nonprotein_solvation_pipeline_counts_waters_without_structural_defaults(tmp_path):
    from mdscope_test_support import solvation_universe
    u = solvation_universe()
    u.residues[0].resname = "LIG"
    top, traj = tmp_path / "solvent.pdb", tmp_path / "solvent.xtc"
    u.atoms.write(str(top))
    with mda.Writer(str(traj), n_atoms=len(u.atoms)) as writer:
        for _ in u.trajectory:
            writer.write(u)
    doc = new_document()
    doc["systems"] = {"ligand": {"topology": str(top), "replicates": [{"label": "rep1", "trajectory": str(traj)}]}}
    doc["project"]["output"] = str(tmp_path / "results")
    doc["selections"].update(protein="invalid syntax", features="invalid syntax")
    doc["analyses"]["solvation"].update(enabled=True, water_selection="resname SOL",
        targets={"ligand": {"oxygen": "name OE1 OE2"}})
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(doc))
    counts = validate_selections(load_config(path), tmp_path)["ligand"]
    assert counts["protein"] is None
    run_one(path, 0)
    saved = pd.read_parquet(tmp_path / "results/ligand/rep1/extensions/solvation/counts.parquet")
    assert saved.sort_values("time_ps").n_solvent.tolist() == [2, 1]
