"""Generate synthetic MDscope projects: the validation oracle and the quickstart.

The molecular geometry is specified in docs/developer/sdc4-fixtures.md.
No production trajectories or force-field resources are required. Outputs are
software fixtures, not molecular dynamics evidence.
"""
from copy import deepcopy
from pathlib import Path
import argparse
import json
import shutil

import MDAnalysis as mda
import numpy as np
import pandas as pd
import yaml

from mdscope.config import DEFAULTS, Replicate
from test_sdc4_geometry_references import _fixture


def _write_topology(pdb):
    u, partial = _fixture()
    u.atoms.chainIDs = ["H" if c == "HS" else c for c in u.atoms.chainIDs]
    # Adjacent identical PDB residue IDs are one residue on re-read.
    u.residues.resids = np.arange(1, len(u.residues)+1)
    u.atoms.write(str(pdb))
    return u, partial


def make_project(tmp_path, chunk):
    pdb, dcd = tmp_path/"mixed.pdb", tmp_path/"mixed.dcd"
    u, partial = _write_topology(pdb)
    initial = u.atoms.positions.copy()
    with mda.Writer(str(dcd), n_atoms=len(u.atoms)) as writer:
        for i in range(8):
            xyz = initial.copy()
            xyz[:12, 0] += .03*i
            u.atoms.positions = xyz
            writer.write(u.atoms)
    cfg = deepcopy(DEFAULTS)
    for spec in cfg["analyses"].values():
        spec["enabled"] = False
    for name in ("alignment", "rmsd", "rmsf", "radius_of_gyration"):
        cfg["analyses"][name]["enabled"] = True
    cfg["selections"]["groups"] = {c: {"selection": f"chainID {c} and name CA"} for c in ("A", "B")}
    cfg["analyses"]["sdc4_mechanisms"] = partial["analyses"]["sdc4_mechanisms"]
    cfg["analyses"]["sdc4_mechanisms"]["hs_chains"] = {"A1": "chainID H"}
    for name in ("bilayer_normal", "lipid_enrichment"):
        cfg["analyses"][name].update(enabled=True, species_selections={"POPC": "resname POPC and name P"})
    cfg["analyses"]["lipid_enrichment"].update(footprint_selection="resindex 5:11 or resindex 17:23", slab_half_width_A=12., shell_edges_A=[0, 20])
    cfg["analyses"]["density_1d"] = {"enabled": True, "runs": {
        "protein": {"selections": {"protein": "protein"}, "axis": "z", "bin_width_A": 5.},
        "lipid": {"selections": {"lipid": "resname POPC"}, "axis": "z", "bin_width_A": 5.},
    }}
    cfg["analyses"]["residue_membrane_distance"] = {"enabled": True, "runs": {
        c: {"selection": f"chainID {c} and name CA", "membrane_selection": "resname POPC"}
        for c in ("A", "B")
    }}
    cfg["analyses"]["solvation"] = {
        "enabled": True, "preset": "kinetics",
        "targets": {"protein": {"backbone_CA": "protein and name CA"}},
        "water_selection": "resname TIP3 and name OH2",
        "search_radius_A": 6., "rdf_bin_width_A": 1.,
    }
    cfg["trajectory"].update(start=0, stop=8, stride=2, chunk_frames=chunk)
    cfg["cache"]["aligned_coordinates"] = False
    cfg["statistics"].update(bootstrap_draws=100, families={"profiles": False, "interactions": False, "solvation_rdfs": False})
    cfg["plots"]["enabled"] = False
    cfg["systems"] = {"test": {"topology": str(pdb), "replicates": [{"trajectory": str(dcd), "label": "rep1"}]}}
    return Replicate("test", 1, pdb, dcd, "rep1"), cfg



def write_project(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    for name in ("mixed.pdb", "mixed.dcd", "validation.yaml", "reference.json", "results"):
        if (directory/name).exists():
            raise FileExistsError(f"Use a fresh validation directory: {directory/name}")
    _, cfg = make_project(directory, 3)
    cfg["project"].update(name="synthetic_validation", output="results")
    cfg["systems"]["test"]["topology"] = "mixed.pdb"
    cfg["systems"]["test"]["replicates"][0]["trajectory"] = "mixed.dcd"
    cfg["progress"]["enabled"] = False
    cfg["analyses"]["bilayer_normal"]["neighbor_radius_A"] = 30.
    (directory/"validation.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    shutil.copyfile(Path(__file__).parent/"data/validation_reference.json", directory/"reference.json")
    return directory/"validation.yaml"


QUICKSTART_CONDITIONS = {"synthetic_reference": 0., "synthetic_displaced": 1.}


def write_quickstart(directory: Path):
    """Two synthetic conditions x three replicas for the documented quickstart.

    Chain B is translated along x away from chain A by the condition offset
    (Angstrom) plus seeded per-frame jitter in [0, 0.2) Angstrom. The displaced
    condition therefore has a larger radius of gyration by construction.
    """
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise FileExistsError(f"Use a fresh quickstart directory: {directory}")
    u, _ = _write_topology(directory/"synthetic.pdb")
    initial = u.atoms.positions.copy()
    rng = np.random.default_rng(2026)
    systems = {}
    for condition, offset in QUICKSTART_CONDITIONS.items():
        replicates = []
        for replica in (1, 2, 3):
            name = f"{condition}_rep{replica}.dcd"
            with mda.Writer(str(directory/name), n_atoms=len(u.atoms)) as writer:
                for jitter in rng.uniform(0., .2, size=20):
                    xyz = initial.copy()
                    xyz[:12, 0] += offset + jitter
                    u.atoms.positions = xyz
                    writer.write(u.atoms)
            replicates.append({"trajectory": name, "label": f"rep{replica}"})
        systems[condition] = {"topology": "synthetic.pdb", "replicates": replicates}
    analyses = {name: {"enabled": False} for name in DEFAULTS["analyses"]}
    for name in ("alignment", "rmsd", "rmsf", "radius_of_gyration"):
        analyses[name]["enabled"] = True
    cfg = {
        "project": {"name": "synthetic_quickstart_fixture", "output": "results"},
        "progress": {"enabled": False},
        "systems": systems,
        "analyses": analyses,
        "plots": {"enabled": False},
    }
    path = directory/"quickstart.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def verify_project(directory: Path):
    reference = json.loads((directory/"reference.json").read_text())
    root = directory/"results/test/rep1"
    for relative, spec in reference["tables"].items():
        table = pd.read_parquet(root/relative)
        if len(table) != spec["rows"]:
            raise AssertionError(f"{relative}: row count {len(table)} != {spec['rows']}")
        for column, expected in spec.get("columns", {}).items():
            np.testing.assert_allclose(table[column], expected, rtol=0,
                atol=spec["atol"], err_msg=f"{relative}: {column}")
        for column, expected in spec.get("sums", {}).items():
            np.testing.assert_allclose(table[column].sum(), expected, rtol=0,
                atol=spec["atol"], err_msg=f"{relative}: sum({column})")
    for relative, spec in reference["arrays"].items():
        np.testing.assert_allclose(np.load(root/relative, allow_pickle=False), spec["values"],
            rtol=0, atol=spec["atol"], err_msg=relative)
    return len(reference["tables"]) + len(reference["arrays"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("create", "verify", "quickstart"))
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    if args.action == "create":
        print(write_project(args.directory.resolve()))
    elif args.action == "quickstart":
        print(write_quickstart(args.directory.resolve()))
    else:
        print(f"Verified {verify_project(args.directory.resolve())} independent reference products")
