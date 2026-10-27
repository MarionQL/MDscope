"""Opt-in smoke test using published MDscope release trajectories.

Run only on the dedicated real-data GitHub Actions workflow. This is a
runtime/selection test, not a physical validation or convergence analysis.
"""
from __future__ import annotations

import gc
import json
import os
import traceback
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import yaml

from mdscope.analysis.runner import run_one

DATA = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data")).resolve()
ROOT = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve()
ROOT.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "carb_sdc4": (
        "carb-sdc4-step7_production.tpr",
        "carb-sdc4-rep1_test.xtc",
    ),
    "nocarb_sdc4": (
        "nocarb-sdc4-step7_production.tpr",
        "nocarb-sdc4-rep1_test.xtc",
    ),
}


def inspect(system: str, topology: Path, trajectory: Path) -> dict:
    print(f"\n=== {system}: topology/trajectory inspection ===", flush=True)
    u = mda.Universe(str(topology), str(trajectory))
    first = u.trajectory[0]
    print(f"atoms={len(u.atoms)} residues={len(u.residues)} frames={len(u.trajectory)} "
          f"first_time_ps={first.time} box={first.dimensions.tolist()}", flush=True)
    assert len(u.atoms) > 0 and len(u.trajectory) > 0
    assert np.all(np.isfinite(first.positions))
    assert np.all(np.asarray(first.dimensions[:3]) > 0)
    protein = u.select_atoms("protein")
    ca = u.select_atoms("protein and name CA")
    print(f"protein_atoms={len(protein)} protein_CA={len(ca)}", flush=True)
    assert len(protein) and len(ca), "MDAnalysis protein/CA selection is empty"

    species = {}
    for name, atom_name in (
        ("POPC", "P"), ("POPE", "P"), ("POPI", "P"),
        ("DOPC", "P"), ("DPPC", "P"), ("PSM", "P"),
        ("CHL1", "O3"),
    ):
        selection = f"resname {name} and name {atom_name}"
        n = len(u.select_atoms(selection))
        if n:
            species[name] = selection
        print(f"lipid_anchor {name}/{atom_name}={n}", flush=True)
    for name in ("GM1a", "BGAL", "CER160", "CER1", "TIP3"):
        n = len(u.select_atoms(f"resname {name}"))
        print(f"resname {name} atoms={n}", flush=True)
    for name in ("GM1a", "BGAL", "CER160"):
        n = len(u.select_atoms(f"segid {name}"))
        print(f"segid {name} atoms={n}", flush=True)
    for name in ("GM1a", "BGAL", "CER160"):
        try:
            n = len(u.select_atoms(f"chainID {name}"))
            print(f"chainID {name} atoms={n}", flush=True)
        except Exception as exc:
            print(f"chainID {name} unavailable: {type(exc).__name__}", flush=True)
    print("non-water residue counts (top 25):", flush=True)
    from collections import Counter
    nonwater = Counter(str(r.resname) for r in u.residues if str(r.resname) != "TIP3")
    print(nonwater.most_common(25), flush=True)
    assert len(species) >= 2, "Insufficient lipid anchors for a bilayer test"
    del u
    gc.collect()
    return species


def build_config(system: str, topology: Path, trajectory: Path, species: dict) -> Path:
    # Restrict all analyses to the first eight frames. We are checking that
    # real TPR+XTC inputs execute, NOT estimating equilibrated observables.
    lipid_types = [name for name in ("POPC", "POPE", "DOPC", "DPPC", "POPI", "PSM") if name in species]
    assert lipid_types
    lipid_sel = "resname " + " ".join(lipid_types[:2])
    cfg = {
        "project": {"name": f"real_smoke_{system}", "output": str(ROOT / system)},
        "trajectory": {"start": 0, "stop": 8, "stride": 1, "chunk_frames": 4},
        "progress": {"enabled": True},
        "systems": {system: {"topology": str(topology),
                             "replicates": [{"trajectory": str(trajectory), "label": "rep1"}]}},
        "selections": {
            "protein": "protein", "alignment": "protein and name CA",
            "features": "protein and name CA", "ligand": None,
        },
        "analyses": {
            "alignment": {"enabled": True},
            "rmsd": {"enabled": True},
            "rmsf": {"enabled": False},
            "radius_of_gyration": {"enabled": True, "selection": "protein and name CA"},
            "density_1d": {"enabled": False},
            "membrane": {"enabled": False},
            "membrane_order": {
                "enabled": True, "axis": "z", "lipid_selection": lipid_sel,
                "per_lipid_type": True, "per_leaflet": True,
                "leaflet_anchor_selection": f"({lipid_sel}) and name P",
            },
            "bilayer_normal": {
                "enabled": True, "axis": "z", "species_selections": species,
                "neighbor_radius_A": 15.0, "min_neighbors": 4,
                "save_lipid_angles": False,
            },
            "lipid_enrichment": {
                "enabled": True, "axis": "z", "species_selections": species,
                "protein_selection": "protein and name CA",
                "inner_radius_A": 0.0, "outer_radius_A": 25.0,
            },
            "lipid_hydrogen_bonds": {
                "enabled": True, "lipid_selection": lipid_sel,
                "exclude_same_residue": True,
            },
            "residue_membrane_distance": {"enabled": False},
            "dssp": {"enabled": False},
            "sasa": {"enabled": False},
            "salt_bridges": {"enabled": False},
            "hydrogen_bonds": {"enabled": False},
            "ligand_contacts": {"enabled": False},
            "lifetimes": {"enabled": False},
        },
        "statistics": {"enabled": False},
        "plots": {"enabled": False},
        "cache": {"aligned_coordinates": False},
    }
    path = ROOT / f"{system}.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def check_outputs(system: str):
    rep = ROOT / system / system / "rep1"
    paths = {
        "streaming": rep / "summary" / "timeseries.parquet",
        "order": rep / "profiles" / "membrane_order" / "timeseries.parquet",
        "normal": rep / "extensions" / "bilayer_normal" / "timeseries.parquet",
        "enrichment": rep / "extensions" / "lipid_enrichment" / "timeseries.parquet",
        "enrichment_counts": rep / "extensions" / "lipid_enrichment" / "counts.parquet",
        "lipid_hbonds": rep / "interactions" / "lipid_hydrogen_bonds" / "count_timeseries.parquet",
    }
    result = {}
    for name, path in paths.items():
        assert path.exists(), f"{system}: missing {name} output {path}"
        if path.is_file():
            df = pd.read_parquet(path)
            assert not df.empty, f"{system}: empty {name} output"
            result[name] = {"rows": len(df), "columns": len(df.columns)}
            if name in {"order", "normal", "enrichment"}:
                metrics = df.select_dtypes(include="number").drop(columns=["time_ps"], errors="ignore")
                finite = int(np.isfinite(metrics.to_numpy(dtype=float)).sum())
                result[name]["finite_metric_values"] = finite
                print(f"NUMERICS {system} {name}: finite={finite} / {metrics.size}", flush=True)
                assert finite > 0, f"{system}: all {name} numeric metrics are non-finite"
                print(f"SAMPLE {system} {name}: {df.iloc[0].to_dict()}", flush=True)
        else:
            result[name] = {"files": len(list(path.rglob("*")))}
        print(f"OUTPUT {system} {name}: {result[name]}", flush=True)
    return result


def main():
    report = {}
    failed = False
    for system, (top_name, xtc_name) in SYSTEMS.items():
        try:
            topology, trajectory = DATA / top_name, DATA / xtc_name
            assert topology.is_file() and trajectory.is_file()
            species = inspect(system, topology, trajectory)
            config = build_config(system, topology, trajectory, species)
            print(f"Running eight-frame MDscope integration: {system}", flush=True)
            run_one(config, 0)
            outputs = check_outputs(system)
            report[system] = {"status": "passed", "species": sorted(species),
                              "outputs": outputs}
        except Exception as exc:
            failed = True
            traceback.print_exc()
            report[system] = {"status": "failed", "error": repr(exc)}
        finally:
            gc.collect()
    (ROOT / "smoke_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
