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

from mdscope.pipeline.runner import run_one
from mdscope.config import Replicate, load_config
from mdscope.observables.membrane_order import run_membrane_order
from mdscope.observables.bilayer_normal import run_bilayer_normal
from mdscope.observables.lipid_enrichment import run_lipid_enrichment
from mdscope.interactions.lipid_hydrogen_bonds import run_lipid_hydrogen_bonds
from mdscope.observables.residue_membrane_distance import run_residue_membrane_distance

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
        ("DOPC", "P"), ("DPPC", "P"), ("PSM", "P"), ("POPS", "P"),
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
    # Glycosphingolipids in these CHARMM-GUI topologies use a CER160
    # ceramide anchor and carbohydrate residues in the same chainID.
    # In the no-headgroup control the same CER160 residues are chainID CER160.
    for name in ("GM1a", "BGAL", "CER160"):
        selection = f"chainID {name} and resname CER160 and name NF"
        n = len(u.select_atoms(selection))
        print(f"ceramide_anchor {name}={n}", flush=True)
        if n:
            species[name] = selection
    print("non-water residue counts (top 25):", flush=True)
    from collections import Counter
    nonwater = Counter(str(r.resname) for r in u.residues if str(r.resname) != "TIP3")
    print(nonwater.most_common(25), flush=True)
    assert len(species) >= 2, "Insufficient lipid anchors for a bilayer test"
    del u
    gc.collect()
    return species


def build_config(system: str, topology: Path, trajectory: Path, species: dict) -> Path:
    # Restrict all analyses to the first forty frames. We are checking that
    # real TPR+XTC inputs execute, NOT estimating equilibrated observables.
    lipid_types = [name for name in ("POPC", "POPE", "DOPC", "DPPC", "POPI", "PSM") if name in species]
    assert lipid_types
    lipid_sel = "resname " + " ".join(lipid_types[:2])
    cfg = {
        "project": {"name": f"real_smoke_{system}", "output": str(ROOT / system)},
        "trajectory": {"start": 0, "stop": 40, "stride": 1, "chunk_frames": 10},
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
                "footprint_selection": "protein and (resid 146:174 or resid 416:444)",
                "shell_edges_A": [0, 5, 10, 20],
            },
            "lipid_hydrogen_bonds": {
                "enabled": True, "lipid_selection": lipid_sel,
                "exclude_same_residue": True,
            },
            "residue_membrane_distance": {
                "enabled": True,
                "selection": "protein and name CA and resid 150:170",
                "expand_residue_atoms": False,
                "membrane_selection": "resname POPC and name P",
                "exclude_hydrogens": True,
                "groups": False,
            },
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
        "enrichment": rep / "extensions" / "lipid_footprint" / "timeseries.parquet",
        "enrichment_counts": rep / "extensions" / "lipid_footprint" / "counts.parquet",
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
    events_path = rep / "interactions" / "lipid_hydrogen_bonds" / "events.parquet"
    assert events_path.is_file(), f"{system}: missing lipid H-bond events table"
    events = pd.read_parquet(events_path)
    result["lipid_hbond_events"] = len(events)
    print(f"H-BOND EVENTS {system}: {len(events)}", flush=True)
    counts = pd.read_parquet(paths["lipid_hbonds"])
    print(f"H-BOND COUNTS {system}: {counts.head().to_dict('records')}", flush=True)
    return result



def compare_fused_to_independent(system, config, topology, trajectory):
    """Scientific parity: fused production outputs vs original three readers.

    The reference computations intentionally reopen the XTC *only in CI*.
    Production never calls them after the fused streaming stage.
    """
    cfg = load_config(config)
    rep = Replicate(system, 1, topology, trajectory, "rep1")
    fused = ROOT / system / system / "rep1"
    independent = ROOT / f"{system}_independent"
    independent.mkdir(parents=True, exist_ok=True)
    checks = {}
    for name, fn, path in (
        ("membrane_order", run_membrane_order, "profiles/membrane_order"),
        ("bilayer_normal", run_bilayer_normal, "extensions/bilayer_normal"),
        ("lipid_enrichment", run_lipid_enrichment, "extensions/lipid_footprint"),
    ):
        fn(rep, cfg, independent)
        left = fused / path
        right = independent / path
        files = ["timeseries.parquet"]
        if name == "membrane_order":
            files += ["site_mapping.parquet", "profile_index.parquet"]
        if name == "lipid_enrichment":
            files += ["counts.parquet", "selection_audit.parquet"]
        for file in files:
            actual = pd.read_parquet(left / file)
            expected = pd.read_parquet(right / file)
            # Exact column order and row identity are part of the public API.
            pd.testing.assert_frame_equal(
                actual, expected, check_dtype=False,
                check_exact=False, rtol=1e-6, atol=1e-6,
                obj=f"{system}/{name}/{file}",
            )
        if name == "membrane_order":
            index = pd.read_parquet(left / "profile_index.parquet")
            for file in index["file"]:
                with np.load(left / file) as a, np.load(right / file) as b:
                    assert set(a.files) == set(b.files)
                    for key in a.files:
                        if a[key].dtype.kind in "fiu":
                            np.testing.assert_allclose(
                                a[key], b[key], rtol=1e-6, atol=1e-6,
                                equal_nan=True, err_msg=f"{system}/{file}/{key}"
                            )
                        else:
                            np.testing.assert_array_equal(a[key], b[key])
        checks[name] = {"status": "identical_within_tolerance", "frames": 40}
        print(f"FUSED PARITY {system}/{name}: PASS", flush=True)
    # The fused MDAnalysis HBA uses the same detector, stepped once per
    # streaming frame. Compare all events, time-series and molecule exclusions
    # against the original independent trajectory traversal.
    run_lipid_hydrogen_bonds(rep, cfg, independent)
    for file in ("events.parquet", "count_timeseries.parquet", "molecule_audit.parquet"):
        actual = pd.read_parquet(fused / "interactions/lipid_hydrogen_bonds" / file)
        expected = pd.read_parquet(independent / "interactions/lipid_hydrogen_bonds" / file)
        pd.testing.assert_frame_equal(
            actual, expected, check_dtype=False,
            check_exact=False, rtol=1e-6, atol=1e-6,
            obj=f"{system}/lipid_hydrogen_bonds/{file}",
        )
    checks["lipid_hydrogen_bonds"] = {"status": "identical_within_tolerance"}
    print(f"FUSED PARITY {system}/lipid_hydrogen_bonds: PASS", flush=True)
    ext_timing = json.loads((fused / "fused_extensions_timings.json").read_text())
    assert "lipid_hydrogen_bonds" in ext_timing["enabled"]
    # Exact per-residue minimum distances and all metadata must be
    # invariant to whether the frames were read independently.
    run_residue_membrane_distance(rep, cfg, independent)
    left = fused / "profiles/residue_membrane_distance/default"
    right = independent / "profiles/residue_membrane_distance/default"
    for filename in ("distance_matrix_A.npy", "time_ps.npy", "membrane_atom_indices.npy"):
        np.testing.assert_allclose(
            np.load(left / filename), np.load(right / filename),
            rtol=1e-6, atol=1e-6, equal_nan=True,
            err_msg=f"{system}/residue_membrane_distance/{filename}",
        )
    for filename in ("median_profile.parquet", "residues.parquet", "selection_audit.parquet"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(left / filename), pd.read_parquet(right / filename),
            check_dtype=False, check_exact=False, rtol=1e-6, atol=1e-6,
            obj=f"{system}/residue_membrane_distance/{filename}",
        )
    checks["residue_membrane_distance"] = {"status": "identical_within_tolerance"}
    print(f"FUSED PARITY {system}/residue_membrane_distance: PASS", flush=True)
    timing = json.loads((fused / "fused_membrane_timings.json").read_text())
    assert timing["shared_anchor_mappings"] == 1
    stages = json.loads((fused / "stage_timings.json").read_text())["stages"]
    assert "streaming" in {row["stage"] for row in stages}
    assert not {"membrane-order", "bilayer_normal", "lipid_enrichment"}.intersection(
        row["stage"] for row in stages
    )
    return {"checks": checks, "shared_anchor_mappings": timing["shared_anchor_mappings"],
            "fused_compute_seconds": timing["compute_seconds_by_analysis"]}

def main():
    report = {}
    failed = False
    for system, (top_name, xtc_name) in SYSTEMS.items():
        try:
            topology, trajectory = DATA / top_name, DATA / xtc_name
            assert topology.is_file() and trajectory.is_file()
            species = inspect(system, topology, trajectory)
            config = build_config(system, topology, trajectory, species)
            print(f"Running forty-frame MDscope integration: {system}", flush=True)
            run_one(config, 0)
            outputs = check_outputs(system)
            parity = compare_fused_to_independent(system, config, topology, trajectory)
            report[system] = {"status": "passed", "species": sorted(species),
                              "outputs": outputs, "fused_parity": parity}
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
