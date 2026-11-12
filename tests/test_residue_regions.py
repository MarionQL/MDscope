"""Regression coverage for residue adjacency and named SASA profiles."""
from copy import deepcopy
import json

import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS
from mdscope.statistics.profiles import analyze_profiles
from mdscope.statistics.residue_regions import scan_saved_residue_regions


def test_saved_scan_requires_chain_adjacency_and_complete_replica_support(tmp_path):
    folder = tmp_path / "profiles"
    folder.mkdir()
    paths = [f"{system}/rep{rep}" for system in ("a", "b") for rep in (1, 2)]
    (tmp_path / "manifest.json").write_text(json.dumps({
        "status": "complete", "families": {"profiles": "enabled"}, "replicate_paths": paths,
    }))
    # A1–3 pass, A4 overlaps, A5 is incomplete, A6–7 change in the
    # opposite direction and B8 cannot extend that run across chains.
    effects = {"A:ALA1": .5, "A:ALA2": .6, "A:ALA3": .7, "A:ALA4": .05,
               "A:ALA5": .5, "A:ALA6": -.2, "A:ALA7": -.2, "B:ALA8": -.2}
    rows = []
    for path in paths:
        system, replica = path.split("/")
        for label, effect in effects.items():
            if path == "b/rep2" and label == "A:ALA5":
                continue
            rows.append({
                "system": system, "replicate": replica, "position_label": label,
                "profile": "dssp_helix_occupancy",
                "value": .2 + .1 * (replica == "rep2") + effect * (system == "b"),
            })
    pd.DataFrame(rows).to_parquet(folder / "replicate_profiles.parquet")
    row, = pd.read_parquet(scan_saved_residue_regions(tmp_path)).to_dict("records")
    assert (row["chain"], row["start_residue"], row["end_residue"]) == ("A", 1, 3)
    assert row["difference_system2_minus_system1"] == pytest.approx(.6)
    assert row["min_replica_separation"] == pytest.approx(.4)


def test_named_sasa_and_burial_reach_profile_statistics_and_region_scan(tmp_path):
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["families"] = {"profiles": {"enabled": True, "bootstrap": False}}
    cfg["analyses"]["sasa"] = {
        "enabled": True, "runs": {"isolated": {}, "environment": {"reference_run": "isolated"}},
    }
    paths = {}
    for system, exposed in (("a", 2.), ("b", 1.)):
        paths[system] = []
        for replica in ("rep1", "rep2"):
            rep = tmp_path / system / replica
            paths[system].append(rep)
            frame = pd.DataFrame({
                "time_ps": np.repeat([0., 10.], 3),
                "topology_residue_index": np.tile([0, 1, 2], 2),
                "residue": ["A:ALA1", "A:VAL2", "A:GLY3"] * 2,
            })
            sasa = rep / "profiles/sasa"
            for name, value in (("isolated", 3.), ("environment", exposed)):
                (sasa / name).mkdir(parents=True)
                frame.assign(sasa_nm2=value).to_parquet(sasa / name / "residue_sasa.parquet")
            frame.assign(buried_sasa_nm2=3. - exposed).to_parquet(
                sasa / "environment/residue_burial.parquet",
            )
    analyze_profiles(cfg, tmp_path, paths)
    dest = tmp_path / "project_summary/statistics/profiles"
    comparisons = pd.read_parquet(dest / "system_comparisons.parquet")
    assert set(comparisons.profile) == {
        "sasa_mean__isolated", "sasa_mean__environment", "sasa_buried_nm2__environment",
    }
    regions = pd.read_parquet(dest / "residue_regions.parquet").set_index("profile")
    assert regions["difference_system2_minus_system1"].to_dict() == {
        "sasa_buried_nm2__environment": 1., "sasa_mean__environment": -1.,
    }
