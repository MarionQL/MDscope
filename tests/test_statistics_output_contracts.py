"""No-data preview base columns must exist in actual statistics writer output."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS
from mdscope.statistics.contracts import statistics_output_contract


@pytest.mark.parametrize("family", ["time_series", "profiles", "interactions", "solvation_rdfs", "diffusion", "coupling"])
def test_statistics_base_contract_matches_writer(tmp_path, monkeypatch, family):
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"].update(bootstrap=False, bootstrap_draws=100)
    cfg["statistics"]["hypothesis_tests"]["enabled"] = False
    paths = {"a": [], "b": []}
    for system in paths:
        for i in range(2):
            rep = tmp_path / system / f"rep{i}"
            (rep / "summary").mkdir(parents=True)
            t = np.arange(30, dtype=float)
            pd.DataFrame({"time_ps": t, "rmsd_A": 1 + .1*t+i, "radius_of_gyration_A": 2+.2*t+i}).to_parquet(rep / "summary/timeseries.parquet")
            paths[system].append(rep)
    if family == "time_series":
        from mdscope.statistics.time_series import analyze_time_series
        analyze_time_series(cfg, tmp_path, paths)
        folder = "time_series"
    elif family == "profiles":
        import mdscope.statistics.profiles as mod
        monkeypatch.setattr(mod, "collect_temporal_profiles", lambda rep, cfg: {
            "profile": {"times": np.arange(30.), "positions": np.arange(2.), "labels": np.array(["one", "two"]),
                        "values": np.arange(60.).reshape(30, 2), "reducer": "mean"}})
        mod.analyze_profiles(cfg, tmp_path, paths)
        folder = "profiles"
    elif family == "interactions":
        from mdscope.statistics.interactions import analyze_interactions
        for reps in paths.values():
            for rep in reps:
                directory = rep / "interactions/salt_bridges"
                directory.mkdir(parents=True)
                keys = {"residue_1": "ASP1", "residue_2": "LYS2", "interaction_class": "interchain"}
                pd.DataFrame([dict(keys, time_ps=float(t)) for t in range(5)]).to_parquet(directory / "presence.parquet")
                pd.DataFrame([dict(keys, start_ps=0., lifetime_ps=5., left_censored=True, completed=True)]).to_parquet(directory / "lifetimes.parquet")
        analyze_interactions(cfg, tmp_path, paths)
        folder = "interactions"
    elif family == "solvation_rdfs":
        from mdscope.statistics.solvation import analyze_solvation_rdfs
        for reps in paths.values():
            for rep in reps:
                directory = rep / "extensions/solvation"
                directory.mkdir(parents=True)
                pd.DataFrame({"solvent": ["water"]*5, "target_class": ["ligand"]*5, "site_class": ["oxygen"]*5,
                    "r_inner_A": np.arange(5.), "r_outer_A": np.arange(5.)+1, "g_bulk": np.ones(5)}).to_parquet(directory / "rdf.parquet")
        analyze_solvation_rdfs(cfg, tmp_path, paths)
        folder = "solvation"
    elif family == "diffusion":
        from mdscope.statistics.diffusion import analyze_multi_origin_diffusion
        cfg["analyses"]["membrane"]["msd"].update(enabled=True, multi_origin=True)
        for reps in paths.values():
            for rep in reps:
                pd.DataFrame({"metric": ["lipid"]*30, "lag_ps": np.arange(30.), "msd_A2": 4*np.arange(30.),
                              "n_time_origins": 100-np.arange(30)}).to_parquet(rep / "summary/multi_origin_msd.parquet")
        analyze_multi_origin_diffusion(cfg, tmp_path, paths)
        folder = "diffusion"
    else:
        from mdscope.statistics.coupling import analyze_coupling
        cfg["analyses"]["coupling_screen"] = {"enabled": True, "include_defaults": False, "block_frames": 3,
            "pairs": [{"name": "pair", "x": "rmsd_A", "y": "radius_of_gyration_A"}]}
        analyze_coupling(cfg, tmp_path, paths)
        folder = "coupling"
    contract = statistics_output_contract(family)
    actual = pd.read_parquet(tmp_path / "project_summary/statistics" / folder / contract["filename"])
    assert not actual.empty
    assert set(contract["columns"]) <= set(actual.columns)
