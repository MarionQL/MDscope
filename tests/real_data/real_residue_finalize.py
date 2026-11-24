"""Bounded release-asset regression; not convergence or production performance."""
from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import pandas as pd
import yaml
from real_trajectory_smoke import DATA, ROOT, SYSTEMS

from mdscope.config import iter_replicates, load_config
from mdscope.core.sequence import residue_label
from mdscope.observables.dssp import _run_dssp
from mdscope.observables.sasa import _run_sasa
from mdscope.pipeline.runner import finalize_project, run_one
from mdscope.statistics.residue_correlations import correlate_replica


def main():
    output = ROOT / "residue_finalize"
    if output.exists():
        raise RuntimeError(f"Use a fresh output directory: {output}")
    ROOT.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path("synd4_full_mechanisms.yaml").read_text())
    cfg["project"] = {"name": "release_residue_finalize", "output": str(output)}
    cfg["systems"] = {
        system: {"topology": str(DATA / top),
                 "replicates": [{"label": "rep1", "trajectory": str(DATA / traj)}]}
        for system, (top, traj) in SYSTEMS.items()
    }
    cfg["trajectory"] = {"start": 0, "stop": 40, "stride": 1, "chunk_frames": 10}
    enabled = {"alignment", "rmsd", "dssp", "sasa", "sdc4_mechanisms"}
    for name in list(cfg["analyses"]):
        if name not in enabled:
            cfg["analyses"][name] = {"enabled": False}
    mechanism = cfg["analyses"]["sdc4_mechanisms"]
    mechanism["local_order"] = {"enabled": False}
    mechanism["ceramide_scaffold_pooling"] = False
    cfg["cache"]["aligned_coordinates"] = False
    cfg["statistics"]["bootstrap_draws"] = 100
    cfg["statistics"]["matched_window_ps"] = None
    cfg["plots"] = {"enabled": False}
    path = ROOT / "residue_finalize.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    resolved = load_config(path)
    diagnostics = {}
    for index, rep in enumerate(iter_replicates(resolved, path)):
        u = mda.Universe(str(rep.topology), str(rep.trajectory))
        protein = u.select_atoms("protein")
        labels = [residue_label(r) for r in protein.residues]
        assert set(protein.chainIDs) == {"PROA", "PROB"}
        assert labels == [residue_label(r.atoms[0]) for r in protein.residues]
        assert len(set(labels)) == len(labels)
        diagnostics[rep.system] = {
            "frames": len(u.trajectory), "chainIDs": sorted(set(protein.chainIDs)),
            "segids": sorted(set(protein.segids)),
            "residue_identities": [{"index": int(r.ix), "label": residue_label(r),
                                    "segid": r.segid} for r in protein.residues],
        }
        u.trajectory.close()
        run_one(path, index)
        repout = output / rep.system / rep.label
        # Independent producers on three source frames must agree with fused
        # output by exact time/index/label. Geometry is checked elsewhere.
        u = mda.Universe(str(rep.topology), str(rep.trajectory))
        independent_cfg = {**resolved, "trajectory": {**resolved["trajectory"], "stop": 3}}
        independent = ROOT / f"independent_{rep.system}"
        _run_sasa(u, independent_cfg, independent)
        _run_dssp(u, independent_cfg, independent)
        u.trajectory.close()
        for relative in ("profiles/dssp/assignments.parquet",
                         "profiles/sasa/isolated/residue_sasa.parquet",
                         "profiles/sasa/environment/residue_sasa.parquet"):
            actual = pd.read_parquet(repout / relative)
            expected = pd.read_parquet(independent / relative)
            keys = ["time_ps", "topology_residue_index", "residue"]
            actual = actual[actual.time_ps.isin(expected.time_ps)].sort_values(keys).reset_index(drop=True)
            expected = expected.sort_values(keys).reset_index(drop=True)
            pd.testing.assert_frame_equal(actual, expected, check_dtype=False, rtol=1e-10, atol=1e-10)
        correlations = correlate_replica(repout, min_frames=30)
        counts = correlations.status.value_counts().to_dict()
        diagnostics[rep.system]["correlation_status"] = counts
        diagnostics[rep.system]["max_paired_frames"] = int(correlations.n_paired_frames.max())
        eligible = correlations.n_paired_frames.ge(30)
        if eligible.any():
            assert correlations.status.eq("ok").any(), counts
        else:
            print(f"INSUFFICIENT SUPPORT: {rep.system}: {counts}", flush=True)
        print(f"CORRELATIONS {rep.system}: {counts}", flush=True)
    # Finalization is an offline contract, even with actual trajectory products.
    original = mda.Universe
    def forbidden(*args, **kwargs):
        raise AssertionError("Finalization reopened a trajectory")
    mda.Universe = forbidden
    try:
        finalize_project(path, cfg=resolved, root=output)
        plotting = {**resolved, "plots": {"enabled": True, "include": [
            "profiles:sasa_mean__isolated:profile", "profiles:dssp_helix_occupancy:profile",
        ]}}
        finalize_project(path, cfg=plotting, root=output, phase="plotting")
    finally:
        mda.Universe = original
    manifest = json.loads((output / "project_summary/statistics/manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["derived"]["residue_correlations"] == "complete"
    assert list((output / "project_summary/figures").rglob("*.png"))
    (ROOT / "residue_finalize_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    print("PASS: short-release identities, fused/independent producers, offline finalize and plotting")


if __name__ == "__main__":
    main()
