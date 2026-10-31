"""Tiny trajectory regression: the streaming coordinator must still write RMSD and RMSF."""
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from mdscope.analysis.config import DEFAULTS,Replicate
from mdscope.analysis.pipeline.streaming import run_streaming_pass

def test_streaming_rmsd_rmsf_after_module_split(tmp_path, monkeypatch):
    monkeypatch.setenv("MDSCOPE_PROFILE", "1")
    pdb=tmp_path/"tiny.pdb"
    xyz=np.array([[10.,10.,10.],[13.,10.,10.],[10.,13.,10.],[10.,10.,13.]])
    lines=["CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1"]
    for i,(x,y,z) in enumerate(xyz,1):
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}           C")
    pdb.write_text("\n".join(lines+["TER","END"])+"\n")
    u=mda.Universe(str(pdb))
    assert len(u.select_atoms("protein and name CA"))==4
    trj=tmp_path/"tiny.dcd"
    with mda.Writer(str(trj),n_atoms=len(u.atoms)) as writer:
        for j in range(8):
            # No whole-body translation; perturb one site to produce nonzero RMSD/RMSF.
            positions=xyz.copy()
            positions[3,2]+=.1*j
            u.atoms.positions=positions
            u.dimensions=np.array([50,50,50,90,90,90])
            writer.write(u.atoms)
    cfg=deepcopy(DEFAULTS)
    cfg["selections"].update(protein="protein",alignment="protein and name CA",features="protein and name CA")
    cfg["trajectory"].update(start=None,stop=None,stride=1,chunk_frames=3)
    cfg["analyses"]["membrane"]["enabled"]=False
    cfg["analyses"]["dssp"]["enabled"]=False
    cfg["analyses"]["sasa"]["enabled"]=False
    cfg["analyses"]["salt_bridges"]["enabled"]=False
    cfg["analyses"]["hydrogen_bonds"]["enabled"]=False
    cfg["analyses"]["rmsd"]["enabled"]=True
    cfg["analyses"]["rmsf"]["enabled"]=True
    cfg["analyses"]["radius_of_gyration"]["enabled"]=True
    cfg["cache"]["aligned_coordinates"]=False
    rep=Replicate("test",1,pdb,trj,"rep1")
    result=run_streaming_pass(rep,cfg,tmp_path/"results","smoke-compat")
    assert result["frames"]==8
    raw=tmp_path/"results"/"test"/"rep1"/"raw"
    tables=sorted(raw.glob("timeseries-part-*.parquet"))
    data=pd.concat([pd.read_parquet(path) for path in tables],ignore_index=True)
    assert len(data)==8
    assert "rmsd_A" in data and np.isfinite(data["rmsd_A"]).all()
    assert data["rmsd_A"].iloc[0]<1e-8
    rmsf=pd.read_parquet(raw/"rmsf.parquet")
    assert len(rmsf)==4 and np.isfinite(rmsf["rmsf_A"]).all()
    profile = json.loads((raw.parent / "streaming_phase_profile.json").read_text())
    assert profile["frames_completed"] == 8
    stages = {row["stage"] for row in profile["timings"]}
    assert {"initialization", "trajectory_decode", "other_frame_work",
            "finalization"}.issubset(stages)
