"""Fault-injection integration tests for checkpoint restart equivalence."""
from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda
from mdscope.analysis.config import DEFAULTS, Replicate
from mdscope.analysis import streaming


def _system(tmp_path):
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
    return rep, cfg


def _snapshot(root):
    raw = root / "test" / "rep1" / "raw"
    timeseries = pd.concat(
        [pd.read_parquet(p) for p in sorted(raw.glob("timeseries-part-*.parquet"))],
        ignore_index=True,
    )
    rmsf = pd.read_parquet(raw / "rmsf.parquet")
    return timeseries, rmsf


@pytest.mark.parametrize("failure_point", ["mid_chunk", "after_chunk_write"])
def test_interrupted_streaming_equals_clean_run(tmp_path, monkeypatch, failure_point):
    rep, cfg = _system(tmp_path)
    reference = tmp_path / "reference"
    interrupted = tmp_path / "interrupted"
    streaming.run_streaming_pass(rep, cfg, reference, "restart-fixture")
    original_flush = streaming._flush
    original_rg = streaming.compute_rg
    counts = {"n": 0}

    if failure_point == "mid_chunk":
        def crash_during_frame(*args, **kwargs):
            counts["n"] += 1
            if counts["n"] == 5:  # chunk of 3 already committed; frame 5 uncommitted
                raise RuntimeError("simulated mid-chunk crash")
            return original_rg(*args, **kwargs)
        monkeypatch.setattr(streaming, "compute_rg", crash_during_frame)
    else:
        def crash_after_write(*args, **kwargs):
            original_flush(*args, **kwargs)
            counts["n"] += 1
            if counts["n"] == 2:  # chunk persisted but checkpoint not advanced
                raise RuntimeError("simulated crash before checkpoint")
        monkeypatch.setattr(streaming, "_flush", crash_after_write)

    with pytest.raises(RuntimeError, match="simulated"):
        streaming.run_streaming_pass(rep, cfg, interrupted, "restart-fixture")
    monkeypatch.setattr(streaming, "_flush", original_flush)
    monkeypatch.setattr(streaming, "compute_rg", original_rg)

    streaming.run_streaming_pass(rep, cfg, interrupted, "restart-fixture")
    actual_t, actual_r = _snapshot(interrupted)
    expected_t, expected_r = _snapshot(reference)
    assert len(actual_t) == len(expected_t) == 8
    assert actual_t["time_ps"].is_monotonic_increasing
    assert actual_t["time_ps"].is_unique
    pd.testing.assert_frame_equal(actual_t, expected_t, check_exact=False, rtol=1e-10, atol=1e-10)
    pd.testing.assert_frame_equal(actual_r, expected_r, check_exact=False, rtol=1e-10, atol=1e-10)


def test_unusable_checkpoint_preserves_existing_output(tmp_path):
    rep, cfg = _system(tmp_path)
    output = tmp_path / "interrupted"
    raw = output / "test" / "rep1" / "raw"
    cache = output / "test" / "rep1" / "cache"
    raw.mkdir(parents=True)
    cache.mkdir(parents=True)
    persisted = raw / "timeseries-part-000000.parquet"
    persisted.write_bytes(b"committed bytes")
    (cache / "stream_checkpoint.npz").write_bytes(b"corrupt checkpoint")
    with pytest.raises(RuntimeError, match="Invalid streaming checkpoint"):
        streaming.run_streaming_pass(rep, cfg, output, "restart-fixture")
    assert persisted.read_bytes() == b"committed bytes"


def test_missing_checkpoint_never_deletes_existing_parts(tmp_path):
    rep, cfg = _system(tmp_path)
    root = tmp_path / "recover"
    raw = root / "test" / "rep1" / "raw"
    raw.mkdir(parents=True)
    part = raw / "timeseries-part-000000.parquet"
    part.write_bytes(b"user output; do not delete")
    with pytest.raises(RuntimeError, match="Unsafe streaming restart"):
        streaming.run_streaming_pass(rep, cfg, root, "restart-fixture")
    assert part.read_bytes() == b"user output; do not delete"


def test_checkpoint_missing_committed_part_preserves_remaining(tmp_path):
    rep, cfg = _system(tmp_path)
    root = tmp_path / "recover"
    streaming.run_streaming_pass(rep, cfg, root, "restart-fixture")
    raw = root / "test" / "rep1" / "raw"
    missing = raw / "timeseries-part-000001.parquet"
    missing.unlink()
    survivor = raw / "timeseries-part-000002.parquet"
    before = survivor.read_bytes()
    with pytest.raises(RuntimeError, match="missing committed timeseries parts"):
        streaming.run_streaming_pass(rep, cfg, root, "restart-fixture")
    assert survivor.read_bytes() == before
