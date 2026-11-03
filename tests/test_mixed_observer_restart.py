"""Mixed raw/fused/specialized observers survive a checkpoint interruption."""
from copy import deepcopy
from types import SimpleNamespace

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS, Replicate
from mdscope.pipeline import streaming, runner
from test_sdc4_geometry_references import _fixture


def _project(tmp_path, chunk):
    u, partial = _fixture()
    u.atoms.chainIDs = ["H" if c == "HS" else c for c in u.atoms.chainIDs]
    # Adjacent identical PDB residue IDs are one residue on re-read.
    u.residues.resids = np.arange(1, len(u.residues)+1)
    pdb, dcd = tmp_path/"mixed.pdb", tmp_path/"mixed.dcd"
    u.atoms.write(str(pdb))
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
    cfg["analyses"]["sdc4_mechanisms"] = partial["analyses"]["sdc4_mechanisms"]
    cfg["analyses"]["sdc4_mechanisms"]["hs_chains"] = {"A1": "chainID H"}
    for name in ("bilayer_normal", "lipid_enrichment"):
        cfg["analyses"][name].update(enabled=True, species_selections={"POPC": "resname POPC and name P"},
                                     protein_selection="protein")
    cfg["analyses"]["density_1d"] = {"enabled": True, "runs": {
        "protein": {"selections": {"protein": "protein"}, "axis": "z", "bin_width_A": 5.},
        "lipid": {"selections": {"lipid": "resname POPC"}, "axis": "z", "bin_width_A": 5.},
    }}
    cfg["trajectory"].update(start=0, stop=8, stride=2, chunk_frames=chunk)
    cfg["cache"]["aligned_coordinates"] = False
    cfg["statistics"].update(bootstrap_draws=100, families={"profiles": False, "interactions": False, "solvation_rdfs": False})
    cfg["plots"]["enabled"] = False
    cfg["systems"] = {"test": {"topology": str(pdb), "replicates": [{"trajectory": str(dcd), "label": "rep1"}]}}
    return Replicate("test", 1, pdb, dcd, "rep1"), cfg


def _final_tables(root):
    return {str(p.relative_to(root)): pd.read_parquet(p) for p in root.rglob("*.parquet")
            if "cache" not in p.parts and "raw" not in p.parts}


@pytest.mark.parametrize("chunk,failure", [(1, "during_frame"), (2, "after_write")])
def test_mixed_restart_and_offline_finalization(tmp_path, monkeypatch, chunk, failure):
    rep, cfg = _project(tmp_path, chunk)
    clean, resumed = tmp_path/"clean", tmp_path/"resumed"
    reference_cfg = deepcopy(cfg)
    reference_cfg["trajectory"]["chunk_frames"] = 3
    streaming.run_streaming_pass(rep, reference_cfg, clean, "reference")
    calls = 0
    target = "compute_rg" if failure == "during_frame" else "_flush"
    original = getattr(streaming, target)
    def interrupted(*a, **k):
        nonlocal calls
        calls += 1
        if failure == "after_write":
            result = original(*a, **k)
        if calls == (3 if failure == "during_frame" else 2):
            raise RuntimeError("injected mixed-observer interruption")
        return result if failure == "after_write" else original(*a, **k)
    monkeypatch.setattr(streaming, target, interrupted)
    with pytest.raises(RuntimeError, match="injected"):
        streaming.run_streaming_pass(rep, cfg, resumed, "restart")
    monkeypatch.setattr(streaming, target, original)
    streaming.run_streaming_pass(rep, cfg, resumed, "restart")

    # Summary consolidation and the actual project finalizer must be offline.
    rep.topology.unlink(); rep.trajectory.unlink()
    def forbidden(*a, **k):
        raise AssertionError("offline finalization attempted to reopen a trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    for root in (clean, resumed):
        runner.summarize_replicate(root/"test/rep1", cfg)
        runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=root,
                                state=SimpleNamespace(config_hash="mixed"), phase="statistics")
    expected, actual = _final_tables(clean), _final_tables(resumed)
    assert expected.keys() == actual.keys()
    assert any("density_1d/protein" in path for path in actual)
    assert any("density_1d/lipid" in path for path in actual)
    assert any("sdc4_mechanisms" in path for path in actual)
    for path in expected:
        pd.testing.assert_frame_equal(actual[path], expected[path], check_exact=False, rtol=1e-9, atol=1e-9)
    table = actual["test/rep1/summary/timeseries.parquet"]
    assert len(table) == 4 and table.time_ps.is_unique and table.time_ps.is_monotonic_increasing
    np.testing.assert_allclose(table.time_ps, [0., 2., 4., 6.], rtol=1e-6, atol=1e-6)
    stats = actual["project_summary/statistics/time_series/replicate_statistics.parquet"]
    assert set(stats.replicate) == {"rep1"}
    relative = "test/rep1/extensions/sdc4_mechanisms/lipid_density_maps.npz"
    with np.load(clean/relative) as left, np.load(resumed/relative) as right:
        assert int(right["nframes"]) == 4
        assert left.files == right.files
        for key in left.files:
            np.testing.assert_array_equal(left[key], right[key])


@pytest.mark.parametrize("times", [[0.,0.,1.], [0.,2.,1.], [0.,np.inf], [0.,np.nan]])
def test_invalid_persisted_times_fail_before_archiving_existing_results(tmp_path, times):
    folder = tmp_path/"test/rep1/summary"
    folder.mkdir(parents=True)
    pd.DataFrame({"time_ps": times, "metric": np.ones(len(times))}).to_parquet(folder/"timeseries.parquet")
    previous = tmp_path/"project_summary/keep.txt"
    previous.parent.mkdir(); previous.write_text("existing valid results")
    cfg = {"systems": {"test": {"topology": "absent.pdb", "replicates": [{"label": "rep1", "trajectory": "absent.dcd"}]}}}
    with pytest.raises(RuntimeError, match="test/rep1: invalid time_ps"):
        runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=tmp_path,
                                state=SimpleNamespace(config_hash="bad"), phase="statistics")
    assert previous.read_text() == "existing valid results"
    assert not (tmp_path/".archives").exists()
