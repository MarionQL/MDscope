"""Regression tests for restart marker and persisted summary integrity (#128)."""

import json

import pandas as pd

from mdscope.pipeline.runner import RunState


def test_runstate_accepts_completed_legacy_marker_regardless_of_hash(tmp_path):
    state = RunState(tmp_path)
    state.marker("stage").write_text(json.dumps({
        "config_hash": "expected", "status": "failed"
    }))
    assert not state.complete("stage")
    state.marker("stage").write_text(json.dumps({
        "config_hash": "other", "status": "complete"
    }))
    assert state.complete("stage")
    state.marker("stage").write_text("broken JSON")
    assert not state.complete("stage")
    state.mark("stage", {"status": "failed"})
    assert state.complete("stage")


def test_summary_marker_requires_inventory_and_valid_parquet(tmp_path):
    state = RunState(tmp_path)
    path = tmp_path / "replica" / "summary" / "timeseries.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0.0, 10.0], "value": [1., 2.]}).to_parquet(path)
    state.mark("summary")
    # Historical markers can reuse readable expected products.
    assert state.complete("summary", required_artifacts=(path,))
    state.mark("summary", {"artifacts": {str(path): path.stat().st_size}})
    assert state.complete("summary", required_artifacts=(path,))

    path.unlink()
    assert not state.complete("summary", required_artifacts=())
    pd.DataFrame({"time_ps": [0.0]}).to_parquet(path)
    assert not state.complete("summary", required_artifacts=(path,))
    state.mark("summary", {"artifacts": {str(path): path.stat().st_size}})
    assert state.complete("summary", required_artifacts=(path,))

    # Preserve the file length to check footer integrity independently of size.
    path.write_bytes(b"X" * path.stat().st_size)
    assert not state.complete("summary", required_artifacts=(path,))


def test_summary_marker_rejects_missing_secondary_artifact(tmp_path):
    state = RunState(tmp_path)
    summary = tmp_path / "summary"
    summary.mkdir()
    whole = summary / "timeseries.parquet"
    grouped = summary / "group_timeseries.parquet"
    pd.DataFrame({"time_ps": [0.0]}).to_parquet(whole)
    pd.DataFrame({"group": ["A"], "time_ps": [0.0]}).to_parquet(grouped)
    state.mark("summary", {"artifacts": {
        str(p): p.stat().st_size for p in (whole, grouped)
    }})
    assert state.complete("summary", required_artifacts=(whole, grouped))
    grouped.unlink()
    assert not state.complete("summary", required_artifacts=(whole,))


def test_finalizer_rejects_missing_standalone_profile_before_archiving(tmp_path):
    from copy import deepcopy
    import pytest
    from mdscope.config import DEFAULTS
    from mdscope.pipeline.runner import finalize_project
    cfg = deepcopy(DEFAULTS)
    for spec in cfg['analyses'].values():
        spec['enabled'] = False
    cfg['analyses']['dssp'].update(enabled=True, backend='multiprocessing')
    cfg['systems'] = {'test': {'topology': 'absent.pdb', 'replicates': [
        {'trajectory': 'absent.xtc', 'label': 'rep1'}]}}
    summary = tmp_path / 'test/rep1/summary'
    summary.mkdir(parents=True)
    pd.DataFrame({'time_ps': [0., 1.], 'value': [1., 2.]}).to_parquet(summary/'timeseries.parquet')
    previous = tmp_path/'project_summary/keep.txt'
    previous.parent.mkdir()
    previous.write_text('preserve')
    with pytest.raises(RuntimeError, match='dssp'):
        finalize_project(tmp_path/'unused.yaml', cfg, tmp_path, RunState(tmp_path), phase='statistics')
    assert previous.read_text() == 'preserve'


def test_completed_run_detects_missing_companion_without_trajectory_read(tmp_path, monkeypatch):
    import MDAnalysis as mda
    import pytest
    from mdscope.pipeline import runner
    from test_mixed_observer_restart import _project
    rep, cfg = _project(tmp_path, 3)
    root = tmp_path/'results'
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    directory = root/'test/rep1/profiles/residue_membrane_distance/A'
    (directory/'time_ps.npy').unlink()
    primary = (directory/'distance_matrix_A.npy').read_bytes()
    def forbidden(*a, **k):
        raise AssertionError('validation reopened a trajectory')
    monkeypatch.setattr(mda, 'Universe', forbidden)
    with pytest.raises(RuntimeError, match='time_ps.npy'):
        runner.run_replicate(rep, cfg, root, state)
    assert (directory/'distance_matrix_A.npy').read_bytes() == primary


import pytest


@pytest.mark.parametrize('damage', ['json', 'npz', 'schema', 'shape'])
def test_finalizer_rejects_invalid_companions_before_mutation(tmp_path, damage):
    import numpy as np
    from mdscope.pipeline import runner
    from test_mixed_observer_restart import _project
    rep, cfg = _project(tmp_path, 3)
    root = tmp_path/'results'
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    base = root/'test/rep1'
    if damage == 'json':
        path = base/'extensions/sdc4_mechanisms/selection_audit.json'
        path.write_text('[]')
    elif damage == 'npz':
        path = base/'extensions/sdc4_mechanisms/lipid_density_maps.npz'
        path.write_bytes(b'not an archive')
    elif damage == 'schema':
        path = base/'extensions/solvation/counts.parquet'
        pd.DataFrame({'wrong': [1]}).to_parquet(path)
    else:
        path = base/'profiles/residue_membrane_distance/A/time_ps.npy'
        np.save(path, [0.])
    previous = root/'project_summary/keep.txt'
    previous.parent.mkdir(); previous.write_text('preserve')
    with pytest.raises(RuntimeError, match='invalid required'):
        runner.finalize_project(tmp_path/'unused.yaml', cfg, root, state, phase='statistics')
    assert previous.read_text() == 'preserve'
