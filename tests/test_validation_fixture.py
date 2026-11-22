"""The documented YAML fixture checks numerical facts through the installed CLI."""
import shutil

import MDAnalysis as mda
import pandas as pd
import pytest
from mdscope_test_support import invoke_cli as invoke
from validation_fixture import verify_project, write_project, write_quickstart


def test_portable_yaml_run_and_independent_references(tmp_path, monkeypatch):
    original = tmp_path/'original'
    write_project(original)
    moved = tmp_path/'moved'
    shutil.move(original, moved)
    invoke(tmp_path, 'run', moved/'validation.yaml')

    def forbidden(*args, **kwargs):
        raise AssertionError('Reference verification must not reopen trajectory inputs')
    monkeypatch.setattr(mda, 'Universe', forbidden)
    assert verify_project(moved) == 7
    # A numerically wrong but perfectly readable output must fail the oracle.
    path = moved/'results/test/rep1/summary/timeseries.parquet'
    table = pd.read_parquet(path)
    table['sdc4_local_thickness_A'] -= 1.
    table.to_parquet(path)
    with pytest.raises(AssertionError, match='sdc4_local_thickness_A'):
        verify_project(moved)


def test_fixture_creation_preserves_existing_inputs(tmp_path):
    sentinel = tmp_path/'mixed.pdb'
    sentinel.write_text('existing input')
    with pytest.raises(FileExistsError, match='fresh validation directory'):
        write_project(tmp_path)
    assert sentinel.read_text() == 'existing input'
    assert not (tmp_path/'mixed.dcd').exists()


def test_quickstart_creation_preserves_existing_inputs(tmp_path):
    sentinel = tmp_path/'notes.txt'
    sentinel.write_text('existing input')
    with pytest.raises(FileExistsError, match='fresh quickstart directory'):
        write_quickstart(tmp_path)
    assert [path.name for path in tmp_path.iterdir()] == ['notes.txt']
