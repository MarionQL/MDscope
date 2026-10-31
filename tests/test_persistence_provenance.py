"""Issue #28: historical imports and crash-safe output contracts."""
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from mdscope import config
from mdscope.core import persistence, provenance


@pytest.mark.parametrize("name", [
    "ensure_dir", "atomic_json", "atomic_text", "atomic_parquet",
])
def test_persistence_legacy_imports_are_identical(name):
    assert getattr(config, name) is getattr(persistence, name)


@pytest.mark.parametrize("name", ["config_hash", "environment_metadata"])
def test_provenance_legacy_imports_are_identical(name):
    assert getattr(config, name) is getattr(provenance, name)


def test_atomic_text_and_json_roundtrip_and_replace(tmp_path):
    target = tmp_path / "nested" / "state.txt"
    config.atomic_text(target, "old")
    config.atomic_text(target, "new")
    assert target.read_text() == "new"
    metadata = tmp_path / "nested" / "data.json"
    config.atomic_json(metadata, {"x": 1})
    config.atomic_json(metadata, {"x": 2})
    assert json.loads(metadata.read_text()) == {"x": 2}
    assert set(target.parent.iterdir()) == {target, metadata}


def test_parquet_roundtrip_and_replace(tmp_path):
    dest = tmp_path / "nested" / "summary.parquet"
    config.atomic_parquet(dest, pd.DataFrame({"time_ps": [1., 2.], "metric": [4, 5]}))
    config.atomic_parquet(dest, pd.DataFrame({"time_ps": [3.], "metric": [6]}))
    pd.testing.assert_frame_equal(pd.read_parquet(dest), pd.DataFrame({"time_ps": [3.], "metric": [6]}))
    assert list(dest.parent.iterdir()) == [dest]


@pytest.mark.parametrize("writer,initial,replacement", [
    (persistence.atomic_text, "previous", "replacement"),
    (persistence.atomic_json, {"value": 1}, {"value": 2}),
    (persistence.atomic_parquet, pd.DataFrame({"x": [1]}), pd.DataFrame({"x": [2]})),
])
def test_failed_replace_preserves_existing_file_and_removes_temporary(
    tmp_path, monkeypatch, writer, initial, replacement,
):
    target = tmp_path / "saved"
    writer(target, initial)
    old = target.read_bytes()

    def broken_replace(source, destination):
        raise OSError("injected replacement failure")

    monkeypatch.setattr(persistence.os, "replace", broken_replace)
    with pytest.raises(OSError, match="injected replacement failure"):
        writer(target, replacement)
    assert target.read_bytes() == old
    assert list(tmp_path.iterdir()) == [target]


def test_config_hash_matches_previous_deterministic_contract():
    obj = {"b": [2, 3], "a": {"value": 1}}
    expected = hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]
    assert config.config_hash(obj) == expected
    assert config.config_hash({"a": {"value": 1}, "b": [2, 3]}) == expected
