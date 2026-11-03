"""Resume safety: committed counters and orphan chunk pruning."""
from pathlib import Path
import pytest
from mdscope.core.checkpoints import (validate_progress, validate_committed_parts, prune_uncommitted_parts,
                                     chunk_artifacts, validate_artifact_inventory, prune_artifact_orphans)


@pytest.mark.parametrize("frames,parts,chunk,total", [
    (0, 0, 3, 8), (3, 1, 3, 8), (6, 2, 3, 8), (8, 3, 3, 8),
])
def test_valid_restart_boundaries(frames, parts, chunk, total):
    validate_progress(frames, parts, chunk, total)


@pytest.mark.parametrize("frames,parts,chunk,total", [
    (-1, 0, 3, 8), (9, 3, 3, 8), (3, 0, 3, 8),
    (3, 2, 3, 8), (0, -1, 3, 8), (3, 1, 0, 8),
])
def test_invalid_restart_boundaries_rejected(frames, parts, chunk, total):
    with pytest.raises(ValueError):
        validate_progress(frames, parts, chunk, total)


def test_prune_only_orphan_numbered_parts(tmp_path):
    names = [
        "timeseries-part-000000.parquet",
        "timeseries-part-000001.parquet",
        "timeseries-part-000002.parquet",
        "aligned-part-000002.npz",
        "aligned-part-000001.npz",
        "unrelated-part-foo.parquet",
        "notes.txt",
        "rmsf.parquet",
    ]
    for name in names:
        (tmp_path / name).write_bytes(b"fixture")
    removed = prune_uncommitted_parts(tmp_path, 2)
    assert {p.name for p in removed} == {
        "timeseries-part-000002.parquet", "aligned-part-000002.npz",
    }
    assert (tmp_path / "timeseries-part-000000.parquet").exists()
    assert (tmp_path / "timeseries-part-000001.parquet").exists()
    assert (tmp_path / "aligned-part-000001.npz").exists()
    assert (tmp_path / "unrelated-part-foo.parquet").exists()
    assert (tmp_path / "notes.txt").exists()


def test_reject_orphan_chunks_without_checkpoint(tmp_path):
    orphan = tmp_path / "timeseries-part-000000.parquet"
    orphan.write_bytes(b"preserved")
    with pytest.raises(ValueError, match="Checkpoint is absent"):
        validate_committed_parts(tmp_path, 0)
    assert orphan.read_bytes() == b"preserved"


def test_reject_missing_committed_chunk(tmp_path):
    (tmp_path / "timeseries-part-000000.parquet").write_bytes(b"first")
    (tmp_path / "timeseries-part-000002.parquet").write_bytes(b"orphan")
    with pytest.raises(ValueError, match="missing committed.*000001"):
        validate_committed_parts(tmp_path, 2)
    assert (tmp_path / "timeseries-part-000002.parquet").exists()


def test_allow_uncommitted_chunk_after_complete_checkpoint(tmp_path):
    for i in (0, 1, 2):
        (tmp_path / f"timeseries-part-{i:06d}.parquet").write_bytes(b"chunk")
    validate_committed_parts(tmp_path, 2)


def test_inventory_covers_owned_formats_without_reading_artifact_contents(tmp_path, monkeypatch):
    paths = [
        "raw/timeseries-part-000000.parquet",
        "cache/aligned-part-000000.npz",
        "cache/fused_membrane/bilayer_normal/part-000000.parquet",
        "cache/fused_extensions/sasa__isolated/part-000000.parquet",
        "cache/fused_distance/run1/part-000000.npz",
        "cache/sdc4_mechanisms/map-000000.npz",
        "cache/solvation/bridges-000000.parquet",
        "cache/lipid_contacts/matrix-000000.parquet",
        "profiles/density_1d/run1/frames/protein__part-000000.npz",
        "profiles/density_1d/run1/frames/123456__part-000000.npz",
    ]
    unrelated = ["raw/unrelated-part-999999.parquet", "cache/unrelated-999999.npz",
                 "profiles/density_1d/run1/123456.parquet"]
    for name in paths + unrelated:
        file = tmp_path/name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes(b"metadata only, not a valid data container")
    def forbidden(*a, **k):
        raise AssertionError("checkpoint inventory read scientific data contents")
    monkeypatch.setattr(Path, "read_bytes", forbidden)
    records = chunk_artifacts(tmp_path, 0)
    assert set(records) == set(paths)
    for name in paths:
        (tmp_path/name.replace("000000", "000001")).write_bytes(b"orphan")
    temporary = tmp_path/"cache/sdc4_mechanisms/map-000001.tmp.npz"
    temporary.write_bytes(b"aborted writer")
    validate_artifact_inventory(tmp_path, records, 1)
    prune_artifact_orphans(tmp_path, 1)
    assert not temporary.exists()
    assert all((tmp_path/name).is_file() for name in paths + unrelated)
    assert not any((tmp_path/name.replace("000000", "000001")).exists() for name in paths)
