"""Resume safety: committed counters and orphan chunk pruning."""
from pathlib import Path
import pytest
from mdscope.analysis.core.checkpoints import validate_progress, prune_uncommitted_parts


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
