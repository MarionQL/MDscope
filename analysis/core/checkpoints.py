"""Minimal checkpoint integrity checks for streamed chunk files."""
from __future__ import annotations
from pathlib import Path
import re

_PART = re.compile(r"-part-(\d{6})\.(?:parquet|npz)$")


def validate_progress(processed: int, next_part: int, chunk: int, total: int) -> None:
    """Reject impossible or non-committed frame/part combinations."""
    if chunk < 1 or total < 0:
        raise ValueError("Invalid frame/chunk bounds")
    if processed < 0 or processed > total or next_part < 0:
        raise ValueError("Checkpoint frame/part counters are out of range")
    expected = (processed + chunk - 1) // chunk
    if next_part != expected:
        raise ValueError(
            f"Checkpoint part counter mismatch: {next_part} parts for "
            f"{processed} frames with chunk size {chunk}"
        )


def prune_uncommitted_parts(folder: Path, next_part: int) -> list[Path]:
    """Delete only known numbered uncommitted part artifacts in *folder*."""
    removed = []
    for path in Path(folder).iterdir():
        match = _PART.search(path.name)
        if match and int(match.group(1)) >= next_part and path.is_file():
            path.unlink()
            removed.append(path)
    return removed
