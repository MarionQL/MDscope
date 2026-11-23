"""Atomic filesystem output helpers used by analysis and reporting.

Preserves the established MDscope on-disk formats and replacement semantics.
"""
from __future__ import annotations

import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def archive_output(source: Path, archive_root: Path, label: str) -> Path:
    """Move previous products aside without deleting scientific results."""
    ensure_dir(archive_root)
    candidate = archive_root / label
    index = 1
    while candidate.exists():
        candidate = archive_root / f"{label}-{index}"
        index += 1
    source.rename(candidate)
    return candidate


@contextmanager
def publish_directory(generated: Path, destination: Path, archive_root: Path, label: str):
    """Publish a generation and roll back if replacement or its commit fails.

    Both directories must be on the same filesystem. The two renames are not
    a cross-process atomic swap; callers must serialize writers.
    """
    if not generated.is_dir():
        raise FileNotFoundError(f"Prepared generation does not exist: {generated}")
    ensure_dir(destination.parent)
    previous = archive_output(destination, archive_root, label) if destination.exists() else None
    published = False
    try:
        generated.rename(destination)
        published = True
        yield destination
    except BaseException:
        # Also restore on cancellation: a dependent commit must succeed before
        # the replacement can outlive this transaction.
        if published:
            destination.rename(generated)
        if previous is not None:
            previous.rename(destination)
        raise


def atomic_json(path: str | Path, obj: Any) -> None:
    # Write to a sibling temporary file first so an interrupted job
    # never leaves a half-written state file.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp = tempfile.mkstemp(prefix=path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(obj, handle, indent=2, default=str)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_text(path: str | Path, text: str) -> None:
    # Atomic replacement matters because resolved configuration is also
    # used to audit/restart long runs. UTF-8 keeps reports (e.g. findings.md)
    # identical across platforms instead of using the locale code page.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp = tempfile.mkstemp(prefix=path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_parquet(path: str | Path, df: pd.DataFrame) -> None:
    # Use a unique temporary file so concurrent/aborted writers cannot
    # collide, then atomically replace the destination.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        try:
            df.to_parquet(tmp, index=False)
        except ImportError as exc:
            raise RuntimeError(
                "Parquet output requires pyarrow (recommended) or "
                "fastparquet."
            ) from exc
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def atomic_npz(path: str | Path, **arrays) -> None:
    """Publish a compressed array archive without exposing partial output."""
    import numpy as np

    path = Path(path)
    ensure_dir(path.parent)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            np.savez_compressed(handle, **arrays)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_parquet_parts(path: str | Path, parts, *, fill_missing=None) -> int:
    """Atomically join saved chunks using bounded Arrow record batches.

    Missing columns retain nulls unless the caller explicitly supplies a value
    (e.g. zero for a measured contact class absent from a particular chunk).
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    path = Path(path)
    parts = list(parts)
    if not parts:
        atomic_parquet(path, pd.DataFrame())
        return 0
    schema = pa.unify_schemas([pq.read_schema(part) for part in parts],
                             promote_options="permissive")
    ensure_dir(path.parent)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    rows = 0
    try:
        with pq.ParquetWriter(temporary, schema) as writer:
            for part in parts:
                for batch in pq.ParquetFile(part).iter_batches(batch_size=65536):
                    columns = [
                        batch.column(batch.schema.get_field_index(field.name)).cast(field.type)
                        if field.name in batch.schema.names else
                        pa.nulls(batch.num_rows, type=field.type) if fill_missing is None else
                        pa.array([fill_missing] * batch.num_rows, type=field.type)
                        for field in schema
                    ]
                    writer.write_batch(pa.RecordBatch.from_arrays(columns, schema=schema))
                    rows += batch.num_rows
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return rows
