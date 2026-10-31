"""Atomic filesystem output helpers used by analysis and reporting.

Preserves the established MDscope on-disk formats and replacement semantics.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any
import pandas as pd

def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


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
    # used to audit/restart long runs.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp = tempfile.mkstemp(prefix=path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
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


