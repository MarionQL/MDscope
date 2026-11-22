"""Stable configuration hashing and runtime environment metadata."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path
from typing import Any


def config_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _source_identity() -> dict[str, Any]:
    """Never confuse the caller's checkout with MDscope's own source tree."""
    root = Path(__file__).resolve().parents[3]
    if not (root / ".git").exists():
        return {"git_commit": None, "git_dirty": None,
                "git_unavailable_reason": "MDscope is installed outside a git checkout"}
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        changes = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            text=True, stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        return {"git_commit": None, "git_dirty": None,
                "git_unavailable_reason": f"Cannot inspect MDscope source checkout: {type(exc).__name__}"}
    return {"git_commit": commit, "git_dirty": bool(changes.strip()),
            "git_unavailable_reason": None}


def _installed_version(distribution: str) -> str | None:
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return None


def environment_metadata() -> dict[str, Any]:
    versions = {
        key: _installed_version(distribution)
        for key, distribution in (
            ("mdscope", "mdscope"), ("mdanalysis", "MDAnalysis"),
            ("numpy", "numpy"), ("pandas", "pandas"), ("pyarrow", "pyarrow"),
            ("pyyaml", "PyYAML"), ("matplotlib", "matplotlib"),
            ("freesasa", "freesasa"),
        )
    }
    return {
        "created_unix": time.time(),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "dependency_versions": versions,
        "mdanalysis": versions["mdanalysis"],  # Retain legacy metadata consumer.
        **_source_identity(),
    }
