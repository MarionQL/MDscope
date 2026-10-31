"""Stable configuration hashing and runtime environment metadata."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from typing import Any

def config_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def environment_metadata() -> dict[str, Any]:
    import MDAnalysis as mda
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        commit = None
    return {
        "created_unix": time.time(),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "mdanalysis": mda.__version__,
        "git_commit": commit,
    }


