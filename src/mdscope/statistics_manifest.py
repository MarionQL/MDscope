"""Read-only statistics manifest facts shared by saved-data consumers.

Publication requires an explicit completed lifecycle. Legacy readers may
accept a missing manifest or missing method metadata, but cannot silently
promote a manifest with no completion status to a completed publication.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class StatisticsManifest:
    path: Path
    state: str
    payload: dict | None = None
    error: str | None = None

    @property
    def status(self):
        return self.payload.get("status") if self.payload is not None else None

    @property
    def methods_present(self):
        return self.payload is not None and "methods" in self.payload

    def family_disabled(self, family: str) -> bool:
        return (self.payload is not None and
                self.payload.get("families", {}).get(family) == "disabled")

    @property
    def complete(self) -> bool:
        return self.state == "complete"


def read_statistics_manifest(path: str | Path) -> StatisticsManifest:
    """Separate parsing/shape/lifecycle facts from consumer admission policy."""
    path = Path(path)
    if not path.is_file():
        return StatisticsManifest(path, "missing")
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError, TypeError) as exc:
        return StatisticsManifest(path, "invalid", error=str(exc))
    if not isinstance(payload, dict):
        return StatisticsManifest(path, "invalid", error="expected a JSON object")
    families = payload.get("families", {})
    if not isinstance(families, dict) or any(value not in ('enabled', 'disabled') for value in families.values()):
        return StatisticsManifest(path, "invalid", error="families must be an object")
    if "methods" in payload and not isinstance(payload["methods"], dict):
        return StatisticsManifest(path, "invalid", error="methods must be an object")
    if "status" in payload and not isinstance(payload['status'], str):
        return StatisticsManifest(path, "invalid", error="status must be a string")
    if "status" not in payload:
        return StatisticsManifest(path, "legacy", payload=payload)
    return StatisticsManifest(
        path, "complete" if payload["status"] == "complete" else "incomplete",
        payload=payload,
    )
