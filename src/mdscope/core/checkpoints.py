"""Minimal checkpoint integrity checks for streamed chunk files."""
from __future__ import annotations
from pathlib import Path
import re

_PART = re.compile(r"-part-(\d{6})\.(?:parquet|npz)$")
_CHUNK = re.compile(r"(?:^|-)(\d{6})\.(?:parquet|npz)$")
_NUMBERED = re.compile(r"(?:^|-)(\d{6})(?:[._-].*)?\.(?:parquet|npz)$")


def numbered_artifacts(repout: Path):
    """Enumerate owned chunk files, including malformed copies/temporary parts."""
    repout = Path(repout)
    for folder in (repout / "raw", repout / "cache", repout / "profiles/density_1d"):
        for path in folder.rglob("*"):
            if folder == repout / "raw" and not path.name.startswith((
                    "timeseries-part-", "group_timeseries-part-", "salt_bridges-part-", "ligand_contacts-part-")):
                continue
            if folder == repout / "cache" and path.parent == folder and not path.name.startswith("aligned-part-"):
                continue
            if folder == repout / "profiles/density_1d" and path.parent.name != "frames":
                continue
            # A configured target may itself be named '123456'. The writer's
            # final part marker, not a numeric label, identifies its chunk.
            match = _NUMBERED.search(path.name.rsplit("part-", 1)[-1])
            if match and path.is_file():
                yield path, int(match.group(1))


def chunk_artifacts(repout: Path, part: int) -> dict:
    """Fingerprint files actually written for one chunk, including sparse events.

    Only filesystem metadata is read, not trajectory/event contents. The
    inventory lives in the existing atomic NPZ checkpoint; no sidecar state.
    """
    records = {}
    for path, index in numbered_artifacts(repout):
        if index == part and _CHUNK.search(path.name):
            stat = path.stat()
            records[path.relative_to(repout).as_posix()] = [stat.st_size, stat.st_mtime_ns]
    return records


def validate_artifact_inventory(repout: Path, records: dict, next_part: int) -> None:
    """Reject missing/replaced committed artifacts before pruning or rewriting.

    Size/mtime fingerprints detect ordinary edits and damage, not bit rot or
    adversarial modifications preserving both metadata fields.
    """
    root = Path(repout).resolve()
    if not isinstance(records, dict) or not records:
        raise ValueError("Invalid committed artifact inventory")
    for relative, expected in records.items():
        path = root / relative
        match = _CHUNK.search(path.name)
        if (not isinstance(expected, list) or len(expected) != 2
                or any(type(v) is not int or v < 0 for v in expected)
                or match is None or int(match.group(1)) >= next_part
                or Path(relative).is_absolute() or ".." in Path(relative).parts
                or not path.resolve().is_relative_to(root)):
            raise ValueError(f"Invalid committed artifact record: {relative}")
        if not path.is_file():
            raise ValueError(f"Missing committed artifact: {relative}")
        stat = path.stat()
        if [stat.st_size, stat.st_mtime_ns] != expected:
            raise ValueError(f"Changed committed artifact: {relative}")
    current = {path.relative_to(root).as_posix() for path, part in numbered_artifacts(root)
               if part < next_part}
    extra = current - set(records)
    if extra:
        raise ValueError("Unrecorded committed artifact: " + ", ".join(sorted(extra)))


def prune_artifact_orphans(repout: Path, next_part: int) -> None:
    """Remove numbered uncommitted chunks only after inventory validation."""
    for path, part in numbered_artifacts(repout):
        if part >= next_part:
            path.unlink()


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
    known_prefixes = ("timeseries-", "group_timeseries-", "salt_bridges-",
                      "ligand_contacts-", "aligned-")
    for path in Path(folder).iterdir():
        if not path.name.startswith(known_prefixes):
            continue
        match = _PART.search(path.name)
        if match and int(match.group(1)) >= next_part and path.is_file():
            path.unlink()
            removed.append(path)
    return removed



def validate_committed_parts(folder: Path, next_part: int) -> None:
    """Require every committed base timeseries chunk before resuming.

    A missing checkpoint with existing chunks is also unsafe: the caller must
    not interpret persisted user output as disposable uncommitted data.
    """
    folder = Path(folder)
    existing = {
        int(match.group(1))
        for path in folder.glob("timeseries-part-*.parquet")
        if (match := _PART.search(path.name)) is not None
    }
    if next_part == 0 and existing:
        raise ValueError(
            "Checkpoint is absent or records zero parts but persisted "
            "timeseries chunks exist; existing output must be reviewed"
        )
    missing = set(range(next_part)) - existing
    if missing:
        raise ValueError(
            "Checkpoint references missing committed timeseries parts: "
            + ", ".join(f"{index:06d}" for index in sorted(missing))
        )
