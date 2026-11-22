"""Write a frame-strided copy of the real GROMACS test data for quick local checks.

    python scripts/stride_real_data.py --stride 5 [--source examples/data]
                                       [--dest examples/data/stride5]

Sources must match tests/data/real_test_data.sha256. Every ``stride``-th frame
(starting with frame 0) is rewritten with its original time and box; TPRs are
copied unchanged. The destination gets ``provenance.json`` and ``SHA256SUMS``.
Strided data are for smoke/parity iteration only: time-dependent quantities at
the coarser spacing are not scientific results, and CI-equivalent evidence still
uses the full files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tests" / "data" / "real_test_data.sha256"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest() -> dict[str, str]:
    entries = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            entries[name.lstrip("*")] = digest
    return entries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stride", type=int, required=True)
    parser.add_argument("--source", type=Path, default=ROOT / "examples" / "data")
    parser.add_argument("--dest", type=Path)
    args = parser.parse_args()
    if args.stride < 2:
        parser.error("--stride must be at least 2")
    dest = args.dest or args.source / f"stride{args.stride}"
    if dest.exists():
        parser.error(f"{dest} already exists; remove it to regenerate")

    import MDAnalysis as mda

    expected = read_manifest()
    for name, digest in expected.items():
        actual = sha256(args.source / name)
        if actual != digest:
            raise SystemExit(f"{args.source / name}: SHA-256 {actual} != manifest {digest}")
    print(f"Verified {len(expected)} source files against {MANIFEST.relative_to(ROOT)}")

    dest.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f".{dest.name}-", dir=dest.parent))
    try:
        systems = {}
        for trajectory in sorted(n for n in expected if n.endswith(".xtc")):
            condition = trajectory.split("-", 1)[0]
            topology = next(n for n in expected if n.startswith(f"{condition}-") and n.endswith(".tpr"))
            shutil.copy2(args.source / topology, work / topology)
            universe = mda.Universe(str(args.source / topology), str(args.source / trajectory))
            times = []
            with mda.Writer(str(work / trajectory), n_atoms=len(universe.atoms)) as writer:
                for ts in universe.trajectory[::args.stride]:
                    writer.write(universe.atoms)
                    times.append(float(ts.time))
            systems[condition] = {
                "topology": topology, "trajectory": trajectory,
                "source_frames": universe.trajectory.n_frames, "frames": len(times),
                "first_time_ps": times[0], "last_time_ps": times[-1],
                "frame_spacing_ps": times[1] - times[0] if len(times) > 1 else None,
            }
            print(f"{trajectory}: {universe.trajectory.n_frames} -> {len(times)} frames "
                  f"({times[0]:.0f}-{times[-1]:.0f} ps)")
        outputs = {p.name: sha256(p) for p in sorted(work.iterdir())}
        (work / "SHA256SUMS").write_text(
            "".join(f"{digest}  {name}\n" for name, digest in outputs.items()), encoding="utf-8")
        (work / "provenance.json").write_text(json.dumps({
            "derived_from": {name: expected[name] for name in sorted(expected)},
            "stride": args.stride, "first_frame": 0,
            "mdanalysis": mda.__version__, "systems": systems, "outputs": outputs,
            "note": "Strided copy for quick local smoke checks; not CI-equivalent evidence.",
        }, indent=2) + "\n", encoding="utf-8")
        work.rename(dest)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    print(f"Wrote {dest}")


if __name__ == "__main__":
    main()
