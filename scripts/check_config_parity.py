"""Compare YAML resolution, config hashes and CLI JSON plans against main.

This is an integration check, not a golden-file mock: the main-branch source is
archived and executed in a distinct Python process using the same inputs.
It requires a Git checkout with the remote main ref available or fetchable.
"""
from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parents[1]
CASES = {
    "minimal": {"systems": {"one": {
        "topology": "example.pdb", "replicates": [{"trajectory": "rep1.xtc", "label": "rep1"}]
    }}},
    "membrane_order": {
        "systems": {"one": {"topology": "example.pdb", "replicates": ["rep1.xtc"]}},
        "analyses": {"membrane_order": {"enabled": True, "tail_preset": "charmm36",
                                      "lipid_selection": "resname POPC",
                                      "leaflet_anchor_selection": "name P"}},
    },
    "density": {
        "systems": {"one": {"topology": "example.pdb", "replicates": ["rep1.xtc"]}},
        "analyses": {"density_1d": {"enabled": True, "selections": {"protein": "protein"}}},
    },
}

COLLECT = r"""
import json
import subprocess
import sys
from mdscope.config import load_config
from mdscope.core.provenance import config_hash
from mdscope.planning import build_plan

results = {}
for name, path in json.load(sys.stdin).items():
    cfg = load_config(path, check_paths=False)
    expected = build_plan(cfg, path, validation_mode="schema_only")
    process = subprocess.run(
        [sys.executable, "-m", "mdscope.cli", "plan", path, "--format", "json"],
        capture_output=True, text=True, check=True,
    )
    actual = json.loads(process.stdout)
    if actual != expected:
        raise AssertionError(f"{name}: CLI and planning API differ")
    results[name] = {"resolved": cfg, "hash": config_hash(cfg), "plan": actual}
print(json.dumps(results, sort_keys=True))
"""


def snapshot(source: Path, inputs: dict[str, str]) -> dict:
    env = {**os.environ, "PYTHONPATH": str(source / "src")}
    process = subprocess.run(
        [sys.executable, "-c", COLLECT],
        input=json.dumps(inputs), text=True, capture_output=True,
        cwd=REPO, env=env, check=False,
    )
    if process.returncode:
        raise RuntimeError(f"{source}: parity probe failed:\n{process.stderr}")
    return json.loads(process.stdout)


def main() -> int:
    maintained = sorted(REPO.glob("synd4_*.yaml"))
    if not maintained:
        raise RuntimeError("No maintained synd4_*.yaml files found")
    with tempfile.TemporaryDirectory(prefix="mdscope-main-parity-") as scratch:
        scratch = Path(scratch)
        fixture_dir = scratch / "fixture_cases"
        fixture_dir.mkdir()
        inputs = {p.name: str(p) for p in maintained}
        for name, raw in CASES.items():
            path = fixture_dir / f"{name}.yaml"
            path.write_text(yaml.safe_dump(raw, sort_keys=True), encoding="utf-8")
            inputs[f"fixture:{name}"] = str(path)

        subprocess.run(
            ["git", "fetch", "--no-tags", "--depth=1", "origin", "main"],
            cwd=REPO, check=True, capture_output=True, text=True,
        )
        archive = subprocess.run(
            ["git", "archive", "FETCH_HEAD"], cwd=REPO,
            check=True, capture_output=True,
        )
        baseline = scratch / "main"
        baseline.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive.stdout), mode="r:") as tar:
            tar.extractall(baseline, filter="data")

        before = snapshot(baseline, inputs)
        after = snapshot(REPO, inputs)
        problems = []
        for name in sorted(inputs):
            for key in ("resolved", "hash", "plan"):
                if before[name][key] != after[name][key]:
                    problems.append(f"{name}: {key} differs from main")
            print(f"{name}: {'MATCH' if not any(p.startswith(name + ':') for p in problems) else 'DIFF'}")
        if problems:
            print("\n".join(problems), file=sys.stderr)
            return 1
        print(f"Parity confirmed for {len(inputs)} configurations against main")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
