"""Compare YAML resolution, config hashes and CLI JSON plans against main.

This is an integration check, not a golden-file mock: the main-branch source is
archived and executed in a distinct Python process using the same inputs.
It requires a Git checkout with the remote main ref available or fetchable.
"""
from __future__ import annotations

import io
import copy
import hashlib
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



def legacy_plan_projection(historical, current):
    """Validate the additive #298 report descriptors without waiving old fields.

    Exact source IDs and legacy destinations stay public. New destinations are
    additive metadata for the already-existing PDF switch, not changed science.
    """
    additions = {
        'consolidated_destination': 'project_summary/figures/statistical_reports/{registered_analysis}/statistics_tables.pdf',
        'csv_companion': 'statistics_tables.csv',
        'selection_semantics': 'Exact source rows within one analysis report; legacy PNG destination applies only with table_reports_pdf: false',
    }
    if isinstance(historical, dict) and isinstance(current, dict):
        result = {}
        for key, value in current.items():
            if key not in historical and key in additions:
                if value != additions[key] or not current.get('id', '').endswith((':bootstrap_table', ':comparison_table')):
                    raise ValueError(f'Unexpected additive reporting metadata: {key}={value!r}')
            else:
                result[key] = legacy_plan_projection(historical.get(key), value)
        return result
    if isinstance(historical, list) and isinstance(current, list) and len(historical) == len(current):
        return [legacy_plan_projection(before, after) for before, after in zip(historical, current)]
    return current

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

        # Old code must consume the *main-branch* YAML, not changed YAML from
        # the PR: new user opt-ins are intentionally unknown to main. Feed the
        # same archived historical config to both versions, at one stable path.
        legacy = scratch / "legacy_configs"
        legacy.mkdir()
        for path in maintained:
            historical = baseline / path.name
            if not historical.is_file():
                raise RuntimeError(f"Missing maintained main-branch config: {path.name}")
            legacy_path = legacy / path.name
            legacy_path.write_bytes(historical.read_bytes())
            inputs[path.name] = str(legacy_path)

        before = snapshot(baseline, inputs)
        after = snapshot(REPO, inputs)

        # Newly introduced disabled-by-default opt-ins may add resolved keys
        # and therefore change the ordinary full config hash. Compare the
        # preexisting scientific contract exactly; do not waive changed values.
        additive_defaults = {
            ("analyses", "sdc4_mechanisms", "interprotomer_contacts"): False,
            ("analyses", "sdc4_mechanisms", "tm_local_tilt"): False,
            ("analyses", "sdc4_mechanisms", "ceramide_scaffold_pooling"): False,
            ("statistics", "residue_correlations"): {
                "enabled": False, "pairs": None, "min_frames": 30,
            },
        }
        problems = []
        for name in sorted(inputs):
            historical = before[name]["resolved"]
            current = copy.deepcopy(after[name]["resolved"])
            for path, expected in additive_defaults.items():
                base_node, current_node = historical, current
                for part in path[:-1]:
                    base_node = base_node.get(part, {})
                    current_node = current_node.get(part, {})
                field = path[-1]
                if field not in base_node and field in current_node:
                    if current_node[field] != expected:
                        problems.append(f"{name}: nondefault new setting {'.'.join(path)}")
                    else:
                        current_node.pop(field)
            if historical != current:
                problems.append(f"{name}: legacy resolved science differs from main")
            # The normalized legacy hash must remain identical to main even
            # though a new project records new disabled options in its hash.
            legacy_hash = hashlib.sha256(
                json.dumps(current, sort_keys=True, default=str).encode()
            ).hexdigest()[:16]
            if before[name]["hash"] != legacy_hash:
                problems.append(f"{name}: legacy config hash differs from main")
            if before[name]["plan"] != legacy_plan_projection(before[name]["plan"], after[name]["plan"]):
                problems.append(f"{name}: legacy plan differs from main")
            print(f"{name}: {'MATCH' if not any(p.startswith(name + ':') for p in problems) else 'DIFF'}")

        # The modified maintained YAML must also resolve via the *new* code's
        # CLI and planning API, although it cannot be read by historical code.
        modern_inputs = {path.name: str(path) for path in maintained}
        snapshot(REPO, modern_inputs)
        if problems:
            print("\n".join(problems), file=sys.stderr)
            return 1
        print(f"Legacy parity confirmed for {len(inputs)} configs; "
              f"{len(modern_inputs)} current YAML files also validated")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
