"""Validate full production YAML selectors against published real TPRs."""
from pathlib import Path
from types import SimpleNamespace
import os
import tempfile

import yaml

from mdscope.cli import cmd_validate

root = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data")).resolve()
original = yaml.safe_load(Path("synd4_full_mechanisms.yaml").read_text())
original["project"]["output"] = str(root / "validate_only")
original["trajectory"]["start"] = 0
original["trajectory"]["stride"] = 1
original["systems"] = {
    "synd4_glycolipids": {
        "label": "Synd4 + Glycolipids",
        "topology": str(root / "carb-sdc4-step7_production.tpr"),
        "replicates": [{
            "label": "rep1",
            "trajectory": str(root / "carb-sdc4-rep1_test.xtc"),
        }],
    },
    "synd4_no_glycolipids": {
        "label": "Synd4 - Glycolipids",
        "topology": str(root / "nocarb-sdc4-step7_production.tpr"),
        "replicates": [{
            "label": "rep1",
            "trajectory": str(root / "nocarb-sdc4-rep1_test.xtc"),
        }],
    },
}
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "validated_real.yaml"
    path.write_text(yaml.safe_dump(original, sort_keys=False))
    cmd_validate(SimpleNamespace(config=str(path)))
print("FULL CONFIG VALIDATION: PASS", flush=True)
