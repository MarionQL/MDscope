"""Execute the command block of docs/quickstart.md so the walkthrough stays current."""
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from mdscope_test_support import mdscope_executable

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT/"docs/quickstart.md"
DESTINATION = "../mdscope-quickstart"
# CI's checkout and editable install stand in for these two documented lines.
SETUP = [
    "git clone https://github.com/MarionQL/MDscope.git && cd MDscope",
    "python -m pip install -e '.[test]'",
]


def test_documented_quickstart_reaches_findings(tmp_path):
    text = DOC.read_text(encoding="utf-8")
    block = re.search(r"```bash\n(.*?)```", text, re.S).group(1).splitlines()
    assert block[:len(SETUP)] == SETUP
    destination = (tmp_path/"quickstart").as_posix()
    for line in block[len(SETUP):]:
        argv = [token.replace(DESTINATION, destination) for token in shlex.split(line)]
        if argv[0] == "python":
            argv[0] = sys.executable
        elif argv[0] == "mdscope":
            argv[0] = mdscope_executable()
        result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, timeout=300, check=False)
        assert result.returncode == 0, f"{line}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    assert "radius_of_gyration_A" in result.stdout

    output = Path(destination)
    documented = re.findall(r"^\| `(results/[^`]+)` \|", text, re.M)
    assert len(documented) == 5
    for relative in documented:
        assert (output/relative).is_file(), relative
    statistics = output/"results/project_summary/statistics"
    report = (statistics/"findings.md").read_text(encoding="utf-8")
    assert "synthetic_displaced minus synthetic_reference" in report

    # Independent closed form: chain separation 3 Å + offset + mean jitter 0.1 Å.
    var_z = np.var([36., 32., 28., 24., 20., 12., 8., 4., 0., -4., -8., -12.])
    expected = np.sqrt(var_z + (4.1/2)**2) - np.sqrt(var_z + (3.1/2)**2)
    findings = pd.read_parquet(statistics/"findings.parquet")
    metric = findings.observable.map(lambda value: json.loads(value)["metric"])
    assert set(metric) == {"radius_of_gyration_A"}
    assert set(findings.estimator) == {"mean", "median"}
    assert (findings.status == "above_zero").all()
    # 60 uniform jitters per condition: sampling error of the separation
    # difference is ~0.01 Å, i.e. ~6e-4 Å in radius of gyration.
    np.testing.assert_allclose(findings.difference_system2_minus_system1, expected, atol=2e-3)
    np.testing.assert_allclose(findings[["p_value", "q_value"]], [[.1, .2]]*2)

    coverage = pd.read_parquet(statistics/"coverage/contrast_coverage.parquet")
    rmsd = coverage[coverage.observable.str.contains('"rmsd_A"', na=False)]
    assert len(rmsd) == 2 and (rmsd.status == "includes_zero").all()
    profiles = coverage[coverage.family == "profiles"]
    assert len(profiles) and (profiles.status == "includes_zero").all()
