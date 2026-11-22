"""Configuration package boundaries and compatibility after #259."""
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

from mdscope import config
from mdscope.config import DEFAULTS, AnalysisRun, Replicate, analysis_runs, resolve_config


def test_package_exports_match_canonical_submodules():
    from mdscope.config.defaults import DEFAULTS as defined_defaults
    from mdscope.config.runs import AnalysisRun as defined_run
    from mdscope.config.runs import Replicate as defined_replica
    from mdscope.config.runs import analysis_runs as defined_runs
    from mdscope.config.validation import validate_config

    assert DEFAULTS is defined_defaults
    assert AnalysisRun is defined_run
    assert Replicate is defined_replica
    assert analysis_runs is defined_runs
    assert config.validate_config is validate_config


def test_resolved_config_stays_stable_and_unmodified():
    original = deepcopy(DEFAULTS)
    raw = {
        "systems": {"a": {"topology": "/tmp/example.pdb", "replicates": [
            {"trajectory": "/tmp/example.xtc", "label": "rep1"}]}},
        "analyses": {"membrane_order": {"enabled": True, "tail_preset": "charmm36",
                                       "lipid_selection": "resname POPC",
                                       "leaflet_anchor_selection": "name P"}},
    }
    first = resolve_config(raw, base=Path("/tmp"))
    again = resolve_config(raw, base=Path("/tmp"))
    from mdscope.core.provenance import config_hash

    assert first == again
    assert config_hash(first) == config_hash(again)
    assert DEFAULTS == original
    assert first["analyses"]["membrane_order"]["tails"]["sn1"]["carbon_pattern"] == r"^C3(\d+)$"


def test_core_profiling_can_import_without_config_cycle():
    code = (
        "import sys\n"
        "from mdscope.core.profiling import PhaseProfiler\n"
        "assert 'mdscope.config' not in sys.modules\n"
        "from mdscope.core.persistence import atomic_parquet, atomic_json\n"
        "import mdscope.config as config\n"
        "assert not hasattr(config, 'atomic_parquet')\n"
        "assert not hasattr(config, 'atomic_json')\n"
    )
    p = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, check=False, timeout=45)
    assert p.returncode == 0, p.stderr
