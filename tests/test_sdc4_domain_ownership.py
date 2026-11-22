"""Domain boundaries after moving system-specific SDC4 execution (#258)."""
import importlib.util
import subprocess
import sys


def test_sdc4_analyses_and_generic_hydration_have_distinct_owners():
    from mdscope.systems.sdc4.mechanisms import SDC4Mechanisms
    from mdscope.systems.sdc4.potassium_bridge_audit import topology_sites
    from mdscope.observables.hydration import HydrationObserver
    from mdscope.contracts.outputs import HYDRATION_COUNT_COLUMNS

    assert callable(topology_sites)
    assert SDC4Mechanisms.__module__ == "mdscope.systems.sdc4.mechanisms"
    assert HydrationObserver.COUNT_COLUMNS is HYDRATION_COUNT_COLUMNS


def test_superseded_import_paths_are_absent():
    code = """
import importlib.util
for old in ("mdscope.mechanisms.sdc4_mechanisms",
            "mdscope.mechanisms.potassium_bridge_audit",
            "mdscope.mechanisms.hs_hydration"):
    assert importlib.util.find_spec(old) is None, old
from mdscope.observables.hydration import HydrationObserver
from mdscope.systems.sdc4.mechanisms import SDC4Mechanisms
"""
    p = subprocess.run([sys.executable, "-c", code], check=False,
                       capture_output=True, text=True, timeout=45)
    assert p.returncode == 0, p.stderr
