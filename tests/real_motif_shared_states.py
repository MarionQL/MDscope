"""Release TPR/XTC parity for shared motif hydration and generic state inputs.

Bounded smoke evidence, not evidence of state convergence or metastability.
Run with MDSCOPE_REAL_DATA and MDSCOPE_REAL_OUTPUT, as the real smoke workflow.
"""
from copy import deepcopy
import json
import os
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd

from mdscope.config import atomic_json
from mdscope.core.motifs import resolve_motifs
from mdscope.observables.solvation import SolvationObserver
from mdscope.statistics.shared_states import fit_shared_states


def run():
    data = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data")).resolve()
    root = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve() / "motif_states"
    sequence = "DDEDVVGPGQESDDFELSGSGDLDDLED"
    selection = "protein and name CA and chainID PROA PROB"
    sites = "name OD1 OD2 OE1 OE2"
    legacy = {"analyses": {"solvation": {"enabled": True, "preset": "standard",
        "water_selection": "resname TIP3 and name OH2",
        "targets": {"motif": {"acidic_patch": {"sequence": sequence,
            "selection": selection, "atom_selection": sites}}}}}}
    shared = deepcopy(legacy)
    shared["selections"] = {"motifs": {"acidic_patch": {"sequence": sequence,
        "selection": selection, "chains": ["PROA", "PROB"], "overlap": True}}}
    shared["analyses"]["solvation"]["targets"]["motif"]["acidic_patch"] = {
        "motif": "acidic_patch", "atom_selection": sites}
    paths, audits = {}, {}
    for system in ("carb", "nocarb"):
        u = mda.Universe(str(data / f"{system}-sdc4-step7_production.tpr"),
                         str(data / f"{system}-sdc4-rep1_test.xtc"))
        audit = resolve_motifs(u, shared["selections"]["motifs"])["acidic_patch"]
        assert set(audit) == {"acidic_patch_PROA_1", "acidic_patch_PROB_1"}
        outputs = [root / mode / system / "rep1" for mode in ("legacy", "shared")]
        observers = [SolvationObserver(u, cfg, out) for cfg, out in zip((legacy, shared), outputs)]
        for before, after in zip(observers[0].targets, observers[1].targets):
            assert before[:3] == after[:3]
            np.testing.assert_array_equal(before[3].indices, after[3].indices)
        assert len(observers[0].targets) == len(observers[1].targets) == 2
        for frame, ts in enumerate(u.trajectory[:20]):
            for observer in observers:
                observer.update(ts, frame, float(ts.time))
        for observer in observers:
            observer.flush(0)
            observer.finalize()
        for product in ("counts", "species_timeseries", "rdf"):
            tables = [pd.read_parquet(out / f"extensions/solvation/{product}.parquet") for out in outputs]
            pd.testing.assert_frame_equal(*tables, check_exact=True)
        paths[system] = [outputs[1]]
        audits[system] = audit
    cfg = {"analyses": {"shared_states": {"enabled": True, "preset": None,
        "features": ["patch_water"], "n_states": 2, "bootstrap_draws": 20,
        "feature_sources": {"patch_water": {"source": "solvation",
            "group": "solvation/water/motif/acidic_patch", "metric": "mean_n_solvent",
            "units": "molecules"}}}}}
    out = fit_shared_states(cfg, root / "shared", paths)
    assignments = pd.read_parquet(out / "frame_assignments.parquet")
    assert len(assignments) == 40
    model = json.loads((out / "model.json").read_text())
    assert model["preset"] is None and model["n_replicates"] == 2
    assert not (out / "state_conditioned_nxip.parquet").exists()
    atomic_json(root / "validation.json", {"systems": audits, "frames_per_system": 20,
        "exact_parity_products": ["counts", "species_timeseries", "rdf"],
        "generic_state_frames": len(assignments),
        "limitation": "Short smoke fixture; no equilibrium, kinetic or causal inference."})
    print("SHARED MOTIF/STATE REAL PARITY: PASS", flush=True)


if __name__ == "__main__":
    run()
