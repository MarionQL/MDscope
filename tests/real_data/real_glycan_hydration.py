"""Release TPR/XTC parity: SDC4 hs_auto vs an equivalent generic glycan definition.

Checks the six bond-verified HS attachments and that the generic consumer of
the fused hydration observer reproduces the SDC4 per-sugar and union counts.
Bounded smoke evidence; run via `scripts/run_real_data_checks.sh smoke`.
"""
import os
from pathlib import Path

import MDAnalysis as mda
import pandas as pd
import yaml

from mdscope.core.persistence import atomic_json
from mdscope.core.glycans import discover_covalent_glycans
from mdscope.observables.hydration import OUTPUT, HydrationObserver

FRAMES = 20
ATTACHMENTS = {"PROA_S39", "PROA_S61", "PROA_S63", "PROB_S309", "PROB_S331", "PROB_S333"}


def run():
    data = Path(os.environ.get("MDSCOPE_REAL_DATA", "examples/data")).resolve()
    root = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve() / "glycan_hydration"
    production = yaml.safe_load(Path("synd4_full_mechanisms.yaml").read_text())["analyses"]["sdc4_mechanisms"]
    hs = {key: production["hs_auto"][key] for key in
          ("glycan_selection", "attachment_selection", "min_glycan_residues", "expected_chains")}
    water = production["water_oxygen_selection"]
    cutoff = float(production["hs_sugar_water_cutoff_A"])
    configs = {
        "sdc4": {"analyses": {"sdc4_mechanisms": {
            "hs_sugar_water_cutoff_A": cutoff, "water_oxygen_selection": water,
            "hs_auto": {"enabled": True, **hs}}}},
        "generic": {"selections": {"glycans": {"heparan_sulfate": hs}},
                    "analyses": {"hydration": {
                        "enabled": True, "cutoff_A": cutoff, "water_oxygen_selection": water,
                        "protein_selection": "protein and resid 1:5", "glycans": ["heparan_sulfate"]}}},
    }
    u = mda.Universe(str(data / "carb-sdc4-step7_production.tpr"),
                     str(data / "carb-sdc4-rep1_test.xtc"))
    _, _, audit = discover_covalent_glycans(u, hs)
    assert set(audit) == ATTACHMENTS, sorted(audit)
    assert {item["protein_atom"] for item in audit.values()} == {"OG"}
    observers = {name: HydrationObserver(u, cfg, root / name) for name, cfg in configs.items()}
    for frame, ts in enumerate(u.trajectory[:FRAMES]):
        for observer in observers.values():
            observer.update(ts, frame, float(ts.time))
    for observer in observers.values():
        observer.flush(0)
        observer.finalize()
    tables = {}
    for name, compartment in (("sdc4", "hs"), ("generic", "heparan_sulfate")):
        directory = root / name / "extensions/hydration"
        tables[name] = {
            product: pd.read_parquet(directory / f"{product}.parquet")
            .query("compartment == @compartment").drop(columns="compartment").reset_index(drop=True)
            for product in ("counts", "aggregate_counts", "target_metadata", "glycan_attachments")}
    for product in tables["sdc4"]:
        pd.testing.assert_frame_equal(tables["sdc4"][product], tables["generic"][product], check_exact=True)
    sugars = pd.read_parquet(root / "sdc4" / OUTPUT)
    assert set(sugars.hs_chain) == ATTACHMENTS and sugars.time_ps.nunique() == FRAMES
    atomic_json(root / "validation.json", {
        "attachments": audit, "frames": FRAMES,
        "sugars": int(len(tables["sdc4"]["target_metadata"])),
        "exact_parity_products": list(tables["sdc4"]),
        "limitation": "Short smoke fixture of one HS topology; no other glycan chemistry."})
    print(f"GLYCAN HYDRATION REAL PARITY: PASS ({len(tables['sdc4']['target_metadata'])} sugars, "
          f"{FRAMES} frames)", flush=True)


if __name__ == "__main__":
    run()
