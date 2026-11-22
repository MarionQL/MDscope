"""Issue #169: YAML glycan definitions reuse fused unique-water hydration.

Non-SDC4 fixture: one protein chain ``A`` with an O-linked (Thr12 OG1) and an
N-linked (Asn40 ND2) linear disaccharide. Expected counts are hand-derived
from the coordinates below (cutoff 1.5 Å, cubic 50 Å box).
"""
from copy import deepcopy
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
import yaml
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.config import resolve_config
from mdscope.core.glycans import discover_covalent_glycans, glycan_sugars
from mdscope.pipeline.fused_extensions import FusedExtensions
from mdscope.registry import fused_output, required_artifacts

ROOT = Path(__file__).resolve().parents[1]
# (name, residue) per atom; residues: (resname, resid)
RESIDUES = [("THR", 12), ("ASN", 40), ("BGALNA", 101), ("BGAL", 102),
            ("BGLCNA", 201), ("BGLCNA", 202), ("TIP3", 301), ("TIP3", 302),
            ("TIP3", 303), ("BGAL", 103)]
ATOMS = [  # name, residue, position
    ("N", 0, [5, 40, 0]), ("OG1", 0, [7, 40, 0]),
    ("N", 1, [20, 40, 0]), ("ND2", 1, [22, 40, 0]),
    ("C1", 2, [9, 5, 0]), ("O6", 2, [10, 0, 0]),
    ("C1", 3, [12, 5, 0]), ("O4", 3, [12, 0, 0]),
    ("C1", 4, [30, 30, 0]), ("O6", 4, [30, 25, 0]),
    ("C1", 5, [1, 30, 0]), ("O4", 5, [0.5, 25, 0]),
    # One water shared by both O-glycan sugars; one across the x boundary.
    ("OH2", 6, [11, 0, 0]), ("OH2", 7, [49.7, 25, 0]), ("OH2", 8, [25, 10, 10]),
    ("C1", 9, [15, 15, 15]),
]
BONDS = [(1, 4), (4, 5), (5, 6), (6, 7), (3, 8), (8, 9), (9, 10), (10, 11)]


def _universe(extra_bonds=()):
    u = mda.Universe.empty(len(ATOMS), n_residues=len(RESIDUES),
                           atom_resindex=[a[1] for a in ATOMS], trajectory=True)
    u.add_TopologyAttr("names", [a[0] for a in ATOMS])
    u.add_TopologyAttr("resnames", [r[0] for r in RESIDUES])
    u.add_TopologyAttr("resids", [r[1] for r in RESIDUES])
    u.add_TopologyAttr("chainIDs", ["A"] * len(ATOMS))
    u.add_bonds(BONDS + list(extra_bonds))
    positions = np.asarray([a[2] for a in ATOMS], np.float32)[None]
    u.load_new(positions, format=MemoryReader,
               dimensions=np.array([[50., 50., 50., 90., 90., 90.]]))
    return u


O_GLYCAN = {"glycan_selection": "resname BGALNA BGAL and not resid 103",
            "attachment_selection": "resname THR and name OG1"}
N_GLYCAN = {"glycan_selection": "resname BGLCNA",
            "attachment_selection": "resname ASN and name ND2"}


def _cfg(glycans=("o_glycan", "n_glycan"), definitions=None):
    return {
        "selections": {"protein": "resname THR ASN",
                       "glycans": definitions or {"o_glycan": O_GLYCAN, "n_glycan": N_GLYCAN}},
        "trajectory": {"start": 0, "stop": 1, "stride": 1, "chunk_frames": 1},
        "time": {"origin_ps": 0., "dt_ps": 10.},
        "analyses": {"hydration": {
            "enabled": True, "cutoff_A": 1.5, "water_oxygen_selection": "name OH2",
            "glycans": list(glycans)}},
    }


def _run(root, cfg, universe=None):
    universe = universe or _universe()
    fused = FusedExtensions(universe, cfg, root)
    fused.prune_from(0)
    fused.update(universe.trajectory[0], 0)
    fused.flush(0)
    fused.finalize()
    return root / "extensions/hydration"


def _aggregate(table, compartment, level, group):
    rows = table.loc[table.compartment.eq(compartment) & table.level.eq(level)
                     & table.group_id.eq(group)]
    assert len(rows) == 1, (compartment, level, group)
    return int(rows.water_count.iloc[0])


def test_non_sdc4_glycoprotein_per_sugar_chain_and_total_hydration(tmp_path):
    out = _run(tmp_path, _cfg())
    meta = pd.read_parquet(out / "target_metadata.parquet")
    glycan = meta.loc[meta.compartment.isin(["o_glycan", "n_glycan"])]
    assert glycan[["compartment", "target_id", "chain", "sugar_position", "resname",
                   "resid", "topology_residue_index", "n_sites"]].values.tolist() == [
        ["o_glycan", "A_T12:1", "A_T12", 1, "BGALNA", 101, 2, 1],
        ["o_glycan", "A_T12:2", "A_T12", 2, "BGAL", 102, 3, 1],
        ["n_glycan", "A_N40:1", "A_N40", 1, "BGLCNA", 201, 4, 1],
        ["n_glycan", "A_N40:2", "A_N40", 2, "BGLCNA", 202, 5, 1],
    ]
    counts = pd.read_parquet(out / "counts.parquet").set_index(["compartment", "target_id"])
    # The shared water counts once for each sugar; the zero is a measured zero.
    assert [int(counts.loc[key, "water_count"]) for key in [
        ("o_glycan", "A_T12:1"), ("o_glycan", "A_T12:2"),
        ("n_glycan", "A_N40:1"), ("n_glycan", "A_N40:2")]] == [1, 1, 0, 1]
    aggregate = pd.read_parquet(out / "aggregate_counts.parquet")
    assert _aggregate(aggregate, "o_glycan", "chain", "A_T12") == 1
    assert _aggregate(aggregate, "o_glycan", "sugar_type", "A_T12:BGAL") == 1
    assert _aggregate(aggregate, "o_glycan", "total", "all") == 1
    assert _aggregate(aggregate, "n_glycan", "sugar_type", "A_N40:BGLCNA") == 1
    assert _aggregate(aggregate, "n_glycan", "total", "all") == 1
    assert _aggregate(aggregate, "protein", "total", "all") == 0
    attachments = pd.read_parquet(out / "glycan_attachments.parquet")
    assert attachments[["compartment", "chain", "protein_chain", "protein_resname",
                        "protein_resid", "protein_atom", "protein_atom_index",
                        "protein_topology_residue_index", "n_sugars",
                        "sugar_types"]].values.tolist() == [
        ["o_glycan", "A_T12", "A", "THR", 12, "OG1", 1, 0, 2, "BGALNA,BGAL"],
        ["n_glycan", "A_N40", "A", "ASN", 40, "ND2", 3, 1, 2, "BGLCNA,BGLCNA"],
    ]
    # Per-chain and total unions are saved scalar identities for statistics.
    from mdscope.statistics.time_series import collect_time_series
    series = collect_time_series(tmp_path)
    assert series[("hydration/n_glycan/chain/A_N40", "unique_water_count")].value.tolist() == [1]
    assert series[("hydration/o_glycan/total/all", "unique_water_per_site")].value.tolist() == [.5]
    entry = fused_output(tmp_path, "hydration")
    required = required_artifacts(entry, _cfg())
    assert out / "glycan_attachments.parquet" in required
    assert all(path.is_file() for path in required)
    # No SDC4 product is written for a generic glycoprotein.
    assert not (tmp_path / "extensions/sdc4_mechanisms").exists()


def test_one_definition_may_span_chemistries_without_merging_chains(tmp_path):
    both = {"glycan_selection": "resname BGALNA BGLCNA or (resname BGAL and not resid 103)",
            "attachment_selection": "(resname THR and name OG1) or (resname ASN and name ND2)",
            "expected_chains": 2}
    out = _run(tmp_path, _cfg(["glycans"], {"glycans": both}))
    aggregate = pd.read_parquet(out / "aggregate_counts.parquet")
    assert _aggregate(aggregate, "glycans", "chain", "A_T12") == 1
    assert _aggregate(aggregate, "glycans", "chain", "A_N40") == 1
    # Same residue name on two chains stays chain-scoped.
    assert set(aggregate.loc[aggregate.level.eq("sugar_type"), "group_id"]) == {
        "A_T12:BGALNA", "A_T12:BGAL", "A_N40:BGLCNA"}
    assert _aggregate(aggregate, "glycans", "total", "all") == 2


def test_generic_definition_reproduces_sdc4_hs_preset(tmp_path):
    from test_fused_hydration import _universe as hs_universe
    hs = {"glycan_selection": "resname BXYL BGAL", "attachment_selection": "name OG"}
    sdc4 = {"trajectory": {"start": 0, "stop": 2, "stride": 1, "chunk_frames": 1},
            "time": {"origin_ps": 10., "dt_ps": 5.},
            "analyses": {"sdc4_mechanisms": {"hs_sugar_water_cutoff_A": 1.4,
                                             "water_oxygen_selection": "name OW",
                                             "hs_auto": hs}}}
    generic = {"selections": {"protein": "resname SER ALA", "glycans": {"heparan": hs}},
               "trajectory": sdc4["trajectory"], "time": sdc4["time"],
               "analyses": {"hydration": {"enabled": True, "cutoff_A": 1.4,
                                          "water_oxygen_selection": "name OW",
                                          "protein_selection": "resname ALA",
                                          "glycans": ["heparan"]}}}
    tables = {}
    for name, cfg in (("sdc4", sdc4), ("generic", generic)):
        universe = hs_universe()
        fused = FusedExtensions(universe, cfg, tmp_path / name)
        for i, ts in enumerate(universe.trajectory):
            fused.update(ts, i)
            fused.flush(i)
        fused.finalize()
        counts = pd.read_parquet(tmp_path / name / "extensions/hydration/counts.parquet")
        aggregate = pd.read_parquet(tmp_path / name / "extensions/hydration/aggregate_counts.parquet")
        compartment = "hs" if name == "sdc4" else "heparan"
        tables[name] = [
            counts.loc[counts.compartment.eq(compartment)].drop(columns="compartment")
            .reset_index(drop=True),
            aggregate.loc[aggregate.compartment.eq(compartment)].drop(columns="compartment")
            .reset_index(drop=True)]
    for before, after in zip(tables["sdc4"], tables["generic"]):
        pd.testing.assert_frame_equal(before, after)
    assert tables["generic"][0].target_id.tolist()[:3] == ["PROA_S39:1", "PROA_S39:2", "PROA_S39:3"]


@pytest.mark.parametrize("change, message", [
    # Omitting the middle sugar leaves a bond to an unselected residue.
    ({"glycan_selection": "resname BGALNA"}, "outside glycan_selection"),
    ({"attachment_selection": "resname THR"}, "multiple attachment atoms"),
    ({"expected_chains": 2}, "Expected 2 covalently attached glycan chains"),
    ({"glycan_selection": "resname XYZ"}, "matched no glycan"),
    ({"attachment_selection": "resname SER and name OG"}, "matched no glycan or attachment"),
    ({"attachment_selection": "resname BGAL and name O4"}, "must not overlap"),
    ({"min_glycan_residues": 3}, "has 2 residues"),
])
def test_missing_or_ambiguous_linkage_fails_explicitly(change, message):
    # The ambiguous case bonds the sugar to a second Thr atom.
    universe = _universe(extra_bonds=[(0, 4)] if change.get("attachment_selection") == "resname THR" else ())
    with pytest.raises(ValueError, match=message):
        discover_covalent_glycans(universe, {**O_GLYCAN, **change})


def test_branched_glycan_and_unattached_fragment_fail():
    branched = _universe(extra_bonds=[(5, 15)])
    spec = {**O_GLYCAN, "glycan_selection": "resname BGALNA BGAL"}
    with pytest.raises(ValueError, match="branches at BGALNA 101"):
        glycan_sugars(branched, spec)
    with pytest.raises(ValueError, match="no covalent protein attachment"):
        discover_covalent_glycans(_universe(), spec)  # resid 103 is unbonded


def test_absent_chains_are_errors_not_zero_hydration(tmp_path):
    # The selected residue exists but has no covalent linkage to the protein.
    cfg = _cfg(["missing"], {"missing": {**O_GLYCAN, "glycan_selection": "resid 103"}})
    with pytest.raises(ValueError, match="no covalent protein attachment"):
        _run(tmp_path, cfg)
    assert not (tmp_path / "extensions/hydration").exists()


def test_overlapping_definitions_are_rejected(tmp_path):
    cfg = _cfg(["a", "b"], {"a": O_GLYCAN, "b": O_GLYCAN})
    with pytest.raises(ValueError, match="target selections overlap"):
        _run(tmp_path, cfg)


def _document(**hydration):
    from mdscope.studio.model import new_document
    doc = new_document()
    doc["selections"]["glycans"] = {"o_glycan": dict(O_GLYCAN)}
    doc["analyses"]["hydration"] = {"enabled": True, "water_oxygen_selection": "name OH2",
                                    **hydration}
    return doc


def test_glycan_yaml_roundtrip_plan_and_omitted_vs_empty(tmp_path):
    from mdscope.planning import build_plan
    from mdscope.studio.model import import_document
    doc = _document(glycans=["o_glycan"])
    resolved = resolve_config(doc, check_paths=False)
    assert resolved["selections"]["glycans"] == {"o_glycan": O_GLYCAN}
    assert import_document(yaml.safe_dump(doc))["selections"]["glycans"] == doc["selections"]["glycans"]
    plan = build_plan(resolved, tmp_path / "analysis.yaml")
    hydration = next(item for item in plan["analyses"] if item["analysis"] == "hydration")
    # Effective settings show the generic defaults: no inherited chain count.
    assert hydration["glycan_definitions"] == {"o_glycan": {
        **O_GLYCAN, "min_glycan_residues": 1, "expected_chains": None}}
    assert any(o["path_within_replica"].endswith("hydration/glycan_attachments.parquet")
               for o in plan["registered_replica_outputs"])
    from mdscope.registry import state_feature_catalog
    assert {"hydration/o_glycan/total/all"} == {
        c["group"] for c in state_feature_catalog(resolved)["candidates"] if c["source"] == "hydration"}
    # Omitted and explicit-empty consumer lists select no glycan targets or products.
    for doc in (_document(), _document(glycans=[])):
        resolved = resolve_config(doc, check_paths=False)
        assert resolved["analyses"]["hydration"]["glycans"] == []
        plan = build_plan(resolved, tmp_path / "analysis.yaml")
        assert not any("glycan_attachments" in o["path_within_replica"] for o in plan["registered_replica_outputs"])
    # Disabled hydration keeps the definition but plans no hydration work.
    disabled = _document(glycans=["o_glycan"])
    disabled["analyses"]["hydration"]["enabled"] = False
    plan = build_plan(resolve_config(disabled, check_paths=False), tmp_path / "analysis.yaml")
    assert not any(item["analysis"] == "hydration" for item in plan["analyses"])


@pytest.mark.parametrize("edit, message", [
    (lambda d: d["analyses"]["hydration"].update(glycans=["unknown"]), "hydration.glycans"),
    (lambda d: d["analyses"]["hydration"].update(glycans=None), "hydration.glycans"),
    (lambda d: d["analyses"]["hydration"].update(glycans=["o_glycan", "o_glycan"]), "hydration.glycans"),
    (lambda d: d["selections"].update(glycans={"hs": O_GLYCAN}), "names must"),
    (lambda d: d["selections"].update(glycans=[]), "selections.glycans must be a mapping"),
    (lambda d: d["selections"]["glycans"]["o_glycan"].update(min_glycan_residues=0), "min_glycan_residues"),
    (lambda d: d["selections"]["glycans"]["o_glycan"].update(expected_chains=0), "expected_chains"),
    (lambda d: d["selections"]["glycans"]["o_glycan"].pop("attachment_selection"), "attachment_selection"),
    (lambda d: d["selections"]["glycans"]["o_glycan"].update(attachment="x"), "Unknown YAML option"),
])
def test_invalid_glycan_choices_fail_at_yaml_boundary(edit, message):
    doc = _document(glycans=["o_glycan"])
    edit(doc)
    with pytest.raises(ValueError, match=message):
        resolve_config(doc, check_paths=False)


def test_maintained_sdc4_yaml_keeps_hs_preset_and_no_generic_glycans(tmp_path):
    from mdscope.config import load_config
    from mdscope.planning import build_plan
    for name in ("synd4_full_mechanisms.yaml", "synd4_full_mechanisms_finalize_fast.yaml",
                 "synd4_baseline.yaml"):
        cfg = load_config(ROOT / name, check_paths=False)
        assert cfg["selections"]["glycans"] == {}
        assert cfg["analyses"]["hydration"]["glycans"] == []
        plan = build_plan(cfg, ROOT / name)
        assert not any("glycan_attachments" in o["path_within_replica"] for o in plan["registered_replica_outputs"])
    hs = load_config(ROOT / "synd4_full_mechanisms.yaml", check_paths=False)["analyses"]["sdc4_mechanisms"]["hs_auto"]
    assert hs == {"enabled": True,
                  "glycan_selection": "resname AGLCNA AIDOA BGLCA BGAL BXYL and chainID PROA PROB",
                  "attachment_selection": "protein and resname SER and name OG",
                  "min_glycan_residues": 3, "expected_chains": 6}
    bad = deepcopy(hs)
    bad["min_glycan_residues"] = 1
    from mdscope.core.glycans import validate_glycan_spec
    with pytest.raises(ValueError, match=">= 2"):
        validate_glycan_spec(bad, "sdc4_mechanisms.hs_auto", min_residues=2)
