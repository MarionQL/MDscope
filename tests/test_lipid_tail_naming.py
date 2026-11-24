"""Lipid-tail naming presets, explicit overrides and topology preflight (issue #165).

Two materially different all-atom POPC naming conventions share one geometry:

* CHARMM36 (as in the SDC4 release TPRs): sn-1 C31..C316 on glycerol C3,
  sn-2 C21..C218 on glycerol C2.
* Maciejewski-Rog (MacRog, OPLS-AA based; sequential atom numbering), POPC
  names from NMRlipids BilayerData ``Molecules/membrane/POPC/
  mappingPOPCMaciejewskiRog1.yaml`` at commit 93ef01f: sn-1 C36 (carbonyl),
  C38..C52; sn-2 C15 (carbonyl), C17..C33; glycerol C34/C13/C12.

Every C-H bond of carbon ``p`` in a chain has a chain-specific known angle to
z, so the analytic -S_CH profile identifies both the position and sn chain.
"""
from __future__ import annotations

import re
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
import yaml
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.config import DEFAULTS, resolve_config
from mdscope.core.provenance import config_hash
from mdscope.errors import ConfigError
from mdscope.observables.membrane_order import (
    _order_sites,
    _tail_spec,
    order_preflight,
    order_setup,
)
from mdscope.pipeline.fused_membrane import FusedMembrane
from mdscope.pipeline.preflight import validate_selections
from mdscope.planning import build_plan
from mdscope.systems.sdc4.mechanisms import _normal_order_setup

CHAIN_LENGTH = {"sn1": 16, "sn2": 18}  # palmitoyl, oleoyl
DOUBLE_BOND = {"sn2": (9, 10)}
MACROG_TAILS = {
    "sn1": {"carbon_names": ["C36"] + [f"C{36 + p}" for p in range(2, 17)],
            "positions": list(range(1, 17))},
    "sn2": {"carbon_names": ["C15"] + [f"C{15 + p}" for p in range(2, 19)],
            "positions": list(range(1, 19))},
}


# Four ectodomain and five GGIVG-containing TM residues per synthetic protomer.
# Lipid-order assertions still select only POPC atoms.
SDC4_DOMAINS = {
    "tm_chains": {
        "A": "protein and chainID A and resid 9:13 and name CA",
        "B": "protein and chainID B and resid 18:22 and name CA",
    },
    "ectodomains": {
        "A": "protein and chainID A and resid 5:8 and name CA",
        "B": "protein and chainID B and resid 14:17 and name CA",
    },
}


def _candidate_sites(universe, spec):
    return _order_sites(universe, spec)[0][0]


def _cos2(chain, position):
    return (0.05 if chain == "sn1" else 0.45) + 0.02 * position


def expected_minus_sch(chain):
    positions = np.arange(2, CHAIN_LENGTH[chain] + 1)
    return positions, -0.5 * (3 * np.array([_cos2(chain, p) for p in positions]) - 1)


def _names(convention, chain, position):
    """Carbon name and its hydrogen names for one acyl carbon."""
    n = CHAIN_LENGTH[chain]
    double = position in DOUBLE_BOND.get(chain, ())
    if convention == "charmm36":
        carbon = f"C{3 if chain == 'sn1' else 2}{position}"
        suffix = ("X", "Y", "Z") if chain == "sn1" else ("R", "S", "T")
        hydrogens = ([f"H{position}1"] if double
                     else [f"H{position}{suffix[0]}", f"H{position}{suffix[1]}"])
        if position == n:
            hydrogens.append(f"H{position}{suffix[2]}")
    else:
        number = (36 if chain == "sn1" else 15) + (position if position > 1 else 0)
        carbon = f"C{number}"
        hydrogens = [f"H{number}1"] if double else [f"H{number}1", f"H{number}2"]
        if position == n:
            hydrogens.append(f"H{number}3")
    return carbon, ([] if position == 1 else hydrogens)


def _glycerol(convention):
    """(sn3 carbon, sn2 carbon, sn1 carbon, phosphate O, P, ester O, carbonyl O per chain)."""
    if convention == "charmm36":
        return {"carbons": (("C1", ["HA", "HB"]), ("C2", ["HS"]), ("C3", ["HX", "HY"])),
                "OP": "O11", "P": "P", "sn2": ("O21", "O22"), "sn1": ("O31", "O32")}
    return {"carbons": (("C12", ["H121", "H122"]), ("C13", ["H131"]), ("C34", ["H341", "H342"])),
            "OP": "O11", "P": "P8", "sn2": ("O14", "O16"), "sn1": ("O35", "O37")}


class _Builder:
    def __init__(self):
        self.names, self.elements, self.resindex, self.xyz, self.bonds = [], [], [], [], []
        self.resnames, self.resids = [], []

    def residue(self, resname):
        self.resnames.append(resname)
        self.resids.append(len(self.resids) + 1)
        return len(self.resnames) - 1

    def atom(self, residue, name, element, xyz, bonded_to=None):
        self.names.append(name)
        self.elements.append(element)
        self.resindex.append(residue)
        self.xyz.append(xyz)
        index = len(self.names) - 1
        if bonded_to is not None:
            self.bonds.append((bonded_to, index))
        return index

    def universe(self, *, frames=5, bonds=True):
        u = mda.Universe.empty(len(self.names), n_residues=len(self.resnames),
                               atom_resindex=self.resindex, trajectory=True)
        u.add_TopologyAttr("names", self.names)
        u.add_TopologyAttr("elements", self.elements)
        u.add_TopologyAttr("types", self.elements)
        u.add_TopologyAttr("resnames", self.resnames)
        u.add_TopologyAttr("resids", self.resids)
        u.add_TopologyAttr("segids", ["MEMB"])
        u.add_TopologyAttr("chainIDs", ["M"] * len(self.names))
        if bonds:
            u.add_TopologyAttr("bonds", self.bonds)
        xyz = np.asarray(self.xyz, np.float32)
        u.load_new(np.repeat(xyz[None], frames, axis=0), format=MemoryReader,
                   dimensions=np.tile([60., 60., 100., 90., 90., 90.], (frames, 1)))
        return u


def _lipid(builder, convention, x, y, z_head, direction, *, hydrogens=True):
    names = _glycerol(convention)
    residue = builder.residue("POPC")
    p = builder.atom(residue, names["P"], "P", [x, y, z_head])
    op = builder.atom(residue, names["OP"], "O", [x, y, z_head + direction * 0.8], p)
    glycerol = []
    previous = op
    for k, (carbon, hs) in enumerate(names["carbons"]):
        c = builder.atom(residue, carbon, "C", [x - 1 + k, y + 1.5, z_head + direction * 1.5], previous)
        for j, h in enumerate(hs):
            if hydrogens:
                builder.atom(residue, h, "H", [x - 1 + k, y + 2.5 + j * 0.3, z_head + direction * 1.5], c)
        glycerol.append(c)
        previous = c
    for chain, offset, glycerol_carbon in (("sn1", -2.0, glycerol[2]), ("sn2", 2.0, glycerol[1])):
        ester, carbonyl = names[chain]
        o = builder.atom(residue, ester, "O", [x + offset, y, z_head + direction * 2.2], glycerol_carbon)
        previous = o
        for position in range(1, CHAIN_LENGTH[chain] + 1):
            carbon, hs = _names(convention, chain, position)
            centre = np.array([x + offset, y, z_head + direction * (2.5 + 1.25 * position)])
            c = builder.atom(residue, carbon, "C", centre, previous)
            if position == 1:
                builder.atom(residue, carbonyl, "O", centre + [0, 1.2, 0], c)
            cos = np.sqrt(_cos2(chain, position))
            sin = np.sqrt(1 - cos * cos)
            for phi, h in zip((0., np.pi, np.pi / 2), hs):
                if hydrogens:
                    vector = 1.09 * np.array([sin * np.cos(phi), sin * np.sin(phi), cos])
                    builder.atom(residue, h, "H", centre + vector, c)
            previous = c
    return residue


def bilayer(convention, *, hydrogens=True, extra=()):
    """Two lipids per leaflet; upper tails point down, lower tails point up."""
    builder = _Builder()
    for x, y, z, direction in ((10, 10, 70, -1), (30, 10, 70, -1),
                               (10, 40, 30, 1), (30, 40, 30, 1)):
        _lipid(builder, convention, x, y, z, direction, hydrogens=hydrogens)
    for add in extra:
        add(builder)
    return builder


def _psm(builder):
    """CHARMM36 sphingomyelin chains as in the SDC4 TPRs: C1S..C18S, C1F (amide
    carbonyl, no hydrogen)..C16F. Every C-H bond is perpendicular to z, so
    -S_CH = 0.5 at each mapped position."""
    residue = builder.residue("PSM")
    p = builder.atom(residue, "P", "P", [50, 25, 70])
    for prefix, count, x in (("S", 18, 48.), ("F", 16, 52.)):
        previous = p
        for k in range(1, count + 1):
            c = builder.atom(residue, f"C{k}{prefix}", "C", [x, 25, 68 - 1.25 * k], previous)
            if not (prefix == "F" and k == 1):
                builder.atom(residue, f"H{k}{prefix}", "H", [x + 1.09, 25, 68 - 1.25 * k], c)
            previous = c


def test_charmm36_maps_sphingolipid_chains_without_pooling_them_into_sn_chains(tmp_path):
    cfg = _resolved({"tail_preset": "charmm36", "lipid_selection": "resname POPC PSM",
                     "leaflet_anchor_selection": "name P", "per_leaflet": False})
    profiles, mapping = _run_order(bilayer("charmm36", extra=(_psm,)).universe(), cfg, tmp_path)
    for tail, first, last in (("sphingosine", 1, 18), ("n_acyl", 2, 16)):
        for profile in (f"lipid_order__{tail}", f"lipid_order__PSM__{tail}"):
            positions, values = profiles[profile]
            np.testing.assert_array_equal(positions, np.arange(first, last + 1), err_msg=profile)
            np.testing.assert_allclose(values, 0.5, atol=1e-5, err_msg=profile)
    # Pooled sn-1/sn-2 stay POPC-only: PSM chains never enter them.
    for chain in ("sn1", "sn2"):
        positions, expected = expected_minus_sch(chain)
        np.testing.assert_array_equal(profiles[f"lipid_order__{chain}"][0], positions)
        np.testing.assert_allclose(profiles[f"lipid_order__{chain}"][1],
                                   np.tile(expected, (5, 1)), atol=1e-5)
    psm = mapping[mapping.resname == "PSM"]
    assert set(psm["tail"]) == {"sphingosine", "n_acyl"}
    assert set(psm.n_lipids) == {1}


def _cholesterol(builder):
    residue = builder.residue("CHL1")
    o = builder.atom(residue, "O3", "O", [50, 50, 30])
    c = builder.atom(residue, "C3", "C", [50, 50, 31.4], o)
    builder.atom(residue, "H3", "H", [50, 51, 31.4], c)


def _resolved(order, *, local_order=None):
    raw = {"systems": {"s": {"topology": "s.pdb", "replicates": [{"trajectory": "s.dcd"}]}},
           "analyses": {"membrane_order": {"enabled": True, **order}}}
    if local_order is not None:
        raw["analyses"]["sdc4_mechanisms"] = {
            "enabled": True, **SDC4_DOMAINS, "lipid_contact_selection": "resname POPC",
            "lipid_anchors": {"POPC": "resname POPC and name P"},
            "local_order": local_order}
    return resolve_config(raw)


def _run_order(universe, cfg, root):
    observer = FusedMembrane(universe, cfg, root, names=("membrane_order",))
    observer.prune_from(0)
    for frame, ts in enumerate(universe.trajectory):
        observer.update(ts, frame)
    observer.flush(0)
    observer.finalize()
    out = root / "profiles" / "membrane_order"
    profiles = {}
    for row in pd.read_parquet(out / "profile_index.parquet").itertuples(index=False):
        with np.load(out / row.file) as data:
            profiles[row.profile] = (data["positions"], data["values"])
    return profiles, pd.read_parquet(out / "site_mapping.parquet")


def test_two_naming_conventions_reproduce_the_same_analytic_sn_profiles(tmp_path):
    charmm = _resolved({"tail_preset": "charmm36", "lipid_selection": "resname POPC",
                        "leaflet_anchor_selection": "name P"})
    macrog = _resolved({"tails": MACROG_TAILS, "lipid_selection": "resname POPC",
                        "leaflet_anchor_selection": "name P8"})
    results = {}
    for name, cfg in (("charmm36", charmm), ("macrog", macrog)):
        results[name] = _run_order(bilayer(name).universe(), cfg, tmp_path / name)
    for name, (profiles, mapping) in results.items():
        for chain in ("sn1", "sn2"):
            positions, expected = expected_minus_sch(chain)
            for profile in (f"lipid_order__{chain}", f"lipid_order__POPC__{chain}",
                            f"lipid_order__{chain}__upper", f"lipid_order__POPC__{chain}__lower"):
                observed_positions, values = profiles[profile]
                np.testing.assert_array_equal(observed_positions, positions, err_msg=f"{name} {profile}")
                np.testing.assert_allclose(values, np.tile(expected, (5, 1)), atol=1e-5,
                                           err_msg=f"{name} {profile}")
        assert set(mapping.n_lipids) == {4}
    identity = ["resname", "tail", "carbon_position", "n_lipids"]
    pd.testing.assert_frame_equal(results["charmm36"][1][identity], results["macrog"][1][identity])


def test_default_tails_label_the_correct_charmm36_chains(tmp_path):
    """Omitted naming and the preset must identify sn-1 as C3x, sn-2 as C2x."""
    cfg = _resolved({"lipid_selection": "resname POPC", "leaflet_anchor_selection": "name P"})
    assert "tail_preset" not in cfg["analyses"]["membrane_order"]
    assert cfg["analyses"]["membrane_order"]["tails"] == DEFAULTS["analyses"]["membrane_order"]["tails"]
    assert cfg["analyses"]["membrane_order"]["tails"] == {
        "sn1": {"carbon_pattern": r"^C3(\d+)$"},
        "sn2": {"carbon_pattern": r"^C2(\d+)$"},
    }
    profiles, mapping = _run_order(bilayer("charmm36").universe(), cfg, tmp_path)
    for chain in ("sn1", "sn2"):
        positions, expected = expected_minus_sch(chain)
        np.testing.assert_array_equal(profiles[f"lipid_order__{chain}"][0], positions)
        np.testing.assert_allclose(profiles[f"lipid_order__{chain}"][1],
                                   np.tile(expected, (5, 1)), atol=1e-5)
        subset = mapping[mapping["tail"] == chain]
        assert set(subset["atom_name"].str.match(r"^C" + ("3" if chain == "sn1" else "2"))) == {True}


def test_legacy_results_require_new_config_hash_and_rerun(tmp_path):
    """Pre-fix stored sn1=C2x and sn2=C3x must not masquerade as corrected data."""
    options = {"lipid_selection": "resname POPC", "leaflet_anchor_selection": "name P"}
    corrected = _resolved(options)
    legacy = _resolved({**options, "tails": {
        "sn1": {"carbon_pattern": r"^C2(\d+)$"},
        "sn2": {"carbon_pattern": r"^C3(\d+)$"},
    }})
    assert config_hash(corrected) != config_hash(legacy)
    universe = bilayer("charmm36").universe()
    new_profiles, _ = _run_order(universe, corrected, tmp_path / "corrected")
    old_profiles, _ = _run_order(bilayer("charmm36").universe(), legacy, tmp_path / "legacy")
    for chain, old_chain in (("sn1", "sn2"), ("sn2", "sn1")):
        np.testing.assert_array_equal(
            new_profiles[f"lipid_order__{chain}"][0],
            old_profiles[f"lipid_order__{old_chain}"][0],
        )
        np.testing.assert_allclose(
            new_profiles[f"lipid_order__{chain}"][1],
            old_profiles[f"lipid_order__{old_chain}"][1],
            atol=1e-5,
        )


def test_default_sdc4_local_order_follows_correct_sn_chain_names():
    universe = bilayer("charmm36").universe()
    local = {"enabled": True, "lipid_selection": "resname POPC", "min_carbons_per_tail": 4}
    naming = _resolved({})["analyses"]["membrane_order"]
    setup = _normal_order_setup(universe, local, naming)
    for site in setup["sites"]:
        assert site.atom_name.startswith("C3" if site.tail == "sn1" else "C2")


@pytest.mark.parametrize("order", [{"tail_preset": "charmm36"}, {}])
def test_charmm_patterns_on_sequentially_numbered_topology_fail_bond_preflight(order):
    """^C2/^C3 also match MacRog C20-C29 and C30-C39: plausible but wrong.

    C29 (sn-2 chain) bonds to C30 in the other configured tail and C36 bonds
    to C38 two positions away; the bond graph rejects the mapping.
    """
    spec = _resolved({**order, "lipid_selection": "resname POPC"})["analyses"]["membrane_order"]
    universe = bilayer("macrog").universe()
    with pytest.raises(ValueError, match="contradicts topology bonds"):
        _candidate_sites(universe, spec)
    with pytest.raises(ValueError, match="contradicts topology bonds"):
        order_setup(universe, spec)
    # Without bonds the same mapping cannot be checked; preflight says so.
    unbonded = order_preflight(bilayer("macrog").universe(bonds=False), spec)
    assert unbonded["chains"] == {"bond_graph": False, "residue_tails": 8, "verified_residue_tails": 0}
    assert unbonded["ch_pairs"] == "distance fallback"


def test_split_tail_and_non_carbon_mappings_are_rejected():
    universe = bilayer("charmm36").universe()
    gap = {"sn1": {"carbon_names": ["C32", "C33", "C35", "C36"]}}
    spec = _resolved({"tail_preset": "charmm36", "tails": gap,
                      "lipid_selection": "resname POPC"})["analyses"]["membrane_order"]
    with pytest.raises(ValueError, match="2 disconnected bonded segments"):
        _candidate_sites(universe, spec)
    oxygen = {"sn1": {"carbon_names": ["O31", "C32", "C33", "C34"]}}
    spec = _resolved({"tail_preset": "charmm36", "tails": oxygen,
                      "lipid_selection": "resname POPC"})["analyses"]["membrane_order"]
    with pytest.raises(ValueError, match="non-carbon atom POPC tail sn1 atom O31"):
        _candidate_sites(universe, spec)
    shared = {"sn2": {"carbon_pattern": r"^C3(\d+)$"}}
    spec = _resolved({"tail_preset": "charmm36", "tails": shared,
                      "lipid_selection": "resname POPC"})["analyses"]["membrane_order"]
    with pytest.raises(ValueError, match="maps POPC atom C31 to both tails sn1 and sn2"):
        _candidate_sites(universe, spec)


def test_united_atom_tails_report_unavailable_order_instead_of_zero():
    spec = _resolved({"tails": MACROG_TAILS, "lipid_selection": "resname POPC"})["analyses"]["membrane_order"]
    with pytest.raises(ValueError, match="no explicit lipid hydrogens"):
        order_preflight(bilayer("macrog", hydrogens=False).universe(), spec)


def test_leaflet_anchor_must_cover_every_order_lipid():
    spec = _resolved({"tail_preset": "charmm36", "lipid_selection": "resname POPC",
                      "leaflet_anchor_selection": "name P and resid 1:3"})["analyses"]["membrane_order"]
    with pytest.raises(ValueError, match="matched no atom in 1 order lipids"):
        order_setup(bilayer("charmm36").universe(), spec)


def _validated(tmp_path, builder, raw_order, *, local_order=None):
    if local_order is not None:
        for x in (0., 10.):
            for position, name in enumerate(
                ("ALA",) * 4 + ("GLY", "GLY", "ILE", "VAL", "GLY")
            ):
                residue = builder.residue(name)
                builder.atom(residue, "CA", "C", [x, 0., 40. - position])
    u = builder.universe(frames=1)
    if local_order is not None:
        u.atoms[-18:].chainIDs = ["A"] * 9 + ["B"] * 9
    topology = tmp_path / "bilayer.pdb"
    with mda.Writer(str(topology), n_atoms=len(u.atoms), bonds="all") as writer:
        writer.write(u.atoms)
    analyses = {name: {"enabled": False} for name in (
        "alignment", "rmsd", "rmsf", "radius_of_gyration", "dssp", "salt_bridges", "lifetimes")}
    analyses["membrane_order"] = {"enabled": True, **raw_order}
    if local_order is not None:
        analyses["sdc4_mechanisms"] = {
            "enabled": True, **SDC4_DOMAINS, "require_motifs": False, "require_gxxxg": True,
            "motif_selection": "name CA",
            "lipid_contact_selection": "resname POPC",
            "lipid_anchors": {"POPC": "resname POPC and name P P8"}, "local_order": local_order}
    raw = {"systems": {"s": {"topology": topology.name,
                             "replicates": [{"trajectory": topology.name}]}},
           "selections": {"protein": "protein" if local_order is not None else "resname POPC"},
           "analyses": analyses}
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(raw))
    cfg = resolve_config(raw, base=tmp_path)
    lines = []
    counts = validate_selections(cfg, tmp_path, emit=lines.append)
    return cfg, path, counts, "\n".join(lines)


def test_validate_shows_counts_and_explicit_exclusions_per_species_and_tail(tmp_path):
    builder = bilayer("charmm36", extra=(_psm, _cholesterol))
    cfg, path, counts, text = _validated(tmp_path, builder, {
        "tail_preset": "charmm36", "lipid_selection": "resname POPC PSM CHL1",
        "leaflet_anchor_selection": "name P",
        "lipid_overrides": {"DOPC": {"sn1": {"carbon_pattern": r"^C3(\d+)$"}}}})
    assert ("membrane_order: tail_preset=charmm36, 5 lipids, C-H pairs from topology bonds, "
            "tail chains bond-verified for 10/10 residue tails") in text
    assert "membrane_order POPC sn1: 4/4 lipids, positions 2-16, 124 C-H bonds" in text
    assert "membrane_order POPC sn2: 4/4 lipids, positions 2-18, 132 C-H bonds" in text
    assert "membrane_order POPC sphingosine: unavailable for all 4 lipids (4 no matching carbons)" in text
    assert "membrane_order PSM sphingosine: 1/1 lipids, positions 1-18, 18 C-H bonds" in text
    assert "membrane_order PSM n_acyl: 1/1 lipids, positions 2-16, 15 C-H bonds" in text
    assert "membrane_order PSM sn1: unavailable for all 1 lipids (1 no matching carbons)" in text
    assert "membrane_order CHL1: 1 residues excluded by exclude_resnames" in text
    assert "membrane_order: lipid_overrides.DOPC matches no selected residue" in text
    assert "membrane_order initial leaflets at topology coordinates: upper=3, lower=2" in text
    audit = counts["s"]["membrane_order"]
    assert audit["multi_residue_molecules"] == 0
    psm = {row["tail"]: row["mapped_residues"] for row in audit["species"] if row["resname"] == "PSM"}
    assert psm == {"sn1": 0, "sn2": 0, "sphingosine": 1, "n_acyl": 1}
    plan = build_plan(cfg, path, validation_mode="topology_aware", topology_counts=counts)
    order_entry = next(item for item in plan["analyses"] if item["analysis"] == "membrane_order")
    assert order_entry["tail_naming"] == {
        "tail_preset": "charmm36",
        "tails": {"sn1": {"carbon_pattern": r"^C3(\d+)$"}, "sn2": {"carbon_pattern": r"^C2(\d+)$"},
                  "sphingosine": {"carbon_pattern": r"^C(\d+)S$"},
                  "n_acyl": {"carbon_pattern": r"^C(\d+)F$"}},
        "lipid_overrides": {"DOPC": {"sn1": {"carbon_pattern": r"^C3(\d+)$"}}}}
    assert plan["topology_counts"]["s"]["membrane_order"]["species"] == audit["species"]


def test_validate_fails_with_system_context_for_wrong_mapping(tmp_path):
    with pytest.raises(ConfigError, match=r"s: membrane_order tail mapping contradicts topology bonds"):
        _validated(tmp_path, bilayer("macrog"), {"tail_preset": "charmm36",
                                                 "lipid_selection": "resname POPC"})


def test_sdc4_local_order_uses_membrane_order_naming(tmp_path):
    universe = bilayer("macrog").universe()
    local = {"enabled": True, "lipid_selection": "resname POPC", "min_carbons_per_tail": 4}
    # The SDC4 observer previously always used the default CHARMM-like patterns.
    with pytest.raises(ValueError, match="contradicts topology bonds"):
        _normal_order_setup(universe, local)
    naming = _resolved({"tails": MACROG_TAILS}, local_order=local)["analyses"]["membrane_order"]
    setup = _normal_order_setup(universe, local, naming)
    assert sorted({(s.tail, s.position) for s in setup["sites"]}) == sorted(
        [("sn1", p) for p in range(2, 17)] + [("sn2", p) for p in range(2, 19)])
    _, _, counts, text = _validated(tmp_path, bilayer("macrog"), {"tails": MACROG_TAILS,
                                    "lipid_selection": "resname POPC"}, local_order=local)
    assert "SDC4 local_order POPC sn2: 4/4 lipids, positions 2-18" in text
    assert "SDC4 TM A: sequence=GGIVG" in text
    assert "SDC4 TM B: sequence=GGIVG" in text
    assert counts["s"]["sdc4_local_order"]["lipids_without_anchor"] == 0


@pytest.mark.parametrize("domain", ["tm_chains", "ectodomains"])
def test_sdc4_preflight_rejects_lipid_residues_as_protein_domains(tmp_path, domain):
    cfg, _, _, _ = _validated(
        tmp_path, bilayer("macrog"),
        {"tails": MACROG_TAILS, "lipid_selection": "resname POPC"},
        local_order={"enabled": True, "lipid_selection": "resname POPC"},
    )
    # The old fixture used POPC residues as fake protein domains.
    cfg["analyses"]["sdc4_mechanisms"][domain]["A"] = "resid 1"
    with pytest.raises(ConfigError, match=rf"SDC4 {domain}\.A includes residues outside"):
        validate_selections(cfg, tmp_path)


def test_tail_preset_resolution_precedence_and_round_trip():
    explicit = {"carbon_names": ["C32", "C33", "C34", "C35"]}
    cfg = _resolved({"tail_preset": "charmm36", "tails": {"sn1": explicit,
                                                          "sn3": {"carbon_pattern": r"^C4(\d+)$"}}})
    order = cfg["analyses"]["membrane_order"]
    assert order["tails"] == {"sn1": explicit, "sn2": {"carbon_pattern": r"^C2(\d+)$"},
                              "sphingosine": {"carbon_pattern": r"^C(\d+)S$"},
                              "n_acyl": {"carbon_pattern": r"^C(\d+)F$"},
                              "sn3": {"carbon_pattern": r"^C4(\d+)$"}}
    assert resolve_config(yaml.safe_load(yaml.safe_dump(cfg))) == cfg
    # Without a preset an explicit tail still merges into the default tail.
    legacy = _resolved({"tails": {"sn1": explicit}})["analyses"]["membrane_order"]["tails"]
    assert legacy["sn1"] == {"carbon_pattern": r"^C3(\d+)$", **explicit}
    assert legacy["sn2"] == {"carbon_pattern": r"^C2(\d+)$"}
    # A per-lipid naming override replaces the shared naming for that lipid only.
    spec = {"tails": {"sn1": {"carbon_names": ["C21"], "positions": [1]}},
            "lipid_overrides": {"DOPC": {"sn1": {"carbon_pattern": r"^C3(\d+)$"}}}}
    assert _tail_spec(spec, "DOPC", "sn1") == {"carbon_pattern": r"^C3(\d+)$"}
    assert _tail_spec(spec, "POPC", "sn1") == spec["tails"]["sn1"]


@pytest.mark.parametrize("order,message", [
    ({"tail_preset": "charmm"}, "tail_preset must be null or one of ['charmm36']"),
    ({"tail_preset": "charmm36", "tails": ["sn1"]}, "membrane_order.tails must be a mapping"),
    ({"lipid_overrides": {"POPC": {"sn3": {"carbon_pattern": "^C4(\\d+)$"}}}},
     "lipid_overrides.POPC.sn3: tail 'sn3' is not defined"),
    ({"tails": {"sn1": {"carbon_pattern": "^C3(\\d+$"}}}, "not a valid regular expression"),
    ({"tails": {"sn1": {"carbon_names": ["C31", "C32"], "positions": [1]}}},
     "tails.sn1.positions must match carbon_names"),
])
def test_invalid_tail_naming_fails_before_topology(order, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        _resolved(order)


def test_disabled_membrane_order_accepts_preset_and_keeps_production_yaml_unchanged():
    disabled = resolve_config({"systems": {"s": {"topology": "s.pdb", "replicates": [{"trajectory": "s.dcd"}]}},
                               "analyses": {"membrane_order": {"tail_preset": "charmm36"}}})
    assert disabled["analyses"]["membrane_order"]["tails"]["sn1"] == {"carbon_pattern": r"^C3(\d+)$"}
    for name in ("synd4_baseline.yaml", "synd4_full_mechanisms.yaml",
                 "synd4_full_mechanisms_finalize_fast.yaml"):
        with open(name) as handle:
            raw = yaml.safe_load(handle)
        order = resolve_config(raw)["analyses"]["membrane_order"]
        expected_preset = name != "synd4_baseline.yaml"
        assert order.get("tail_preset") == ("charmm36" if expected_preset else None)
        assert {key: order["tails"][key] for key in ("sn1", "sn2")} == DEFAULTS["analyses"]["membrane_order"]["tails"]
        if expected_preset:
            assert set(order["tails"]) == {"sn1", "sn2", "sphingosine", "n_acyl"}
        else:
            assert set(order["tails"]) == {"sn1", "sn2"}
        assert order["lipid_overrides"] == {}


def test_documented_user_defined_example_is_the_tested_mapping():
    text = Path("docs/configuration.md").read_text(encoding="utf-8")
    section = text[text.index("## Lipid-tail naming for -S_CH"):]
    block = section[section.index("```yaml") + 7:section.index("```", section.index("```yaml") + 7)]
    order = yaml.safe_load(block)["analyses"]["membrane_order"]
    assert order["tails"] == MACROG_TAILS


def test_production_order_statistics_feed_six_column_analysis_report(tmp_path, monkeypatch):
    """Real order writer -> profile inference -> offline reporting, no relabeling."""
    from mdscope.reporting.statistical_tables import collect_statistical_reports
    from mdscope.statistics.profiles import analyze_profiles
    cfg = _resolved({'tail_preset': 'charmm36', 'lipid_selection': 'resname POPC PSM',
                     'leaflet_anchor_selection': 'name P'})
    cfg['statistics'].update(bootstrap_draws=8, hypothesis_tests={'enabled': False})
    cfg['progress'] = {'enabled': False}
    reps = {}
    for system, count in [('A', 2), ('B', 1), ('C', 1)]:
        reps[system] = []
        for i in range(count):
            rep = tmp_path / system / f'rep{i}'
            _run_order(bilayer('charmm36', extra=(_psm,)).universe(), cfg, rep)
            reps[system].append(rep)
    analyze_profiles(cfg, tmp_path, reps)
    wanted = [f'lipid_order__POPC__{tail}{leaf}' for tail in ('sn1', 'sn2')
              for leaf in ('__lower', '__upper', '')]
    cfg['plots'].update(enabled=True, bootstrap_tables=False,
                        include=[f'profiles:{profile}:comparison_table' for profile in wanted])
    # Plotting/collection cannot reopen trajectories or rerun inference.
    monkeypatch.setattr(mda, 'Universe', lambda *a, **k: pytest.fail('reporting opened a trajectory'))
    rows = collect_statistical_reports(cfg, tmp_path)
    assert set(rows.profile) == set(wanted)
    assert rows.analysis.unique().tolist() == ['membrane_order']
    assert set(zip(rows.system_1, rows.system_2)) == {('A', 'B'), ('A', 'C'), ('B', 'C')}
    assert set(rows.condition_1_n_replicates) == {1, 2}
    assert len(rows) > 40
    for tail in ('sn1', 'sn2'):
        positions, expected = expected_minus_sch(tail)
        actual = rows.loc[(rows.profile == f'lipid_order__POPC__{tail}') &
                          (rows.system_1 == 'A') & (rows.system_2 == 'B')].sort_values('position')
        np.testing.assert_allclose(actual.condition_1_mean, expected, atol=1e-5)
