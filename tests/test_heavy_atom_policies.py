"""Heavy-atom policy regression tests for alternate topology conventions."""
import numpy as np
import pytest
import MDAnalysis as mda

from mdscope.core.atom_selection import heavy_atom_mask, select_heavy_atoms
from mdscope.biophysical_common import heavy as biophysical_heavy
from mdscope.mechanisms.mechanistic_common import heavy as mechanistic_heavy


def atoms(with_elements):
    u = mda.Universe.empty(5, trajectory=True)
    u.add_TopologyAttr("names", ["C1", "H1", "1H2", "HG", "CA"])
    if with_elements:
        u.add_TopologyAttr("elements", ["C", "H", "H", "H", "C"])
    return u.atoms


def test_without_elements_preserves_legacy_policies():
    ag = atoms(False)
    assert list(biophysical_heavy(ag).names) == ["C1", "1H2", "CA"]
    assert list(select_heavy_atoms(ag, policy="extended").names) == ["C1", "1H2", "CA"]
    assert list(mechanistic_heavy(ag).names) == ["C1", "CA"]


def test_with_elements_preserves_legacy_policies():
    ag = atoms(True)
    assert list(biophysical_heavy(ag).names) == ["C1", "CA"]
    assert list(select_heavy_atoms(ag, policy="extended").names) == ["C1", "1H2", "CA"]
    assert list(mechanistic_heavy(ag).names) == ["C1", "CA"]


def test_element_metadata_takes_precedence_only_for_mechanistic():
    u = mda.Universe.empty(2, trajectory=True)
    u.add_TopologyAttr("names", ["HG", "1H"])
    u.add_TopologyAttr("elements", ["C", "C"])
    assert list(mechanistic_heavy(u.atoms).names) == ["HG", "1H"]
    assert list(select_heavy_atoms(u.atoms, policy="extended").names) == ["1H"]


def test_unknown_policy_fails_explicitly():
    with pytest.raises(ValueError, match="Unknown heavy-atom policy"):
        heavy_atom_mask(atoms(False), policy="guess")
