"""Explicit heavy-atom selection policies for heterogeneous MD topologies.

These policies deliberately preserve three existing, non-equivalent rules.
No chemistry rule is silently substituted for another.
"""
from __future__ import annotations
import re
import numpy as np
from MDAnalysis.exceptions import NoDataError


def heavy_atom_mask(atoms, *, policy):
    """Return a boolean mask under an explicitly named historical policy.

    'biophysical': reject names starting H or atoms with element H.
    'extended': reject names starting H, with no element lookup.
    'mechanistic': prefer element metadata; if absent reject names matching
    optional leading digits followed by H.
    """
    names = np.asarray(atoms.names, dtype=str)
    if policy == "extended":
        return ~np.char.startswith(np.char.upper(names), "H")
    if policy == "biophysical":
        mask = ~np.char.startswith(np.char.upper(names), "H")
        # Preserve the original per-atom fallback for unavailable elements.
        for i, atom in enumerate(atoms):
            try:
                element = str(getattr(atom, "element", "")).upper()
            except (AttributeError, ValueError, NotImplementedError, NoDataError):
                element = ""
            if element == "H":
                mask[i] = False
        return mask
    if policy == "mechanistic":
        try:
            elements = np.asarray(atoms.elements, dtype=str)
            return np.char.upper(elements) != "H"
        except (AttributeError, ValueError, NotImplementedError, NoDataError):
            return np.array([not re.match(r"^([0-9]+)?H", str(n), re.I)
                             for n in names], dtype=bool)
    raise ValueError(f"Unknown heavy-atom policy: {policy!r}")


def select_heavy_atoms(atoms, *, policy):
    return atoms[heavy_atom_mask(atoms, policy=policy)]
