"""Radius of gyration from an MDAnalysis AtomGroup (Å)."""
from __future__ import annotations
import numpy as np
def compute_rg(atoms)->float:
    return float(atoms.radius_of_gyration()) if len(atoms) else float("nan")
