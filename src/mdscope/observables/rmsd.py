"""RMSD calculation, independent of RMSF and structural observable collection."""
from __future__ import annotations
import numpy as np
from MDAnalysis.analysis.rms import rmsd as _rmsd
def compute_rmsd(aligned_positions: np.ndarray, reference_positions: np.ndarray)->float:
    """Å RMSD after caller-specified structural fitting (no further superposition)."""
    if aligned_positions.shape!=reference_positions.shape:
        raise ValueError("RMSD selection size differs from reference")
    if len(aligned_positions)==0:return float("nan")
    return float(_rmsd(aligned_positions,reference_positions,center=False,superposition=False))
