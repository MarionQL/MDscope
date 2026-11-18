"""Shared solvation normalization and persisted episode-column contract.

Imported by both the source-frame solvation observer and offline
postprocessors. Keep the established RDF reference-density formula and
episode column ordering unchanged.
"""
import numpy as np


EPISODE_COLUMNS = ["solvent", "target_id", "solvent_id", "start_frame", "end_frame",
                   "duration_ps", "left_censored", "right_censored",
                   "sampled_duration_ps", "duration_convention"]


def rdf_from_counts(counts, edges, density_exposure):
    """Observed pairs / expected pairs at the supplied reference density exposure."""
    shell_volume = (4 * np.pi / 3) * np.diff(np.asarray(edges, float) ** 3)
    if density_exposure <= 0:
        return np.full(len(counts), np.nan)
    return np.asarray(counts, float) / (shell_volume * density_exposure)
