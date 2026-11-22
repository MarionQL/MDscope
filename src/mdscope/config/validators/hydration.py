"""Validation rules for hydration configuration."""
import numpy as np


def validate_hydration(cfg):
    from ...core.glycans import validate_glycan_spec, validate_glycans
    glycans = cfg["selections"].get("glycans", {})
    validate_glycans(glycans)
    hydration = cfg["analyses"].get("hydration", {})
    solvation = cfg["analyses"].get("solvation", {})
    references = hydration.get("glycans", [])
    if (not isinstance(references, list) or len(set(map(str, references))) != len(references)
            or any(not isinstance(name, str) or name not in glycans for name in references)):
        raise ValueError("hydration.glycans must list distinct selections.glycans names ([] selects none)")
    if hydration.get("enabled"):
        if not (hydration.get("water_oxygen_selection")
                or cfg["analyses"].get("sdc4_mechanisms", {}).get("water_oxygen_selection")
                or solvation.get("water_selection")):
            raise ValueError("hydration requires a water oxygen selection")
        if (not np.isfinite(float(hydration.get("cutoff_A", 3.5)))
                or float(hydration.get("cutoff_A", 3.5)) <= 0):
            raise ValueError("hydration.cutoff_A must be finite and positive")
