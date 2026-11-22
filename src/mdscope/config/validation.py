"""Resolved YAML validation, dispatched in the established rule order."""
from pathlib import Path

from .validators.inputs import validate_inputs
from .validators.structure import validate_structure
from .validators.membrane import validate_membrane
from .validators.solvation import validate_solvation
from .validators.hydration import validate_hydration
from .validators.sdc4 import validate_sdc4
from .validators.statistics import validate_shared_states
from .validators.plots import validate_plots
from .validators.statistics import validate_statistics
from .validators.residue_response import residue_response_settings as residue_response_settings


def validate_config(cfg: dict, base: Path = Path("."), *, check_paths: bool = True) -> None:
    """Validate input paths, analysis selections and inference settings without mutating cfg."""
    validate_inputs(cfg, base, check_paths=check_paths)
    validate_structure(cfg)
    validate_membrane(cfg)
    validate_solvation(cfg)
    validate_hydration(cfg)
    validate_sdc4(cfg)
    validate_shared_states(cfg)
    validate_plots(cfg)
    validate_statistics(cfg)
