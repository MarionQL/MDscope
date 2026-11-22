"""Compatibility imports for the former SDC4-only state-model location."""
from ..statistics.shared_states import (
                                         _write_matched_populations as _write_matched_populations,
)
from ..statistics.shared_states import fit_shared_states as fit_shared_states

__all__ = ["fit_shared_states"]
