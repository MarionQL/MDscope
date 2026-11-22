"""Residue-response options for saved statistical profiles."""
import math
from copy import deepcopy

from ..defaults import DEFAULTS

def residue_response_settings(cfg):
    """Resolve and validate the saved-data exploratory screen, independently of UI."""
    supplied = cfg.get("statistics", {}).get("residue_response", {})
    if not isinstance(supplied, dict):
        raise ValueError("statistics.residue_response must be a mapping")
    defaults = DEFAULTS["statistics"]["residue_response"]
    unknown = set(supplied) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown statistics.residue_response options: {sorted(unknown)}")
    settings = {**deepcopy(defaults), **deepcopy(supplied)}
    if type(settings["enabled"]) is not bool:
        raise ValueError("statistics.residue_response.enabled must be boolean")
    for key in ("profiles", "chains"):
        values = settings[key]
        if values is not None and (not isinstance(values, list) or
                any(not isinstance(v, str) or not v.strip() for v in values) or
                len(set(values)) != len(values)):
            raise ValueError(f"statistics.residue_response.{key} must be null or distinct nonempty strings")
    from ...contracts.metric_selection import RESIDUE_RESPONSE_PROFILE
    unsupported = [name for name in settings["profiles"] or []
                   if not RESIDUE_RESPONSE_PROFILE.match(name)]
    if unsupported:
        raise ValueError(f"Unsupported statistics.residue_response profiles: {unsupported}")
    length = settings["min_length"]
    if type(length) is not int or length < 2:
        raise ValueError("statistics.residue_response.min_length must be an integer >= 2")
    cosine = settings["minimum_cosine"]
    if (isinstance(cosine, bool) or not isinstance(cosine, (int, float)) or
            not math.isfinite(cosine) or not 0 <= cosine <= 1):
        raise ValueError("statistics.residue_response.minimum_cosine must be finite between 0 and 1")
    return settings
