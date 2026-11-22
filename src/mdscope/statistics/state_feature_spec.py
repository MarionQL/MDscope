"""Lightweight validation shared by YAML authoring and saved-state execution."""
import math

from ..output_paths import SAVED_SCALAR_SOURCE_FILES

SDC4_FEATURES = [
    "sdc4_ectodomain_A_height_A", "sdc4_ectodomain_B_height_A",
    "sdc4_ectodomain_A_tilt_deg", "sdc4_ectodomain_B_tilt_deg",
    "sdc4_tm_crossing_angle_deg", "sdc4_local_thickness_A",
]
STATE_SOURCES = frozenset(SAVED_SCALAR_SOURCE_FILES)


def validate_state_features(spec):
    preset = spec.get("preset", "sdc4")
    if preset not in (None, "sdc4"):
        raise ValueError("shared_states.preset must be sdc4 or null")
    if spec.get("alignment", "exact") != "exact":
        raise ValueError("shared_states.alignment supports exact time grids only")
    if spec.get("scaling", "standardize") != "standardize":
        raise ValueError("shared_states.scaling supports standardize only")
    if type(spec.get("seed", 2026)) is not int or spec.get("seed", 2026) < 0:
        raise ValueError("shared_states.seed must be a nonnegative integer")
    features = spec.get("features")
    if features is not None and (not isinstance(features, list) or not features
            or any(not isinstance(name, str) or not name or name == "time_ps" for name in features)
            or len(set(features)) != len(features)):
        raise ValueError("shared_states.features must be a nonempty list of distinct names or null for the sdc4 preset")
    sources = spec.get("feature_sources", {})
    if not isinstance(sources, dict):
        raise ValueError("shared_states.feature_sources must be a mapping")
    if preset == "sdc4" and sources:
        raise ValueError("Set shared_states.preset: null to use feature_sources")
    if preset is None and (not features or set(sources) != set(features)):
        raise ValueError("shared_states.features and feature_sources must define the same nonempty names")
    for name, source in sources.items():
        if not isinstance(source, dict) or source.get("source") not in STATE_SOURCES:
            raise ValueError(f"shared_states feature {name}: unsupported saved-product source")
        for field in ("metric", "units"):
            if not isinstance(source.get(field), str) or not source[field].strip():
                raise ValueError(f"shared_states feature {name}: {field} must be nonempty")
        if source.get("group") is not None and (not isinstance(source["group"], str) or not source["group"]):
            raise ValueError(f"shared_states feature {name}: group must be null or nonempty")
        prefix = {"solvation": "solvation/", "hydration": "hydration/",
                  "solvation_bridges": "solvation_bridges/", "lipid_contacts": "lipid_contacts/"}.get(source["source"])
        if prefix and not (source.get("group") or "").startswith(prefix):
            raise ValueError(f"shared_states feature {name}: group must start with {prefix}")
    for field, minimum in (("n_states", 2), ("bootstrap_draws", 20), ("training_frames_per_replicate", 2)):
        value = spec.get(field, {"n_states": 3, "bootstrap_draws": 300, "training_frames_per_replicate": 2000}[field])
        if type(value) is not int or value < minimum:
            raise ValueError(f"shared_states.{field} must be an integer >= {minimum}")
    if spec.get("training_frames_per_replicate", 2000) < spec.get("n_states", 3):
        raise ValueError("shared_states.training_frames_per_replicate must be >= n_states")
    fraction = spec.get("min_complete_fraction", .9)
    if type(fraction) not in (int, float) or not math.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError("shared_states.min_complete_fraction must be finite in (0, 1]")
    block = spec.get("block_frames")
    if block is not None and (type(block) is not int or block < 2):
        raise ValueError("shared_states.block_frames must be null or an integer >= 2")
    if type(spec.get("seed", 2026)) is not int or spec.get("seed", 2026) < 0:
        raise ValueError("shared_states.seed must be a nonnegative integer")
