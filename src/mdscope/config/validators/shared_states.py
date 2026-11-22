"""Validation rules for shared states configuration."""


def validate_shared_states(cfg):
    from ...contracts.state_feature_spec import validate_state_features
    mechanisms = cfg["analyses"].get("sdc4_mechanisms", {})
    shared = cfg["analyses"].get("shared_states", {})
    validate_state_features(shared)
    if shared.get("enabled") and shared.get("preset", "sdc4") == "sdc4" and not mechanisms.get("enabled"):
        raise ValueError("shared_states preset sdc4 requires sdc4_mechanisms.enabled")
