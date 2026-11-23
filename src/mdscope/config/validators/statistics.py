"""Cohesive statistics, shared-state and residue-response configuration validation."""
import math
from copy import deepcopy

from ..defaults import DEFAULTS
from ..runs import replicate_specs
from .common import _validate_time_window


def validate_statistics(cfg):
    stats = cfg.get("statistics", {})
    residue_response_settings(cfg)
    residue_correlations_settings(cfg)
    from ...contracts.statistics import validate_metric_selection
    validate_metric_selection(stats)
    for field in ("block", "hypothesis_tests"):
        if stats.get(field) is not None and not isinstance(stats[field], dict):
            raise ValueError(f"statistics.{field} must be a mapping")
    if int(stats.get("bootstrap_draws", 2000)) < 100:
        raise ValueError("statistics.bootstrap_draws must be >= 100")
    if not 0 < float(stats.get("ci", 0.95)) < 1:
        raise ValueError("statistics.ci must be between 0 and 1")
    for key in ("interaction_include_kinds", "interaction_exclude_kinds"):
        value = stats.get(key)
        if value is not None and (not isinstance(value, list) or any(not isinstance(x, str) for x in value)):
            raise ValueError(f"statistics.{key} must be a list of kinds")
    if int(stats.get("interaction_progress_every", 1000)) < 1:
        raise ValueError("statistics.interaction_progress_every must be positive")
    _validate_time_window(stats.get("scalar_window_ps"), "statistics.scalar_window_ps")
    if stats.get("matched_window_ps") != "auto":
        _validate_time_window(stats.get("matched_window_ps"), "statistics.matched_window_ps")
    from ...contracts.statistics import STATISTICS_FAMILIES
    families = stats.get("families", {})
    if not isinstance(families, dict) or set(families) - set(STATISTICS_FAMILIES):
        raise ValueError("statistics.families has unknown families or is not a mapping")
    for name, option in families.items():
        if isinstance(option, bool):
            continue
        if not isinstance(option, dict) or set(option) - (
                {"enabled", "bootstrap", "bootstrap_draws", "ci", "hypothesis_tests", "hierarchy"}
                | ({"adaptive_hierarchy"} if name == "time_series" else set())):
            raise ValueError(f"statistics.families.{name}: unsupported options")
        if "adaptive_hierarchy" in option:
            adaptive = option["adaptive_hierarchy"]
            counts = {"min_draws", "max_draws", "batch_size", "stable_checks"}
            tolerances = {"absolute_tolerance", "relative_tolerance"}
            if (not isinstance(adaptive, dict) or not counts <= adaptive.keys()
                    or adaptive.keys() - counts - tolerances
                    or any(type(adaptive[k]) is not int or adaptive[k] < 1 for k in counts)
                    or not 100 <= adaptive["min_draws"] <= adaptive["max_draws"]
                    or any(type(adaptive.get(k, 0)) not in (int, float)
                           or not math.isfinite(adaptive.get(k, 0)) or adaptive.get(k, 0) < 0
                           for k in tolerances)
                    or not any(adaptive.get(k, 0) > 0 for k in tolerances)):
                raise ValueError("statistics.families.time_series.adaptive_hierarchy requires valid "
                                 "draw bounds, batch_size, stable_checks and a positive tolerance")
        for key in ("enabled", "bootstrap", "hypothesis_tests"):
            if key in option and not isinstance(option[key], bool):
                raise ValueError(f"statistics.families.{name}.{key} must be boolean")
        if "bootstrap_draws" in option and int(option["bootstrap_draws"]) < 100:
            raise ValueError(f"statistics.families.{name}.bootstrap_draws must be >= 100")
        if "ci" in option and not 0 < float(option["ci"]) < 1:
            raise ValueError(f"statistics.families.{name}.ci must be between 0 and 1")
        hierarchy = option.get("hierarchy", {})
        if (not isinstance(hierarchy, dict) or set(hierarchy) - {"enabled", "block_within_replicate", "resample_replicates"}
                or any(not isinstance(v, bool) for v in hierarchy.values())):
            raise ValueError(f"statistics.families.{name}.hierarchy must contain boolean controls")
    if not isinstance(stats.get("bootstrap", True), bool):
        raise ValueError("statistics.bootstrap must be boolean")
    hierarchy = stats.get("hierarchy", {})
    if (not isinstance(hierarchy, dict) or set(hierarchy) - {"enabled", "block_within_replicate", "resample_replicates"}
            or any(not isinstance(v, bool) for v in hierarchy.values())):
        raise ValueError("statistics.hierarchy must contain boolean controls")
    hypotheses = stats.get("hypothesis_families", {})
    if not isinstance(hypotheses, dict):
        raise ValueError("statistics.hypothesis_families must map names to metric lists")
    seen = set()
    for name, metrics in hypotheses.items():
        if not isinstance(name, str) or not name.strip() or name == "exploratory" or not isinstance(metrics, list) or not metrics or any(not isinstance(x, str) for x in metrics):
            raise ValueError("hypothesis families require named nonempty metric lists; exploratory is reserved")
        if seen.intersection(metrics) or len(set(metrics)) != len(metrics):
            raise ValueError("A metric cannot belong to multiple hypothesis families")
        seen.update(metrics)
    block = stats.get("block", {}) or {}
    if str(block.get("method", "moving")) != "moving":
        raise ValueError("statistics.block.method currently supports only 'moving'")
    length = block.get("length", "auto")
    if length != "auto" and int(length) < 1:
        raise ValueError("statistics.block.length must be 'auto' or an integer >= 1")
    if float(block.get("tau_multiplier", 2.0)) <= 0:
        raise ValueError("statistics.block.tau_multiplier must be > 0")
    max_acf_lag = block.get("max_acf_lag")
    if max_acf_lag is not None and int(max_acf_lag) < 1:
        raise ValueError("statistics.block.max_acf_lag must be null or >= 1")
    tests = stats.get("hypothesis_tests", {}) or {}
    tests_active = tests.get("enabled", True) or any(
        isinstance(option, dict) and option.get("enabled", True) and option.get("hypothesis_tests", False)
        for option in families.values())
    if tests_active:
        design = str(tests.get("design", "independent"))
        if design not in {"independent", "paired"}:
            raise ValueError(
                "statistics.hypothesis_tests.design must be independent or "
                "paired"
            )
        if design == "paired" and len(cfg.get("systems", {})) > 1:
            label_sets = []
            for system in cfg["systems"]:
                labels = {
                    spec["id"]
                    for spec in replicate_specs(cfg, system)
                }
                label_sets.append((system, labels))
            reference_system, reference_labels = label_sets[0]
            for system, labels in label_sets[1:]:
                if labels != reference_labels:
                    raise ValueError(
                        "paired hypothesis tests require identical replicate "
                        f"labels across systems; {reference_system!r} and "
                        f"{system!r} differ"
                    )
        if str(tests.get("method", "permutation")) != "permutation":
            raise ValueError(
                "statistics.hypothesis_tests.method currently supports only permutation"
            )
        if int(tests.get("permutation_draws", 10000)) < 100:
            raise ValueError(
                "statistics.hypothesis_tests.permutation_draws must be >= 100"
            )
        if int(tests.get("exact_max_assignments", 100000)) < 1:
            raise ValueError(
                "statistics.hypothesis_tests.exact_max_assignments must be >= 1"
            )
        if str(tests.get("multiple_testing", "fdr_bh")) not in {"fdr_bh", "none"}:
            raise ValueError("statistics.hypothesis_tests.multiple_testing must be fdr_bh or none")


def validate_shared_states(cfg):
    from ...contracts.state_feature_spec import validate_state_features
    mechanisms = cfg["analyses"].get("sdc4_mechanisms", {})
    shared = cfg["analyses"].get("shared_states", {})
    validate_state_features(shared)
    if shared.get("enabled") and shared.get("preset", "sdc4") == "sdc4" and not mechanisms.get("enabled"):
        raise ValueError("shared_states preset sdc4 requires sdc4_mechanisms.enabled")

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
    from ...contracts.statistics import RESIDUE_RESPONSE_PROFILE
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

def residue_correlations_settings(cfg):
    """Validate explicitly descriptive per-replica residue associations."""
    supplied = cfg.get("statistics", {}).get("residue_correlations", {})
    if not isinstance(supplied, dict):
        raise ValueError("statistics.residue_correlations must be a mapping")
    defaults = DEFAULTS["statistics"]["residue_correlations"]
    unknown = set(supplied) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown statistics.residue_correlations options: {sorted(unknown)}")
    settings = {**deepcopy(defaults), **deepcopy(supplied)}
    if type(settings["enabled"]) is not bool:
        raise ValueError("statistics.residue_correlations.enabled must be boolean")
    if type(settings["min_frames"]) is not int or settings["min_frames"] < 3:
        raise ValueError("statistics.residue_correlations.min_frames must be an integer >= 3")
    pairs = settings["pairs"]
    if pairs is not None:
        from ...contracts.statistics import RESIDUE_CORRELATION_METRICS
        if (not isinstance(pairs, list) or any(
                not isinstance(pair, (list, tuple)) or len(pair) != 2
                or any(not isinstance(name, str) or name not in RESIDUE_CORRELATION_METRICS for name in pair)
                or pair[0] == pair[1] for pair in pairs)
                or len(set(tuple(pair) for pair in pairs)) != len(pairs)):
            raise ValueError("statistics.residue_correlations.pairs must contain distinct supported metric pairs")
    return settings
