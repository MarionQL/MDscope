"""Dependency-free statistics metric selection and eligible profile contracts."""
import re

# Residue-resolved protein observables with comparable saved labels. Shared by
# YAML validation and the saved-data screen so eligibility has one definition.
RESIDUE_RESPONSE_PROFILE = re.compile(
    r"^(?:rmsf$|dssp_(?:helix|sheet|coil)_occupancy$|"
    r"sasa_mean(?:__|$)|sasa_buried_nm2__|"
    r"residue_occlusion__|residue_water_coordination__all$|"
    r"residue_membrane_distance__(?!.*__group_))"
)


# Every statistics family and its directory under project_summary/statistics/.
# The one family list for YAML validation, the dispatcher and read-only planning.
STATISTICS_FAMILIES = {
    "time_series": "time_series",
    "profiles": "profiles",
    "interactions": "interactions",
    "solvation_rdfs": "solvation",
    "diffusion": "diffusion",
    "coupling": "coupling",
    "lipid_contacts": "lipid_contacts",
}


SIMPLE_FAMILY_KEYS = {
    "profiles": ("profile",),
    "solvation_rdfs": ("solvent", "target_class", "site_class"),
    "diffusion": ("metric",),
    "coupling": ("mechanism",),
    "lipid_contacts": ("level", "leaflet", "species_a", "species_b", "metric"),
}


# All families have descriptive point summaries. Whether to compute additional
# inferential products can be configured independently for each selected metric.
# The default remains True, so older selectors keep their original behavior.
# Capabilities reflect writer outputs, not GUI assumptions.
ESTIMATOR_FLAGS = {
    "time_series": frozenset({"ci", "hypothesis_tests"}),
    "profiles": frozenset({"ci", "hypothesis_tests"}),
    "interactions": frozenset({"ci", "hypothesis_tests"}),
    "solvation_rdfs": frozenset({"ci"}),
    "diffusion": frozenset({"ci"}),
    "coupling": frozenset({"ci"}),
    "lipid_contacts": frozenset({"ci"}),
}


def _valid_estimator_flags(family, entry):
    for key in ("inference", "ci", "hypothesis_tests"):
        if key in entry and not isinstance(entry[key], bool):
            return False
    if any(key in entry and key not in ESTIMATOR_FLAGS[family]
           for key in ("ci", "hypothesis_tests")):
        return False
    if entry.get("inference") is False and (
        entry.get("ci") is True or entry.get("hypothesis_tests") is True
    ):
        return False
    return True


INTERACTION_METRICS = frozenset({
    "occupancy", "event_rate_per_ns", "mean_lifetime_ps",
    "median_lifetime_ps", "total_lifetime_ps",
})


def validate_metric_selection(stats):
    selection = stats.get("metric_selection")
    if selection is None:
        return
    if not isinstance(selection, dict) or set(selection) - {"time_series", "interactions", *SIMPLE_FAMILY_KEYS}:
        raise ValueError("statistics.metric_selection contains unsupported family")
    for family, entries in selection.items():
        if not isinstance(entries, list):
            raise ValueError(f"statistics.metric_selection.{family} must be a list")
        seen = set()
        if family in SIMPLE_FAMILY_KEYS:
            keys = SIMPLE_FAMILY_KEYS[family]
            for entry in entries:
                if (not isinstance(entry, dict) or set(entry) - (set(keys) | {"inference", "ci", "hypothesis_tests"})
                        or any(key not in entry for key in keys)
                        or not _valid_estimator_flags(family, entry)
                        or any(not isinstance(entry.get(key), str) or not entry[key].strip()
                               for key in keys)):
                    raise ValueError(
                        f"statistics.metric_selection.{family} entries require "
                        f"nonempty {', '.join(keys)}"
                    )
                identity = tuple(entry[key] for key in keys)
                if identity in seen:
                    raise ValueError(f"Duplicate statistics.metric_selection.{family}: {identity!r}")
                seen.add(identity)
            continue
        if family == "interactions":
            for entry in entries:
                if (not isinstance(entry, dict) or set(entry) - {"kind", "interaction", "metric", "inference", "ci", "hypothesis_tests"}
                        or not _valid_estimator_flags(family, entry)
                        or not isinstance(entry.get("kind"), str) or not entry["kind"].strip()
                        or ("interaction" in entry and
                            (not isinstance(entry["interaction"], str) or not entry["interaction"].strip()))
                        or ("metric" in entry and
                            (not isinstance(entry["metric"], str)
                             or entry["metric"] not in INTERACTION_METRICS))):
                    raise ValueError(
                        "statistics.metric_selection.interactions entries require a "
                        "nonempty kind, optional exact interaction and a supported metric"
                    )
                key = (entry["kind"], entry.get("interaction"), entry.get("metric"))
                if key in seen:
                    raise ValueError(f"Duplicate statistics.metric_selection.interactions: {key!r}")
                seen.add(key)
            continue
        for entry in entries:
            if (not isinstance(entry, dict) or set(entry) - {"metric", "group", "inference", "ci", "hypothesis_tests"}
                    or not _valid_estimator_flags(family, entry)
                    or not isinstance(entry.get("metric"), str)
                    or not entry["metric"].strip()
                    or (entry.get("group") is not None
                        and (not isinstance(entry["group"], str) or not entry["group"].strip()))):
                raise ValueError(
                    "statistics.metric_selection.time_series entries require a nonempty "
                    "metric and optional nonempty group"
                )
            key = (entry.get("group"), entry["metric"])
            if key in seen:
                raise ValueError(f"Duplicate statistics.metric_selection.time_series: {key!r}")
            seen.add(key)


