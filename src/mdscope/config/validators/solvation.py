"""Solvation, motif and glycan hydration validation."""
import numpy as np


def validate_motif_definitions(cfg):
    from ...core.motifs import validate_motifs
    validate_motifs(cfg["selections"].get("motifs", {}))
    local_motifs = cfg["analyses"].get("local_environment", {}).get("motifs")
    if local_motifs is not None:
        if not isinstance(local_motifs, dict):
            raise ValueError("local_environment.motifs must map labels to a regex or {motif: name}")
        for label, value in local_motifs.items():
            if isinstance(value, dict):
                if set(value) != {"motif"} or value["motif"] not in cfg["selections"].get("motifs", {}):
                    raise ValueError(f"local_environment motif {label}: reference must be "
                                     "{motif: <selections.motifs name>}")
            else:
                validate_motifs({label: {"pattern": value}})


def validate_solvation_targets(cfg):
    solvation = cfg["analyses"].get("solvation", {})
    if solvation.get("enabled"):
        if solvation.get("preset", "standard") not in ("standard", "detailed", "kinetics"):
            raise ValueError("solvation.preset must be standard, detailed or kinetics")
        if not isinstance(solvation.get("targets"), dict) or not solvation["targets"]:
            raise ValueError("solvation.targets must define explicit chemical-site selections")
        if solvation.get("ion_selection") and solvation.get("ion_selections"):
            raise ValueError("solvation: specify ion_selection OR ion_selections")
        if not (solvation.get("water_selection") or solvation.get("ion_selection")
                or solvation.get("ion_selections")):
            raise ValueError("solvation requires water_selection or ion selections")
        for kind, groups in solvation["targets"].items():
            if not isinstance(groups, dict) or not groups:
                raise ValueError(f"solvation.targets.{kind} must map labels to sites")
            for label, selection in groups.items():
                if isinstance(selection, dict) and "motif" in selection:
                    if (not isinstance(selection["motif"], str)
                            or selection["motif"] not in cfg["selections"].get("motifs", {})):
                        raise ValueError(f"solvation target {kind}/{label}: unknown motif reference")
                    if any(key in selection for key in ("sequence", "selection", "group_by")):
                        raise ValueError("solvation motif references cannot override sequence, selection or group_by")
                    continue
                if isinstance(selection, dict) and "sequence" in selection:
                    import re
                    sequence = selection["sequence"]
                    if not isinstance(sequence, str) or not re.fullmatch(
                        "[ACDEFGHIKLMNPQRSTVWY]+", sequence.upper()
                    ):
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}.sequence must be an exact protein sequence")
                    if "selection" in selection and (
                        not isinstance(selection["selection"], str)
                        or not selection["selection"].strip()
                    ):
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}.selection must be nonempty")
                    if selection.get("group_by", "selection") != "selection":
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}: sequence targets are grouped by motif")
                    continue
                actual = selection if isinstance(selection, str) else (
                    selection.get("selection") if isinstance(selection, dict) else None)
                if not isinstance(actual, str) or not actual.strip():
                    raise ValueError(f"solvation.targets.{kind}.{label} needs a selection")
        if solvation.get("ion_selections") and not isinstance(solvation["ion_selections"], dict):
            raise ValueError("solvation.ion_selections must be a mapping")
        radius = float(solvation.get("search_radius_A", 8.0))
        width = float(solvation.get("rdf_bin_width_A", .2))
        if not (np.isfinite(radius) and np.isfinite(width) and 0 < width <= radius):
            raise ValueError("solvation: invalid search radius or RDF bin width")
        if any(not 0 < float(value) <= radius
               for value in (solvation.get("shell_cutoffs_A") or {}).values()):
            raise ValueError("solvation shell cutoffs must be positive and <= search radius")


def validate_solvation(cfg):
    validate_motif_definitions(cfg)
    validate_solvation_targets(cfg)

def validate_hydration(cfg):
    from ...core.glycans import validate_glycans
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
