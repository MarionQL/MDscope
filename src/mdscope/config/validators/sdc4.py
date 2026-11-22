"""Validation rules for sdc4 configuration."""
import numpy as np
from ...core.glycans import validate_glycan_spec
from ..runs import analysis_runs, _validate_path_component, analysis_enabled


def validate_sdc4(cfg):
    mechanisms = cfg["analyses"].get("sdc4_mechanisms", {})
    sugar_cutoff = mechanisms.get("hs_sugar_water_cutoff_A")
    if sugar_cutoff is not None and (
            not mechanisms.get("enabled") or not (mechanisms.get("hs_auto") or {}).get("enabled")
            or not mechanisms.get("water_oxygen_selection")
            or not np.isfinite(float(sugar_cutoff)) or float(sugar_cutoff) <= 0):
        raise ValueError("sdc4_mechanisms.hs_sugar_water_cutoff_A requires enabled SDC4 and hs_auto, water_oxygen_selection and a positive cutoff")
    water_cutoff = mechanisms.get("residue_water_cutoff_A")
    if water_cutoff is not None and (not mechanisms.get("residue_occlusion")
            or not mechanisms.get("water_oxygen_selection")
            or not np.isfinite(float(water_cutoff)) or float(water_cutoff) <= 0):
        raise ValueError("sdc4_mechanisms.residue_water_cutoff_A requires residue_occlusion, water_oxygen_selection and a positive cutoff")
    if not isinstance(mechanisms.get("residue_occlusion", False), bool):
        raise ValueError("sdc4_mechanisms.residue_occlusion must be a boolean")
    surface_sites = mechanisms.get("residue_surface_sites", {})
    if not isinstance(surface_sites, dict) or any(
            not isinstance(name, str) or not name.strip()
            or not isinstance(selection, str) or not selection.strip()
            for name, selection in surface_sites.items()):
        raise ValueError("sdc4_mechanisms.residue_surface_sites must map site names to selections")
    environment_run = mechanisms.get("environment_sasa_run")
    if mechanisms.get("enabled") and environment_run is not None:
        sasa_runs = {run.run_id: run.spec for run in analysis_runs(cfg, "sasa")}
        reference = mechanisms.get("sasa_run", "default")
        if (not isinstance(environment_run, str) or environment_run not in sasa_runs
                or sasa_runs[environment_run].get("reference_run") != reference):
            raise ValueError("sdc4_mechanisms.environment_sasa_run must reference its sasa_run")
    if mechanisms.get("enabled"):
        tm = mechanisms.get("tm_chains") or {}
        ecto = mechanisms.get("ectodomains") or {}
        hs = mechanisms.get("hs_chains") or {}
        sites = mechanisms.get("hs_attachment_sites") or {}
        if not isinstance(tm, dict) or len(tm) != 2 or set(tm) != set(ecto):
            raise ValueError("sdc4_mechanisms needs two matching tm_chains and ectodomains")
        hs_auto = mechanisms.get("hs_auto") or {}
        if hs_auto.get("enabled"):
            if hs or sites:
                raise ValueError("Use hs_auto OR explicit hs_chains/hs_attachment_sites, not both")
            validate_glycan_spec(hs_auto, "sdc4_mechanisms.hs_auto", min_residues=2)
        elif not isinstance(hs, dict) or not isinstance(sites, dict) or set(hs) != set(sites):
            raise ValueError("sdc4_mechanisms.hs_chains and hs_attachment_sites must match")
        for group in (tm, ecto, hs, sites):
            for key, selection in group.items():
                _validate_path_component(str(key), "sdc4 mechanism group")
                if not isinstance(selection, str) or not selection.strip():
                    raise ValueError("SDC4 group selections must be nonempty strings")
        if not mechanisms.get("lipid_anchors") or not mechanisms.get("lipid_contact_selection"):
            raise ValueError("sdc4_mechanisms requires lipid_anchors and lipid_contact_selection")
        if not mechanisms.get("surface_anchor_species") or not isinstance(mechanisms["surface_anchor_species"], list):
            raise ValueError("sdc4_mechanisms.surface_anchor_species must be a nonempty list")
        if (hs or hs_auto.get("enabled")) and not mechanisms.get("water_oxygen_selection"):
            raise ValueError("HS hydration requires water_oxygen_selection")
        if str(mechanisms.get("axis", "z")).lower() != "z":
            raise ValueError("sdc4_mechanisms supports z-normal membranes only")
        if float(mechanisms.get("map_extent_A", 40)) <= 0 or float(mechanisms.get("map_bin_A", 2)) <= 0:
            raise ValueError("sdc4_mechanisms map extent and bin size must be positive")
        if int(mechanisms.get("local_min_lipids", 8)) < 8:
            raise ValueError("sdc4_mechanisms.local_min_lipids must be >= 8")
        if mechanisms.get("local_order", {}).get("enabled") and not mechanisms["local_order"].get("lipid_selection"):
            raise ValueError("sdc4_mechanisms.local_order.lipid_selection is required")

def validate_potassium_bridge_audit(cfg):
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        analyses = cfg["analyses"]
        solvation = analyses["solvation"]
        motif = solvation["targets"].get("motif", {}).get("acidic_patch", {})
        hs = analyses["sdc4_mechanisms"]["hs_auto"]
        if not (solvation.get("water_selection") and isinstance(motif, dict)
                and motif.get("sequence") and motif.get("selection")
                and hs.get("glycan_selection") and hs.get("attachment_selection")):
            raise ValueError("potassium_bridge_audit requires solvation water and acidic_patch "
                             "selectors and sdc4_mechanisms.hs_auto selections")
