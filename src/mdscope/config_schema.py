"""Reject unknown YAML option names without restricting user-defined identifiers.

This checks *names* before the existing semantic and topology validation. Templates
come from config.DEFAULTS; the additions below describe supported optional
options and mapping keys whose names are supplied by the user (not the program).
"""

from __future__ import annotations

from copy import deepcopy
from difflib import get_close_matches


def _extend(mapping, **fields):
    mapping.update(fields)
    return mapping


def option_schema(defaults):
    schema = deepcopy(defaults)
    schema["systems"] = {
        "*": {
            "label": None,
            "topology": None,
            "replicates": [{
                "trajectory": None,
                "label": None,
                "plot_label": None,
                "msd_trajectory": None,
            }],
        }
    }
    schema["selections"]["groups"] = {
        "*": {"selection": None, "alignment_selection": None,
              "label": None, "linestyle": None}
    }

    analyses = schema["analyses"]
    for name in analyses:
        analyses[name]["label"] = None

    # Valid analyses not represented in DEFAULTS have their own executors.
    analyses["lipid_contacts"] = {
        "enabled": None, "label": None, "axis": None,
        "allow_absent_species": None, "selection": None,
        "species_selections": {"*": None},
        "cutoff_A": None, "batch_atoms": None,
    }
    analyses["local_environment"] = {
        "enabled": None, "label": None,
        "protein_selection": None, "lipid_selection": None,
        "headgroup_selection": None, "glycan_selection": None,
        "ion_selection": None, "anchor_selection": None,
        "contact_cutoff_A": None, "motifs": {"*": None},
    }
    analyses["coupling_screen"] = {
        "enabled": None, "label": None, "include_defaults": None,
        "pairs": None, "block_frames": None,
    }

    for name in ("rmsd", "radius_of_gyration"):
        analyses[name]["groups"] = None

    for name in ("density_1d", "residue_membrane_distance", "sasa"):
        analyses[name]["runs"] = {"*": None}  # Each run is checked below.

    density = analyses["density_1d"]
    _extend(density, membrane_center_selection=None, half_width_A=None)
    density["selections"] = {
        "*": {"selection": None, "label": None, "mode": None, "groups": None}
    }
    analyses["membrane"]["msd"]["multi_origin_window_ps"] = None

    tail = {"carbon_pattern": None, "carbon_names": None, "positions": None}
    order = analyses["membrane_order"]
    order["tails"] = {"*": tail}
    order["lipid_overrides"] = {"*": {"*": tail}}

    for name in ("bilayer_normal", "lipid_enrichment"):
        analyses[name]["species_selections"] = {"*": None}

    _extend(analyses["hydrogen_bonds"],
            between=None, donors_sel=None, hydrogens_sel=None,
            acceptors_sel=None, restrict_search_to_protein=None)
    _extend(analyses["dssp"], heavyatom_names=None, hydrogen_name=None)

    sasa = analyses["sasa"]
    _extend(sasa, selection=None, target_selection=None, occluder_selection=None,
            atom_sasa_selection=None, atom_sasa_sequences=None,
            atom_sasa_expected_matches=None, reference_run=None,
            algorithm=None, periodic=None, probe_radius_A=None,
            n_points=None, n_slices=None, n_threads=None)
    sasa["element_radii_A"] = {"*": None}

    mechanism = analyses["sdc4_mechanisms"]
    _extend(mechanism, sasa_run=None, environment_sasa_run=None,
            motif_lipid_residence=None, accessible_area_enabled=None,
            accessible_area_grid_A=None, protein_footprint_exclusion_A=None,
            protein_footprint_slab_A=None)
    for field in ("tm_chains", "ectodomains", "hs_chains",
                  "hs_attachment_sites", "lipid_anchors", "ion_selections",
                  "motifs", "residue_surface_sites", "lipid_moiety_selections"):
        mechanism[field] = {"*": None}
    _extend(mechanism["local_order"],
            by_species_leaflet=None, min_carbons_per_tail=None)

    solvation = analyses["solvation"]
    _extend(solvation, sasa_run=None, neighbor_method=None,
            grid_width_A=None, grid_spacing_A=None)
    solvation["hs_auto"] = deepcopy(mechanism["hs_auto"])
    solvation["ion_selections"] = {"*": None}
    solvation["shell_cutoffs_A"] = {"*": None}
    solvation["targets"] = {
        "*": {"*": {
            "selection": None, "sequence": None,
            "atom_selection": None, "group_by": None,
        }}
    }
    solvation["aqueous_reference"] = {
        "membrane_selection": None, "axis": None, "half_width_A": None,
    }
    solvation["grid_frames"] = {"*": None}
    solvation["options"] = dict.fromkeys((
        "events", "bridges", "water_hbond_bridges",
        "track_target_classes", "bridge_target_pairs", "bridge_spatial_ion",
        "hbond_da_cutoff_A", "hbond_angle_deg",
    ))

    stats = schema["statistics"]
    _extend(stats, bootstrap=None, scalar_window_ps=None,
            matched_window_ps=None, interaction_min_occupancy=None,
            interaction_include_kinds=None, interaction_exclude_kinds=None,
            interaction_progress_every=None, adaptive_hierarchy=None,
            metric_selection={"time_series": [{"metric": None, "group": None, "inference": None, "ci": None, "hypothesis_tests": None}],
                              "interactions": [{"kind": None, "interaction": None, "metric": None, "inference": None, "ci": None, "hypothesis_tests": None}],
                              "profiles": [{"profile": None, "inference": None, "ci": None, "hypothesis_tests": None}],
                              "solvation_rdfs": [{"solvent": None, "target_class": None,
                                                  "site_class": None, "inference": None, "ci": None}],
                              "diffusion": [{"metric": None, "inference": None, "ci": None}],
                              "coupling": [{"mechanism": None, "inference": None, "ci": None}]})
    stats["hierarchy"]["enabled"] = None
    family = {
        "enabled": None, "bootstrap": None, "bootstrap_draws": None,
        "ci": None, "hypothesis_tests": None,
        "hierarchy": {"enabled": None, "block_within_replicate": None,
                      "resample_replicates": None},
        "adaptive_hierarchy": {
            "min_draws": None, "max_draws": None, "batch_size": None,
            "stable_checks": None, "absolute_tolerance": None,
            "relative_tolerance": None,
        },
    }
    stats["families"] = {"*": family}
    stats["adaptive_hierarchy"] = deepcopy(family["adaptive_hierarchy"])
    stats["hypothesis_families"] = {"*": None}

    plots = schema["plots"]
    for name in ("solvation", "sdc4", "shared_states", "diffusion",
                 "lipid_enrichment", "lipid_contacts", "other"):
        plots[name] = None
    plots["include"] = None
    plots["options"] = {"*": {"dpi": None}}
    plots["families"] = {"*": None}  # Existing validator checks family names.
    plots["interaction_exclude_kinds"] = None

    return schema


def _check_keys(value, allowed, path):
    if not isinstance(value, dict) or not isinstance(allowed, dict):
        return
    for key, nested in value.items():
        if key in allowed:
            child = allowed[key]
        elif "*" in allowed:
            child = allowed["*"]
        else:
            full_path = f"{path}.{key}" if path else str(key)
            closest = get_close_matches(str(key), [str(k) for k in allowed], n=1)
            hint = (f"; did you mean {path + '.' if path else ''}{closest[0]}?"
                    if closest else "")
            raise ValueError(f"Unknown YAML option {full_path}{hint}")
        child_path = f"{path}.{key}" if path else str(key)
        if path.startswith("analyses.") and path.count(".") == 1 and key == "runs":
            # Run names are user-defined but their options match the parent
            # analysis, excluding recursive runs. Booleans remain supported.
            if isinstance(nested, dict):
                run_options = {k: v for k, v in allowed.items() if k != "runs"}
                for run_id, run_spec in nested.items():
                    _check_keys(run_spec, run_options, f"{child_path}.{run_id}")
        elif isinstance(nested, list) and isinstance(child, list) and child:
            for i, element in enumerate(nested):
                _check_keys(element, child[0], f"{child_path}[{i}]")
        else:
            _check_keys(nested, child, child_path)


def validate_known_keys(raw, defaults):
    """Reject misspelled options; existing config validation checks their values."""
    if not isinstance(raw, dict):
        return
    _check_keys(raw, option_schema(defaults), "")
