"""Canonical config defaults, presets and deep merge semantics."""
from __future__ import annotations

from copy import deepcopy

DEFAULT_INSTANCE_ID = "default"
MULTI_RUN_ANALYSES = {"density_1d", "residue_membrane_distance", "sasa"}
# Analyses whose summaries accept ``analyses.<name>.groups`` (named selection groups).
GROUP_ANALYSES = ("rmsd", "radius_of_gyration")

# Opt-in acyl-tail naming presets (``membrane_order.tail_preset``); never
# inferred from residue names. CHARMM36 glycerophospholipids: glycerol C3
# carries the sn-1 chain (C31...), C2 the sn-2 chain (C21...). CHARMM36
# sphingolipids (PSM, CER*) have no sn positions: the sphingoid base is
# C1S...C18S from the headgroup carbon, the amide-linked acyl chain C1F
# (carbonyl)...C16F. They are separate tails, never pooled with sn-1/sn-2.
TAIL_PRESETS = {
    "charmm36": {
        "sn1": {"carbon_pattern": r"^C3(\d+)$"},
        "sn2": {"carbon_pattern": r"^C2(\d+)$"},
        "sphingosine": {"carbon_pattern": r"^C(\d+)S$"},
        "n_acyl": {"carbon_pattern": r"^C(\d+)F$"},
    },
}

DEFAULTS = {
    "project": {"name": "mdscope_project", "output": "mdscope_results"},
    "progress": {"enabled": True},
    "trajectory": {"start": None, "stop": None, "stride": 1, "chunk_frames": 1000},
    "time": {"dt_ps": None, "origin_ps": 0.0},
    "selections": {
        "protein": "protein",
        "alignment": "protein and name CA",
        "features": "protein and name CA",
        "ligand": None,
        "groups": {}, "motifs": {}, "glycans": {},
    },
    "analyses": {
        "alignment": {"enabled": True, "selection": None},
        "rmsd": {"enabled": True},
        "rmsf": {"enabled": True},
        "radius_of_gyration": {"enabled": True, "selection": None},
        "density_1d": {
            "enabled": False,
            "axis": "z",
            "bin_width_A": 1.0,
            "center_selection": None,
            "selections": {},
        },
        "membrane": {
            "enabled": False,
            "axis": "z",
            "headgroup_selection": None,
            "area_per_lipid": {"enabled": True},
            "thickness": {"enabled": True},
            "msd":{
                "enabled": False,
                "selection": None,
                "include_3d": False,
                "remove_drift": False,
                "slope_start_fraction": 0.2,
                "slope_end_fraction": 0.8,
                "multi_origin": False,
                "multi_origin_fit_diagnostics": False,
                "multi_origin_max_lag_frames": None,
                "multi_origin_particle_batch": 64
            },
        },
        "membrane_order": {
            "enabled": False,
            "axis": "z",
            "lipid_selection": None,
            "exclude_resnames": ["CHOL", "CHL1"],
            "hydrogen_cutoff_A": 1.25,
            "min_carbons_per_tail": 4,
            "per_lipid_type": True,
            "per_leaflet": True,
            "leaflet_anchor_selection": None,
            "tails": {
                "sn1": {"carbon_pattern": r"^C3(\d+)$"},
                "sn2": {"carbon_pattern": r"^C2(\d+)$"},
            },
            "lipid_overrides": {},
        },
        "residue_membrane_distance": {
            "enabled": False,
            "selection": None,
            "lipid_resnames":[
                "POPC",
                "POPE",
                "POPI",
                "DOPC",
                "DPPC",
                "PSM",
                "CHL1",
                "BGAL",
                "GM1a",
                "CER1"
            ],
            "membrane_selection": None,
            "exclude_hydrogens": True,
            "expand_residue_atoms": True,
            "initial_search_cutoff_A": 6.0,
            "search_growth_factor": 2.0,
            "fallback_block_atoms": 4096,
            "contact_cutoff_A": 5.0,
            "heatmap_max_A": None,
            "groups": False,
        },
        "dssp": {
            "enabled": True,
            "selection": None,
            "guess_hydrogens": True,
            "backend": "serial",
            "n_workers": None,
            "n_parts": None,
        },
        "bilayer_normal": {"enabled": False, "axis": "z", "species_selections": {},
                           "allow_absent_species": False,
                           "neighbor_radius_A": 15.0, "min_neighbors": 5,
                           "save_lipid_angles": False, "angle_chunk_frames": 100},
        "lipid_enrichment": {"enabled": False, "axis": "z", "species_selections": {},
                             "allow_absent_species": False,
                             "footprint_selection": None, "shell_edges_A": [0, 5, 10, 20],
                             "exclusion_radius_A": 3., "slab_half_width_A": 8.,
                             "area_grid_size": 64, "contact_cutoff_A": 5.},
        "lipid_hydrogen_bonds": {"enabled": False, "lipid_selection": None,
                                 "donors_sel": None, "hydrogens_sel": None, "acceptors_sel": None,
                                 "d_h_cutoff_A": 1.2, "d_a_cutoff_A": 3.0,
                                 "angle_cutoff_deg": 150.0,"exclude_same_residue": True,
                                 "exclude_same_molecule": True,
                                 "update_selections": False},
        "sdc4_mechanisms": {
            "enabled": False, "axis": "z",
            "protein_selection": None, "residue_occlusion": False, "residue_water_cutoff_A": None,
            "hs_sugar_water_cutoff_A": None,
            "tm_chains": {}, "ectodomains": {},
            "hs_chains": {}, "hs_attachment_sites": {},
            "hs_auto": {"enabled": False, "glycan_selection": None,
                        "attachment_selection": None,
                        "min_glycan_residues": 3, "expected_chains": None},
            "lipid_anchors": {}, "allow_absent_species": False,
            "surface_anchor_species": ["POPC", "POPE", "POPI", "POPS", "DOPC", "DPPC", "PSM"],
            "lipid_contact_selection": None,
            "water_oxygen_selection": None, "ion_selections": {},
            "motif_selection": None, "motifs": {"NxIP": "N.IP"},
            "require_motifs": True, "require_gxxxg": False,
            "map_extent_A": 40.0, "map_bin_A": 2.0,
            "local_radius_A": 20.0, "local_min_lipids": 8,
            "contact_cutoff_A": 5.0, "hydration_cutoff_A": 3.5,
            "ion_cutoff_A": 4.0,
            "local_order": {"enabled": False, "lipid_selection": None},
        },
        "potassium_bridge_audit": {"enabled": False},
        "shared_states": {
            "enabled": False, "features": None, "preset": "sdc4",
            "feature_sources": {}, "alignment": "exact", "scaling": "standardize",
            "n_states": 3, "seed": 2026, "bootstrap_draws": 300,
            "training_frames_per_replicate": 2000,
            "min_complete_fraction": 0.90, "block_frames": None,
        },
        "salt_bridges": {"enabled": True, "cutoff_A": 4.0},
        "ligand_contacts": {"enabled": False, "cutoff_A": 4.0, "save_atom_events": True},
        "hydrogen_bonds": {
            "enabled": False,
            "mode": "protein_protein",
            "d_h_cutoff_A": 1.2,
            "d_a_cutoff_A": 3.0,
            "angle_cutoff_deg": 150.0,
            "update_selections": False,
        },
        "sasa": {"enabled": False},
        "solvation": {"enabled": False, "preset": "standard",
                      "targets": {}, "water_selection": None, "ion_selection": None,
                      "ion_selections": {},
                      "search_radius_A": 8.0, "rdf_bin_width_A": 0.2,
                      "shell_cutoffs_A": {"water": 3.5, "ions": 4.0},
                      "options": {}},
        "hydration": {"enabled": False, "water_oxygen_selection": None,
                      "protein_selection": None, "membrane_selection": None,
                      "glycans": [], "cutoff_A": 3.5},
        "lifetimes": {"enabled": True, "gap_factor": 1.5},
    },
    "statistics": {
        "enabled": True,
        "residue_response": {"enabled": True, "profiles": None, "chains": None,
                             "min_length": 3, "minimum_cosine": 0.7},
        "bootstrap_draws": 2000,
        "ci": 0.95,
        "seed": 2026,
        "block": {
            "method": "moving",
            "length": "auto",
            "tau_multiplier": 2.0,
            "max_acf_lag": None,
        },
        "hierarchy": {
            "block_within_replicate": True,
            "resample_replicates": True,
        },
        "hypothesis_tests": {
            "enabled": True,
            "design": "independent",
            "method": "permutation",
            "permutation_draws": 10000,
            "exact_max_assignments": 100000,
            "multiple_testing": "fdr_bh",
        },
    },
    "cache": {"aligned_coordinates": True, "dtype": "float32", "keep": True},
    "plots": {
        "enabled": False,
        "include": None,
        "options": {},
        "min_occupancy": 0.05,
        "time_series": True,
        "time_series_rolling_window_ns": None,
        "distributions": False,
        "autocorrelation": False,
        "bootstrap_tables": False,
        "table_reports_pdf": False,
        "profiles": True,
        "interactions": True,
        "replicate_lightness_range": [0.28, 0.72],
        "font_sizes": {"base": 13, "title": 16, "suptitle": 18, "legend": 10, "ticks": 11},
        "dpi": 250,
    },
}


def merge(base: dict, update: dict) -> dict:
    out = deepcopy(base)
    for key, value in update.items():
        out[key] = merge(
            out[key],
            value
        ) if isinstance(value, dict) and isinstance(out.get(key), dict) else deepcopy(value)
    return out


