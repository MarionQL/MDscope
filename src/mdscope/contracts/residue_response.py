"""Supported profile identities for residue-response screening."""
import re

# Residue-resolved protein observables with comparable saved labels. Shared by
# YAML validation and the saved-data screen so eligibility has one definition.
RESIDUE_RESPONSE_PROFILE = re.compile(
    r"^(?:rmsf$|dssp_(?:helix|sheet|coil)_occupancy$|"
    r"sasa_mean(?:__|$)|sasa_buried_nm2__|"
    r"residue_occlusion__|residue_water_coordination__all$|"
    r"residue_membrane_distance__(?!.*__group_))"
)

