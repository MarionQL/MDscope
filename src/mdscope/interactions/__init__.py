"""Interaction-table, contact and hydrogen-bond API.

Explicit named exports prevent unrelated pandas/NumPy/persistence imports from
escaping via this package. Scientific implementations live in their own modules.
"""
from .interaction_tables import (
    atom_label,
    chain_label,
    count_timeseries,
    event_lifetimes,
    interaction_key_columns,
    interaction_tables,
    lifetime_summary,
    load_parts,
    postprocess_interactions,
    residue_label,
    sample_dt,
    summarize_interaction_kind,
)
from .salt_bridges import initialize_salt_bridges, salt_bridge_events
from .ligand_contacts import (
    initialize_ligand_contacts,
    _minimum_distance,
    ligand_contact_events,
)
from .hydrogen_bonds import run_hydrogen_bonds

__all__ = [
    "atom_label",
    "chain_label",
    "count_timeseries",
    "event_lifetimes",
    "interaction_key_columns",
    "interaction_tables",
    "lifetime_summary",
    "load_parts",
    "postprocess_interactions",
    "residue_label",
    "sample_dt",
    "summarize_interaction_kind",
    "initialize_salt_bridges",
    "salt_bridge_events",
    "initialize_ligand_contacts",
    "ligand_contact_events",
    "run_hydrogen_bonds",
]
