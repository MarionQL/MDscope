"""Lipid-molecule residence from sparse motif contact events."""
from ..core.residence import summarize_contact_residence


def summarize_lipid_residence(events, frame_times, gap_frames=0):
    """Sampled lipid episodes; gaps do not imply continuous physical binding."""
    return summarize_contact_residence(
        events, frame_times, identity_columns=("motif",),
        species_column="species", molecule_column="molecule_id", gap_frames=gap_frames,
    )


def summarize_intermittent_residence(events, frame_times, gap_frames=1):
    """Gap-tolerant episodes, reported separately from continuous ones."""
    return summarize_lipid_residence(events, frame_times, gap_frames=gap_frames)
