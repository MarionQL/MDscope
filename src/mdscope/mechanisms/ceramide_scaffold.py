"""Chemistry-matched, per-frame ceramide scaffold contact pooling.

Pool GM1a+BGAL molecule numerators/denominators *before* forming ratios.
CER160 is its own cohort. Glycan-only contacts are deliberately excluded.
"""
import numpy as np
import pandas as pd

COLUMNS = (
    "time_ps", "region", "leaflet", "moiety", "cohort",
    "available_molecules", "locally_available_molecules",
    "contacted_molecules", "locally_contacted_molecules",
    "contact_fraction", "matched_local_contact_fraction",
)


def pool_ceramide_scaffolds(specificity: pd.DataFrame):
    """Return comparable ceramide_head/tail fractions without fictitious zeros."""
    needed = {
        "time_ps", "region", "leaflet", "species", "moiety",
        "available_molecules", "locally_available_molecules",
        "contacted_molecules", "locally_contacted_molecules",
    }
    if not needed.issubset(specificity):
        raise ValueError(f"Ceramide pooling requires {sorted(needed - set(specificity))}")
    selected = specificity.loc[
        specificity.species.isin(["GM1a", "BGAL", "CER160"]) &
        specificity.moiety.isin(["ceramide_head", "ceramide_tail"]),
        list(needed),
    ].copy()
    if selected.empty:
        return pd.DataFrame(columns=COLUMNS)
    key = ["time_ps", "region", "leaflet", "species", "moiety"]
    if selected.duplicated(key).any():
        raise ValueError("Duplicate ceramide specificity molecule counts")
    selected["cohort"] = np.where(
        selected.species == "CER160", "CER160", "GM1a+BGAL"
    )
    numeric = ["available_molecules", "locally_available_molecules",
               "contacted_molecules", "locally_contacted_molecules"]
    if not np.isfinite(selected[numeric].to_numpy(float)).all() or (
        selected[numeric].to_numpy(float) < 0
    ).any():
        raise ValueError("Ceramide scaffold molecule counts must be finite and nonnegative")
    pooled = selected.groupby(
        ["time_ps", "region", "leaflet", "moiety", "cohort"],
        sort=True, as_index=False,
    )[numeric].sum()
    for numerator, denominator, field in (
        ("contacted_molecules", "available_molecules", "contact_fraction"),
        ("locally_contacted_molecules", "locally_available_molecules",
         "matched_local_contact_fraction"),
    ):
        n = pooled[numerator].to_numpy(float)
        d = pooled[denominator].to_numpy(float)
        if np.any(n > d):
            raise ValueError("Ceramide contacts exceed available molecules")
        pooled[field] = np.divide(
            n, d, out=np.full(len(n), np.nan), where=d > 0
        )
    return pooled[list(COLUMNS)]
