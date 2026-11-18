"""Names of statistical families shared by execution and read-only planning.

The list contains identifiers only; importing it must not load estimators.
"""

STATISTICS_FAMILIES = (
    "time_series", "profiles", "interactions", "solvation_rdfs",
    "diffusion", "coupling",
)
