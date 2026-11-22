"""Smallest attainable permutation p-value for a replicate design.

Lightweight on purpose: schema-only planning imports this without loading the
statistics estimators. It mirrors ``statistics.common._permutation``, whose
two-sided p is the fraction of assignments with ``|T| >= |T_observed|``.
"""
from __future__ import annotations

import math


def permutation_floor(n_1: int, n_2: int, *, n_matched: int | None, tests: dict) -> tuple[float, str]:
    """Return ``(value, kind)`` with kind ``exact``, ``monte_carlo_resolution`` or ``unknown``.

    Exact enumeration: the most extreme assignment always ties with its mirror
    (relabelled groups, or every sign flipped) when groups are equal-sized or
    paired, giving ``2/total``. With unequal unpaired groups the mirror is not
    a valid assignment, so the best case is ``1/total``. Monte Carlo uses the
    ``(extreme + 1) / (draws + 1)`` estimator, whose resolution is not an exact
    design floor.
    """
    tests = tests or {}
    if tests.get("design", "independent") == "paired":
        if not n_matched:
            return float("nan"), "unknown"
        total, mirrored = 2 ** n_matched, True
    else:
        if not n_1 or not n_2:
            return float("nan"), "unknown"
        total, mirrored = math.comb(n_1 + n_2, n_1), n_1 == n_2
    if total < 2:
        return float("nan"), "unknown"
    if total <= int(tests.get("exact_max_assignments", 100000)):
        return (2 if mirrored else 1) / total, "exact"
    return 1 / (int(tests.get("permutation_draws", 10000)) + 1), "monte_carlo_resolution"
