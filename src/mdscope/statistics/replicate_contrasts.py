"""Intervals for estimands computed once per independent replica."""
from hashlib import sha256
from itertools import combinations

import numpy as np
import pandas as pd


def replicate_mean_contrasts(values, keys, system_repouts, stats):
    """Bootstrap complete replica estimates for each system pair.

    These intervals propagate between-replica variability. Their input is
    one value per independent trajectory, not a temporal bootstrap sample.
    """
    if values.duplicated([*keys, "system", "replicate"]).any():
        raise ValueError("Duplicate replica estimates for a contrast")
    expected = {system: {path.name for path in paths} for system, paths in system_repouts.items()}
    paired = stats.get("hypothesis_tests", {}).get("design") == "paired"
    hierarchy = stats.get("hierarchy", {})
    draws = int(stats.get("bootstrap_draws", 2000))
    level = float(stats.get("ci", .95))
    seed = int(stats.get("seed", 2026))
    bootstrap = (draws > 0 and stats.get("bootstrap", True)
                 and hierarchy.get("enabled", True) and hierarchy.get("resample_replicates", True))
    rows = []
    resamples = {}
    for identity, group in values.groupby(keys, sort=True):
        identity = identity if isinstance(identity, tuple) else (identity,)
        by_system = {
            system: table.loc[np.isfinite(table.value)].set_index("replicate").value
            for system, table in group.groupby("system", sort=False)
        }
        for first, second in combinations(expected, 2):
            a = by_system.get(first, pd.Series(dtype=float)).sort_index()
            b = by_system.get(second, pd.Series(dtype=float)).sort_index()
            if a.empty and b.empty and not group.system.isin([first, second]).any():
                continue
            row = dict(zip(keys, identity), system_1=first, system_2=second,
                       n_replicates_1=len(a), n_replicates_2=len(b),
                       expected_replicates_1=len(expected[first]), expected_replicates_2=len(expected[second]),
                       difference_system2_minus_system1=np.nan, difference_ci_low=np.nan,
                       difference_ci_high=np.nan, ci_level=level, bootstrap_draws=0,
                       ci_method="replicate_percentile", ci_status=None)
            if paired:
                common = a.index.intersection(b.index)
                left, right = a.loc[common].to_numpy(), b.loc[common].to_numpy()
            else:
                left, right = a.to_numpy(), b.to_numpy()
            if len(left) and len(right):
                row["difference_system2_minus_system1"] = float(right.mean() - left.mean())
            if set(a.index) != expected[first] or set(b.index) != expected[second]:
                row["ci_status"] = "incomplete_replica_coverage"
            elif min(len(left), len(right)) < 2 or (paired and not a.index.equals(b.index)):
                row["ci_status"] = "insufficient_replicates"
            elif not bootstrap:
                row["ci_status"] = "bootstrap_disabled"
            else:
                # Use the same replica resamples for each observable with
                # these system/replica keys so correlated outputs stay aligned.
                key = (first, second, tuple(a.index), tuple(b.index))
                if key not in resamples:
                    salt = int.from_bytes(sha256(repr(key).encode()).digest()[:8], "little")
                    rng = np.random.default_rng(np.random.SeedSequence([seed, salt]))
                    i = rng.integers(len(left), size=(draws, len(left)))
                    j = i if paired else rng.integers(len(right), size=(draws, len(right)))
                    resamples[key] = i, j
                i, j = resamples[key]
                boot = right[j].mean(axis=1) - left[i].mean(axis=1)
                alpha = (1 - level) / 2
                row["difference_ci_low"], row["difference_ci_high"] = np.quantile(boot, [alpha, 1 - alpha])
                row["bootstrap_draws"] = draws
            rows.append(row)
    return pd.DataFrame(rows)
