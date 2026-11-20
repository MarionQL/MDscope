"""Leaflet-resolved lipid contact counts, lipid classes and bootstrap inference.

Uses the *saved* molecule-unique pair counts; a frame, not an atom contact,
is the temporal observation. Statistical replication is by independent runs.
"""
from __future__ import annotations

from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import statistics_family_config
from ..core.persistence import atomic_parquet, ensure_dir
from .common import (_block_bootstrap, _ci, _hierarchical,
                     _hierarchical_difference, _temporal, _regular_times,
                     selected_rng)

METRICS = ("pair_count_undirected", "neighbors_per_molecule", "enrichment")
IDENTITY = ("level", "leaflet", "species_a", "species_b", "metric")
FRAME_ID = ("time_ps", "leaflet")
REQUIRED = {*FRAME_ID, "species_a", "species_b", "molecules_a", "molecules_b",
            "pair_count_undirected", "total_contact_pairs", "neighbors_per_molecule",
            "enrichment"}


def _validate_species(table: pd.DataFrame) -> None:
    missing = REQUIRED - set(table.columns)
    if missing:
        raise ValueError(f"lipid contact matrix missing fields: {sorted(missing)}")
    if table.empty:
        raise ValueError("lipid contact matrix contains no frames")
    if table.duplicated([*FRAME_ID, "species_a", "species_b"]).any():
        raise ValueError("duplicate lipid-contact matrix cell in a frame")
    if not np.isfinite(table.time_ps.to_numpy(float)).all():
        raise ValueError("lipid contact frame times must be finite")
    for leaflet, data in table.groupby("leaflet"):
        times = np.sort(data.time_ps.unique())
        if not _regular_times(times):
            raise ValueError(f"lipid contacts: irregular frame times in {leaflet}")
        groups = data.groupby("time_ps", sort=False)
        expected = len(data.species_a.unique()) ** 2
        if groups.size().ne(expected).any():
            raise ValueError(f"incomplete species contact matrix in {leaflet}")
        for field in ("total_contact_pairs",):
            if groups[field].nunique().gt(1).any():
                raise ValueError(f"inconsistent {field} within lipid frame")
        by_a = data.groupby([*FRAME_ID, "species_a"]).molecules_a.nunique()
        if by_a.gt(1).any():
            raise ValueError("inconsistent molecule population within lipid frame")
    if (table.pair_count_undirected < 0).any():
        raise ValueError("negative lipid pair count")


def _ordered_contacts(table: pd.DataFrame) -> pd.DataFrame:
    """One record per unordered species pair; verify AB and BA agree."""
    left = table.species_a.astype(str).to_numpy()
    right = table.species_b.astype(str).to_numpy()
    if (table.loc[:, ["species_a", "species_b"]].isna().any(axis=None)):
        raise ValueError("missing lipid species identity")
    symmetric = table.assign(
        first=np.minimum(left, right), second=np.maximum(left, right),
    )
    checks = symmetric.groupby([*FRAME_ID, "first", "second"])
    if checks.pair_count_undirected.nunique(dropna=False).gt(1).any():
        raise ValueError("AB/BA contact counts disagree")
    return (symmetric.loc[left <= right,
                          [*FRAME_ID, "first", "second", "pair_count_undirected"]]
            .rename(columns={"first": "species_a", "second": "species_b"}))


def _classes(table: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Aggregate *unique* pairs and recompute class denominators by frame.

    Distinct A-B species pairs belonging to the same class contribute once.
    A-A also contributes once, not twice. Then re-expand into a directed
    square matrix only for display and neighbors-per-reference-molecule.
    """
    species = set(table.species_a.astype(str)) | set(table.species_b.astype(str))
    missing = species - set(mapping)
    if missing:
        raise ValueError(f"lipid_contacts.classes missing species {sorted(missing)}")
    undirected = _ordered_contacts(table)
    undirected["class_a"] = undirected.species_a.map(mapping)
    undirected["class_b"] = undirected.species_b.map(mapping)
    aa = undirected.class_a.to_numpy(str)
    bb = undirected.class_b.to_numpy(str)
    undirected["class_a"], undirected["class_b"] = np.minimum(aa, bb), np.maximum(aa, bb)
    counts = (undirected.groupby([*FRAME_ID, "class_a", "class_b"], as_index=False)
              .pair_count_undirected.sum())
    population = (table.loc[:, [*FRAME_ID, "species_a", "molecules_a"]]
                  .drop_duplicates([*FRAME_ID, "species_a"]))
    population["class"] = population.species_a.map(mapping)
    population = (population.groupby([*FRAME_ID, "class"], as_index=False)
                  .molecules_a.sum())
    total_contacts = (table.groupby(list(FRAME_ID), as_index=False)
                      .total_contact_pairs.first())
    total_n = population.groupby(list(FRAME_ID), as_index=False).molecules_a.sum()
    total_n = total_n.rename(columns={"molecules_a": "total_molecules"})
    classes = sorted(set(mapping[s] for s in species))
    frames = total_contacts.merge(total_n, on=list(FRAME_ID), validate="one_to_one")
    desired = pd.DataFrame(list(product(classes, classes)),
                           columns=["species_a", "species_b"])
    complete = frames.merge(desired, how="cross")
    pop_a = population.rename(columns={"class": "species_a", "molecules_a": "molecules_a"})
    pop_b = population.rename(columns={"class": "species_b", "molecules_a": "molecules_b"})
    complete = complete.merge(pop_a, on=[*FRAME_ID, "species_a"], validate="many_to_one")
    complete = complete.merge(pop_b, on=[*FRAME_ID, "species_b"], validate="many_to_one")
    complete["class_a"] = np.minimum(complete.species_a, complete.species_b)
    complete["class_b"] = np.maximum(complete.species_a, complete.species_b)
    complete = complete.merge(counts, on=[*FRAME_ID, "class_a", "class_b"],
                              how="left", validate="many_to_one")
    complete.pair_count_undirected = complete.pair_count_undirected.fillna(0)
    same = complete.species_a.eq(complete.species_b)
    a, b = complete.molecules_a.to_numpy(float), complete.molecules_b.to_numpy(float)
    possible = np.where(same, a * (a - 1) / 2, a * b)
    n = complete.total_molecules.to_numpy(float)
    possible_total = n * (n - 1) / 2
    expected = np.divide(
        complete.total_contact_pairs.to_numpy(float) * possible, possible_total,
        out=np.full(len(complete), np.nan), where=possible_total > 0,
    )
    cnt = complete.pair_count_undirected.to_numpy(float)
    complete["neighbors_per_molecule"] = np.divide(
        cnt * np.where(same, 2, 1), a,
        out=np.full(len(complete), np.nan), where=a > 0,
    )
    complete["enrichment"] = np.divide(
        cnt, expected, out=np.full(len(complete), np.nan),
        where=expected > 0,
    )
    # Absent class is not observed zero. Keep count missing as well.
    complete.loc[(a == 0) | (b == 0), list(METRICS)] = np.nan
    return complete.loc[:, [*FRAME_ID, "species_a", "species_b", *METRICS]]


def contact_frames(table: pd.DataFrame, classes: dict[str, list[str]]) -> pd.DataFrame:
    """Return symmetric species/class matrices with no duplicated counting."""
    _validate_species(table)
    mapping = {}
    for name, members in classes.items():
        if not isinstance(members, list) or not members:
            raise ValueError(f"lipid class {name!r} must contain species names")
        for species in members:
            if species in mapping:
                raise ValueError(f"lipid species {species!r} belongs to multiple classes")
            mapping[species] = str(name)
    data = table.loc[:, [*FRAME_ID, "species_a", "species_b", *METRICS]].copy()
    # A species with zero molecules has no pair-count estimand.
    missing = (table.molecules_a.to_numpy(float) <= 0) | (table.molecules_b.to_numpy(float) <= 0)
    data.loc[missing, list(METRICS)] = np.nan
    data.insert(0, "level", "species")
    class_data = _classes(table, mapping)
    class_data.insert(0, "level", "class")
    return pd.concat([data, class_data], ignore_index=True)


def _prepare_replicate(frame: pd.DataFrame, cfg: dict, draws: int, seed: int, identity):
    """Frame-complete features share a conservative leaflet block length."""
    times = np.sort(frame.time_ps.unique())
    if not _regular_times(times):
        raise ValueError(f"irregular lipid-contact times: {identity}")
    long = frame.melt(id_vars=[*FRAME_ID, "level", "species_a", "species_b"],
                      value_vars=list(METRICS), var_name="metric", value_name="value")
    values = long.pivot(index="time_ps", columns=list(IDENTITY), values="value")
    values = values.reindex(times)
    if values.index.has_duplicates or not len(values):
        raise ValueError("invalid lipid contact frame grid")
    matrix = values.to_numpy(float)
    valid = np.isfinite(matrix).all(axis=0)
    # Determine autocorrelation from each valid channel. Sharing the maximum
    # block size within one leaflet is deliberately conservative and keeps
    # resampling coherent across contact classes.
    block = 1
    for col in np.flatnonzero(valid):
        block = max(block, int(_temporal(matrix[:, col], times, cfg)["block_length_frames"]))
    rng = selected_rng(seed, "lipid_contacts", *identity)
    boot = (np.full((draws, matrix.shape[1]), np.nan)
            if draws else np.empty((0, matrix.shape[1])))
    if draws and valid.any():
        boot[:, valid] = _block_bootstrap(matrix[:, valid], block, draws, rng)
    points = np.where(valid, matrix.mean(axis=0), np.nan)
    return list(values.columns), points, boot, valid, block


def analyze_lipid_contacts(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    cfg = statistics_family_config(cfg, "lipid_contacts")
    if cfg is None or not cfg.get("analyses", {}).get("lipid_contacts", {}).get("enabled"):
        return
    settings = cfg["statistics"]
    draws = int(settings.get("bootstrap_draws", 2000))
    level = float(settings.get("ci", .95))
    seed = int(settings.get("seed", 2026))
    resample = settings.get("hierarchy", {}).get("resample_replicates", True)
    classes = cfg["analyses"]["lipid_contacts"].get("classes") or {}
    if not classes:
        raise ValueError("analyses.lipid_contacts.classes is required for class-level statistics")
    outdir = ensure_dir(Path(root) / "project_summary/statistics/lipid_contacts")
    prepared = {}
    rep_rows = []
    for system, reps in system_repouts.items():
        for repout in reps:
            source = Path(repout) / "extensions/lipid_contacts/matrix.parquet"
            if not source.is_file():
                raise FileNotFoundError(f"Missing lipid contact matrix: {source}")
            frame = contact_frames(pd.read_parquet(source), classes)
            for (group_level, leaflet), sub in frame.groupby(["level", "leaflet"], sort=True):
                identity = (system, Path(repout).name, group_level, leaflet)
                keys, points, boot, valid, block = _prepare_replicate(
                    sub, cfg, draws, seed, identity,
                )
                prepared[identity] = dict(keys=keys, point=points,
                                          boot=boot, valid=valid)
                for i, (_, _, a, b, metric) in enumerate(keys):
                    low, high = _ci(boot[:, i], level)
                    rep_rows.append(dict(system=system, replicate=Path(repout).name,
                                         level=group_level, leaflet=leaflet,
                                         species_a=a, species_b=b, metric=metric,
                                         mean=float(points[i]), mean_ci_low=low,
                                         mean_ci_high=high, block_length_frames=block,
                                         n_frames=int(len(sub.time_ps.unique())),
                                         bootstrap_draws=draws, ci_level=level))
    if not prepared:
        return
    systems = {}
    rows = []
    for system, reps in system_repouts.items():
        for group_level, leaflet in sorted({(k[2], k[3]) for k in prepared if k[0] == system}):
            keys = sorted({feature for rep in reps
                           for feature in prepared[(system, Path(rep).name, group_level, leaflet)]["keys"]})
            vectors = []
            sample = []
            complete = np.ones(len(keys), bool)
            for rep in reps:
                item = prepared[(system, Path(rep).name, group_level, leaflet)]
                mapping = {key: i for i, key in enumerate(item["keys"])}
                idx = [mapping.get(k) for k in keys]
                complete &= np.array([i is not None and item["valid"][i] for i in idx])
                vectors.append(np.array([item["point"][i] if i is not None else np.nan
                                         for i in idx]))
                sample.append(np.column_stack([
                    item["boot"][:, i] if i is not None else np.full(draws, np.nan)
                    for i in idx
                ]))
            point = np.mean(np.stack(vectors), axis=0)
            point[~complete] = np.nan
            rng = selected_rng(seed, "lipid_contacts", system, group_level, leaflet, "hierarchy")
            boot = _hierarchical(sample, draws, rng, "mean", resample)
            if len(boot):
                boot[:, ~complete] = np.nan
            systems[system, group_level, leaflet] = (keys, point, boot, complete)
            for i, (_, _, a, b, metric) in enumerate(keys):
                low, high = _ci(boot[:, i], level)
                rows.append(dict(system=system, level=group_level, leaflet=leaflet,
                                 species_a=a, species_b=b, metric=metric,
                                 mean=float(point[i]), mean_ci_low=low,
                                 mean_ci_high=high, n_replicates=len(reps),
                                 bootstrap_draws=draws, ci_level=level))
    contrast = []
    for first, second in combinations(system_repouts, 2):
        for group_level, leaflet in sorted({(key[1], key[2]) for key in systems if key[0] == first}):
            if (second, group_level, leaflet) not in systems:
                continue
            ka, ma, ba, va = systems[first, group_level, leaflet]
            kb, mb, bb, vb = systems[second, group_level, leaflet]
            map_a, map_b = {key:i for i,key in enumerate(ka)}, {key:i for i,key in enumerate(kb)}
            for key in sorted(set(ka) | set(kb)):
                ia, ib = map_a.get(key), map_b.get(key)
                available = (ia is not None and ib is not None and va[ia] and vb[ib])
                diff = float(mb[ib] - ma[ia]) if available else np.nan
                if available and draws:
                    if settings.get("hypothesis_tests", {}).get("design") == "paired":
                        # Preserve deliberately matched replica labels.
                        a = {}
                        b = {}
                        for rep in system_repouts[first]:
                            name = Path(rep).name
                            item = prepared[first, name, group_level, leaflet]
                            index = {k:i for i,k in enumerate(item["keys"])}.get(key)
                            if index is not None:
                                a[name] = item["boot"][:, index][:, None]
                        for rep in system_repouts[second]:
                            name = Path(rep).name
                            item = prepared[second, name, group_level, leaflet]
                            index = {k:i for i,k in enumerate(item["keys"])}.get(key)
                            if index is not None:
                                b[name] = item["boot"][:, index][:, None]
                        delta = _hierarchical_difference(
                            a, b, draws, selected_rng(seed, "lipid_contacts", first,
                                                       second, group_level, leaflet, *key),
                            "mean", resample,
                        ).ravel()
                    else:
                        delta = bb[:, ib] - ba[:, ia]
                    low, high = _ci(delta, level)
                    pgt = float(np.mean(delta > 0)) if np.isfinite(delta).all() else np.nan
                else:
                    low, high, pgt = np.nan, np.nan, np.nan
                _, _, a, b, metric = key
                contrast.append(dict(system_1=first, system_2=second,
                                     level=group_level, leaflet=leaflet,
                                     species_a=a, species_b=b, metric=metric,
                                     difference_system2_minus_system1=diff,
                                     difference_ci_low=low, difference_ci_high=high,
                                     probability_difference_gt_0=pgt,
                                     ci_excludes_zero=bool(np.isfinite(low) and
                                                           np.isfinite(high) and
                                                           (low > 0 or high < 0)),
                                     bootstrap_draws=draws, ci_level=level))
    atomic_parquet(outdir / "replicate_statistics.parquet", pd.DataFrame(rep_rows))
    atomic_parquet(outdir / "system_statistics.parquet", pd.DataFrame(rows))
    atomic_parquet(outdir / "system_comparisons.parquet", pd.DataFrame(contrast))
