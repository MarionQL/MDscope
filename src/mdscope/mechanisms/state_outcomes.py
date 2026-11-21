"""SDC4-only descriptive outcomes downstream of shared geometric states."""
import numpy as np
import pandas as pd
from ..core.persistence import atomic_parquet


class SDC4StateOutcomes:
    """Accumulate preset outcomes without changing generic feature training."""

    def __init__(self):
        self.conditioned = []
        self.lipid_state_rows = []
        self.nxip_lipid_rows = []
        self.nxip_hs_rows = []
        self.nxip_membrane_rows = []
        self.nxip_shift_null_rows = []
        self.lipid_replica_rows = []
        self.nxip_conditioned = []
        self.nxip_joint = []

    def add(self, df, *, system, replicate, repout, labels, times, n_states, temporal_valid):
        k = n_states
        # Conditional NxIP outcomes: define states exclusively from the
        # configured geometric features; do not include NxIP as a predictor.
        # This table is descriptive within each independent replicate.
        for chain in ("A", "B"):
            prefix = f"sdc4_motif_NxIP_PRO{chain}_"
            endpoints = [
                col for col in df.columns if col.startswith(prefix)
                and (col.endswith("_steric_occlusion_fraction")
                     or col.endswith("_sasa_nm2")
                     or col.endswith("_glycan_contacts")
                     or col.endswith("_protein_neighbor_residues"))
            ]
            for endpoint in endpoints:
                values = df[endpoint].to_numpy(float)
                for state in range(k):
                    mask = (labels == state) & np.isfinite(values)
                    other = (labels != state) & np.isfinite(values)
                    self.nxip_conditioned.append({
                        "system": system, "replicate": replicate,
                        "chain": chain, "state": state, "endpoint": endpoint,
                        "frames_in_state": int(mask.sum()),
                        "frames_outside_state": int(other.sum()),
                        "mean_in_state": float(values[mask].mean()) if mask.any() else np.nan,
                        "mean_outside_state": float(values[other].mean()) if other.any() else np.nan,
                        "within_replicate_difference": (
                            float(values[mask].mean() - values[other].mean())
                            if mask.any() and other.any() else np.nan
                        ),
                    })
        # Threshold sensitivity, not a data-fitted "binding state".
        # 0.25/0.50/0.75 are prespecified occluded-heavy-atom fractions.
        for chain in ("A", "B"):
            endpoint = f"sdc4_motif_NxIP_PRO{chain}_1_steric_occlusion_fraction"
            if endpoint not in df:
                continue
            values = df[endpoint].to_numpy(float)
            finite = np.isfinite(values)
            for threshold in (.25, .50, .75):
                exposed = values <= threshold
                for state in range(k):
                    mask = finite & (labels == state)
                    n = int(mask.sum())
                    self.nxip_joint.append({
                        "system": system, "replicate": replicate,
                        "chain": chain, "state": state, "threshold": threshold,
                        "frames_in_state": n,
                        "frames_exposed_in_state": int((mask & exposed).sum()),
                        "exposed_fraction_given_state": (
                            float((mask & exposed).sum() / n) if n else np.nan
                        ),
                    })
        # Circular-shift control preserves each NxIP trace's temporal
        # autocorrelation but breaks exact framewise alignment to states.
        # Exploratory only: nonstationary trajectories can invalidate this null.
        for chain in ("A", "B"):
            endpoint = f"sdc4_motif_NxIP_PRO{chain}_1_steric_occlusion_fraction"
            if endpoint not in df:
                continue
            y = df[endpoint].to_numpy(float)
            if len(y) < 4 or not temporal_valid:
                continue
            shifts = np.unique(np.linspace(
                1, len(y) - 1, min(100, len(y) - 1), dtype=int
            ))
            for state in range(k):
                state_mask = labels == state
                mask = np.isfinite(y)
                if not np.any(state_mask & mask) or not np.any(~state_mask & mask):
                    continue
                observed = float(y[state_mask & mask].mean() -
                                 y[~state_mask & mask].mean())
                null = []
                for shift in shifts:
                    shifted = np.roll(y, int(shift))
                    valid = np.isfinite(shifted)
                    if (state_mask & valid).any() and (~state_mask & valid).any():
                        null.append(float(
                            shifted[state_mask & valid].mean() -
                            shifted[~state_mask & valid].mean()
                        ))
                if not null:
                    continue
                lo, hi = np.quantile(null, (.025, .975))
                self.nxip_shift_null_rows.append({
                    "system": system, "replicate": replicate,
                    "chain": chain, "state": state,
                    "observed_difference": observed,
                    "shift_null_ci_low": float(lo),
                    "shift_null_ci_high": float(hi),
                    "n_shifts": len(null),
                })
        # Within-replica NxIP exposure conditional on ectodomain–membrane
        # contact: no frame-level p-values and no fitted cutoffs.
        for chain in ("A", "B"):
            nxip_name = f"sdc4_motif_NxIP_PRO{chain}_1_steric_occlusion_fraction"
            contact_name = f"sdc4_ectodomain_{chain}_membrane_contact_state"
            if nxip_name not in df or contact_name not in df:
                continue
            nxip = df[nxip_name].to_numpy(float)
            contacts = df[contact_name].to_numpy(float)
            available = np.isfinite(nxip) & np.isfinite(contacts)
            for contact_value in (0, 1):
                mask = available & (contacts == contact_value)
                self.nxip_membrane_rows.append({
                    "system": system, "replicate": replicate, "chain": chain,
                    "membrane_contact_state": contact_value,
                    "frames": int(mask.sum()),
                    "mean_nxip_occlusion": float(nxip[mask].mean()) if mask.any() else np.nan,
                    "nxip_exposed_fraction_threshold_0_5": (
                        float(np.mean(nxip[mask] <= .5)) if mask.any() else np.nan
                    ),
                })
        # HS-self.conditioned NxIP occlusion uses molecule-resolved HS contacts
        # already computed during the same trajectory pass.
        for chain in ("A", "B"):
            endpoint = f"sdc4_motif_NxIP_PRO{chain}_1_steric_occlusion_fraction"
            hs_prefix = f"sdc4_motif_NxIP_PRO{chain}_1_hs_"
            if endpoint not in df:
                continue
            for hs_col in (col for col in df.columns
                           if col.startswith(hs_prefix)
                           and col.endswith("_contact_state")):
                occlusion = df[endpoint].to_numpy(float)
                contact = df[hs_col].to_numpy(float)
                valid = np.isfinite(occlusion) & np.isfinite(contact)
                yes = valid & (contact > 0)
                no = valid & (contact == 0)
                self.nxip_hs_rows.append({
                    "system": system, "replicate": replicate,
                    "chain": chain,
                    "hs_chain": hs_col[len(hs_prefix):-len("_contact_state")],
                    "frames_contact": int(yes.sum()),
                    "frames_no_contact": int(no.sum()),
                    "occlusion_contact_mean": float(occlusion[yes].mean()) if yes.any() else np.nan,
                    "occlusion_no_contact_mean": float(occlusion[no].mean()) if no.any() else np.nan,
                })
        # A frame-aligned downstream outcome; lipid contacts do not enter
        # geometric state training or assignments.
        lipid_path = repout / "extensions/sdc4_mechanisms/lipid_region_specificity.parquet"
        if lipid_path.is_file():
            lipid_data = pd.read_parquet(lipid_path)
            if not lipid_data.empty:
                merged = lipid_data.merge(
                    pd.DataFrame({"time_ps": times, "state": labels}),
                    on="time_ps", how="inner", validate="many_to_one",
                )
                summary = lipid_data.loc[lipid_data.moiety.eq("any")].groupby(
                    ["region", "species", "leaflet"], sort=True
                )
                for (region, species, leaflet), item in summary:
                    self.lipid_replica_rows.append({
                        "system": system, "replicate": replicate,
                        "region": region, "species": species, "leaflet": leaflet,
                        "frames_available": int(item.contact_fraction.notna().sum()),
                        "mean_contact_fraction": float(item.contact_fraction.mean()),
                        "mean_composition_enrichment": float(item.composition_enrichment.mean()),
                        "mean_local_enrichment": float(item.local_composition_enrichment.mean()),
                        "mean_matched_local_contact_fraction": float(item.get(
                            "matched_local_contact_fraction", pd.Series(dtype=float)).mean()),
                        "mean_matched_local_enrichment": float(item.get(
                            "matched_local_composition_enrichment", pd.Series(dtype=float)).mean()),
                        "frames_matched_local_available": int(item.get(
                            "matched_local_contact_fraction", pd.Series(dtype=float)).notna().sum()),
                        "frames_matched_local_enrichment": int(item.get(
                            "matched_local_composition_enrichment", pd.Series(dtype=float)).notna().sum()),
                    })
                # Issue #101: contact-self.conditioned NxIP exposure, computed
                # within replica and matched exactly by time. Species absent
                # in a frame are excluded, not mislabelled noncontact.
                lip = lipid_data.loc[
                    lipid_data.moiety.eq("any")
                    & lipid_data.region.str.startswith("motif_NxIP_")
                ].copy()
                for (region, species, leaflet), group in lip.groupby(
                    ["region", "species", "leaflet"], sort=True
                ):
                    chain = "A" if "_PROA_" in region else "B" if "_PROB_" in region else None
                    if chain is None:
                        continue
                    endpoint = f"sdc4_motif_NxIP_PRO{chain}_1_steric_occlusion_fraction"
                    if endpoint not in df:
                        continue
                    outcome = df[["time_ps", endpoint]].merge(
                        group[["time_ps", "contacted_molecules", "available_molecules"]],
                        on="time_ps", how="inner", validate="one_to_one",
                    )
                    outcome = outcome.loc[
                        outcome.available_molecules.gt(0)
                        & outcome[endpoint].notna()
                    ]
                    contact = outcome.contacted_molecules.gt(0)
                    self.nxip_lipid_rows.append({
                        "system": system, "replicate": replicate,
                        "chain": chain, "species": species, "leaflet": leaflet,
                        "frames_contact": int(contact.sum()),
                        "frames_no_contact": int((~contact).sum()),
                        "occlusion_contact_mean": (
                            float(outcome.loc[contact, endpoint].mean())
                            if contact.any() else np.nan
                        ),
                        "occlusion_no_contact_mean": (
                            float(outcome.loc[~contact, endpoint].mean())
                            if (~contact).any() else np.nan
                        ),
                    })
                for keys, subgroup in merged.groupby(
                    ["region", "species", "leaflet", "moiety", "state"],
                    sort=True,
                ):
                    valid = subgroup[subgroup["available_molecules"] > 0]
                    self.lipid_state_rows.append({
                        "system": system, "replicate": replicate,
                        "region": keys[0], "species": keys[1],
                        "leaflet": keys[2], "moiety": keys[3],
                        "state": int(keys[4]),
                        "frames_in_state": int(len(subgroup)),
                        "frames_with_species": int(len(valid)),
                        "mean_contact_fraction": (
                            float(valid["contact_fraction"].mean())
                            if len(valid) else np.nan
                        ),
                    })
        hs_names = sorted({
            name[len("sdc4_hs_"):-len("_hydration_water_per_heavy_atom")]
            for name in df.columns if name.startswith("sdc4_hs_")
            and name.endswith("_hydration_water_per_heavy_atom")
        })
        for hs in hs_names:
            prefix = f"sdc4_hs_{hs}_"
            if prefix + "membrane_contact_state" not in df:
                continue
            for state in range(k):
                for contact in (0, 1):
                    mask = (labels == state) & (
                        df[prefix + "membrane_contact_state"].to_numpy(int) == contact
                    )
                    if not np.any(mask):
                        continue
                    row = {
                        "system": system, "replicate": replicate, "hs_chain": hs,
                        "state": state, "membrane_contact_state": contact,
                        "frames": int(mask.sum()),
                        "hydration_water_per_heavy_atom_mean": float(
                            df.loc[mask, prefix + "hydration_water_per_heavy_atom"].mean()
                        ),
                    }
                    for col in df.columns:
                        if col.startswith(prefix + "ion_") and col.endswith("_count"):
                            row[col[len(prefix):] + "_mean"] = float(df.loc[mask, col].mean())
                    self.conditioned.append(row)

    def write(self, out, rng, draws):
        nxip_df = pd.DataFrame(self.nxip_conditioned, columns=[
            "system", "replicate", "chain", "state", "endpoint",
            "frames_in_state", "frames_outside_state", "mean_in_state",
            "mean_outside_state", "within_replicate_difference",
        ])
        atomic_parquet(out / "state_conditioned_nxip.parquet", nxip_df)
        atomic_parquet(
            out / "nxip_state_exposure_threshold_sensitivity.parquet",
            pd.DataFrame(self.nxip_joint, columns=[
                "system", "replicate", "chain", "state", "threshold",
                "frames_in_state", "frames_exposed_in_state",
                "exposed_fraction_given_state",
            ]),
        )
        # A condition contrast is computed only from independent replicate-level
        # contrasts. With three trajectories per condition these intervals are
        # exploratory; no frame-based p-values or significance stars are emitted.
        comparisons = []
        if not nxip_df.empty:
            systems_for_comparison = list(dict.fromkeys(nxip_df["system"]))
            if len(systems_for_comparison) == 2:
                for (chain, state, endpoint), group in nxip_df.groupby(
                    ["chain", "state", "endpoint"], sort=True
                ):
                    arrays = [
                        group.loc[group["system"] == system,
                                  "within_replicate_difference"].dropna().to_numpy(float)
                        for system in systems_for_comparison
                    ]
                    if min(map(len, arrays)) < 2:
                        continue
                    difference = float(np.mean(arrays[1]) - np.mean(arrays[0]))
                    indices0 = rng.integers(0, len(arrays[0]), (draws, len(arrays[0])))
                    indices1 = rng.integers(0, len(arrays[1]), (draws, len(arrays[1])))
                    boot = arrays[1][indices1].mean(axis=1) - arrays[0][indices0].mean(axis=1)
                    lo, hi = np.quantile(boot, [.025, .975])
                    comparisons.append({
                        "chain": chain, "state": int(state), "endpoint": endpoint,
                        "system_0": systems_for_comparison[0],
                        "system_1": systems_for_comparison[1],
                        "n_replicates_0": len(arrays[0]),
                        "n_replicates_1": len(arrays[1]),
                        "difference_system1_minus_system0": difference,
                        "replicate_bootstrap_ci_low": float(lo),
                        "replicate_bootstrap_ci_high": float(hi),
                        "inference": "exploratory; no multiplicity-adjusted test",
                    })
        atomic_parquet(out / "state_conditioned_nxip_comparisons.parquet",
                       pd.DataFrame(comparisons, columns=[
                           "chain", "state", "endpoint", "system_0", "system_1",
                           "n_replicates_0", "n_replicates_1",
                           "difference_system1_minus_system0",
                           "replicate_bootstrap_ci_low", "replicate_bootstrap_ci_high",
                           "inference",
                       ]))
        atomic_parquet(
            out / "state_conditioned_lipid_specificity.parquet",
            pd.DataFrame(self.lipid_state_rows, columns=[
                "system", "replicate", "region", "species", "leaflet", "moiety",
                "state", "frames_in_state", "frames_with_species",
                "mean_contact_fraction",
            ]),
        )
        lipid_replica_df = pd.DataFrame(self.lipid_replica_rows, columns=[
            "system", "replicate", "region", "species", "leaflet",
            "frames_available", "mean_contact_fraction",
            "mean_composition_enrichment", "mean_local_enrichment",
            "mean_matched_local_contact_fraction", "mean_matched_local_enrichment",
            "frames_matched_local_available", "frames_matched_local_enrichment",
        ])
        atomic_parquet(out / "lipid_specificity_replicate_statistics.parquet",
                       lipid_replica_df)
        # Equally weight independent trajectories, never individual frames.
        aggregate = []
        for (system, region, species, leaflet), subset in lipid_replica_df.groupby(
            ["system", "region", "species", "leaflet"], sort=True
        ):
            for measure in ("mean_contact_fraction", "mean_composition_enrichment",
                            "mean_local_enrichment", "mean_matched_local_contact_fraction",
                            "mean_matched_local_enrichment"):
                vals = subset[measure].dropna().to_numpy(float)
                if not len(vals):
                    continue
                lo, hi = (np.nan, np.nan)
                if len(vals) >= 2:
                    sample = rng.integers(0, len(vals), size=(draws, len(vals)))
                    lo, hi = np.quantile(vals[sample].mean(axis=1), [.025, .975])
                aggregate.append({
                    "system": system, "region": region, "species": species,
                    "leaflet": leaflet, "measure": measure,
                    "n_replicates": len(vals),
                    "replicate_mean": float(vals.mean()),
                    "replicate_sd": float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
                    "replicate_bootstrap_ci_low": float(lo),
                    "replicate_bootstrap_ci_high": float(hi),
                })
        atomic_parquet(out / "lipid_specificity_system_statistics.parquet",
                       pd.DataFrame(aggregate, columns=[
                           "system", "region", "species", "leaflet", "measure",
                           "n_replicates", "replicate_mean", "replicate_sd",
                           "replicate_bootstrap_ci_low", "replicate_bootstrap_ci_high",
                       ]))
        atomic_parquet(
            out / "nxip_state_circular_shift_null.parquet",
            pd.DataFrame(self.nxip_shift_null_rows, columns=[
                "system", "replicate", "chain", "state",
                "observed_difference", "shift_null_ci_low",
                "shift_null_ci_high", "n_shifts",
            ]),
        )
        atomic_parquet(
            out / "nxip_membrane_contact_conditioning.parquet",
            pd.DataFrame(self.nxip_membrane_rows, columns=[
                "system", "replicate", "chain", "membrane_contact_state",
                "frames", "mean_nxip_occlusion",
                "nxip_exposed_fraction_threshold_0_5",
            ]),
        )
        atomic_parquet(
            out / "nxip_hs_contact_conditioning.parquet",
            pd.DataFrame(self.nxip_hs_rows, columns=[
                "system", "replicate", "chain", "hs_chain",
                "frames_contact", "frames_no_contact",
                "occlusion_contact_mean", "occlusion_no_contact_mean",
            ]),
        )
        atomic_parquet(
            out / "nxip_lipid_contact_conditioning.parquet",
            pd.DataFrame(self.nxip_lipid_rows, columns=[
                "system", "replicate", "chain", "species", "leaflet",
                "frames_contact", "frames_no_contact",
                "occlusion_contact_mean", "occlusion_no_contact_mean",
            ]),
        )
        atomic_parquet(out / "state_conditioned_hs_solvation.parquet", pd.DataFrame(self.conditioned))
