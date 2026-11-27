"""Authoritative persisted-output contracts for MDscope analyses.

Planning and execution both consume these paths, companion artifacts and
schemas. Executable replica tasks live in pipeline.replica_tasks.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import analysis_enabled, analysis_runs
from .contracts.outputs import (
    FUSED_OUTPUTS,
    GLYCAN_ATTACHMENTS_FILE,
    HS_SUGAR_HYDRATION_PATH,
    OUTPUT_DESCRIPTORS,
    RESIDUE_OCCLUSION_FILE,
)
from .contracts.outputs import RUN_SCOPED_OUTPUTS as RUN_SCOPED_OUTPUTS
from .contracts.outputs import (
    artifact_schema as artifact_schema,
)


# Static schemas live in contracts.outputs. Missing column/shape metadata
# means the producer's output is dynamic or not yet catalogued, not empty.
@dataclass(frozen=True)
class FusedOutput:
    name: str
    path: Path
    run_id: str | None = None

    @property
    def checkpoint_key(self) -> str:
        """Unambiguous current-format identifier, independent of old markers."""
        return self.name if self.run_id is None else f"{self.name}-{self.run_id}"


def fused_output(repout: Path, name: str, *, run_id: str | None = None) -> FusedOutput:
    """Return the expected artifact for a named analysis instance.

    A run named ``default`` maps to the unscoped product when one exists.
    """
    descriptor = OUTPUT_DESCRIPTORS.get(name)
    if run_id == "default" and descriptor is not None and descriptor.primary_path:
        run_id = None
    if run_id is None:
        if descriptor is None or descriptor.primary_path is None:
            raise ValueError(f"{name} requires a run_id")
        return FusedOutput(name, Path(repout) / descriptor.primary_path)
    if descriptor is None or descriptor.named_run is None:
        raise ValueError(f"{name} does not yet support independently named runs")
    from .config import _validate_path_component
    run_id = _validate_path_component(str(run_id), f"{name} run ID")
    directory, filename = descriptor.named_run
    return FusedOutput(name, Path(repout) / directory / run_id / filename, run_id)


def fused_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Plan enabled fused outputs without creating readers or opening trajectories."""
    # Share production scheduling predicates without importing heavy producers.
    from .config.runs import (
        distance_enabled_names as distance_enabled,
    )
    from .config.runs import (
        extensions_enabled_names as extensions_enabled,
    )
    from .config.runs import (
        membrane_enabled_names as membrane_enabled,
    )
    from .config.runs import (
        sdc4_enabled,
    )
    entries = [fused_output(repout, name) for name in membrane_enabled(cfg)]
    extensions = extensions_enabled(cfg)
    entries.extend(fused_output(repout, name) for name in extensions
                   if name != "sasa")
    if "sasa" in extensions:
        entries.extend(fused_output(repout, "sasa", run_id=run.run_id)
                       for run in analysis_runs(cfg, "sasa"))
    entries.extend(fused_output(repout, "residue_membrane_distance", run_id=run_id)
                   for run_id in distance_enabled(cfg))
    if sdc4_enabled(cfg):
        entries.append(fused_output(repout, "sdc4_mechanisms"))
    return tuple(entries)



def analysis_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """The same persisted contracts apply to fused and standalone scheduling."""
    entries = list(fused_output_plan(cfg, repout))
    names = {entry.name for entry in entries}
    if cfg['analyses'].get('dssp', {}).get('enabled') and 'dssp' not in names:
        entries.append(fused_output(repout, 'dssp'))
    if cfg['analyses'].get('sasa', {}).get('enabled') and 'sasa' not in names:
        entries.extend(fused_output(repout, 'sasa', run_id=run.run_id)
                       for run in analysis_runs(cfg, 'sasa'))
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        entries.append(fused_output(repout, "potassium_bridge_audit"))
    return tuple(entries)


def saved_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Readiness contracts, including products written by the core stream.

    Paths come from the same contracts used by authoring and statistics; this
    does not add work to the executable fused schedule.
    """
    entries = list(analysis_output_plan(cfg, repout))
    names = {entry.name for entry in entries}
    for contract in scalar_contracts(cfg) + structural_profile_contracts(cfg):
        name = contract["analysis"]
        if name not in names:
            entries.append(FusedOutput(name, Path(repout) / contract["path_within_replica"]))
            names.add(name)
    return tuple(entries)


def core_saved_contracts(cfg, name):
    """Core stream companions and columns, from the public output contracts."""
    if name in FUSED_OUTPUTS:
        return []
    return [item for item in scalar_contracts(cfg) + structural_profile_contracts(cfg)
            if item["analysis"] == name]


def required_artifacts(entry: FusedOutput, cfg: dict) -> tuple[Path, ...]:
    """Required companions, without opening trajectories or inventing a manifest.

    Sparse event products that are legitimately absent after no events are not
    required. Large numerical tables are checked by the caller using metadata.
    """
    directory, name = entry.path.parent, entry.name
    core = core_saved_contracts(cfg, name)
    if core:
        # Core primaries are one directory below the replica root.
        return tuple(dict.fromkeys(entry.path.parent.parent / item["path_within_replica"]
                                   for item in core))
    descriptor = OUTPUT_DESCRIPTORS.get(name)
    companions = descriptor.companions if descriptor is not None else ()
    paths = [entry.path, *(directory / filename for filename in companions)]
    spec = cfg['analyses'].get(name, {})
    if name == 'sasa':
        spec = next(run.spec for run in analysis_runs(cfg, name)
                    if run.run_id == (entry.run_id or 'default'))
        if spec.get('atom_sasa_selection'):
            paths.append(directory/'atom_sasa.parquet')
        if spec.get('reference_run') is not None:
            paths.extend(directory/f for f in ('residue_burial.parquet', 'residue_burial_summary.parquet'))
    if name == 'sdc4_mechanisms' and spec.get('residue_occlusion', False):
        paths.append(directory / RESIDUE_OCCLUSION_FILE)
    if name == 'sdc4_mechanisms' and spec.get('interprotomer_contacts', False):
        paths.extend(directory/f for f in (
            'interprotomer_pair_frames.parquet', 'interprotomer_residue_frames.parquet',
            'interprotomer_pair_summary.parquet', 'interprotomer_residue_summary.parquet'))
    if name == 'sdc4_mechanisms' and spec.get('ceramide_scaffold_pooling', False):
        paths.append(directory/'ceramide_scaffold_pooled.parquet')
    if name == 'sdc4_mechanisms' and spec.get('motif_lipid_residence', False):
        paths.extend(directory/f for f in ('motif_lipid_contact_events.parquet',
            'motif_lipid_residence.parquet', 'motif_lipid_residence_intermittent.parquet',
            'lipid_region_specificity.parquet', 'lipid_residue_contacts.parquet'))
    if name == 'hydration' and spec.get('enabled') and spec.get('glycans'):
        # Optional for HS-only runs, so pre-existing SDC4 outputs stay reusable.
        paths.append(directory / GLYCAN_ATTACHMENTS_FILE)
    if name == 'hydration' and cfg['analyses'].get('sdc4_mechanisms', {}).get('hs_sugar_water_cutoff_A') is not None:
        paths.append(Path(entry.path).parents[2] / HS_SUGAR_HYDRATION_PATH)
    if name == 'solvation':
        options = spec.get('options') or {}
        bridges = options.get('bridges', spec.get('preset', 'standard') != 'standard')
        if options.get('events', spec.get('preset', 'standard') != 'standard') or bridges or options.get('water_hbond_bridges'):
            paths.extend(directory/f for f in ('episodes.parquet', 'survival.parquet'))
        if options.get('bridge_target_pairs'):
            paths.append(directory/'bridge_pair_summary.parquet')
        if options.get("bridge_spatial_ion") is not None:
            paths.extend(directory/name for name in ("bridge_spatial_reference.parquet",
                                                      "bridge_pair_timeseries.parquet"))
        if spec.get('aqueous_reference'):
            paths.append(directory/'aqueous_reference.parquet')
        if spec.get('grid_frames'):
            paths.append(directory/'water_density_3d.parquet')
    return tuple(paths)


STRUCTURAL_SCALAR_METRICS = {"rmsd": "rmsd_A",
                             "radius_of_gyration": "radius_of_gyration_A"}


def structural_scalar_contracts(cfg: dict, *, include_figures=True) -> list[dict]:
    """Configured scalar identities from the shared streaming summary writer.

    Group RMSD measures the intersection with selections.features. Rg uses
    its own selection intersection. The global scalar remains written even
    when a groups list is supplied; lists filter grouped rows only.
    """
    from .config import analysis_groups
    from .contracts.plot_validation import describe_plot, safe_plot_stem

    groups = analysis_groups(cfg)
    result = []
    for analysis, metric in STRUCTURAL_SCALAR_METRICS.items():
        spec = cfg["analyses"].get(analysis, {})
        if not spec.get("enabled", False):
            continue
        selected = spec.get("groups")
        ids = [None] + [key for key in groups if selected is None or key in selected]
        for group in ids:
            stem = safe_plot_stem(metric + (f"__group_{group}" if group else ""))
            measurement = (cfg["selections"]["features"] if analysis == "rmsd"
                           else spec.get("selection") or cfg["selections"]["protein"])
            fit = (groups[group].get("alignment_selection") if group else None)
            fitted = analysis == "rmsd" and cfg["analyses"]["alignment"].get("enabled", True)
            fit = ((fit or cfg["analyses"]["alignment"].get("selection")
                    or cfg["selections"]["alignment"]) if fitted else None)
            result.append({
                "analysis": analysis, "metric": metric, "group": group,
                "label": groups[group]["label"] if group else "Global",
                "columns": ["time_ps", *(["group"] if group else []), metric],
                "columns_status": "verified_subset_of_shared_table",
                "units": {"time_ps": "ps", metric: "angstrom"},
                "row_grain": "one analyzed frame per replicate" + (" per group" if group else ""),
                "path_within_replica": "summary/" + (
                    "group_timeseries.parquet" if group else "timeseries.parquet"),
                "measurement_selection": measurement,
                "group_selection": groups[group]["selection"] if group else None,
                "selection_semantics": "intersection with group_selection" if group else "direct",
                "alignment_selection": fit,
                "coordinate_source": "fitted" if fitted else "source",
                "statistics_selector": {"metric": metric, "group": group},
                "figures": [describe_plot(f"time_series:{stem}:{variant}", table_reports_pdf=cfg.get("plots", {}).get("table_reports_pdf", False), analysis=analysis)
                            for variant in ("timeseries", "distribution", "autocorrelation",
                                            "bootstrap_table", "comparison_table")] if include_figures else [],
            })
    return result


def structural_profile_contracts(cfg):
    """Streaming RMSF products; separate from frame-wise scalar tables."""
    from .config import analysis_groups
    if not cfg["analyses"]["rmsf"].get("enabled", True):
        return []
    base = ["atom_index", "topology_residue_index", "residue", "atom_name", "rmsf_A"]
    result = []
    for grouped in ([False, True] if analysis_groups(cfg) else [False]):
        result.append({
            "analysis": "rmsf", "run_id": "default", "role": "conditional_summary",
            "format": "Parquet", "columns_status": "verified_fixed",
            "columns": ["group", *base] if grouped else [*base, "equivalent_bfactor_A2"],
            "path_within_replica": "raw/group_rmsf.parquet" if grouped else "raw/rmsf.parquet",
            "row_grain": "one selected feature atom per replicate" + (" per group" if grouped else ""),
            "units": {"rmsf_A": "angstrom", **({} if grouped else {"equivalent_bfactor_A2": "angstrom^2"})},
            "coordinate_source": "fitted" if cfg["analyses"]["alignment"].get("enabled", True) else "source",
            "statistics_requirement": "Temporal RMSF inference additionally requires saved aligned-coordinate cache; a static RMSF table alone is insufficient.",
        })
    return result


def scalar_contracts(cfg, *, include_figures=True):
    """Deterministic pre-run scalar identities, including membrane and DSSP."""
    from .contracts.plot_validation import describe_plot
    result = structural_scalar_contracts(cfg, include_figures=include_figures)
    specs = []
    membrane = cfg["analyses"]["membrane"]
    if membrane.get("enabled"):
        if membrane.get("area_per_lipid", {}).get("enabled", True):
            specs.extend(("membrane", f"area_per_lipid{suffix}_A2", "angstrom^2")
                         for suffix in ("_upper", "_lower", ""))
        if membrane.get("thickness", {}).get("enabled", True):
            specs.append(("membrane", "membrane_thickness_A", "angstrom"))
    if cfg["analyses"]["dssp"].get("enabled"):
        specs.extend(("dssp", f"{kind}_fraction", "fraction")
                     for kind in ("helix", "sheet", "coil"))
    if cfg["analyses"].get("bilayer_normal", {}).get("enabled"):
        specs += [("bilayer_normal", "bilayer_normal_valid_fraction", "fraction"),
                  ("bilayer_normal", "bilayer_normal_leaflet_changes", "lipids")]
        for leaflet in ("upper", "lower"):
            specs.extend(("bilayer_normal", f"bilayer_normal_{leaflet}_{stat}_deg", "degrees")
                         for stat in ("mean", "sd"))
    if cfg["analyses"].get("lipid_contacts", {}).get("enabled"):
        specs.extend(("lipid_contacts", metric, "pairs" if "pair" in metric else "lipids") for metric in (
            "lipid_contact_pair_count", "lipid_contact_cross_leaflet_pair_count", "lipid_contact_leaflet_changes"))
    for analysis, metric, unit in specs:
        result.append({
            "analysis": analysis, "metric": metric, "group": None, "label": metric,
            "columns": ["time_ps", metric], "columns_status": "verified_subset_of_shared_table",
            "units": {"time_ps": "ps", metric: unit},
            "row_grain": "one analyzed frame per replicate", "path_within_replica": (FUSED_OUTPUTS[analysis] if analysis in {"lipid_contacts", "bilayer_normal"} else "summary/timeseries.parquet"),
            "measurement_selection": (membrane.get("headgroup_selection") if analysis == "membrane"
                                      else cfg["selections"]["protein"] if analysis == "dssp"
                                      else cfg["analyses"][analysis].get("selection") or cfg["analyses"][analysis].get("species_selections")),
            "group_selection": None, "selection_semantics": "producer-defined population; see scientific contracts",
            "alignment_selection": None, "coordinate_source": "source",
            "statistics_selector": {"metric": metric, "group": None},
            "figures": [describe_plot(f"time_series:{metric}:{variant}", table_reports_pdf=cfg.get("plots", {}).get("table_reports_pdf", False), analysis=analysis)
                        for variant in ("timeseries", "distribution", "autocorrelation", "bootstrap_table", "comparison_table")] if include_figures else [],
        })
    return result


def configured_profile_variables(cfg, *, include_optional=False):
    """Pre-run identities consumed by temporal_rmsf/structure_profiles readers."""
    from .config import analysis_groups
    variables = []
    def add(analysis, profile, unit, prerequisite, run_id="default"):
        variables.append(dict(analysis=analysis, profile=profile, run_id=run_id,
                              units=unit, prerequisite=prerequisite,
                              statistics_selector={"profile": profile}))
    if cfg["analyses"]["rmsf"].get("enabled", True):
        add("rmsf", "rmsf", "angstrom", "Saved aligned-coordinate cache and compatible atom identities required")
        for group in analysis_groups(cfg):
            add("rmsf", f"rmsf__group_{group}", "angstrom", "Group feature intersection and saved coordinate cache required")
    for run in analysis_runs(cfg, "sasa"):
        suffix = "" if run.run_id == "default" else f"__{run.run_id}"
        add("sasa", "sasa_mean" + suffix, "nm^2", "Per-frame residue SASA with compatible residue identities", run.run_id)
        if run.run_id != "default" and run.spec.get("reference_run"):
            add("sasa", "sasa_buried_nm2" + suffix, "nm^2", "Compatible isolated/reference-run burial tables required", run.run_id)
    if cfg["analyses"]["dssp"].get("enabled", False):
        for kind in ("helix", "sheet", "coil"):
            add("dssp", f"dssp_{kind}_occupancy", "fraction", "Per-frame residue assignments required")
    for run in analysis_runs(cfg, "density_1d"):
        for key, raw in run.spec.get("selections", {}).items():
            target = {"selection": raw} if isinstance(raw, str) else raw
            unit = "amu/A^3" if target.get("mode", "number") == "mass" else "1/A^3"
            add("density_1d", f"density_1d__{run.run_id}__{key}", unit,
                "Compatible physical bin grids and per-frame density cache required", run.run_id)
            if target.get("groups"):
                for group in analysis_groups(cfg):
                    add("density_1d", f"density_1d__{run.run_id}__{key}__{group}", unit,
                        "Conditional: target/group intersection must be nonempty; compatible grids required", run.run_id)
    sdc4 = cfg["analyses"].get("sdc4_mechanisms", {})
    if sdc4.get("enabled") and sdc4.get("interprotomer_contacts"):
        for metric, unit in (
            ("interprotomer_contact_occupancy", "fraction"),
            ("interprotomer_partner_count", "residues"),
            ("interprotomer_atom_contact_fraction", "fraction"),
        ):
            add("sdc4_mechanisms", metric, unit,
                "Saved frame-resolved full-protomer interface required")
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        names = [("mean_min_distance_A", "angstrom"), ("median_min_distance_A", "angstrom")]
        if include_optional or run.spec.get("contact_cutoff_A") is not None:
            names.append(("fraction_frames_within_cutoff", "fraction"))
        prefixes = [f"residue_membrane_distance__{run.run_id}"]
        if run.spec.get("groups"):
            prefixes.extend(f"residue_membrane_distance__{run.run_id}__group_{group}" for group in analysis_groups(cfg))
        for prefix in prefixes:
            for metric, unit in names:
                add("residue_membrane_distance", f"{prefix}__{metric}", unit,
                    "Compatible residue identities and per-frame matrix; grouped outputs require nonempty intersections", run.run_id)
    return variables


def density_output_contracts(cfg):
    """Mean density tables saved by DensityProfileAccumulator.finalize."""
    from .config import analysis_groups, profile_dir
    outputs = []
    for run in analysis_runs(cfg, "density_1d"):
        for key, raw in run.spec.get("selections", {}).items():
            target = {"selection": raw} if isinstance(raw, str) else raw
            unit = "amu/A^3" if target.get("mode", "number") == "mass" else "1/A^3"
            outputs.append(dict(analysis="density_1d", run_id=run.run_id, role="conditional_summary",
                format="Parquet", columns_status="verified_fixed", row_grain="one spatial bin per target per replicate",
                columns=[f"{run.spec.get('axis', 'z')}_A", "density", "analysis_instance", "selection_key",
                         "selection_label", "mode", "density_unit", "n_frames"],
                units={"density": unit}, path_within_replica=str(profile_dir(Path(), "density_1d", run.run_id) / f"{key}.parquet")))
            if target.get("groups"):
                for group in analysis_groups(cfg):
                    output = dict(outputs[-1])
                    output["path_within_replica"] = str(profile_dir(Path(), "density_1d", run.run_id) / f"{key}__{group}.parquet")
                    output["row_grain"] = "conditional on nonempty target/group intersection; one spatial bin per replicate"
                    outputs.append(output)
    return outputs


def state_feature_catalog(cfg):
    """Available static identities plus explicit saved-only discovery contracts.

    Topology-dependent molecule/bridge identities are templates until the common
    scalar reader discovers saved rows; they are never advertised as measured.
    """
    from .contracts.outputs import SAVED_SCALAR_SOURCE_FILES
    candidates = []
    for scalar in scalar_contracts(cfg):
        if scalar["path_within_replica"].startswith("summary/"):
            candidates.append({"source": "summary", "group": scalar["group"],
                "metric": scalar["metric"], "units": scalar["units"][scalar["metric"]],
                "analysis": scalar["analysis"], "availability": "conditional_on_saved_data"})
    from .contracts.statistics import GROUP_VARIABLE_UNITS
    solvation = cfg["analyses"].get("solvation", {})
    if solvation.get("enabled"):
        solvents = (["water"] if solvation.get("water_selection") else [])
        solvents += ["ions/" + name for name in (solvation.get("ion_selections") or {})]
        if solvation.get("ion_selection"):
            solvents.append("ions")
        for kind, targets in (solvation.get("targets") or {}).items():
            for label in targets:
                for solvent in solvents:
                    for metric in ("mean_n_solvent", "hydrated_fraction"):
                        units = GROUP_VARIABLE_UNITS[metric]
                        candidates.append({"source": "solvation",
                            "group": f"solvation/{solvent}/{kind}/{label}", "metric": metric,
                            "units": units, "analysis": "solvation",
                            "availability": "conditional_on_saved_data"})
    hydration = cfg["analyses"].get("hydration", {})
    for name in (hydration.get("glycans") or []) if hydration.get("enabled") else []:
        for metric in ("unique_water_count", "unique_water_per_site"):
            units = GROUP_VARIABLE_UNITS[metric]
            candidates.append({"source": "hydration", "group": f"hydration/{name}/total/all",
                               "metric": metric, "units": units, "analysis": "hydration",
                               "availability": "conditional_on_saved_data"})
    sources = {
        "summary": {"group": "null or shared group", "metrics": "numerical saved columns",
                    "prerequisite": "Saved frame-wise global/group summary with selected numerical columns"},
        "solvation": {"group": "solvation/<solvent>/<target_class>/<site_class>",
                      "metrics": ["mean_n_solvent", "hydrated_fraction"],
                      "prerequisite": "Saved per-frame solvent/site aggregates; requested motif targets must match topology"},
        "hydration": {"group": "hydration/<compartment>/<level>/<label>",
                      "metrics": ["unique_water_count", "unique_water_per_site"],
                      "prerequisite": "Saved aggregate hydration counts with topology-resolved group identities"},
        "solvation_bridges": {"group": "solvation_bridges/<JSON [solvent,target_a,target_b]>",
                      "metrics": ["bridge_entity_count", "bridge_present_fraction", "observed_minus_spatial_bridge_fraction"],
                      "prerequisite": "Saved bridge-pair traces; spatial contrast requires its optional reference"},
        "lipid_contacts": {"group": "lipid_contacts/<leaflet>/<species_a>/<species_b>",
                      "metrics": ["neighbors_per_molecule", "enrichment"],
                      "prerequisite": "Saved species/leaflet matrix; unavailable enrichment remains nonfinite"},
    }
    for name, source in sources.items():
        source["paths_within_replica"] = list(SAVED_SCALAR_SOURCE_FILES[name])
        source["identity_status"] = "template_until_saved_product_validation"
    return {"candidates": candidates, "sources": sources,
            "time_requirement": "Exact common time grid within each replicate; complete feature vectors only"}


def statistical_variable_owner(cfg, family, identity):
    """Resolve saved dynamic variables against the producer-owned catalog.

    Unknown identities deliberately remain unmapped; a display label is never
    evidence of analysis ownership. Configured run IDs are matched verbatim.
    """
    from .config.defaults import DEFAULTS, merge
    from .contracts.statistics import (
        DYNAMIC_VARIABLE_OWNERS,
        GROUP_VARIABLE_UNITS,
        INTERACTION_METRICS,
    )
    config = merge(DEFAULTS, cfg)
    variable = str(identity.get('profile' if family == 'profiles' else 'metric', ''))
    if family == 'interactions':
        kind = identity.get('interaction_type')
        return (kind, {'occupancy': 'fraction', 'event_rate_per_ns': '1/ns'}.get(
            variable, 'ps' if variable.endswith('_ps') else 'unspecified')) if kind in config['analyses'] and variable in INTERACTION_METRICS else ('unmapped', 'unspecified')
    descriptors = (configured_profile_variables(config, include_optional=True)
                   if family == 'profiles' else scalar_contracts(config, include_figures=False))
    key = 'profile' if family == 'profiles' else 'metric'
    for item in descriptors:
        if item[key] == variable:
            unit = item['units'] if family == 'profiles' else item['units'][variable]
            return item['analysis'], unit
    group = str(identity.get('group') or '')
    if family == 'time_series' and '/' in group:
        source = group.split('/', 1)[0]
        analysis = 'solvation' if source == 'solvation_bridges' else source
        if source in {'solvation', 'hydration', 'solvation_bridges', 'lipid_contacts'} and variable in GROUP_VARIABLE_UNITS:
            return analysis, GROUP_VARIABLE_UNITS.get(variable, 'unspecified')
    for owner_family, prefix, analysis, unit in DYNAMIC_VARIABLE_OWNERS:
        if family == owner_family and variable.startswith(prefix):
            return analysis, unit
    return 'unmapped', 'unspecified'


def figure_descriptor(cfg, identifier):
    """Join plot presentation with exact producer-owned variable identities."""
    from .config.defaults import DEFAULTS, merge
    from .contracts.plot_validation import describe_plot, safe_plot_stem
    config = merge(DEFAULTS, cfg)
    from .config import normalize_analysis_booleans
    normalize_analysis_booleans(config)
    parts = identifier.split(":")
    # Ownership does not change when the user disables a producer. Eligibility
    # is checked separately by planning and execution.
    ownership_config = config
    for spec in ownership_config['analyses'].values():
        if isinstance(spec, dict):
            spec['enabled'] = True
    owner = None
    if parts[0] in {"time_series", "profiles"} and len(parts) == 3:
        profiles = parts[0] == "profiles"
        descriptors = (configured_profile_variables(ownership_config, include_optional=True)
                       if profiles else scalar_contracts(ownership_config, include_figures=False))
        matches = []
        for item in descriptors:
            variable = item["profile" if profiles else "metric"]
            group = item.get("group")
            stem = safe_plot_stem(variable + (f"__group_{group}" if group else ""))
            if stem == parts[1]:
                matches.append(item["analysis"])
        if len(matches) > 1:
            raise ValueError(f"Ambiguous plot ownership for {identifier}")
        if matches:
            owner = matches[0]
        else:
            # The existing statistics owner handles explicitly contracted dynamic
            # prefixes; unknown identities stay visible as unmapped.
            owner, _ = statistical_variable_owner(ownership_config, parts[0], {
                "profile" if profiles else "metric": parts[1]})
    descriptor = describe_plot(identifier, analysis=owner,
                               table_reports_pdf=config["plots"].get("table_reports_pdf", False))
    if parts[0] == 'overview' and parts[-1] == 'distributions':
        sources = sorted({"{system}/{replicate}/" + item['path_within_replica']
                          for item in scalar_contracts(config, include_figures=False)
                          if item['analysis'] == parts[1]})
        if parts[1] in FUSED_OUTPUTS:
            sources.append('{system}/{replicate}/' + FUSED_OUTPUTS[parts[1]])
        descriptor['source_candidates'] = descriptor['prerequisites'] = sorted(set(sources))
        descriptor['units_note'] = 'Metric-specific raw units; one frame distribution and mean per replica'
    return descriptor
