"""Reusable topology/selection preflight; no trajectory processing or outputs."""
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.exceptions import SelectionError

from ..config import analysis_runs
from ..core.alignment import resolved_selection_counts
from ..core.glycans import discover_covalent_glycans, glycan_definition, glycan_sugars
from ..core.sequence import sequence_matches
from ..errors import ConfigError
from ..observables.residue_membrane_distance import _without_hydrogens
from ..observables.sasa import import_freesasa, sasa_instances


def require_optional_backends(cfg):
    """Fail before topology/trajectory work if an enabled analysis lacks its backend.

    Environment check only: schema-only planning deliberately does not call it.
    """
    if analysis_runs(cfg, "sasa"):
        import_freesasa()


def _emit_order_audit(emit, label, audit):
    """Per-species/tail counts and exclusions from ``order_preflight``."""
    chains = audit["chains"]
    chain_text = (f"tail chains bond-verified for {chains['verified_residue_tails']}/"
                  f"{chains['residue_tails']} residue tails" if chains["bond_graph"]
                  else "no topology bonds: tail chains unverified")
    emit(f"    {label}: tail_preset={audit['tail_preset'] or 'none'}, "
         f"{audit['lipids']} lipids, C-H pairs from {audit['ch_pairs']}, {chain_text}")
    minimum = audit["min_carbons_per_tail"]
    for row in audit["species"]:
        text = (f"{row['mapped_residues']}/{row['residues']} lipids, positions "
                f"{row['positions']}, {row['ch_bonds']} C-H bonds" if row["mapped_residues"]
                else f"unavailable for all {row['residues']} lipids")
        excluded = [f"{row[key]} {reason}" for key, reason in (
            ("unmatched_residues", "no matching carbons"),
            ("below_minimum_residues", f"fewer than {minimum} carbons"),
            ("no_ch_residues", "no C-H pairs")) if row[key]]
        if row["mapped_residues"] and excluded:
            text += "; excluded " + ", ".join(excluded)
        elif excluded:
            text += " (" + ", ".join(excluded) + ")"
        emit(f"    {label} {row['resname']} {row['tail']}: {text}")
    for resname, count in sorted(audit["excluded_resnames"].items()):
        emit(f"    {label} {resname}: {count} residues excluded by exclude_resnames")
    for resname in audit["unused_overrides"]:
        emit(f"    {label}: lipid_overrides.{resname} matches no selected residue")
    if audit.get("multi_residue_molecules"):
        emit(f"    {label}: WARNING {audit['multi_residue_molecules']} bonded molecules "
             "contain several order residues; species and leaflet counts are per residue")
    if audit.get("leaflets"):
        emit(f"    {label} initial leaflets at topology coordinates: "
             f"upper={audit['leaflets']['upper']}, lower={audit['leaflets']['lower']}")


def validate_selections(cfg, base=Path("."), *, systems=None, emit=None):
    """Check topology against configured selections; mismatches raise ``ConfigError``.

    Preflight reads topologies only, so its ``ValueError``/``SelectionError``
    failures describe user inputs rather than numerical faults.
    """
    try:
        return _validate_selections(cfg, base, systems=systems, emit=emit)
    except ConfigError:  # already contextual; keep its original cause
        raise
    except SelectionError as exc:
        raise ConfigError(f"Invalid MDAnalysis selection: {exc}") from exc
    except ValueError as exc:
        raise ConfigError(str(exc)) from exc


def _validate_selections(cfg, base, *, systems, emit):
    emit = emit or (lambda *args, **kwargs: None)
    require_optional_backends(cfg)
    base = Path(base)
    resolved = {}
    for system, spec in cfg["systems"].items():
        if systems is not None and system not in systems:
            continue
        topology = Path(spec["topology"])
        topology = topology if topology.is_absolute() else (base / topology).resolve()
        try:
            universe = mda.Universe(str(topology))
        except (OSError, ValueError) as exc:
            raise ConfigError(
                f"System {system!r}: cannot read topology {topology}: {exc}. "
                "Verify the topology file and format before running analyses."
            ) from exc
        try:
            counts = resolved_selection_counts(universe, cfg)
        except (SelectionError, ValueError, AttributeError) as exc:
            raise ConfigError(
                f"System {system!r}: cannot resolve protein, feature, alignment "
                f"or named group selections in {topology}: {exc}. "
                "Check MDAnalysis selection syntax and topology attributes."
            ) from exc
        from ..core.motifs import resolve_motifs
        counts["motifs"] = resolve_motifs(universe, cfg.get("selections", {}).get("motifs", {}))
        for name, matches in counts["motifs"].items():
            emit(f"    motif {name}: {matches}")
        glycans = cfg.get("selections", {}).get("glycans", {})
        hydration = cfg["analyses"].get("hydration", {})
        consumed = set(hydration.get("glycans") or []) if hydration.get("enabled") else set()
        counts["glycans"] = {}
        for name in glycans:
            definition = glycan_definition(glycans, name)
            audit = discover_covalent_glycans(universe, definition)[2]
            if name in consumed:  # hydration also needs linear sugar positions
                sugars = glycan_sugars(universe, definition)[0]
                for chain, item in audit.items():
                    item["sugars"] = sugars.loc[sugars.chain.eq(chain), "sugar_type"].tolist()
            counts["glycans"][name] = audit
            for chain, item in audit.items():
                emit(f"    glycan {name} {chain}: {item['glycan_residues']} residues "
                     f"{item['glycan_resnames']} linked to {item['protein_chain']}:"
                     f"{item['protein_resname']}{item['protein_resid']} {item['protein_atom']}")
        if cfg["analyses"].get("solvation", {}).get("enabled"):
            from ..observables.solvation import _targets
            solvation = cfg["analyses"]["solvation"]
            hs_auto = solvation.get("hs_auto") or cfg["analyses"].get("sdc4_mechanisms", {}).get("hs_auto") or {}
            _targets(universe, {**solvation, "hs_auto": hs_auto},
                     motifs=cfg.get("selections", {}).get("motifs", {}))
        resolved[system] = counts
        emit(
            f"  {system}: protein={counts['protein']}, "
            f"alignment={counts['alignment']}, features={counts['features']}"
        )
        for role, consumers in counts["required_by"].items():
            emit(f"    {role}: required by {', '.join(consumers)}" if consumers
                 else f"    {role}: not required by enabled analyses")
        if ((counts["required_by"]["protein"] and counts["protein"] == 0)
                or (counts["required_by"]["features"] and counts["features"] == 0)):
            raise ValueError(
                f"System {system!r} has a zero-atom protein/feature selection "
                f"(protein={counts['protein']}, features={counts['features']}) "
                f"in {topology}. Check selections.protein and selections.features."
            )
        if (
            cfg["analyses"]["alignment"].get("enabled", True)
            and counts["alignment"] < 3
        ):
            raise ValueError(
                f"System {system!r}: alignment requires at least three atoms "
                f"(matched {counts['alignment']}) in {topology}. "
                "Check selections.alignment."
            )
        for group_id, group in counts["groups"].items():
            emit(
                f"    group {group_id}: group={group['group_atoms']}, "
                f"features={group['feature_atoms']}, "
                f"alignment={group['alignment_atoms']}"
            )
            if not group["group_atoms"] or (counts["required_by"]["features"] and not group["feature_atoms"]):
                raise ValueError(
                    f"System {system!r}, group {group_id!r} resolves to zero atoms "
                    f"(group={group['group_atoms']}, features={group['feature_atoms']}) "
                    f"in {topology}. Check selections.groups."
                )
            if (
                cfg["analyses"]["alignment"].get("enabled", True)
                and group["alignment_atoms"] < 3
            ):
                raise ValueError(
                    f"System {system!r}, group {group_id!r}: alignment requires at least three atoms "
                    f"(matched {group['alignment_atoms']}) in {topology}. "
                    "Check the group's alignment_selection."
                )
        if cfg["analyses"].get("sasa", {}).get("enabled", False):
            sasa_instances(universe, cfg, validate_only=True)
        for run in analysis_runs(cfg, "density_1d"):
            reference = run.spec.get("membrane_center_selection")
            if reference:
                atoms = universe.select_atoms(reference)
                if len(atoms) < 4 or len(atoms) != len(atoms.residues):
                    raise ValueError(f"{system}: membrane density needs at least four single-atom lipid anchors")
            for label, selection in (run.spec.get("selections") or {}).items():
                selection = selection if isinstance(selection, str) else selection["selection"]
                if not len(universe.select_atoms(selection)):
                    raise ValueError(f"{system}: density_1d[{run.run_id}].{label} selection matched no atoms")
        membrane = cfg["analyses"].get("membrane", {})
        if membrane.get("enabled") and not len(universe.select_atoms(membrane["headgroup_selection"])):
            raise ValueError(f"{system}: membrane headgroup_selection matched no atoms")
        if cfg["analyses"]["membrane_order"].get("enabled"):
            from ..observables.membrane_order import order_preflight
            try:
                counts["membrane_order"] = order_preflight(universe, cfg["analyses"]["membrane_order"])
            except ValueError as exc:
                raise ValueError(f"{system}: {exc}") from exc
            _emit_order_audit(emit, "membrane_order", counts["membrane_order"])
        mechanism = cfg.get("analyses", {}).get("sdc4_mechanisms", {})
        if mechanism.get("enabled", False):
            for kind in ("tm_chains", "ectodomains"):
                for label, selection in mechanism[kind].items():
                    atoms = universe.select_atoms(selection)
                    if len(atoms) < 4:
                        raise ValueError(
                            f"{system}: SDC4 {kind}.{label} matched only {len(atoms)} atoms"
                        )
                    emit(f"    SDC4 {kind}.{label}: {len(atoms)} atoms")
            ca = universe.select_atoms(
                mechanism.get("motif_selection") or
                cfg["selections"]["protein"] + " and name CA"
            )
            motifs = sequence_matches(ca, mechanism.get("motifs") or {"NxIP": "N.IP"})
            if mechanism.get("require_motifs", True) and not motifs:
                raise ValueError(f"{system}: no NxIP motif found in topology sequence")
            emit(f"    SDC4 NxIP motifs: {list(motifs)}")
            lipid_owners = set()
            for name, selection in mechanism["lipid_anchors"].items():
                atoms = universe.select_atoms(selection)
                if not len(atoms):
                    if mechanism.get("allow_absent_species", False):
                        emit(f"    SDC4 lipid {name}: absent (permitted)")
                        continue
                    raise ValueError(f"{system}: missing SDC4 lipid {name}")
                ids = set(map(int, atoms.resindices))
                if ids & lipid_owners:
                    raise ValueError(f"{system}: overlapping SDC4 lipid anchor selections")
                lipid_owners.update(ids)
                emit(f"    SDC4 lipid {name}: {len(ids)} molecules")
            local_order = mechanism.get("local_order") or {}
            if local_order.get("enabled"):
                from ..systems.sdc4.mechanisms import local_order_options
                from ..observables.membrane_order import order_preflight
                try:
                    audit = order_preflight(
                        universe, local_order_options(local_order, cfg["analyses"]["membrane_order"]),
                        leaflets=False, anchor_residues=lipid_owners)
                except ValueError as exc:
                    raise ValueError(f"{system}: SDC4 local_order: {exc}") from exc
                counts["sdc4_local_order"] = audit
                _emit_order_audit(emit, "SDC4 local_order", audit)
                if audit["lipids_without_anchor"] == audit["lipids"]:
                    raise ValueError(
                        f"{system}: no SDC4 local_order lipid is an SDC4 lipid anchor residue; "
                        "local -S_CH would always be unavailable")
                if audit["lipids_without_anchor"]:
                    emit(f"    SDC4 local_order: {audit['lipids_without_anchor']} order lipids "
                         "have no SDC4 lipid anchor and enter only the bulk mean")
            for name, selection in mechanism.get("ion_selections", {}).items():
                count = len(universe.select_atoms(selection))
                if not count:
                    raise ValueError(f"{system}: ion selection {name} matched zero atoms")
                emit(f"    SDC4 ion {name}: {count} atoms")
            if (mechanism.get("hs_auto") or {}).get("enabled"):
                hs, sites, audit = discover_covalent_glycans(
                    universe, mechanism["hs_auto"]
                )
                for name, item in audit.items():
                    emit(
                        f"    HS {name}: {item['glycan_residues']} glycan residues "
                        f"covalently linked to {item['protein_chain']}:"
                        f"{item['protein_resid']} {item['protein_atom']}"
                    )
            elif mechanism.get("hs_chains"):
                for name, selection in mechanism["hs_chains"].items():
                    if not len(universe.select_atoms(selection)):
                        raise ValueError(f"{system}: HS chain {name} selection empty")
            else:
                emit("    WARNING: SDC4 has no mapped HS chains; HS metrics will be absent")
        for kind in ("bilayer_normal","lipid_enrichment"):
            extension=cfg.get("analyses",{}).get(kind,{})
            if not extension.get("enabled",False):
                continue
            lipid_residues=set()
            matched=0
            for species,selection in extension.get("species_selections",{}).items():
                atoms=universe.select_atoms(str(selection))
                if not len(atoms):
                    if extension.get("allow_absent_species",False):
                        emit(f"    {kind}: {species} absent (permitted)")
                        continue
                    raise ValueError(f"{system}: {kind} species {species} matched no atoms")
                ids=set(map(int,atoms.resindices))
                if lipid_residues.intersection(ids):
                    raise ValueError(f"{system}: {kind} lipid anchor residue selections overlap")
                lipid_residues.update(ids)
                matched+=len(ids)
                emit(f"    {kind}: {species} {len(ids)} lipid anchors")
            if matched<4:
                raise ValueError(f"{system}: {kind} requires >=4 lipid anchors")
            if kind=="lipid_enrichment":
                protein=universe.select_atoms(str(
                    extension["footprint_selection"]))
                if not len(protein):
                    raise ValueError(f"{system}: lipid enrichment protein selection matched no atoms")
        hb=cfg.get("analyses",{}).get("lipid_hydrogen_bonds",{})
        if hb.get("enabled",False):
            lipid=universe.select_atoms(str(hb["lipid_selection"]))
            if not len(lipid):
                raise ValueError(f"{system}: lipid H-bond selection matched no atoms")
            emit(f"    lipid_hydrogen_bonds: {len(lipid)} lipid atoms")
        for run in analysis_runs(cfg, "residue_membrane_distance"):
            query = universe.select_atoms(str(run.spec.get("selection") or cfg["selections"]["protein"]))
            if not len(query):
                raise ValueError(f"{system}: residue-membrane {run.run_id} query selection matched no atoms")
            membrane_selection = run.spec.get("membrane_selection") or "resname " + " ".join(map(
                str,
                run.spec.get("lipid_resnames", [])
            ))
            membrane = _without_hydrogens(
                universe.select_atoms(str(membrane_selection)),
                bool(run.spec.get("exclude_hydrogens", True))
            )
            if not len(membrane):
                raise ValueError(
                    f"System {system!r}, residue-membrane run {run.run_id!r} "
                    "matched zero membrane atoms"
                )
            emit(
                f"    residue-membrane {run.run_id}: "
                f"{len(membrane)} membrane atoms"
            )
    return resolved
