"""Reusable topology/selection preflight; no trajectory processing or outputs."""
from pathlib import Path
import MDAnalysis as mda
from ..config import analysis_runs
from ..observables.residue_membrane_distance import _without_hydrogens
from ..mechanisms.hs_topology import discover_covalent_hs
from ..mechanisms.mechanistic_common import sequence_matches
from .alignment import resolved_selection_counts
from ..observables.sasa import sasa_instances


def validate_selections(cfg, base=Path("."), *, systems=None, emit=None):
    emit = emit or (lambda *args, **kwargs: None)
    base = Path(base)
    for system, spec in cfg["systems"].items():
        if systems is not None and system not in systems:
            continue
        topology = Path(spec["topology"])
        topology = topology if topology.is_absolute() else (base / topology).resolve()
        universe = mda.Universe(str(topology))
        counts = resolved_selection_counts(universe, cfg)
        emit(
            f"  {system}: protein={counts['protein']}, "
            f"alignment={counts['alignment']}, features={counts['features']}"
        )
        if counts["protein"] == 0 or counts["features"] == 0:
            raise ValueError(
                f"System {system!r} has a zero-atom protein/feature selection"
            )
        if (
            cfg["analyses"]["alignment"].get("enabled", True)
            and counts["alignment"] < 3
        ):
            raise ValueError(
                f"System {system!r}: alignment requires at least three atoms"
            )
        for group_id, group in counts["groups"].items():
            emit(
                f"    group {group_id}: group={group['group_atoms']}, "
                f"features={group['feature_atoms']}, "
                f"alignment={group['alignment_atoms']}"
            )
            if not group["group_atoms"] or not group["feature_atoms"]:
                raise ValueError(
                    f"System {system!r}, group {group_id!r} resolves to zero atoms"
                )
            if (
                cfg["analyses"]["alignment"].get("enabled", True)
                and group["alignment_atoms"] < 3
            ):
                raise ValueError(
                    f"System {system!r}, group {group_id!r}: alignment requires at least three atoms"
                )
        if cfg["analyses"].get("sasa", {}).get("enabled", False):
            sasa_instances(universe, cfg, validate_only=True)
        for run in analysis_runs(cfg, "density_1d"):
            for label, selection in (run.spec.get("selections") or {}).items():
                selection = selection if isinstance(selection, str) else selection["selection"]
                if not len(universe.select_atoms(selection)):
                    raise ValueError(f"{system}: density_1d[{run.run_id}].{label} selection matched no atoms")
        membrane = cfg["analyses"].get("membrane", {})
        if membrane.get("enabled") and not len(universe.select_atoms(membrane["headgroup_selection"])):
            raise ValueError(f"{system}: membrane headgroup_selection matched no atoms")
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
            for name, selection in mechanism.get("ion_selections", {}).items():
                count = len(universe.select_atoms(selection))
                if not count:
                    raise ValueError(f"{system}: ion selection {name} matched zero atoms")
                emit(f"    SDC4 ion {name}: {count} atoms")
            if (mechanism.get("hs_auto") or {}).get("enabled"):
                hs, sites, audit = discover_covalent_hs(
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
                    extension.get("protein_selection") or cfg["selections"]["protein"]))
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

