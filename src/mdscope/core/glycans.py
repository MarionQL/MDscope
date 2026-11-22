"""Covalent glycan-chain discovery from topology bonds, never proximity.

Glycan residues may share chainIDs with their protein (CHARMM-GUI syndecan
systems do), so chainID cannot separate chains. Each bonded glycan-only
component must bond to exactly one allowed protein attachment atom; any other
bond leaving the glycan selection is an error. Residue names and bonds come
from the force-field topology (TPR/PSF): bare PDB files cannot define chains.
"""
from __future__ import annotations

import re

import numpy as np

GLYCAN_DEFAULTS = {"glycan_selection": None, "attachment_selection": None,
                   "min_glycan_residues": 1, "expected_chains": None}
# Hydration compartments; "hs" is the SDC4 hs_auto preset.
RESERVED_GLYCAN_NAMES = ("protein", "membrane", "hs")


def _positive_int(value, minimum):
    return type(value) is int and value >= minimum


def validate_glycan_spec(spec, label, *, min_residues=1):
    """Validate one discovery definition; ``label`` names it in errors."""
    if not isinstance(spec, dict):
        raise ValueError(f"{label} must be a mapping")
    for key in ("glycan_selection", "attachment_selection"):
        if not isinstance(spec.get(key), str) or not spec[key].strip():
            raise ValueError(f"{label} requires a nonempty {key}")
    if not _positive_int(spec.get("min_glycan_residues", min_residues), min_residues):
        raise ValueError(f"{label}.min_glycan_residues must be an integer >= {min_residues}")
    expected = spec.get("expected_chains")
    if expected is not None and not _positive_int(expected, 1):
        raise ValueError(f"{label}.expected_chains must be null or a positive integer")


def validate_glycans(definitions):
    if not isinstance(definitions, dict):
        raise ValueError("selections.glycans must be a mapping ({} means no glycans)")
    for name, spec in definitions.items():
        if (not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", name)
                or name in RESERVED_GLYCAN_NAMES):
            raise ValueError(
                f"glycan {name!r}: names must be letters, digits, '_', '.' or '-' "
                f"and not {', '.join(RESERVED_GLYCAN_NAMES)}")
        validate_glycan_spec(spec, f"selections.glycans.{name}")


def glycan_definition(definitions, name):
    """Effective definition with generic defaults (no hardcoded chain count)."""
    if name not in definitions:
        raise ValueError(f"Unknown glycan definition {name!r}")
    return {**GLYCAN_DEFAULTS, **definitions[name]}


def _attachment_label(atom, chain):
    from .sequence import AA3
    code = AA3.get(str(atom.resname).upper(), str(atom.resname))
    return re.sub(r"[^A-Za-z0-9_]+", "_", f"{chain}_{code}{int(atom.resid)}")


def discover_covalent_glycans(u, spec):
    """Return (glycan AtomGroups, attachment AtomGroups, audit) keyed by attachment.

    Labels are ``<chain>_<one-letter residue><resid>`` (``PROA_S39``). Raises on
    missing bonds, empty selections, unattached or undersized components,
    bonds to atoms outside both selections, ambiguous attachments and an
    unexpected chain count.
    """
    import MDAnalysis as mda
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    from .sequence import sequence_chain_label

    glysel, att_sel = spec.get("glycan_selection"), spec.get("attachment_selection")
    if not glysel or not att_sel:
        raise ValueError("Glycan discovery requires glycan_selection and attachment_selection")
    gly = u.select_atoms(str(glysel))
    attachment = u.select_atoms(str(att_sel))
    if not len(gly) or not len(attachment):
        raise ValueError("Glycan discovery selections matched no glycan or attachment atoms")
    if np.intersect1d(gly.indices, attachment.indices).size:
        raise ValueError("Glycan attachment_selection must not overlap glycan_selection")
    try:
        bonds = np.asarray(u.bonds.indices, np.int32)
    except (AttributeError, mda.exceptions.NoDataError) as exc:
        raise ValueError("Glycan discovery requires covalent topology bonds (TPR/PSF)") from exc
    if bonds.ndim != 2 or bonds.shape[1] != 2:
        raise ValueError("Invalid bond topology for glycan discovery")
    gly_local = np.full(len(u.atoms), -1, np.int32)
    gly_local[gly.indices] = np.arange(len(gly), dtype=np.int32)
    ga, gb = gly_local[bonds[:, 0]], gly_local[bonds[:, 1]]
    internal = (ga >= 0) & (gb >= 0)
    graph = coo_matrix((
        np.ones(int(internal.sum()) * 2, np.int8),
        (np.r_[ga[internal], gb[internal]], np.r_[gb[internal], ga[internal]]),
    ), shape=(len(gly), len(gly))).tocsr()
    n_components, membership = connected_components(graph, directed=False)
    allowed = np.zeros(len(u.atoms), bool)
    allowed[attachment.indices] = True
    # Each external bond has exactly one glycan endpoint.
    crossing = (ga >= 0) != (gb >= 0)
    gly_end = np.where(ga >= 0, bonds[:, 0], bonds[:, 1])[crossing]
    other_end = np.where(ga >= 0, bonds[:, 1], bonds[:, 0])[crossing]
    stray = ~allowed[other_end]
    if stray.any():
        atoms = [f"{a.resname}{a.resid}:{a.name}" for a in u.atoms[np.unique(other_end[stray])[:5]]]
        raise ValueError(
            "Glycan atoms bond to atoms outside glycan_selection and attachment_selection "
            f"({', '.join(atoms)}); include every residue of each chain and only "
            "validated attachment chemistry")
    owner = membership[gly_local[gly_end]]
    out_gly, out_sites, audit = {}, {}, {}
    min_res = int(spec.get("min_glycan_residues", GLYCAN_DEFAULTS["min_glycan_residues"]))
    for component in range(n_components):
        ag = gly[membership == component]
        if len(ag.residues) < min_res:
            # Dropping a small fragment could misrepresent the chain count.
            raise ValueError(
                f"Glycan component {component} has {len(ag.residues)} residues "
                f"(<{min_res}); verify glycan_selection")
        sites = u.atoms[np.unique(other_end[owner == component])]
        if not len(sites):
            raise ValueError(
                f"Glycan component {component} has no covalent protein attachment; "
                "check glycan/attachment selections and topology bonds")
        if len(sites.residues) != 1:
            raise ValueError(f"Glycan component {component} attaches to multiple protein residues")
        if len(sites) != 1:
            raise ValueError(
                f"Glycan component {component} bonds to multiple attachment atoms "
                f"{sites.names.tolist()}; the linkage is ambiguous")
        site = sites[0]
        chain = sequence_chain_label(site)
        name = _attachment_label(site, chain)
        if name in out_gly:
            raise ValueError(f"Multiple glycan components map to attachment {name}")
        out_gly[name] = ag
        out_sites[name] = sites
        audit[name] = {
            "protein_chain": chain,
            "protein_resid": int(site.resid),
            "protein_resname": str(site.resname),
            "protein_atom": str(site.name),
            "glycan_atoms": len(ag),
            "glycan_residues": len(ag.residues),
            "glycan_resnames": sorted(set(ag.resnames.tolist())),
        }
    expected = spec.get("expected_chains")
    if expected is not None and len(out_gly) != int(expected):
        raise ValueError(
            f"Expected {expected} covalently attached glycan chains; "
            f"found {len(out_gly)}: {list(out_gly)}")
    if not out_gly:
        raise ValueError("No covalently attached glycan chains found")
    return out_gly, out_sites, audit


def glycan_sugars(u, spec):
    """Number each linear chain's sugars from the attachment sugar outward.

    Positions follow inter-residue bonds, independent of residue numbering.
    Branched chains fail: they have no unique attachment-relative position.
    Hydration sites are every atom named ``O*`` within each sugar residue.
    Returns (sugar table, oxygen AtomGroup, per-oxygen sugar row, attachments).
    """
    import pandas as pd
    chains, attachments, _ = discover_covalent_glycans(u, spec)
    resindices = u.atoms.resindices
    left, right = u.bonds.indices.T
    edges = resindices[u.bonds.indices]
    edges = edges[edges[:, 0] != edges[:, 1]]
    rows, atoms, owners = [], [], []
    for chain, group in sorted(chains.items()):
        residues = set(group.resindices)
        attachment = attachments[chain].indices
        roots = set(resindices[right[np.isin(left, attachment)]]) | set(
            resindices[left[np.isin(right, attachment)]])
        roots &= residues
        if len(roots) != 1:
            raise ValueError(f"{chain}: expected one attachment sugar")
        neighbors = {i: set() for i in residues}
        for a, b in edges[np.isin(edges, list(residues)).all(axis=1)]:
            neighbors[a].add(b)
            neighbors[b].add(a)
        current = roots.pop()
        visited = set()
        while current not in visited:
            visited.add(current)
            residue = u.residues[current]
            sugar = group[group.resindices == current]
            oxygen = sugar[np.char.startswith(sugar.names.astype(str), "O")]
            if not len(oxygen):
                raise ValueError(f"{chain}: no oxygen atoms in {residue.resname} {residue.resid}")
            owners.extend([len(rows)] * len(oxygen))
            atoms.extend(oxygen.indices)
            rows.append((chain, len(visited), residue.resname, int(residue.resid),
                         int(current), len(oxygen)))
            onward = neighbors[current] - visited
            if len(neighbors[current]) > (1 if len(visited) == 1 else 2):
                raise ValueError(
                    f"{chain}: glycan branches at {residue.resname} {residue.resid}; "
                    "only linear chains have attachment-relative sugar positions")
            if not onward:
                break
            current = onward.pop()
        if visited != residues:
            raise ValueError(f"{chain}: incomplete glycan traversal")
    columns = ["chain", "sugar_position", "sugar_type", "sugar_resid",
               "topology_residue_index", "oxygen_count"]
    return (pd.DataFrame(rows, columns=columns), u.atoms[atoms],
            np.asarray(owners), attachments)
