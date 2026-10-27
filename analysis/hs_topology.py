"""Discover covalently attached HS chains using topology bonds, never proximity.

For CHARMM-GUI syndecan systems the sugar residues may share PROA/PROB
chainIDs with the protein. ChainID alone is therefore insufficient to
distinguish the three separate HS chains on each protomer. We isolate the
glycan-only bonded components and map each to its unique covalent protein
attachment atom (e.g. Ser OG). This works for any residue numbering.
"""
from __future__ import annotations

import re
import MDAnalysis as mda
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def _label(atom):
    chain = str(atom.chainID) if hasattr(atom, "chainID") else "protein"
    label = f"{chain}_S{int(atom.resid)}"
    return re.sub(r"[^A-Za-z0-9_]+", "_", label)


def discover_covalent_hs(u, spec):
    """Return (glycan AtomGroups, attachment AtomGroups, audit).

    Requires all glycan atoms, the protein covalent attachment atoms, and
    TPR/PSF bond connectivity. Raises on ambiguous/duplicate attachments,
    unmatched components, or unexpected chain counts.
    """
    glysel = spec.get("glycan_selection")
    att_sel = spec.get("attachment_selection")
    if not glysel or not att_sel:
        raise ValueError("hs_auto requires glycan_selection and attachment_selection")
    gly = u.select_atoms(str(glysel))
    attachment = u.select_atoms(str(att_sel))
    if not len(gly) or not len(attachment):
        raise ValueError("HS auto-discovery selections matched no glycan or attachment atoms")
    try:
        bonds = np.asarray(u.bonds.indices, np.int32)
    except (AttributeError, mda.exceptions.NoDataError) as exc:
        raise ValueError("HS auto-discovery requires covalent topology bonds (TPR/PSF)") from exc
    if bonds.ndim != 2 or bonds.shape[1] != 2:
        raise ValueError("Invalid bond topology for HS discovery")
    n = len(u.atoms)
    gly_local = np.full(n, -1, np.int32)
    gly_local[gly.indices] = np.arange(len(gly), dtype=np.int32)
    ga, gb = gly_local[bonds[:, 0]], gly_local[bonds[:, 1]]
    internal = (ga >= 0) & (gb >= 0)
    graph = coo_matrix((
        np.ones(int(internal.sum()) * 2, np.int8),
        (np.r_[ga[internal], gb[internal]],
         np.r_[gb[internal], ga[internal]]),
    ), shape=(len(gly), len(gly))).tocsr()
    n_components, membership = connected_components(graph, directed=False)
    allowed_attachment = np.zeros(n, bool)
    allowed_attachment[attachment.indices] = True
    # Exactly one endpoint in glycan, other in allowed protein attachment.
    ab = bonds[(ga >= 0) & allowed_attachment[bonds[:, 1]]]
    ba = bonds[(gb >= 0) & allowed_attachment[bonds[:, 0]]]
    gly_indices = np.r_[ab[:, 0], ba[:, 1]]
    att_indices = np.r_[ab[:, 1], ba[:, 0]]
    out_gly, out_sites, audit = {}, {}, {}
    min_res = int(spec.get("min_glycan_residues", 3))
    if min_res < 2:
        raise ValueError("hs_auto.min_glycan_residues must be >= 2")
    for component in range(n_components):
        ids = np.flatnonzero(membership == component)
        ag = gly[ids]
        if len(ag.residues) < min_res:
            # Reject orphan small glycan fragments: silently dropping them
            # could misrepresent the number of attached chains.
            raise ValueError(
                f"Unassigned glycan component {component} has "
                f"{len(ag.residues)} residues (<{min_res}); verify glycan_selection"
            )
        matched = membership[gly_local[gly_indices]] == component
        bonded_att = np.unique(att_indices[matched])
        if not len(bonded_att):
            raise ValueError(
                f"Glycan component {component} has no covalent protein attachment; "
                "check glycan/attachment selections and topology bonds"
            )
        protein_sites = u.atoms[bonded_att]
        residue_ids = np.unique(protein_sites.resindices)
        if len(residue_ids) != 1:
            raise ValueError(
                f"Glycan component {component} attaches to multiple protein residues"
            )
        site = protein_sites[0]
        name = _label(site)
        if name in out_gly:
            raise ValueError(f"Multiple glycan components map to attachment {name}")
        out_gly[name] = ag
        out_sites[name] = protein_sites
        audit[name] = {
            "protein_chain": str(site.chainID),
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
            f"Expected {expected} covalently attached HS chains; "
            f"found {len(out_gly)}: {list(out_gly)}"
        )
    if not out_gly:
        raise ValueError("No covalently attached HS chains found")
    return out_gly, out_sites, audit
