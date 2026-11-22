"""Shared sequence motif definitions and topology-only membership resolution."""
import re


def validate_motifs(definitions):
    if not isinstance(definitions, dict):
        raise ValueError("selections.motifs must be a mapping ({} means no motifs)")
    for name, spec in definitions.items():
        if not isinstance(name, str) or not name or "/" in name:
            raise ValueError("Motif names must be nonempty path components")
        if not isinstance(spec, dict) or ("sequence" in spec) == ("pattern" in spec):
            raise ValueError(f"motif {name}: specify exactly one sequence or pattern")
        sequence = spec.get("sequence")
        pattern = spec.get("pattern")
        if sequence is not None:
            if not isinstance(sequence, str) or not re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+", sequence):
                raise ValueError(f"motif {name}: sequence must use uppercase standard amino acids")
        elif not isinstance(pattern, str) or not pattern:
            raise ValueError(f"motif {name}: pattern must be nonempty")
        else:
            try:
                compiled = re.compile(pattern)
            except re.error as exc:
                raise ValueError(f"motif {name}: invalid pattern: {exc}") from exc
            if compiled.match("") is not None:
                raise ValueError(f"motif {name}: pattern must consume residues")
        selection = spec.get("selection", "protein")
        if not isinstance(selection, str) or not selection.strip():
            raise ValueError(f"motif {name}: selection must be nonempty")
        chains = spec.get("chains")
        if chains is not None and (not isinstance(chains, list)
                or any(not isinstance(chain, str) or not chain for chain in chains)
                or len(set(chains)) != len(chains)):
            raise ValueError(f"motif {name}: chains must be null or distinct chain labels; [] selects none")
        if type(spec.get("overlap", False)) is not bool:
            raise ValueError(f"motif {name}: overlap must be boolean")
        if spec.get("missing", "error") not in ("error", "allow"):
            raise ValueError(f"motif {name}: missing must be error or allow")


def resolve_motifs(universe, definitions):
    """Return ordered topology residue indices and identity audit for each match.

    Chain labels use chainID, falling back to segid exactly as the legacy matcher.
    Unknown residues, selection gaps and segment boundaries never bridge a hit.
    """
    from .sequence import sequence_chain_label, sequence_matches
    validate_motifs(definitions)
    resolved = {}
    for name, spec in definitions.items():
        ca = universe.select_atoms(spec.get("selection", "protein")).select_atoms("name CA")
        chains = spec.get("chains")
        if chains is not None:
            ca = ca[[(sequence_chain_label(atom) in chains) for atom in ca]]
        pattern = re.escape(spec["sequence"]) if "sequence" in spec else spec["pattern"]
        matches = sequence_matches(ca, {name: pattern}, allow_overlaps=spec.get("overlap", False))
        if not matches and spec.get("missing", "error") == "error":
            raise ValueError(f"motif {name}: no match found in topology")
        resolved[name] = {}
        for identity, indices in matches.items():
            residues = universe.residues[indices]
            atom = ca[ca.resindices == indices[0]][0]
            resolved[name][identity] = {
                "resindices": [int(i) for i in indices],
                "resids": [int(i) for i in residues.resids],
                "resnames": residues.resnames.tolist(),
                "segindex": int(atom.segindex), "segid": str(getattr(atom, "segid", "")),
                "chain": sequence_chain_label(atom),
            }
    return resolved
