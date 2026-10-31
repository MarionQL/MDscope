# MDscope analysis package

See [the current architecture](../docs/ARCHITECTURE.md) and
[repository rules](../docs/REPOSITORY_RULES.md).

Current canonical packages:
- `core/`: geometry, alignment, atom selection, trajectory sources and checkpoints
- `observables/`: RMSD, RMSF, MSD, radius of gyration, bilayer normals and lipid enrichment
- `reporting/`: PDF/ZIP report, scientific figure rendering and coupling-screen reporting
- `extended/`: numerical support for extended membrane observers

- `interactions/`: hydrogen bonds, salt bridges, ligand contacts and interaction tables
- `mechanisms/`: SDC4, HS topology, local environment and shared states

Remaining flat implementation modules will migrate to `statistics/` and `pipeline/` in later changes.
`cli.py`, `config.py` and `registry.py` remain top-level entry points.

Membrane observers consume physical coordinates before protein-only fitting.
For raw MSD, a separate `msd_trajectory` source can be synchronized with
the primary trajectory. Do not rotate or center the underlying shared frames.
Old Python module aliases are **not** supported.
