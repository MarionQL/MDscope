# Package map

| Path | Responsibility |
|---|---|
| `core/` | Geometry, selections, alignment, coordinate sources, checkpoints |
| `observables/` | Structure, membrane, lipid, SASA, DSSP, MSD and density metrics |
| `interactions/` | H bonds, salt bridges, ligand/lipid contacts |
| `mechanisms/` | SDC4, HS topology, local environments, shared states |
| `statistics/` | Autocorrelation, replicate bootstrap, contrasts |
| `reporting/` | Saved findings and optional figures |
| `pipeline/` | Fused streaming, restart, orchestration and finalization |
| `cli.py`, `config.py`, `registry.py` | CLI, configuration and task registration |

**Coordinate rule:** Membrane observables use source coordinates **before** protein fitting; optional `msd_trajectory` is separately synchronized. Do not rotate or recenter the shared input frames. Old internal Python aliases are unsupported.

[Architecture](../../docs/developer/architecture.md) · [Scientific contracts](../../docs/developer/scientific-contracts.md) · [Contributing](../../docs/developer/contributing.md)
