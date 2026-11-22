# Package map

| Path | Responsibility |
|---|---|
| `config/` | YAML loading, defaults, named runs and ordered per-analysis validators (`validators/`) |
| `core/` | Geometry, selections, alignment, coordinate sources, checkpoints |
| `observables/` | Structure, membrane, lipid, SASA, DSSP, MSD and density metrics |
| `interactions/` | H bonds, salt bridges, ligand/lipid contacts |
| `mechanisms/` | SDC4, HS topology, local environments, shared states |
| `statistics/` | Autocorrelation, replicate bootstrap, contrasts |
| `reporting/` | Saved findings and optional figures |
| `pipeline/` | Fused streaming, executable replica tasks, restart and finalization |
| `cli.py`, `planning.py`, `studio/` | User interfaces, planning and Studio |
| `contracts/` | Per-analysis output descriptors, static paths, companions and shared figure/state/statistics option schemas |
| `output_paths.py`, `registry.py` | Backward-compatible path exports, configuration-dependent companion selection and output planning |

**Coordinate rule:** Membrane observables use source coordinates **before** protein fitting; optional `msd_trajectory` is separately synchronized. Do not rotate or recenter the shared input frames. Old internal Python aliases are unsupported.

[Architecture](../../docs/developer/architecture.md) · [Scientific contracts](../../docs/developer/scientific-contracts.md) · [Contributing](../../docs/developer/contributing.md)
