# MDscope

MDscope analyzes replicated molecular dynamics simulations and compares conditions with replicate-aware bootstrap confidence intervals. It covers protein structure, membrane properties, molecular contacts, and solvation.

## Quick start

From the repository checkout, install MDscope with Python 3.11 or later:

```bash
python -m pip install -e .
mdscope init --output analysis.yaml
mdscope validate analysis.yaml
mdscope run analysis.yaml
```

Use `python -m pip install -e '.[sasa]'` to enable SASA. See the [baseline YAML](synd4_baseline.yaml) for an example and the [workflow guide](docs/workflows.md) for Slurm, restarting runs, and recomputing statistics from saved outputs.

Start with `project_summary/statistics/findings.md` inside the configured output directory. It summarizes bootstrap findings and comparison coverage. The underlying numerical tables are stored alongside it.

## Limitations

- Choose the production frame range explicitly. MDscope does not detect equilibration, and block-bootstrap intervals depend on adequate sampling and approximate stationarity.
- With three independent replicas per condition, an exact unpaired permutation test has only 20 possible assignments. Interpret effect sizes and intervals alongside individual replica results.
- Membrane measurements depend on their geometric assumptions. Membrane-relative density follows translation along a configured global normal and can still broaden with tilt or undulations. Lipid footprint enrichment requires an orthorhombic cell.
- Interaction results depend on topology, atom selections, and frame spacing. Completed-episode lifetimes exclude contacts that cross the observation boundaries.

See the [documentation index](docs/README.md) for configuration, measurement definitions, and validation details.
