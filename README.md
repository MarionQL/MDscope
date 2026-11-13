# MDscope

Replicate-aware molecular dynamics (MD) analysis for protein structure, membranes, contacts, solvation, and condition comparisons.

## Start

Requires Python **3.11+**.

```bash
python -m pip install -e .
# Optional FreeSASA support:
python -m pip install -e '.[sasa]'

mdscope init --output analysis.yaml
mdscope validate analysis.yaml
mdscope run analysis.yaml
```

- Example configuration: [`synd4_baseline.yaml`](synd4_baseline.yaml); full SDC4 setup: [`synd4_full_mechanisms.yaml`](synd4_full_mechanisms.yaml).
- **First result to read:** `<output>/project_summary/statistics/findings.md`; inspect the adjacent Parquet tables and individual replicates before drawing conclusions.
- Restart, Slurm, offline statistics: [Workflows](docs/workflows.md).

## Interpret cautiously

| Risk | Required check |
|---|---|
| Equilibration and sampling | Select production frames explicitly; check stationarity and replica agreement |
| Statistical power | Three replicas per condition give only 20 exact unpaired assignments; prioritize effect size, CIs and individual replicas |
| Membrane geometry | Check leaflet/normal assumptions; membrane-relative density does not remove tilt or undulation |
| Molecular identities | Validate TPR bonds, lipid/HS selections and frame spacing |
| Contact lifetimes | Boundary-censored episodes and sparse sampling cannot establish continuous residence |

[Documentation index](docs/README.md) · [Scientific analyses](docs/analyses.md) · [Validation limits](docs/developer/scientific-audit.md) · [MIT license](LICENSE)
