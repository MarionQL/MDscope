# MDscope full SDC4 analysis + shareable evidence report

## Configuration

Use **`synd4_full_mechanisms.yaml`**. It preserves the **same two systems,
topology paths and all six replicate XTC paths** as
`synd4_baseline.yaml`. The output is separate:
`results_memb_full_mechanisms/` (relative to the YAML directory), so the
existing baseline outputs are not archived or overwritten.

The new YAML enables:

- Original RMSD, RMSF, Rg, membrane geometry, density, residue-membrane
  distances, DSSP, protein H-bonds and salt bridges
- Global lipid C-H order, local bilayer-normal deviations, leaflet-normalized
  lipid enrichment and intermolecular lipid-lipid hydrogen bonds
- Protein-centered GM1a/BGAL/cholesterol density maps
- Ectodomain height, tilt and membrane-contact states
- Sequence-derived NxIP occlusion and direct motif SASA (protein-only)
- **Six individual covalently attached HS chains** with orientation, direct
  glycan-lipid contacts, hydration and K/Cl coordination
- TM crossing, generic helix-helix interface and per-chain GXXXG contacts
- Local membrane thickness, quadratic curvature and C-H order
- One shared conformational-state model across all six replicates
- Hierarchical statistics, replicate-aware comparisons and figures
- A PDF evidence report and a compact review ZIP

### HS chains: actual TPR bond audit

The original two GROMACS TPRs were inspected using the bond topology,
not inferred from spatial proximity. In **both** systems, the following
covalent attachments were found:

| Protomer | Protein attachment residues |
|---|---|
| PROA | SER 39 OG, SER 61 OG, SER 63 OG |
| PROB | SER 309 OG, SER 331 OG, SER 333 OG |

Each protomer has three HS chains. **All three share the PROA or PROB
chainID with the protein**; selecting `chainID PROA` alone would incorrectly
merge three glycans and protein. `hs_auto` instead isolates connected
components of glycan-only topology bonds and requires each component to
have exactly one covalent attachment site. It checks that six chains were
found and fails rather than silently omitting a chain.

This audit applies to the two topologies referenced by the current YAML.
**Rebuilt mature-protein systems must be validated again**, because their
residue numbering and glycosylation may differ.

## Running

On Tinkercliffs or another host with the original topology/trajectory paths:

```bash
python -m pip install reportlab pillow freesasa
python -m mdscope.cli validate synd4_full_mechanisms.yaml
python -m mdscope.cli run synd4_full_mechanisms.yaml
```

Use the existing Slurm/array runner for individual replicates if preferred.
After all six replicates have completed, run:

```bash
python -m mdscope.cli finalize synd4_full_mechanisms.yaml
```

**Note:** `finalize` recomputes and archives the previous project-level
summary, whereas `report` below only regenerates the PDF/ZIP from
existing tables. The full `run` command already finalizes automatically.

## Files to send back to ChatGPT

Once all analyses and project-level finalization finish:

```text
results_memb_full_mechanisms/
  project_summary/
    report/
      MDscope_evidence_report.pdf
      MDscope_review_bundle.zip
```

Upload **`MDscope_evidence_report.pdf`** to request an interpretation
of the condition differences, state populations, HS/membrane mechanisms,
and their statistical support.

For deeper independent verification and new plots, upload
**`MDscope_review_bundle.zip`** too. The ZIP contains:

- Full-precision statistical tables in CSV format
- Shared-state model, centroids, assignments, populations and diagnostics
- All replicate HS selection/attachment audits
- Replicate lipid density/enrichment NPZ maps
- Conditioned HS hydration/ion summary tables
- Up to 250 evenly spaced frame-level rows per replicate (configurable)
- Resolved YAML, environment/run metadata and the PDF itself

It **does not** contain XTC/TPR trajectories, large caches or complete
per-frame data tables, so it should be far smaller than the original MD
outputs. The PDF can be reviewed on its own; the ZIP enables a more
rigorous re-analysis of the underlying statistics.

To regenerate **only** the report after finalization, without opening
trajectories or archiving the project summary:

```bash
python -m mdscope.cli report synd4_full_mechanisms.yaml
```

## Scientific and statistical limitations

The PDF is a reproducible descriptive evidence report, **not an LLM-written
results/discussion section**. It displays condition means, hierarchical
confidence intervals, permutation p-values and FDR q-values, replicate
consistency, shared-state centroids, population/transition diagnostics,
selection audits, lipid maps and methodology caveats.

With three independent replicates per condition, exact permutation
tests have limited p-value resolution. Shared-state populations and
first/second-half divergence help diagnose sampling, but do not establish
equilibrium or causal relationships.

NxIP occlusion is a steric-contact proxy; direct motif SASA is calculated
with the configured FreeSASA atom selection and does not account for
glycan/lipid occlusion unless those atoms are included. HS hydration and
ion contacts are distance-cutoff observables, not electrostatic binding
free energies. Lipid density maps are not corrected for protein-excluded
area. Curvature is a local quadratic height-field approximation.

The selected TM GXXXG motif was detected in PROA but **not** PROB in the
audited older topology. The analysis therefore also reports generic
helix-helix and each available motif-to-opposite-TM contacts instead of
mislabeling a nonexistent second GXXXG motif.

**Before a long production analysis:** run `validate` and verify that
the selected atoms, six HS attachments, NxIP motifs, K/Cl ion names, lipid
anchor counts and time window are as intended. This is especially
important if the topology was rebuilt.
