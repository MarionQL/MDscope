# Execution, HPC workflows, and reporting

This guide consolidates the authoritative technical details formerly scattered across topic and implementation-history files. Sections retain their scientific caveats and validation limitations; consult the source code for the executable contract.


---

### MDscope trajectory pass consolidation

MDscope's regular protein/membrane streaming pass now also computes the three
lipid observables that previously opened the same XTC separately:

- **Lipid C–H order** (`membrane_order`), including original tail/site
  mapping, per-leaflet order and temporal profile NPZ outputs.
- **Local bilayer normals** (`bilayer_normal`), including optional per-lipid
  angles and dynamic leaflet assignment.
- **Protein-centered lipid enrichment** (`lipid_enrichment`), including raw
  counts, selection audits and PBC-aware protein centers.

The observer runs **before any coordinate alignment**, uses the original
topology-based selections, and shares anchor coordinates and leaflet
classification when the normal/enrichment selection specifications match.
Outputs retain their existing paths and column schemas; downstream statistics
and plots do not need to change.

### Number of XTC traversals

The main raw-coordinate frame loop now computes all of these analyses
**without reopening or re-iterating the XTC**:

- RMSD, RMSF, Rg, density, bilayer thickness/APL/MSD, salt bridges and ligand contacts
- Lipid C–H order, bilayer normals and protein-centered lipid enrichment
- Protein–protein/protein–ligand H-bonds and intermolecular lipid H-bonds
- Serial DSSP secondary-structure assignments and FreeSASA residue SASA
- Local-environment protein/lipid/HS/ion contacts and motif distances
- All named residue-to-membrane minimum-distance analyses

For the two published test systems, the initial baseline required **five**
passes, PR #11 reduced that to **two**, and this revision reduces production
to **one** (main streaming, three lipid observables, lipid H-bonds and
residue-to-membrane distances). The independent reference readers in CI are
deliberate and are *not* production passes.

For configurations with no regular protein streaming, all optional observers
share a single observer-only fallback traversal. A non-serial DSSP backend
remains an intentional separate analysis, because its multiprocessing
execution is controlled by MDAnalysis. The full production profile also runs multi-origin MSD and the potassium
bridge audit as separately restartable passes. Extra reads may occur when
rebuilding missing/corrupted outputs.

### Hydrogen bonds and DSSP

HydrogenBondAnalysis and DSSP are not reimplemented with new physical
approximations. The shared frame loop calls the original MDAnalysis
`_single_frame()` detector on each raw frame, then uses the original
MDscope event serialization, filtering and statistics. In particular,
multi-residue glycolipid whole-molecule exclusions are preserved.

**Compatibility:** `_single_frame`, `_prepare` and `_conclude` are
private MDAnalysis interfaces. The CI test suite pins MDAnalysis 2.9.0 for
the real-data parity tests; changing MDAnalysis versions requires rerunning
these tests before production use.

### Residue-to-membrane distance

Each configured instance uses the original exact hybrid minimum-distance
search. Frames are stored as compact, atomic per-chunk NPZ matrices. After
the trajectory completes, the existing profile writer rebuilds the original
float32 memory-mapped distance matrix, residue metadata, and quantile
summaries from these chunks **without iterating the XTC again**.

### Checkpoint and crash recovery

All observer families write atomic part-numbered parquet or NPZ files under
`<replicate>/cache/fused_membrane/`,
`<replicate>/cache/fused_extensions/` and
`<replicate>/cache/fused_distance/`.

A part is persisted before the main streaming checkpoint commits. On
restart, uncommitted part files are discarded and only the interrupted
chunk is replayed. Once all frames have been processed, the observer
finalizers reconstruct the original public output files.

Do **not** delete these intermediate chunks until the streaming completion
marker is written. `fused_membrane_timings.json`,
`fused_extensions_timings.json`, `fused_distance_timings.json` and
`stage_timings.json` report the timing breakdown.

### Scientific correctness

- The lipid-order code uses the same fixed topology C–H bond pairs and
  initial-frame leaflet assignment as the original implementation.
- Bilayer normals and enrichment update leaflet assignment on every frame
  and use minimum-image PBC calculations.
- Protein centers are calculated using circular means in fractional box
  coordinates; coordinates are not modified in the streaming Universe.
- Lipid H-bond detection now shares the main frame loop and retains its bonded
  whole-molecule exclusion for multi-residue glycolipids.

Real-data CI compares the fused outputs to the original independent readers
on both published TPR/XTC systems, including every per-frame table and
membrane-order profile array. This is a numerical equivalence check, **not**
an assessment of sampling convergence or statistical significance.

---

### MDscope full SDC4 analysis + shareable evidence report

### Configuration

Use **`synd4_full_mechanisms.yaml`**. It preserves the **same two systems,
topology paths and all six replicate XTC paths** as
`synd4_baseline.yaml`. The output is separate:
`results_memb_full_analysis/` (relative to the YAML directory), so the
existing baseline outputs are not archived or overwritten.

The new YAML enables:

- Original RMSD, RMSF, Rg, membrane geometry, density, residue-membrane
  distances, DSSP, protein H-bonds and salt bridges
- Global lipid C-H order, local bilayer-normal deviations, leaflet-normalized
  lipid enrichment and intermolecular lipid-lipid hydrogen bonds
- Lipid neighbor matrices, protein-excluded lipid area, and
  protein-centered density maps for all configured lipid species
- Single- and multi-origin MSD with multi-origin fit diagnostics
- Isolated and environment SASA, including membrane/HS burial
- Protein, HS and lipid solvation counts and RDFs, aqueous-reference RDFs,
  and two carboxylate-centered 3D hydration maps
- Acidic-patch/HS bridge occupancy, water hydrogen-bond bridges and K spatial
  references, plus a separate topology, RDF peak, shared-water and steric audit
- Local environment contacts and the descriptive coupling screen
- Ectodomain height, tilt and membrane-contact states
- Sequence-derived NxIP occlusion and direct motif SASA (protein-only)
- **Six individual covalently attached HS chains** with orientation, direct
  glycan-lipid contacts, hydration and K/Cl coordination
- TM crossing, generic helix-helix interface and per-chain GXXXG contacts
- Local membrane thickness, quadratic curvature and C-H order
- One shared conformational-state model across all six replicates
- Hierarchical statistics, replicate-aware comparisons and figures
- A PDF evidence report and a review ZIP containing the full audit tables

Solvent identity tracking is limited to the acidic patches and HS anionic
sites. Other charged protein sites and lipid phosphates contribute counts
and RDFs. Environment SASA uses whole, centered coordinates because the
dimer exceeds the compact-target limit of periodic SASA. Interpret this
burial metric only while the target is clear of the periodic faces.
Ligand contacts require a designated ligand. The finalization-fast
variant consumes the same outputs with atom-identity lipid H-bond statistics
omitted and a relaxed shared-state completeness threshold.

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

### Running

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

### Files to send back to ChatGPT

Once all analyses and project-level finalization finish:

```text
results_memb_full_analysis/
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

### Scientific and statistical limitations

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
