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

### Full SDC4 analysis

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
- Local environment contacts and matched-window coupling comparisons
- Ectodomain height, tilt and membrane-contact states
- Sequence-derived NxIP occlusion and direct motif SASA (protein-only)
- **Six individual covalently attached HS chains** with orientation, direct
  glycan-lipid contacts, hydration and K/Cl coordination
- TM crossing, generic helix-helix interface and per-chain GXXXG contacts
- Local membrane thickness, quadratic curvature and C-H order
- One shared conformational-state model across all six replicates
- Hierarchical statistics, replicate-aware comparisons and a findings report

Solvent identity tracking is limited to the acidic patches and HS anionic
sites. Other charged protein sites and lipid phosphates contribute counts
and RDFs. Environment SASA uses whole, centered coordinates because the
dimer exceeds the compact-target limit of periodic SASA. Interpret this
burial metric only while the target is clear of the periodic faces.
Ligand contacts require a designated ligand. The finalization-fast
variant targets the earlier saved replica products. It omits the later
shared-phosphate distance profile and additional motif atom SASA inputs,
skips atom-identity lipid H-bond statistics, and uses a lower shared-state
completeness threshold.

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
python -m pip install freesasa
python -m mdscope.cli validate synd4_full_mechanisms.yaml
python -m mdscope.cli run synd4_full_mechanisms.yaml
```

Use the existing Slurm/array runner for individual replicates if preferred.
After all six replicates have completed, run:

```bash
python -m mdscope.cli finalize synd4_full_mechanisms.yaml
```

`finalize` and `statistics` recompute derived numerical outputs and
archive the previous project summary. `run` also performs this phase
after processing replicas.

### Lipid footprint and residue hydration

With the updated SDC4 YAML, the normal replicate workflow adds the
`extensions/lipid_footprint/` products and the residue water columns in
`extensions/sdc4_mechanisms/residue_occlusion.parquet`. Completed replicas
reuse their other persisted analyses. When residue occlusion is already
available, the hydration pass reuses its saved contact fractions. Regenerate
project statistics after all replicas finish so the new footprint metrics
replace the former centroid-based enrichment in the consolidated results.

### Membrane-relative density

The SDC4 configurations use the named `membrane_relative` run. At each
frame, MDscope separates the shared phospholipid phosphate anchors into
two periodic leaflets and centers the density bins at the midpoint of
their circular means. The two leaflets have equal weight even when their
populations differ. Positions are measured along the configured planar
membrane-normal axis, z for SDC4. The method follows bilayer translation,
including periodic crossings, and preserves leaflet asymmetry. It uses a
fixed physical grid spanning ±100 Å for these systems. Each frame's slab
volume is its cross-sectional area times bin width.

For completed replicas, run MDscope with the updated YAML to add the new
density profiles. Only the missing named runs require coordinates. The
previous box-referenced profiles are retained under their original ID.
Statistics can rebin old profiles to shared slab edges, but that operation
cannot recover the per-frame membrane midplane. This global profile also
retains broadening from membrane undulations and does not fit a local
surface or a tilted plane.

Project statistics write `statistics/coverage/contrast_coverage.parquet`
under `project_summary`. The audit records valid intervals and missing
comparisons by observable and system pair for time series, profiles,
interactions, RDF shells, diffusion and coupling. An interval excluding
zero is distinct from a multiple-testing-adjusted permutation result.

### Discovering contiguous protein changes

Profile statistics write `project_summary/statistics/profiles/residue_regions.parquet`.
To scan an older completed statistics export directly, run:

```bash
mdscope scan-residue-regions /path/to/project_summary/statistics
```

The scan reads saved replica profiles and their statistics manifest. It
compares the DSSP, RMSF, SASA, and residue–membrane-distance profiles already
present in the export. Rerun the statistics phase to include named SASA runs
and environmental buried area from saved per-replica SASA files.

A region contains at least three consecutively numbered residues on the
same labeled chain. At each residue every expected replica must have a
finite value and the two replica ranges must be strictly separated in the
same direction throughout the region. Each condition needs at least two
replicas. Gaps, ambiguous residue numbering, overlap and direction changes
break regions. Within a profile and system pair, the table ranks regions by
the absolute mean residue contrast and reports the smallest separation of
the two replica ranges in the region. Positive contrasts mean system 2 is
higher. SASA contrasts are per-residue mean differences in nm².

The table is exploratory. Absence from it can reflect overlapping replica
ranges or incomplete data. Region boundaries and ranks are descriptive
and do not carry confidence intervals. Pointwise confidence intervals
remain in the profile comparison table.

RDF shell contrasts use the shell-volume-weighted mean of the saved
normalized RDF up to the configured water or ion cutoff. Diffusion
contrasts refit all replica curves on the same supported physical lags,
with the fit range bounded by the shortest saved curve. These two
families resample independent replica estimates. They require complete
replica coverage and at least two replicas per condition (or two matched
pairs). Diffusion intervals describe the effective MSD slopes over those
lags. Check the per-replica scaling diagnostics and source trajectory
windows before interpreting them as a diffusivity change.

For local upper height, upper curvature and local thickness, statistics
export the fraction of frames with a usable estimate. Conditional
geometry intervals resample the original regular frame grid with the
missing entries preserved. The conditional interval is unavailable when
the observed frames do not support it.

The coupling screen compares replicate Pearson coefficients over the
shared regular frame grid. It keeps missing observations at their original
times and excludes incomplete blocks. Its confidence interval describes
the difference in mean replicate correlation and is exploratory.
Set `statistics.matched_window_ps` to `auto` or explicit `[start, stop]`
bounds to keep a second scalar statistics pass in
`statistics/matched_time_series/`, alongside the default full-window
tables. Explicit bounds must be inside every saved replica time range.
The shared-state analysis compares full- and matched-window state
populations using the same state assignments. Use
`statistics.scalar_window_ps` to restrict the primary scalar pass itself. For a matched trajectory-time MSD
estimate, set `analyses.membrane.msd.multi_origin_window_ps` to those
bounds when generating the replica MSD curves. Fitted curves record their
actual sampled start and end times.

The full SDC4 production configuration includes a common phospholipid
phosphate distance profile alongside the composition-dependent nearest
lipid profiles. It also writes distinct HS residue and chain contact
counts, nonoverlapping near/far chain order by lipid species and leaflet,
and isolated/environment acidic-patch SASA and acidic-oxygen SASA.
The configured Ser138 site has signed C-alpha height relative to a fitted
local phosphate plane and a separate exact nearest-lipid distance,
species and topology-local molecule identifier. Configure additional sites
with `analyses.sdc4_mechanisms.residue_surface_sites`. Height is available
when enough nearby anchors constrain the fit. Scalar statistics derive HS
occupancy and contact-conditional counts from the saved series and export
censored HS contact episodes from complete, regularly sampled traces.
The new trajectory observables require refreshed replicate products.
To apply updated statistics to legacy results first, use the configuration
that generated their existing replicate products.

### Reviewing saved findings

The numerical phase writes the review under `project_summary/statistics/`.

| File | Contents |
| --- | --- |
| `findings.md` | Coverage, missing comparisons and representative interval exclusions by family and condition pair |
| `findings.parquet` | All finite between-condition bootstrap intervals that exclude zero, with exact values and source paths |
| `output_inventory.parquet` | Saved numerical products and Parquet row counts, including replica outputs when present |
| `coverage/contrast_coverage.parquet` | Every registered contrast and its interval status |

The review includes available shared-state population and conditioned
NxIP contrasts. Descriptive outputs appear in the inventory, which excludes
working caches and trajectories. When sharing an export, include
`statistics/` and its sibling `shared_states/` directory if available.

```bash
mdscope report /path/to/project_summary/statistics
```

For diagnostic figures, set `plots.enabled: true` and run
`mdscope plotting analysis.yaml`. Both commands use saved results.

With three independent replicas per condition, exact permutation tests have
limited p-value resolution. Inspect temporal stability, shared-state
populations and replicate agreement before interpreting the contrasts.
Lipid density maps retain their documented coordinate and area references.

### Per-residue occlusion by contact partner

Set `analyses.sdc4_mechanisms.residue_occlusion: true` to enable collection. Each frame records the fraction of each
protein residue's heavy atoms within `contact_cutoff_A` of nonlocal protein,
HS, or lipid heavy atoms, plus the union across partners. The protein
fraction excludes the same residue and residue numbers within two positions
on the same chain. Interchain contacts count even when residue numbers
match. The HS and lipid selections are the same as in the other SDC4
mechanism outputs. The partner fractions can overlap. The union counts
each contacted atom once. HS contacts may include atoms near the covalent
attachment site. These are distance-based contact fractions.

Run `mdscope run analysis.yaml` or use the usual per-replica command.
New replicas collect occlusion during streaming. For completed replicas,
MDscope reuses the existing outputs and reads the trajectory only if the
occlusion table is missing. It checks the saved SDC4 frame count and
timestamps before adding `extensions/sdc4_mechanisms/residue_occlusion.parquet`.
Use the same selections and sampling as the original run. The statistics
phase uses that table for per-residue comparisons and the region scan.

The selected TM GXXXG motif was detected in PROA but **not** PROB in the
audited older topology. The analysis therefore also reports generic
helix-helix and each available motif-to-opposite-TM contacts instead of
mislabeling a nonexistent second GXXXG motif.

**Before a long production analysis:** run `validate` and verify that
the selected atoms, six HS attachments, NxIP motifs, K/Cl ion names, lipid
anchor counts and time window are as intended. This is especially
important if the topology was rebuilt.
