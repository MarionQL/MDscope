# MDscope trajectory pass consolidation

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

## Number of XTC traversals

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
execution is controlled by MDAnalysis. Extra reads may also occur when
rebuilding missing/corrupted outputs.

## Hydrogen bonds and DSSP

HydrogenBondAnalysis and DSSP are not reimplemented with new physical
approximations. The shared frame loop calls the original MDAnalysis
`_single_frame()` detector on each raw frame, then uses the original
MDscope event serialization, filtering and statistics. In particular,
multi-residue glycolipid whole-molecule exclusions are preserved.

**Compatibility:** `_single_frame`, `_prepare` and `_conclude` are
private MDAnalysis interfaces. The CI test suite pins MDAnalysis 2.9.0 for
the real-data parity tests; changing MDAnalysis versions requires rerunning
these tests before production use.

## Residue-to-membrane distance

Each configured instance uses the original exact hybrid minimum-distance
search. Frames are stored as compact, atomic per-chunk NPZ matrices. After
the trajectory completes, the existing profile writer rebuilds the original
float32 memory-mapped distance matrix, residue metadata, and quantile
summaries from these chunks **without iterating the XTC again**.

## Checkpoint and crash recovery

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

## Scientific correctness

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
