# MDscope trajectory pass consolidation

MDscope's regular protein/membrane streaming pass now also computes the three
lipid observables that previously opened the same XTC separately:

- **Lipid C–H order** (\`membrane_order\`), including original tail/site
  mapping, per-leaflet order and temporal profile NPZ outputs.
- **Local bilayer normals** (\`bilayer_normal\`), including optional per-lipid
  angles and dynamic leaflet assignment.
- **Protein-centered lipid enrichment** (\`lipid_enrichment\`), including raw
  counts, selection audits and PBC-aware protein centers.

The observer runs **before any coordinate alignment**, uses the original
topology-based selections, and shares anchor coordinates and leaflet
classification when the normal/enrichment selection specifications match.
Outputs retain their existing paths and column schemas; downstream statistics
and plots do not need to change.

## Number of XTC traversals

For the tested configuration with RMSD/Rg, the three lipid observables and
lipid–lipid hydrogen bonds:

| Mode | Before | Now |
|---|---:|---:|
| Main streaming (RMSD/Rg/etc.) | 1 | 1 |
| Lipid C–H order | 1 | 0 additional |
| Bilayer normals | 1 | 0 additional |
| Lipid enrichment | 1 | 0 additional |
| Lipid H-bonds (MDAnalysis HBA) | 1 | 1 |
| **Total** | **5** | **2** |

If the regular streaming pass is disabled and only the three lipid observables
are requested, a shared geometry-only pass runs once. No unnecessary protein
streaming pass is forced.

Other optional analyses can add passes: DSSP, SASA, residue-to-membrane
distance, local environment and protein hydrogen bonds still have their own
readers. The physical validation test also deliberately traverses trajectories
separately and is **not** part of normal MDscope execution.

## Checkpoint and crash recovery

Fused observers write part-numbered parquet and NPZ files under
\`<replicate>/cache/fused_membrane/\`. A part is written before the streaming
checkpoint commits. On restart, any files at or beyond the next uncommitted
part number are discarded; committed parts are retained. When the trajectory
completes, MDscope assembles the original output files.

Do **not** delete the \`cache/fused_membrane\` intermediates before the
streaming completion marker has been committed; they are necessary for safe
finalization after a crash. Files can be removed manually after successful
completion if the output archive is backed up. Finalized outputs live in the
usual \`profiles/membrane_order\` and \`extensions\` directories.

\`fused_membrane_timings.json\` reports per-analysis compute times (which
exclude shared trajectory decoding), chunk-write/finalization time, and the
number of unique lipid anchor mappings. The parent \`stage_timings.json\`
reports the total streaming wall time.

## Scientific correctness

- The lipid-order code uses the same fixed topology C–H bond pairs and
  initial-frame leaflet assignment as the original implementation.
- Bilayer normals and enrichment update leaflet assignment on every frame
  and use minimum-image PBC calculations.
- Protein centers are calculated using circular means in fractional box
  coordinates; coordinates are not modified in the streaming Universe.
- The lipid H-bond pass remains separate and retains its bonded whole-molecule
  exclusion for multi-residue glycolipids.

Real-data CI compares the fused outputs to the original independent readers
on both published TPR/XTC systems, including every per-frame table and
membrane-order profile array. This is a numerical equivalence check, **not**
an assessment of sampling convergence or statistical significance.
