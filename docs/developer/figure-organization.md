# Figure organization audit — issue #171

## Contract

- Scientific owner and presentation tier are independent of storage/statistics family.
- `contracts/plot_validation.py`: IDs, prerequisites, tiers, legacy and canonical destinations.
- `registry.figure_descriptor`: exact configured variable ownership, reusing the statistics/table owner. Unrecognized identities are `unmapped`; ambiguous normalized identities fail before publication.
- `reporting/plot_ids.py`: reversible filenames, presentation text and disk saving. Existing IDs survive migration.
- `reporting/templates.py`: paginated, unit-explicit heatmaps shared by interaction, lipid-contact and standardized residue-response renderers. `overview.py`: replica violins and saved-effect intervals; existing profile/trace drawing paths remain shared production renderers.
- `figure_generation`: prepare, validate and publish with rollback; obsolete selected pages retire only after successful rendering. Unselected user files survive. Per-analysis tables share `<analysis>/secondary/` with detail plots.

The generated [migration matrix](figure-migration.csv) audits maintained YAML identities plus synthetic saved-data examples. It is review evidence, not a metadata registry. Unknown dynamic identities remain visible instead of acquiring an owner from a guessed biological label. Diffusion/shared-state PNGs move into the hierarchy; their numerical tables stay in place. Historical replica SDC4 PNGs remain rejected legacy products. PDFs retain exact source-row/page identities.

## Reproduce visual review

```bash
python scripts/plot_preview.py --output /tmp/mdscope-full
python scripts/plot_preview.py --output /tmp/mdscope-primary --preset primary
python scripts/plot_preview.py --output /tmp/mdscope-generic-primary --preset primary --generic
python scripts/plot_preview.py --output /tmp/mdscope-generic-full --generic
python scripts/preview_gallery.py --preview /tmp/mdscope-full --output /tmp/mdscope-full/index.html
python scripts/figure_gallery.py --base /tmp/mdscope-before --head /tmp/mdscope-full --output /tmp/mdscope-review/index.html
```

- Previews use production renderers and deterministic saved Parquet/NPZ fixtures; each image says **SYNTHETIC EXAMPLE — NOT MD RESULTS**. Different revisions match by stable ID even when files move.
- Primary views preserve zeros, signed effects and finite replica support. Units/estimands get separate facets. Diamonds on effect plots mean finite ordered CI excludes zero, independently of FDR; distributions show one violin and mean per replica, never frame-pooled inferential CIs.
- Large interaction/species/residue heatmaps and saved table views paginate without ranking away rows. Solvation RDF, trace and survival pages preserve every supported group. Configured occupancy filtering and documented zero-only footprint suppression remain explicit legacy display policies.
- Review includes named groups, uneven/missing replicas, invalid/NaN bounds, all-zero matrices, 61×61 species matrices, source/config staleness, migration, interrupted publication and page-count shrinkage.
- This validates saved-data rendering and organization. Full six-replica production readability/convergence remains an ARC review; no synthetic fixture establishes biological significance.

## Evidence

Measured counts, timings, executed tests and CI artifact links are recorded in the PR and linked from plotting QA #172. The unchanged default fixture renders 51 PNGs before and after migration. The primary SDC4 fixture renders 18 PNGs (51 full); the generic fixture renders 12 primary PNGs (33 full). Primary/full SDC4 previews took 8.90/21.65 s, with child peak RSS about 316/387 MiB in this Python 3.12 environment. These include fixture construction and synthetic stamping, not production trajectory costs. Later fixture additions can change these counts; use `preview_manifest.json` as the exact run inventory.
