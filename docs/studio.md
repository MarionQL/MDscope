# MDscope Studio

Optional local YAML authoring; free, no account or results required.

```bash
pip install -e '.[studio]'  # from a checkout
mdscope-studio
```

Open the localhost URL printed by Streamlit (also usable through VS Code port
forwarding). Save downloaded YAML in the configuration directory set in Studio.
Relative topology, trajectory and output paths resolve from that directory.

1. Set project, condition/system IDs, replica labels/paths and sampling settings
   with guided controls. Add/remove conditions and replicas; the last of each is
   protected. Advanced YAML remains available and preserves supported options.
2. Enter global selections and any reusable MDAnalysis groups. Fitting atoms can
   differ from measurement atoms. Use **Preview selections from topology** to
   inspect actual per-condition matches without loading trajectories. Group RMSD
   measures the intersection with `selections.features`; group Rg intersects its
   own selection. Group filters do not suppress the global scalar.
3. Enable analyses individually. Edit supported backend options, inspect saved
   schema/path contracts, then choose descriptive figures. Everything starts off.
4. Optionally enable statistics and select exact variables and estimators.
   Structural scalar choices are guided; other families accept exact backend
   selectors in YAML. Inference uses independent replicas, with backend coverage
   checks deciding availability.
5. Review and download. Path checks and topology preflight are explicit buttons;
   topology-required mode blocks download until the current document passes.
   Schema-only export remains marked unverified.

No trajectories, saved-result inputs, calculations or plot jobs are opened by the
wizard. Schema placeholders are **ILLUSTRATIVE — NO MD RESULTS**. The figure cards
show the production plot ID and destination; they do not compute example plots.

## Current coverage

| Capability | Status |
|---|---|
| No-data RMSD/Rg schema, selection provenance and exact figure IDs | Guided |
| Figure and inference defaults off; inactive scalar figures removed | Implemented |
| Named groups: add, duplicate, rename, delete | Implemented |
| Actual topology counts, atom preview, segment/residue labels | Explicit preview in Selections and export preflight |
| Supported analysis options | Backend settings + advanced YAML |
| Named runs (SASA, density, residue distance) | Add/delete editor with inherited per-run option panels |
| All registered artifact schemas | Displayed; uncatalogued shapes explicitly unknown |
| RMSD/Rg statistics variables and separate figure selector | Guided |
| Other statistics families | Advanced exact selectors; no fabricated IDs |
| Figure eligibility | Shared catalog; descriptive vs statistics-dependent, explicit dynamic identities |
| Statistics schema previews | Writer-verified base columns, row grain and placeholder rows for six families |
| Option help | Backend field help, selection/reference semantics and unit suffixes |
| Post-run figure editor | Deferred |

Imported options and legacy all-variable/all-figure switches are retained until
explicitly replaced by exact selectors. Explicit figure choices override conflicting
legacy family/variant switches in exported YAML; unrelated switches are preserved.
YAML comments/formatting are not retained. Review export
before replacing a configuration. Protein/features/alignment defaults are parsed
only when enabled analyses consume them. Lipid-only membrane/density and explicit
nonprotein solvation targets can pass preflight and execute without a protein.
Enabling structural analyses restores their selection requirements; this does
not extend protein-specific algorithms to unsupported chemistry. Topology success does not
verify trajectory atom ordering, chemistry for every observer, or convergence.

For wide tables, only guaranteed columns are shown; species/target metric columns
depend on topology and configuration. Array previews give axes and keys, without
inventing dimensions. JSON audit objects retain explicitly variable metadata keys.
Dynamic statistics identities can be entered in the exact-selector tables;
unavailable saved identities are rejected by the statistics backend at execution.
