# MDscope Studio

Optional local YAML authoring; free, no account or results required.

```bash
pip install -e '.[studio]'  # from a checkout
mdscope-studio
```

Open the localhost URL printed by Streamlit (also usable through VS Code port
forwarding). Save downloaded YAML in the configuration directory set in Studio.
Relative topology, trajectory and output paths resolve from that directory.

1. Set project, condition/system IDs, replica paths and sampling settings.
   Reusable measurement groups now live here. Topology parsing starts in a
   background thread, and the panel shows per-system atom/residue counts,
   SegIDs, chain IDs, residue names, and group matches as results arrive.
   The rest of Studio remains editable during this work. No trajectories load.
2. Configure global selections and individual analyses. Group-specific structural
   fitting selections are edited under RMSD, not project setup. Group RMSD
   measures the intersection with `selections.features`; group Rg intersects
   its own selection. Group filters do not suppress the global scalar.
3. Optionally configure statistics and its hierarchical metric selector.
4. Choose optional figures through the hierarchical figure selector.
5. Review and export. Explicit production topology preflight is still required
   when topology-verified export is enabled; background metadata does not
   replace full validation.

The wizard does not open trajectories or run analyses, statistics, or plot jobs.
Saved products are inspected read-only **only** when explicitly requested through
**Inspect existing results**; Studio also reads topologies automatically in a worker for lightweight metadata and group validation.
Schema placeholders are **ILLUSTRATIVE — NO MD RESULTS**. Figure entries show
production plot IDs and destinations without computing example plots.

## Current coverage

| Capability | Status |
|---|---|
| No-data RMSD/Rg schema, selection provenance and exact figure IDs | Guided |
| Figure and inference defaults off; inactive scalar figures removed | Implemented |
| Named groups: add, duplicate, rename, delete | Implemented |
| Actual topology counts, atom preview, segment/residue labels | Automatic, nonblocking metadata/group inspection in Project & systems; explicit export preflight |
| Supported analysis options | Backend settings + advanced YAML |
| Named runs (SASA, density, residue distance) | Add/delete editor with inherited per-run option panels |
| All registered artifact schemas | Displayed; uncatalogued shapes explicitly unknown |
| RMSD/Rg statistics variables and separate figure selector | Guided |
| Other statistics families | Advanced exact selectors; no fabricated IDs |
| Figure eligibility | Shared catalog; descriptive vs statistics-dependent, explicit dynamic identities |
| Statistics schema previews | Writer-verified base columns, row grain and placeholder rows for six families |
| Option help | Backend field help, selection/reference semantics and unit suffixes |
| Post-run figure editor | Deferred |

Imported options and legacy all-variable/all-figure switches are retained on
no-op export until the user edits their hierarchy; there is no separate
Legacy/Exact interface. Explicit figure choices override conflicting
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
