# MDscope Studio

Optional local YAML authoring; free, no account or results required.

```bash
pip install -e '.[studio]'  # from a checkout
mdscope-studio
```

Open the localhost URL printed by Streamlit (also usable through VS Code port
forwarding). Save downloaded YAML in the configuration directory set in Studio.
Relative topology, trajectory and output paths resolve from that directory.

1. Set project, condition/system IDs, replica labels/paths and sampling settings.
2. Enter global selections and any reusable MDAnalysis groups. Fitting atoms can
   differ from measurement atoms. Group RMSD measures the intersection with
   `selections.features`; group Rg intersects its own selection. Group filters
   do not suppress the global scalar.
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
| Actual topology counts, atom preview, segment/residue labels | Production preflight |
| All supported analysis options / named runs | Backend settings + advanced YAML |
| All registered artifact schemas | Displayed; uncatalogued shapes explicitly unknown |
| RMSD/Rg statistics variables and separate figure selector | Guided |
| Other statistics families | Advanced exact selectors; no fabricated IDs |
| All-analysis guided figure eligibility and complete schema glossary | Remaining #173 work |
| Post-run figure editor | Deferred |

Imported options and legacy all-variable/all-figure switches are retained until
explicitly replaced by exact selectors. YAML comments/formatting are not retained. Review export
before replacing a configuration. The existing backend still requires nonempty
protein/features even for some protein-free configurations (#164); Studio reports
that preflight error rather than claiming validation. Topology success does not
verify trajectory atom ordering, chemistry for every observer, or convergence.
