# MDscope Studio

Optional local YAML authoring; free, no account or results required.

```bash
pip install -e '.[studio]'
mdscope-studio
```

- Open Streamlit's localhost URL; VS Code port forwarding also works.
- Save exported YAML in the configuration directory chosen in Studio. Relative input/output paths resolve there.
- Studio reads topology and explicitly requested saved results; it runs no trajectories, inference or production rendering.
- Previews are **ILLUSTRATIVE — NO MD RESULTS**. Figure entries show production IDs/destinations, not computed example plots.

## Author a configuration

1. **Project & systems:** project, conditions, replicas, paths, sampling and reusable groups.
2. **Selections & analyses:** global/analysis settings and named runs. Group fitting belongs under RMSD. Group RMSD intersects `selections.features`; group Rg intersects its own selection. Group filters retain the global scalar.
3. **Statistics:** choose inference and exact metric identities.
4. **Figures:** select bootstrap/replicate summaries, condition comparisons and autocorrelation independently, or edit exact IDs.
5. **Review & export:** run explicit production preflight when topology-verified export is enabled. Background metadata does not replace it.

### Topology inspection

- Nonblocking background parsing shows atom/residue counts, SegIDs, chain IDs, residue names, previews and group matches.
- Group edits reuse cached topology; changed topology path/file revision invalidates it. No XTC is loaded.
- Protein/features/alignment defaults are checked only when consumed by enabled analyses. Supported lipid-only and explicit nonprotein solvation configurations require no protein.
- Structural analyses restore their selection requirements. Successful preflight does not establish trajectory atom order, every observer's chemistry or convergence.

## Editor coverage

| Capability | Interface |
|---|---|
| Project settings, groups, RMSD/Rg, selection provenance | Guided; groups support add/duplicate/rename/delete |
| Supported analysis options | Backend controls plus Advanced YAML |
| Named SASA/density/residue-distance runs | Add/delete, inherited option panels |
| RMSD/Rg inference | Guided variables; other families use exact selectors |
| Figures | Shared catalog, exact IDs, descriptive/statistical prerequisites |
| Artifact schemas | Registered columns, row grain and axes; unknown shapes remain unknown |
| Statistics previews | Writer-verified base schemas for six families |
| Option help | Backend semantics, references and units |
| Post-run figure editor | Deferred |

## Import and export

- New inference/figure defaults are off; inactive scalar figures are removed.
- Imported legacy all-variable/all-figure switches survive no-op export until the corresponding hierarchy is edited; no separate Legacy/Exact mode.
- Explicit figures override conflicting legacy family/variant switches; unrelated switches survive.
- YAML comments/formatting do not survive export. Review before replacing an existing configuration.
- Wide previews show guaranteed columns only; topology-dependent dimensions/metrics and variable JSON keys are not invented.
- Exact dynamic identities are checked against saved data during backend execution.

## Inspect existing results

- Inspection is explicit and read-only, using execution's persisted-input validator.
- Statistics/Figures show replica coverage and validation reasons. Scheduled work is distinct from verified products; incomplete data, scientific-setting mismatches and unknown provenance cannot authorize reuse.
- Draft/results-path changes invalidate inspection; refresh it. Saved-only inputs do not enable trajectory analyses in exported YAML.
- Verified raw RMSF supports static figures. Temporal RMSF inference additionally needs matching feature-atom metadata and the complete aligned-coordinate cache; finalization rechecks it before replacing results.

Settings and reuse rules: [Configuration](configuration.md#existing-results-and-readiness).
