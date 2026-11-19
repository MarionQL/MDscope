# Studio delivery (#173)

Studio authors configuration. It never executes trajectories, inference or
production rendering. New documents opt out of analyses, statistics and figures.

| Responsibility | Owner |
|---|---|
| Configuration schema/defaults/validation | `config`, `config_schema` |
| Numerical output schemas and variable identities | `registry`, existing producers |
| Figure IDs, variants, renderer destinations | `reporting.plot_ids` |
| Authoring figure catalog, active/inactive metadata | `planning.authoring_catalog` |
| Group/run operations, draft commits, exact export reconciliation | `studio.model` |
| Draft input widgets and cached read-only previews | `studio.app` |
| Topology preflight | `core.validation` |

## Accepted authoring behavior

- Widgets edit a detached draft. The model commits copies; reruns commit before
  restarting. Incomplete drafts remain editable; invalid exports remain blocked.
- Shared catalog supplies UI and export figure candidates. Exact export filters
  inactive catalogued selectors/figures without mutating editable choices.
- Group rename/delete migrates catalogued scalar/profile identities and imported
  and in-session figure choices, including density/distance grouped profiles.
- Run deletion protects references and removes catalogued selectors/figures.
  Disabling runs retains editable choices but excludes them from active export.
- Exact mode rejects unowned imported figures with an actionable explanation.
  Legacy mode preserves their configuration; unknown ownership is not fabricated.
- Preflight belongs to the current exported document. No XTC reads or result writes.

## Remaining acceptance work

| Item | Status |
|---|---|
| RMSD/Rg, RMSF, SASA, DSSP, density and residue-distance guided authoring | Implemented; topology/data eligibility conditional |
| Interaction/solvation descriptive figures and exact family selector tables | Implemented; dynamic identities explicitly unverified |
| Complete other producer contracts and statistical schema previews | Remaining |
| Comprehensive chemistry/protein-free preconditions (#164) | Remaining |
| Option glossary and full representative UI acceptance review | Remaining |
| Final current-head validation and #173 closure assessment | Remaining |
| Post-run figure editor | Explicitly deferred |

Use the check policy for validation; record current executed evidence in the PR,
not a growing chronological list in this document.
