# Studio delivery (#173)

Studio authors configuration. It never executes trajectories, inference or
production rendering. New documents opt out of analyses, statistics and figures.

| Responsibility | Owner |
|---|---|
| Configuration schema/defaults/validation | `config`, `config_schema` |
| Numerical output schemas and variable identities | `registry`, existing producers |
| Statistics preview base schemas | `statistics_contracts`, verified against actual writers |
| Figure IDs, variants, renderer destinations | `reporting.plot_ids` |
| Authoring figure catalog, active/inactive metadata | `planning.authoring_catalog` |
| Group/run operations, draft commits, exact export reconciliation | `studio.model` |
| Draft input widgets and cached read-only previews | `studio.app` |
| Guided project/system/replica/sampling inputs | `studio.project_editor`; system mutations in `studio.model` |
| Topology preflight | `core.validation` |
| Saved output verification and offline finalization | `pipeline.persisted_inputs`; execution gate in `pipeline.runner`, Studio adapter in `studio.model` |

## Accepted authoring behavior

- Widgets edit a detached draft. The model commits copies; reruns commit before
  restarting. Incomplete drafts remain editable; invalid exports remain blocked.
- Shared catalog supplies UI and export figure candidates. Exact export filters
  inactive catalogued selectors/figures without mutating editable choices.
- Group rename/delete migrates catalogued scalar/profile identities and imported
  and in-session figure choices, including density/distance grouped profiles.
- Run deletion protects references and removes catalogued selectors/figures.
  Disabling runs retains editable choices but excludes them from active export.
- The single hierarchical selector preserves imported legacy YAML on no-op export.
  When an exact-ID export is requested, unowned imported figure IDs are rejected
  with an actionable explanation; no separate UI mode exists.
- Preflight belongs to the current exported document. No XTC reads or result writes.

## Acceptance coverage

| Item | Status |
|---|---|
| RMSD/Rg, RMSF, SASA, DSSP, density and residue-distance guided authoring | Implemented; topology/data eligibility conditional |
| Interaction/solvation descriptive figures and exact family selector tables | Implemented; dynamic identities explicitly unverified |
| Registered producer contracts | Fixed/base/dynamic columns, array axes/keys, row semantics; JSON keys explicitly variable |
| Statistical schema previews | Six families verified against nonempty writer products; conditional provenance explicit |
| Protein-free prerequisites (#164) | Shared role resolver; real lipid-only membrane/density and ligand/water execution tests; structural negatives |
| Option help and representative UI review | Shared field help; AppTest covers opt-in figures/statistics, groups, imports and named runs |
| Final acceptance evidence / merge gates | See PR #200; CI matrix and applicable real-data checks must pass before merge |
| Post-run figure editor | Explicitly deferred |

Use the check policy for validation; record current executed evidence in the PR,
not a growing chronological list in this document.

## Review limits

- Schema previews do not estimate results. Wide chemistry-dependent columns and
  saved interaction identities remain explicit dynamic contracts.
- Topology preflight is not a universal chemistry audit or trajectory-order check.
  Protein-specific analyses retain their own validated scope.
- Automated Streamlit interaction checks cover behavior; they do not establish
  visual layout quality on every browser/display. Short real-data smoke/parity
  checks establish integration, not statistical convergence.

## Remaining release work

- Review the local UI on representative browser/display sizes. AppTest covers
  behavior, not visual layout.
- Recheck latest-head automatic CI after follow-ups; offline finalization and
  topology preflight hash implementation changes require their own regression
  review even when physical estimators are unchanged.
- Review/merge PR #200 and assess #173 closure. Post-run figure editing stays
  outside this release.
