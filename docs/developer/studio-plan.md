# Studio delivery (#173)

Studio authors configuration; it never starts simulation analysis, inference or
production plotting. New documents explicitly disable every analysis, statistics
family and figure. Imported documents retain their options.

| Boundary | Owner |
|---|---|
| Defaults, supported keys, semantic validation | `config`, `config_schema` |
| Analysis outputs, scalar variable identity | `registry`, exposed by `planning` |
| Figure IDs, sources, destination | `reporting.plot_ids` |
| Metric selection and estimator capabilities | `statistics_selection` |
| Topology parsing and chemistry checks | `core.validation` |
| Editable document, active figure plan | `studio.model` |
| Local optional widgets | `studio.app` |

Delivery order: shared in-memory config resolution → no-data RMSD slice → all
analysis settings and registered outputs → groups/topology preflight → exact
statistics selection → review/export and UI integration tests.

Risk checks:
- Preserve exact selection text, group/run order and imported options.
- A disabled analysis cannot contribute selected figures; retain choices in the
  editable session so re-enabling restores them.
- Never silently enable a statistics family or select all metrics.
- Mark uncatalogued dynamic schemas/IDs unknown; no fabricated column names.
- Validation is tied to the current document; edits invalidate prior preflight.
- Resolve relative paths from an explicit configuration directory. No XTC read,
  result overwrite or output-folder creation during authoring.

Deferred gaps must remain visible in Studio and the PR. Full #173 closure requires
all supported producer schemas/variable identities and their figure eligibility,
plus topology/chemistry coverage beyond structural selections. Advanced YAML
authoring does not imply that every option has a guided widget.

## Delivery progress

- [x] Guided exact identity/estimator tables for all backend statistics families.
- [x] Registered artifact cards show illustrative columns, schema status, units
  and resolved per-replica destinations without saved inputs.
- [x] DSSP assignment columns, code/time units and residue/frame row grain.
- [ ] Complete remaining producer schemas and column definitions.
- [ ] Derive nonstructural figure eligibility and statistics previews from those contracts.
- [ ] Resolve protein-free preconditions (#164) and complete acceptance review.

Validation for this increment: 34 Studio/planner tests passed on Python 3.13.
Previous PR head passed Python checks, figure review and real GROMACS smoke;
new changes require their own CI results.
- [x] H-bond event schemas verified against producer fields, including lipid molecule/species identity.
- [x] Interaction/solvation descriptive figure candidates use production plot IDs.
- [x] Exact profile/interaction selectors offer separate statistical figure choices.
  Saved-data availability remains conditional; unspecified interaction metrics do
  not create guessed statistical plot IDs.
- [x] RMSF global/group atom-profile schemas and source/fitted provenance;
  temporal-inference cache prerequisite visible. Global schema checked against
  actual synthetic pipeline output (not inferred from RMSD tables).
- [x] Membrane APL/thickness and DSSP helix/sheet/coil fractions now share the
  scalar authoring flow: units, schema/path, exact descriptive figures and
  exact-variable statistics/figure selection. Suboption disabling removes the
  corresponding planned metric. Scientific producers remain unchanged.
- [x] Check policy incorporated from merged PR #201; Studio test dependency retained.
- [x] Guided opt-in profile selectors for RMSF/groups, named SASA/burial runs and
  DSSP occupancies, with producer-matching identities, units and prerequisites.
  Exact selectors feed the existing statistical figure picker. Imported
  inference:false remains disabled when editing guided estimator switches.
- [x] Named-run editor for the three backend-supported multi-run analyses:
  add/delete controls and per-run option panels; backend validates IDs, inherited
  settings and references. SASA UI test covers run creation → saved path → exact
  profile selector → opt-in inference with no automatic plots.
- [x] Named-run lifecycle: disabled run selectors excluded from exact export but
  retained in editable state; deletion removes catalogued SASA profile selectors
  and imported figure IDs. Referenced runs are protected from deletion, including
  SDC4 SASA references. Last-run deletion requires explicitly leaving named mode.

## Density and distance authoring milestone

- [x] Density named-run/target profile identities, number/mass units, global and
  conditional group table paths and column previews.
- [x] Distance run mean/median and optional cutoff-fraction profile selectors,
  conditional grouped profiles and separate descriptive per-condition heatmaps.
- [x] Both families feed existing exact statistics selectors/figure choices;
  disabled runs and deleted distance heatmaps are removed from active export.
- [x] Distance `groups` uses the backend boolean contract; only RMSD/Rg use the
  shared-group multiselect. No backend `groups` list is invented.
- [x] Density schema/identity tested against actual saved producer output;
  distance figure deletion and boolean UI contracts covered by integration tests.

Dynamic dimensions, empty topology-dependent intersections, sampling eligibility
and unregistered optional outputs remain explicitly conditional. This completes
this authoring milestone, not all acceptance requirements of #173.
