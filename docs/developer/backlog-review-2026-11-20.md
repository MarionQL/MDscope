# Backlog review — 2026-11-20

Reviewed all **15 open issues** before creating four focused follow-ups. Evidence:
current issue bodies/comments, merged PR metadata, and source at `dc06153`
(the production source reviewed in PR #214). Issue state alone is not proof that
work remains, and a merged preliminary PR is not proof of complete acceptance.
No issues were closed without evidence that all remaining criteria were met.

## Architecture work is now tracked

| Recommendation | Issue | Decision |
|---|---|---|
| PDF/PNG report publication, identity and provenance | [#215](https://github.com/MarionQL/MDscope/issues/215) | New bounded reporting milestone |
| Shared manifest facts and explicit legacy policies | [#216](https://github.com/MarionQL/MDscope/issues/216) | New bounded consumer-policy milestone |
| Saved-result inspection and contextual readiness | [#202](https://github.com/MarionQL/MDscope/issues/202) | Extend existing issue; no duplicate |
| SDC4 numerical export ownership outside rendering | [#217](https://github.com/MarionQL/MDscope/issues/217) | New bounded reporting milestone |

These implement the recommendations in [the architecture review](architecture-review-2026-11-20.md).
The broader incremental ownership umbrella remains #166.

## Priority recommendations

1. **#218 (P1): YAML control over residue-response discovery.** Newly found
   configuration gap: profile statistics currently call the exploratory screen
   automatically; minimum length 3, similarity 0.7 and eligible metric classes
   are source-defined. Figure selection is not numerical-stage selection.
   Give users validated YAML control and document backward-compatible defaults.
2. **#202 (P1): saved-result readiness.** Some marker/config/schema validation
   already exists. Align inspection facts and explicit reuse/finalization
   policies, including time-grid/checkpoint guarantees and missing coverage.
3. **#173 (P1): acceptance audit of merged Studio.** Close stale release gates
   only with representative layout and applicable real-data evidence. PR #200
   is merged; continuing to tell agents to keep it open is obsolete.
4. **#215 then #216 (P2): report visibility and completion consistency.** These
   are reproduced reporting/policy gaps with bounded saved-data verification;
   not claims of incorrect estimators. #217 is an independent smaller follow-up.
5. **#103 (P1 scientific evidence):** retain its priority when interpreting
   acidic-patch/K+/HS mechanisms. Distance-based bridges and temporal references
   do not complete the requested chemistry, spatial control or cooperativity
   validation. This is a scientific interpretation gate, not a reproduced
   blocking runtime defect.

The new YAML-control issue is a bounded slice of #203, not a second region
inference implementation. Choosing its exact schema is implementation work;
this review adds no unsupported keys to user configurations.

## Existing issues: disposition and remaining work

| Issue | Assessment | Action / remaining outcome |
|---|---|---|
| #103 · ion bridging | Partial scientific scope; still relevant | Keep open for chemistry, spatial controls and real-system evidence; do not equate existing bridge counts with cooperativity |
| #157 · distance-resolved membrane hydration | Relevant new measurement | Normalize P2 title; add YAML geometry/bin/species/leaflet control and retained per-replica denominators |
| #163 · generic shared states | Still coupled to SDC4 | Existing feature list is YAML, but loader/validation are SDC4-specific; require YAML feature sources and preserved preset parity |
| #165 · lipid naming presets | Partly supported already | Tails and per-lipid overrides already exist; focus optional verified presets and topology/chemistry preflight |
| #166 · architecture umbrella | Stale implementation status | PRs #184–#190 and recovery work landed; current graph has no SCCs. Replace broad cleanup with focused #215/#216/#202/#217 work |
| #167 · generic motifs | Partial, not complete | PR #204 hardens shared matching; generic YAML definition, aggregation/catalog and CLI integration remain |
| #168 · helix/domain geometry | Still relevant | Generic selected domains/pairs/references must be YAML-addressable; preserve SDC4 kernels and outputs |
| #169 · generic glycan hydration | Still relevant | Reuse fused counting; configure glycan/attachment chemistry, preserve SDC4 attachment identities |
| #170 · membrane references | Still relevant | YAML-supported validated reference models; metadata and unsupported-geometry checks, not new unvalidated physics |
| #171 · generic plot templates | Partly supported already | Exact plot IDs, generic renderers and real-renderer previews exist; reuse those, add only demonstrated template/preset needs |
| #173 · Studio | Stale pre-merge language | PR #200 merged; keep only remaining acceptance/closure work current; no separate Legacy/Exact UI |
| #181 · parallelism | Defer without new bottleneck evidence | Profiling #193 and interaction optimization #197 landed; job-array/worker/RSS/I/O comparison remains, no speculative concurrency |
| #182 · release/data separation | Blocked external prerequisite | Proposed target repo returned 404 to this connection; keep original assets intact until replacement assets/access/consumer checks succeed |
| #202 · saved-result eligibility | Partial implementation; high priority | Recognize current marker/config/schema validation; strengthen shared time/identity/readiness evidence and exact coverage |
| #203 · residue-response regions | Preliminary implementation landed | PR #205 is exploratory; topology identity, adjacency, coverage and selection-aware inference remain; #218 owns YAML controls |

Authoritative current-status notes were added to stale issues without erasing
historical evidence. None of #167/#173/#203 was closed solely because a related
PR merged. New architecture issues #215–#217 retain the existing public YAML,
output and scientific contracts.

## YAML acceptance rule

All relevant generalization issues now require the
[YAML generalization contract](../configuration.md#yaml-contract-for-generalization):
new scientific choices work from validated YAML/CLI, explicit overrides remain,
Studio round-trips the same choices, and the maintained SDC4 workflow stays
reproducible. Kernel-only extractions need no gratuitous configuration knobs.

Issue-specific requirements include feature source selection (#163), verified
lipid presets/overrides (#165), motif scope/policy (#167), chosen helix pairs and
reference (#168), glycan attachment chemistry (#169), supported membrane
references (#170) and optional presentation presets (#171). #157 and #203 also
carry YAML acceptance criteria. Existing scientific invariants remain mandatory.

## Executed validation

- All three maintained files passed `load_config(path, check_paths=False)` and
  resolved-config YAML serialize/reload equality: `synd4_baseline.yaml`,
  `synd4_full_mechanisms.yaml`, `synd4_full_mechanisms_finalize_fast.yaml`.
- Parsed all six YAML blocks in `docs/configuration.md`; corrected its malformed
  rolling-average example and kept the PDF-switch explanation with PDF docs.
- Documentation whitespace/link checks; no runtime implementation changed.

These checks are **schema/semantic and serialization evidence**, not topology,
TPR/XTC identity, real-system execution or numerical parity. No pytest suite or
heavy workflows were run for this documentation/backlog update. The new #218
configuration controls are not implemented by this PR.
