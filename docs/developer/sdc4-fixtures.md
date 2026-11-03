# SDC4 scientific reference fixtures

`tests/test_sdc4_geometry_references.py` defines a 45-atom synthetic system
(one selected atom per residue), without force-field files or production data.
The constructor is also reused by the mixed-observer restart test; expected
values are independent of MDscope geometry helpers.

The box is 100 Å on each side. Chain B precedes A in topology order. B uses
residue numbers 301–312 and A uses 71–82; lipid/HS residue numbers repeat.
Each protein has five ectodomain CA sites at z = 36, 32, 28, 24, 20 Å and seven
TM CA sites at z = 12, 8, 4, 0, −4, −8, −12 Å. The TM sequence is GGIVGAA.
A is at x=0 and B at x=3 Å (contact case) or x=8 Å (separated case); y=0.
Each leaflet has nine POPC P anchors on x,y ∈ {−10,0,10} Å at z=±16 Å.
HS_A has two sites at z=26,30 Å, attached to the A site at z=24 Å.
One water oxygen is at z=29 Å. All HS/water x,y coordinates are zero.

| Observable | Independent expected answer | Tolerance / boundary |
|---|---|---|
| Parallel TM crossing | 0 degrees | 1e−6 degrees; collapsed or tied PCA directions unavailable |
| TM residue-pair contacts at ≤5 Å | 19 at x separation 3 Å; 0 at 8 Å | Exact counts: seven same-height pairs plus twelve neighboring-height pairs |
| GGIVG–GGIVG atom/residue pairs | 13 close; 0 far | Exact: five same-height plus eight neighboring-height pairs |
| Each GGIVG to other complete TM | 14 close; 0 far | Exact: five same-height, four preceding and five following pairs |
| Ectodomain height | 28−16 = 12 Å | 1e−5 Å, relative to the outward leaflet |
| Ectodomain tilt | 0 degrees | 1e−6 degrees, relative to the z membrane normal |
| Local thickness | 16−(−16) = 32 Å | 1e−5 Å; insufficient anchors give NaN |
| Flat leaflet curvature | 0 Å⁻¹ | 1e−8 Å⁻¹; rank-deficient samples give NaN |
| Quadratic leaflet curvature | a+b = 0.03 Å⁻¹ | 1e−9 Å⁻¹ in `test_sdc4_mechanisms.py`; this is the small-slope convention |
| Attachment-to-distal HS extension | 30−24 = 6 Å | 1e−5 Å; orientation cosine 1 |
| Equal-weight two-site HS Rg | 2 Å | 1e−6 Å |
| HS / attachment lipid contacts | 0 | Nearest upper anchor is 10 / 8 Å away |
| Complete NxIP SASA | 1+2+3+4 = 10 nm² | Exact synthetic persisted areas; missing/nonfinite residue makes total unavailable |

Additional persisted references cover A ectodomain tilt 45 degrees at the same
12 Å height and a B TM axis tilted 60 degrees in the yz plane. Relative to
the downward A axis and +x interhelix separation, the crossing angle is
60 degrees and handedness is −1 (angle tolerance 1e−5 degrees). Shifting the
attachment and HS downward by 10 Å gives attachment z=14 and HS z=16/20 Å:
one attachment–lipid pair, two HS–lipid atom pairs, one contacted lipid and a
membrane-contact state of 1, with extension still 6 Å. The original fixture
provides the corresponding no-contact boundary.

Collapsed TM axes persist unavailable angle/handedness; a collapsed HS has
measured zero extension but unavailable orientation/shape anisotropy. Missing
optional HS mappings omit their columns and retain empty audit counts, while
required empty mappings raise. Fewer than four ectodomain CA sites and
triclinic height-field geometry fail explicitly. Together with sparse/rank-
deficient surface tests, mapping ambiguity/missing-motif tests, and the
composition/chemical-molecule references below, these provide positive and
boundary/failure cases for every major observable named in #71.

Every contact/height/shape reference is checked both in the observer's direct
row and after Parquet finalization. The contact cutoff includes equality.
Variants independently shift atoms by box vectors, translate the entire
system across z=0, or rotate it 90 degrees about z. These operations preserve
the listed scalar results. An arbitrary tilt of the bilayer is **not** an
invariance for these fixed-z observables. The SDC4 observer rejects triclinic
boxes; its use of minimum-image distances does not make its xy maps general
triclinic or curved-surface measurements. Centers/axes require compact groups.

Missing optional GGIVG mappings now produce NaN in affected contact columns
and `unavailable_missing_motif` in `selection_audit.json`. A validated chain
with no contacts still produces zero. Multiple GGIVG matches raise an error
instead of being combined. This is an intentional scientific correction;
regenerate affected per-replicate products in a fresh output directory before
recomputing statistics. Previously saved legacy products are not rewritten.

Existing independent references in `test_scientific_reference_oracles.py`
cover composition denominators, projected area exclusion, chemical atom
subsets, stride and censoring. `test_lipid_contact_matrices.py` covers bonded
multi-residue glycolipids and finite-population mixing. These serve different
contracts from the single-site contact oracle here.

## Mixed-observer restart and offline finalization

`test_mixed_observer_restart.py` writes eight PDB/DCD frames, samples stride 2,
and compares a clean chunk-size-3 run against chunk sizes 1 and 2 interrupted
during a frame, after a raw chunk write, or after every observer flush but
before checkpoint replacement. RMSD/RMSF/Rg, chain-specific groups,
bilayer normals, lipid enrichment, two named density runs, two named
residue–membrane distance runs, solvation kinetics and SDC4 share the pass. Final
Parquet tables, schemas and statistics agree within rtol=atol=1e−9; lipid map
and residue–membrane distance arrays agree exactly. The solvation
`shell_stability.parquet` diagnostic is intentionally chunk-dependent: it
estimates cutoff stability over processing windows. Its chunk count is checked
against the configured chunk size; it is not asserted equal across window sizes.
Underlying counts, RDFs, episodes and statistics retain the parity requirement.
DCD times agree with [0,2,4,6] ps within 1e−6 because
the format stores a floating-point time conversion. Replica IDs are asserted.
Inputs are deleted and Universe construction forbidden before consolidation
and project statistics finalization. Invalid saved time axes fail before any
previous project summary is archived. Unequal replica lengths preserve equal
replica weighting, and unavailable species stay NaN in persisted raw summaries.
This does not substitute for six-replica production QA.

### Committed-artifact and replica integrity

The existing atomic `cache/stream_checkpoint.npz` now includes an inventory
version, `(system, replicate label)` identity, and relative file paths with
size/mtime-nanosecond fingerprints. Entries are added for the actual files
written at each committed chunk, including sparse contact-event files,
aligned coordinates, named density/distance runs, fused extensions, solvation,
lipid matrices and SDC4 maps. The checkpoint and scientific output formats
remain NPZ/Parquet; there is no additional sidecar or dependency.

Resume validates identity and fingerprints before constructing observers or
pruning uncommitted files. Missing/modified artifacts and inserted duplicate
committed chunks raise an actionable error and preserve existing file bytes.
An incomplete or legacy checkpoint without this inventory is rejected rather
than silently claiming complete observer coverage. Integer counter types and
agreement between processed frames and global/group/density accumulator
counts are also checked before pruning. Once validation succeeds, all owned
numbered uncommitted artifacts, including density frames and temporary map
files, are discarded for replay; unrelated raw files are preserved.

This is metadata freshness verification, **not a content checksum**. It cannot
detect bit rot or edits preserving both file size and mtime. Recording and
checking inventories enumerate filenames and read filesystem metadata, without
scanning large scientific tables or rereading trajectories. Inventory size
scales with the number of chunk files, not atom/event count. Serialized
inventories are sorted; numerical estimators and RNG draws are unchanged.

Migration: completed legacy results still support offline finalization and
reporting. For an interrupted legacy run, preserve its directory and start
fresh output with the current version; a checkpoint lacking the inventory
cannot establish which specialized files were committed. Moving/copying a
checkpoint to another replica label is rejected even with a shared project
configuration hash. No existing ARC job or saved product is modified here.
Offline finalization also checks checkpoint identity when available and checks
that summary row count agrees with its processed-frame counter, before archiving
an earlier project summary. Cached trajectories are never reopened. Legacy
summaries without a checkpoint remain readable; this does not invent missing
replica provenance for them.

## SASA and statistical references

Real FreeSASA tests additionally use unequal carbon/oxygen and carbon/phosphorus
probe-inflated sphere intersections. With R=3.1 Å, r=(element radius+1.4) Å
and d=4 Å, the occluded cap height is
`h = R − (d² + R² − r²)/(2d)`. Target SASA is
`(4πR² − 2πRh)/100` nm², tolerance 0.006 nm² at 20,000 Shrake–Rupley points.
Only the target residue is emitted; independent occluder wrapping preserves
the result and an invalid periodic box raises. Existing tests retain sphere,
equal-radius cap, burial sign, identity joins and fake-library streaming parity.

`test_publication_references.py` enumerates all 16 equally likely two-stage
replicate/temporal draws of [0,4] and [10,14], rather than only constant inputs.
A separate nonconstant paired oracle uses opposite dictionary insertion order
and explicitly scripted pair/draw indices. These exact arithmetic checks do
not depend on a particular NumPy random sequence. Seed repeatability remains
tested within the installed version; bitwise equivalence across library
versions is not promised. Existing exact 3-vs-3 tests document p≥0.1 for a
two-sided independent permutation test, regardless of the number of frames.

## Remaining production acceptance

This fixture advances #36, #69–#73, #98–#99 and #100/#101 without claiming that
the six original ARC trajectories were rerun. Their topology inventories,
corrected GGIVG outputs, matched windows, convergence, chemical/null sensitivity
and final report review remain required by the corresponding science issues.
Do not close those issues solely because these synthetic tests pass.
