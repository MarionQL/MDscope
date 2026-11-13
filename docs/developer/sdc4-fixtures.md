# SDC4 scientific reference fixtures

These are **synthetic geometry and restart tests**, not force-field simulations or convergence evidence.

## Reference geometry

`tests/test_sdc4_geometry_references.py`: **45 atoms**, cubic **100 Å** box; chain B precedes A in topology order and repeats some lipid/HS residue numbers.

| Part | Coordinates / selection |
|---|---|
| Protein A/B | Ectodomain CA: z=36,32,28,24,20 Å; TM CA: z=12,8,4,0,−4,−8,−12 Å |
| TM sequence | GGIVGAA; x_A=0, x_B=3 Å (contact) or 8 Å (separated) |
| Lipid anchors | Nine POPC P atoms per leaflet; x,y ∈ {−10,0,10} Å; z=±16 Å |
| HS_A | z=26,30 Å; protein attachment at z=24 Å |
| Water | One oxygen at z=29 Å |

### Independent expectations

| Observable | Reference | Boundary |
|---|---|---|
| TM crossing | 0° | Collapsed/tied PCA axes → unavailable |
| TM residue contact pairs ≤5 Å | 19 close; 0 separated | Cutoff equality included |
| GGIVG–GGIVG | 13 close; 0 separated | Missing optional motif → NaN, **not** 0 |
| GGIVG to complete opposite TM | 14 close; 0 separated | Multiple motif matches → error |
| Ectodomain height / tilt | 12 Å / 0° | Fewer than 4 CA sites → error |
| Local thickness / flat curvature | 32 Å / 0 Å⁻¹ | Sparse/degenerate surface → NaN |
| Quadratic curvature | 0.03 Å⁻¹ | Small-slope proxy |
| HS extension / two-site Rg | 6 Å / 2 Å | Collapsed chain: zero extension, unavailable orientation |
| Baseline HS/attachment lipid contacts | 0 | Lowered variant: 2 HS atom pairs, 1 attachment pair |
| Four-residue NxIP SASA | 10 nm² (=1+2+3+4) | Missing/nonfinite site → unavailable |

- Tilted A ectodomain: **45°**, same height. Tilted B TM axis: **60°** crossing and **−1** handedness under defined reference axes.
- Shift HS/attachment down 10 Å: one contacting lipid, membrane-contact state 1; extension remains 6 Å.
- Lattice wrapping, translation and 90° rotation about z preserve applicable scalars. Tilting a membrane relative to the *fixed z reference* does **not**.
- Validate both direct observer rows and saved Parquet outputs. Unsupported SDC4 triclinic local geometry must fail explicitly.

## Restart and finalization

`tests/test_mixed_observer_restart.py`: eight synthetic PDB/DCD frames sampled at stride 2, comparing clean chunk-size-3 execution against interrupted chunk-size-1/2 runs.

| Injected failure | Must preserve |
|---|---|
| During a frame | Committed data; replay uncommitted work |
| After writing a part | Detect/prune only uncommitted parts |
| Before checkpoint replacement | Consistent checkpoint/part boundary |
| Missing/corrupt committed part; copied checkpoint | Fail without overwriting existing outputs |
| Offline finalize after input deletion | Never reopen trajectory; retain replica identity and tables |
| Unequal replica lengths / absent species | Equal replica weighting, missing values not zero |

Tests compare RMSD/RMSF/Rg, density, normals, lipid contacts, SDC4, solvation and multiple named profile runs. Final Parquet agreement uses rtol/atol 1e−9; map/matrix arrays are exact. Solvation `shell_stability.parquet` is *chunk-dependent* and not expected to match under changed chunk boundaries.

### Checkpoint contract

- `cache/stream_checkpoint.npz` records inventory version, (system, replica) ID, sorted relative file paths, size and mtime-ns fingerprints.
- Verify before pruning or constructing restart observers. Remove only owned uncommitted chunks; never unrelated data.
- Incomplete legacy checkpoints lack proof of committed specialized outputs: start fresh. **Completed** legacy results can still be finalized offline.
- These fingerprints are metadata, **not** content hashes; matching size+mtime does not certify bytes.

## Analytical references

- Real FreeSASA heterogeneous spheres/caps: for R=3.1 Å, d=4 Å and r=element radius + 1.4 Å, `h = R − (d²+R²−r²)/(2d)`; exposed area `(4πR²−2πRh)/100` nm². Tolerance **0.006 nm²** at 20,000 surface points.
- `test_publication_references.py`: exactly **16** two-stage temporal/replica outcomes, scripted paired-label permutations, BH/permutation and numerical RNG references; unpinned cross-version RNG bitstreams are not promised.
- Exact two-sided independent 3+3 permutation: smallest p-value **0.10**.
- Additional independent chemistry, lipid PBC/area and residence oracles: `test_scientific_reference_oracles.py`, `test_lipid_contact_matrices.py`.

**Production gate still open:** all six ARC trajectories need chain/motif identities, chemical selections, matched windows, cutoffs, convergence and replicate agreement checked. [Release checklist](release-validation.md).
