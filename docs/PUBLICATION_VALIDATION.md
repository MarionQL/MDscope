# Scientific numerical validation ledger

These deterministic tests are **independent numerical references**, not merely
round-trip execution checks. Run with:

```bash
pytest -q tests/test_publication_references.py
```

The independent reference suite covers direct, finite-sample-normalized FFT
autocorrelation; a zero-crossing automatic block-length case; exact independent
label permutations; paired sign-flip enumeration; Benjamini–Hochberg correction
with missing values and family separation; constant-data block and paired
hierarchical bootstrap controls; and direct, quadratic multi-origin lateral
MSD. Inputs are synthetic and contain no private simulation data.

Bootstrap uses NumPy Generator streams with a fixed seed; repeatability within
a fixed software environment is tested. Numerical results, rather than a
bitwise sequence of draws, should be compared across NumPy/Python versions
unless an RNG bit-generator and version are explicitly pinned.

Independent simulations (replicates), not autocorrelated frames, determine
between-condition test units. Paired tests preserve matching replicate labels.
Units: MSD coordinates are unwrapped Å; lateral MSD is Å²; ACF and corrected
p-values are dimensionless; block lengths are measured in sampled frames.

Outstanding before a publication-quality release: a redistributable YAML
end-to-end reference dataset with audited output tables and tolerances, mixed
observer interrupted-checkpoint restart comparison, profile-axis alignment
and missing-observation reference cases, representative real-GROMACS numerical
comparisons, benchmark resource budget, and final documented license/citation.
Passing these unit tests does not certify the remaining requirements.
