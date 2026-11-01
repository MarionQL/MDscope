# Solvation analysis for the original syndecan-4 CHARMM-GUI models

The repository's `synd4_full_mechanisms.yaml` uses old-build topology paths
on Tinkercliffs. **Those files are not present in this GitHub repository.**
Consequently, atom selections below are candidates inferred from the existing
configuration, not independently validated on the TPRs. Do not enable a full
multi-microsecond run until the preflight succeeds for **both conditions**.

## Minimal common-target configuration

Add this under `analyses:` in a copy of the existing YAML. It is deliberately
disabled, and no costly kinetic preset is enabled by default.

```yaml
  solvation:
    enabled: false
    preset: standard
    water_selection: "resname TIP3 and name OH2"
    ion_selections:
      K: "resname POT"
      Cl: "resname CLA"
    targets:
      protein:
        acidic_sidechain:
          selection: "protein and ((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))"
          group_by: residue
        basic_sidechain:
          selection: "protein and ((resname LYS and name NZ) or (resname ARG and name NE NH1 NH2))"
          group_by: residue
      hs:
        glycan_oxygen:
          selection: "(resname AGLCNA AIDOA BGLCA BGAL BXYL) and (chainID PROA PROB) and name O*"
          group_by: hs_chain
      lipid:
        phospholipid_phosphate:
          selection: "resname POPC POPE POPI POPS DPPC DOPC PSM and name O11 O12 O13 O14"
          group_by: residue
    search_radius_A: 8.0
    rdf_bin_width_A: 0.2
    shell_cutoffs_A: {water: 3.5, ions: 4.0}
```

The `hs_chain` grouping reuses MDscope's bond-based
`discover_covalent_hs` and, when present, the existing
`sdc4_mechanisms.hs_auto` settings. It identifies six attached HS chains
as separate glycan-only bonded components (rather than merging them through
the shared PROA/PROB protein backbone). This works only if glycan atom
selections cover the **complete glycan bond graph** and bond connectivity
is present in the TPR/PSF. The example `name O*` defines measured hydration
sites, not the glycan component discovery selection.

The LYS/ARG selections above are appropriate for *water proximity*, but
ion coordination requires **chemical interpretation by ion sign**.
Cl⁻ near the protein basic groups and K⁺ near acidic groups are not
interchangeable. Sulfate-specific versus generic glycan oxygen groups should
be defined using atom names audited from the actual HS topology.
Likewise, this example's phosphate oxygen names must be checked against
the actual CHARMM lipid topology.

GM1a and BGAL headgroups are missing in the nonglycolipid condition.
Define separate **condition-specific** solvation target selections for those
species or analyze them individually. Do not silently treat a missing
species as having zero hydration. GM1a's glycan residues can span multiple
topology residues; do not assume one residue represents an entire lipid.

## Preflight and real observer benchmark on Tinkercliffs

Run the tool from an editable MDscope installation or with its package on
`PYTHONPATH`. Use a new YAML file whose `solvation.enabled` is true.
You may temporarily switch off other analyses in that YAML to avoid
unrelated validation requirements.

```bash
python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --output solvation_benchmark_standard.json

python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset detailed --output solvation_benchmark_detailed.json
```

Repeat for `synd4_no_glycolipids`. Compare walltime, peak RSS,
per-frame processing rate, and output sizes. The benchmark invokes the
actual solvation observer (including chunk flush and finalization), but
not the rest of MDscope's complete analysis pipeline. Memory numbers
include the loaded Universe and dependencies.

**Scientific limitations:** RDF normalization uses bulk water density, not
the solvent-accessible volume around a site. Analyze saved SASA with water
counts and RDFs, and evaluate cutoff stability across trajectory chunks
and independent replicates. Unresolved RDF minima should not be accepted
as automatic shell cutoffs. With an XTC saved every 100 ps and a trajectory
stride of 10, adjacent analyzed frames would be 1 ns apart; water-contact
episodes then quantify sampled persistence, not picosecond water residence.
Check actual `ts.time` rather than assume the interval.

## Validation criteria

1. The preflight outputs the expected topology atom names, nonzero target
   counts and per-solvent counts for each condition; no unexpected overlaps.
2. Run with standard mode and inspect RDFs and available SASA for burial.
3. Check candidate minimum stability for each chemical class, with matched
   cutoff criteria between conditions.
4. Only enable detailed or kinetics modes on selected target sites once
   output-size and physical frame-interval limitations are understood.
5. Compare per-replicate results; do not use individual frames or hydration
   episodes as independent biological replicates.
