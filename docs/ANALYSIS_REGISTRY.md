# Analysis registration and persisted output contracts

MDscope deliberately has different execution families: the fast streaming
reader, fused observers, specialized analyses, and postprocessing. The shared
catalog in `analysis/registry.py` describes persisted fused outputs and legacy
checkpoint marker names; it **does not** replace optimized fused execution.

To add a fused observable:

1. Implement the observer in the appropriate fused module and add its enable
   predicate to that module's `enabled_names` / `enabled` function.
2. Register its exact final output relative path, legacy checkpoint name and
   execution family in `FUSED_OUTPUTS`. For an analysis supporting multiple
   independently named runs, extend `fused_output` with explicit run IDs.
3. Ensure the observer finalizes that output and the runner checks its
   existence through `fused_output_plan`. A missing final output must not
   count as a successful completed analysis.
4. Add a row to `tests/test_fused_output_catalog.py` and an end-to-end
   test for both streaming and fused-only execution if applicable.

Existing exceptions are intentional: serial DSSP/SASA may be fused; alternate
DSSP backends remain specialized. The `profiles` completion marker is shared
for historical compatibility. Residue-distance run IDs remain distinct.

The catalog currently represents *fused persisted outputs*. Specialized tasks
continue using `AnalysisTask` and their own numerical outputs. Replacing
independent analysis algorithms with a generic plugin framework is explicitly
out of scope.
