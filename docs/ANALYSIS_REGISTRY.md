# Analysis instances and persisted results

MDscope separates **analysis algorithms**, **named instances**, and the **trajectory pass**.
Multiple named instances of an enabled analysis must share the same trajectory traversal,
but retain independent selections, accumulators, output paths, and completion records.

## Current support

- `density_1d` and `residue_membrane_distance` support named `runs` in YAML.
- Other observables currently have a single configuration. Do not assume adding
  `runs` to YAML enables multi-instance execution until that observer supports it.
- Fused observers share a trajectory reader; specialized analyses execute independently
  when their backend or coordinate requirements demand it.

## Adding or extending an analysis

1. Implement the numerical observer in its owning module and its enable logic.
2. If supporting named instances, use `analysis_runs(cfg, name)`, initialize each
   instance once, and update each instance within the **same frame loop**.
3. Give every instance a unique run ID and output directory. Register the expected
   finalized artifact in `FUSED_OUTPUTS`/`fused_output`; do not invent legacy aliases.
4. Validate finalized artifacts before recording completion. Resume/checkpoint keys
   must distinguish instances; a missing output is incomplete even if a marker exists.
5. Update existing contract and streaming/fused-only integration tests where applicable.

`FUSED_OUTPUTS` is a path catalog, **not** a numerical plugin engine or pass scheduler.
There is no requirement to preserve historical internal checkpoint keys. A changed key
can cause work to be repeated; never silently reuse old incomplete results.

Do not create a generic plugin framework merely to register analysis instances.
