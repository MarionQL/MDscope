# Legacy module cleanup inventory

The initial cleanup is intentionally conservative: thin re-export modules are
public import paths and **cannot be proven unused merely because the canonical
implementation moved**.

| Compatibility module | Canonical implementation | Current decision |
| --- | --- | --- |
| `analysis/alignment.py` | `analysis/core/alignment.py` | Retain until production/import consumers are migrated and external compatibility break is approved |
| `analysis/lipid_geometry.py` | `analysis/core/lipid_geometry.py` | Retain pending import audit; removal changes historical public API |
| `analysis/msd.py` | `analysis/observables/msd.py` | Retain pending callers audit |
| `analysis/bilayer_normal.py` | `analysis/observables/bilayer_normal.py` | Retain pending callers audit |
| `analysis/lipid_enrichment.py` | `analysis/observables/lipid_enrichment.py` | Retain pending callers audit |
| `analysis/profiling.py` | `analysis/core/profiling.py` | Retain pending callers audit |

## Cleanup completed

Deleted `tests/test_core_contract.py` after moving its two substantive
numerical regression tests into `tests/test_core_geometry_contracts.py`.
The removed historical-import identity assertions tested only shim forwarding,
not scientific behavior. Canonical core imports are now exercised directly.

## Follow-up criteria for production-file removal

Scan imports in `analysis/`, all scripts, examples, tests and documented
Python API usage; check dynamic/importlib imports and external HPC launchers.
Delete only explicitly deprecated modules with zero repository consumers,
document the import migration, and run both Python matrices and real-trajectory
smoke checks. Keep compatibility adapters where external callers still rely on
the old API. No production module has been removed in this first audit.


## Internal consumer migration (31 October 2026)

Moved imports in ten production files to canonical `core` and
`observables` implementations:

- Runner, streaming, membrane and trajectory-source MSD imports
- Streaming/timeseries alignment and profiling imports
- Fused membrane, registry and mechanistic lipid geometry/observables
- Bilayer normal and lipid enrichment observables' geometry imports

These are import-only changes. Legacy shim files are **intentionally retained**
to preserve downstream user scripts, historical imports and external HPC
workflows; zero in-repository imports is not proof of zero external consumers.

The two rigid-alignment and periodic-center numerical tests were relocated
before deleting the redundant historical import-identity test file.

No maintained CLI, Slurm launcher or public module was moved or removed.
A future breaking release can retire shims only after a deprecation period
and documented import migration. Scientific results should remain unchanged.
