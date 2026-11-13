"""Shared-trajectory observer tests including crash/replay semantics."""
import json

import pandas as pd

from mdscope.pipeline.fused_membrane import enabled_names
from mdscope_test_support import (
    membrane_universe as _universe,
    membrane_config as _cfg,
    execute_membrane as _execute,
)


def test_shared_anchor_and_checkpoint_replay(tmp_path):
    cfg = _cfg()
    assert enabled_names(cfg) == ("bilayer_normal", "lipid_enrichment")
    root_a, root_b = tmp_path / "clean", tmp_path / "resumed"
    result_a = _execute(_universe(), cfg, root_a)
    result_b = _execute(_universe(), cfg, root_b, resume=True)
    assert result_a["bilayer_normal"]["frames"] == 5
    assert result_b["lipid_enrichment"]["frames"] == 5
    for path in (
        "extensions/bilayer_normal/timeseries.parquet",
        "extensions/lipid_footprint/timeseries.parquet",
        "extensions/lipid_footprint/counts.parquet",
        "extensions/lipid_footprint/selection_audit.parquet",
    ):
        pd.testing.assert_frame_equal(
            pd.read_parquet(root_a / path), pd.read_parquet(root_b / path),
            check_dtype=False,
        )
    for root in (root_a, root_b):
        timing = json.loads((root / "fused_membrane_timings.json").read_text())
        assert timing["shared_anchor_mappings"] == 1
        angles = sorted((root / "extensions/bilayer_normal/per_lipid").glob("*.parquet"))
        assert len(angles) == 3
        assert sum(len(pd.read_parquet(path)) for path in angles) == 5 * 8
