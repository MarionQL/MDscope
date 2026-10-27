"""Labels shared across structural profile observables."""
from __future__ import annotations
def _residue_label(res) -> str:
    chain = getattr(res, "chainID", "") or getattr(res, "segid", "") or "?"
    return f"{chain}:{res.resname}{res.resid}"
