"""Regression tests for selective interaction table loading."""
import pandas as pd

from mdscope.interactions import interaction_tables


def test_excluded_kind_is_not_read(tmp_path, monkeypatch):
    base = tmp_path / "interactions"
    for kind in ("hydrogen_bonds", "lipid_hydrogen_bonds"):
        folder = base / kind
        folder.mkdir(parents=True)
        for name in ("presence", "occupancy", "lifetimes", "summary"):
            (folder / f"{name}.parquet").touch()

    reads = []

    def fake_read(path):
        reads.append(path)
        return pd.DataFrame({"time_ps": [1.0]})

    monkeypatch.setattr(pd, "read_parquet", fake_read)
    tables = interaction_tables(
        tmp_path,
        exclude_kinds={"lipid_hydrogen_bonds"},
        table_names=("presence", "lifetimes"),
    )

    assert list(tables) == ["hydrogen_bonds"]
    assert set(tables["hydrogen_bonds"]) == {"presence", "lifetimes"}
    assert {p.name for p in reads} == {"presence.parquet", "lifetimes.parquet"}
    assert all("lipid_hydrogen_bonds" not in str(p) for p in reads)


def test_default_loads_all_existing_tables(tmp_path, monkeypatch):
    folder = tmp_path / "interactions" / "salt_bridges"
    folder.mkdir(parents=True)
    for name in ("presence", "occupancy", "lifetimes", "summary"):
        (folder / f"{name}.parquet").touch()
    monkeypatch.setattr(pd, "read_parquet", lambda path: pd.DataFrame())
    assert set(interaction_tables(tmp_path)["salt_bridges"]) == {
        "presence", "occupancy", "lifetimes", "summary"
    }
