#!/usr/bin/env python3
"""Small readable CSV snapshots of saved Parquet inputs for figure previews.

These are *samples*, not replacements for the original Parquet measurements.
Only a bounded first batch is read; full tables stay in Parquet.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pyarrow.parquet as pq


def export_preview_tables(root: Path, *, sample_rows: int = 5) -> list[dict]:
    root = Path(root)
    if not 1 <= sample_rows <= 25:
        raise ValueError("sample_rows must be between 1 and 25")
    output = root / "table_previews"
    if output.exists():
        raise FileExistsError(f"Refusing to replace table previews in {output}")
    output.mkdir(parents=True)
    records = []
    for source in sorted(root.rglob("*.parquet")):
        if output in source.parents:
            continue
        relative = source.relative_to(root)
        reader = pq.ParquetFile(source)
        schema = reader.schema_arrow
        names = schema.names
        nrows = reader.metadata.num_rows
        sample = next(reader.iter_batches(batch_size=sample_rows), None)
        sample_table = (sample.to_pandas() if sample is not None else
                        __import__("pandas").DataFrame(columns=names))
        sample_path = output / relative.with_suffix(".sample.csv")
        columns_path = output / relative.with_suffix(".columns.csv")
        sample_path.parent.mkdir(parents=True, exist_ok=True)
        sample_table.head(sample_rows).to_csv(sample_path, index=False)
        with columns_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("column", "dtype"))
            writer.writerows((field.name, str(field.type)) for field in schema)
        records.append({
            "source_parquet": relative.as_posix(),
            "n_rows": nrows,
            "n_columns": len(names),
            "column_names": " | ".join(names),
            "sample_rows": min(nrows, sample_rows),
            "sample_csv": sample_path.relative_to(root).as_posix(),
            "columns_csv": columns_path.relative_to(root).as_posix(),
        })
    with (output / "table_inventory.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "source_parquet", "n_rows", "n_columns", "column_names",
            "sample_rows", "sample_csv", "columns_csv",
        ])
        writer.writeheader()
        writer.writerows(records)
    (output / "README.txt").write_text(
        "SYNTHETIC EXAMPLE DATA — NOT MD RESULTS\n"
        "table_inventory.csv: one record for every Parquet source with its size, "
        "columns and links to readable samples.\n"
        "*.columns.csv: all column names and Arrow types (no data omitted).\n"
        "*.sample.csv: at most the first %d rows; do not confuse the sample "
        "with a full numerical export.\n"
        "The original Parquet files remain the authoritative numerical tables.\n"
        % sample_rows
    )
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--sample-rows", type=int, default=5)
    args = parser.parse_args()
    records = export_preview_tables(args.root, sample_rows=args.sample_rows)
    print(json.dumps({"tables": len(records),
                      "inventory": str(args.root / "table_previews/table_inventory.csv")}))


if __name__ == "__main__":
    main()
