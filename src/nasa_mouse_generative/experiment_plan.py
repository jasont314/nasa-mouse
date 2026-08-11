"""Expand staged benchmark YAML into a concrete, auditable run table."""

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import pandas as pd
import yaml

from .models import MODEL_REGISTRY


def expand_matrix(payload: dict) -> pd.DataFrame:
    rows = []
    for phase in payload.get("phases", []):
        axes = phase.get("axes", {})
        keys = list(axes)
        values = [value if isinstance(value, list) else [value] for value in axes.values()]
        for combination in product(*values):
            row = {"phase": phase["name"], "purpose": phase.get("purpose", "")}
            row.update(dict(zip(keys, combination)))
            model = MODEL_REGISTRY.get(str(row.get("model", "")))
            row["native_expression_generator"] = bool(
                model and model.supports_expression_generation
            )
            row["status"] = "planned"
            if row.get("task") == "conditional_generation" and not row[
                "native_expression_generator"
            ]:
                row["status"] = "capability_blocked"
            rows.append(row)
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> Path:
    payload = yaml.safe_load(Path(args.matrix).read_text(encoding="utf-8")) or {}
    table = expand_matrix(payload)
    if not table.empty and "status" in table.columns:
        ordered = [column for column in table.columns if column != "status"]
        table = table[[*ordered, "status"]]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    print(f"wrote {len(table)} planned rows to {output}")
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default="configs/generative/benchmark/experiment_matrix.yaml")
    parser.add_argument(
        "--output", default="outputs/generative/benchmark/summary/experiment_plan.tsv"
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
