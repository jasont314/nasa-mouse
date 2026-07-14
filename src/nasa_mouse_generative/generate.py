"""Generate conditioned synthetic expression from a trained generator run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .adapters import load_adapter
from .conditioning import CategoryEncoder
from .config import load_config
from .preprocessing import FittedPreprocessor


def _overrides(values: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Expected COVARIATE=VALUE, got {value!r}")
        key, item = value.split("=", 1)
        result[key] = item
    return result


def _default_profile(run_dir: Path, covariates: tuple[str, ...]) -> dict[str, str]:
    table = pd.read_csv(run_dir / "train_obs.tsv.gz", sep="\t")
    columns = [covariate for covariate in covariates if covariate in table]
    if not columns:
        return {}
    combinations = (
        table[columns]
        .fillna("__missing__")
        .astype(str)
        .value_counts(dropna=False)
        .reset_index(name="profiles")
        .sort_values(["profiles", *columns], ascending=[False, *([True] * len(columns))])
    )
    row = combinations.iloc[0]
    return {covariate: str(row[covariate]) for covariate in columns}


def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir)
    config = load_config(run_dir / "resolved_config.yaml")
    adapter = load_adapter(
        run_dir, device_spec=args.device or config.execution.device
    )
    if not adapter.supports_generation:
        raise SystemExit(f"{adapter.adapter_id} is a representation model, not a generator")
    encoder = CategoryEncoder.load(run_dir / "categorical_encoder.json")
    preprocessor = FittedPreprocessor.load(run_dir)
    profile = _default_profile(run_dir, encoder.covariates)
    profile.update(_overrides(args.set))
    if args.condition:
        profile["condition"] = args.condition
    for covariate, value in profile.items():
        if covariate in encoder.vocabularies and value not in encoder.vocabularies[covariate]:
            raise ValueError(
                f"Unknown {covariate}={value!r}; choose a value recorded in "
                f"{run_dir / 'categorical_encoder.json'}"
            )
    profiles = [dict(profile) for _ in range(int(args.n))]
    categories = encoder.encode_profiles(profiles)
    transformed = adapter.generate(categories, seed=args.seed or config.training.seed)
    studies = [profile.get("study", "__unknown__")] * len(transformed)
    normalized = preprocessor.inverse_transform(transformed, studies)
    output_dir = Path(args.output_dir) if args.output_dir else run_dir / "generated"
    output_dir.mkdir(parents=True, exist_ok=True)
    condition = profile.get("condition", "unconditioned")
    stem = args.name or f"{condition}_{len(transformed)}"
    path = output_dir / f"{stem}.npz"
    np.savez_compressed(path, transformed=transformed, normalized=normalized)
    summary = {
        "run_dir": str(run_dir),
        "adapter_id": adapter.adapter_id,
        "samples": len(transformed),
        "genes": len(adapter.genes),
        "normalized_units": preprocessor.output_units,
        "conditioning_profile": profile,
        "outputs": {
            "matrix": str(path),
            "genes": str(run_dir / "genes.tsv"),
        },
    }
    summary_path = output_dir / f"{stem}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--name", default="")
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--condition", choices=["flight", "ground_control"], default="")
    parser.add_argument("--set", action="append", default=[], metavar="COVARIATE=VALUE")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
