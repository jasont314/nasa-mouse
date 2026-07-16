"""Evaluate a trained run on an accession-held-out validation or test split."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .adapters import load_adapter
from .adapters.base import ModelAdapter
from .conditioning import CategoryEncoder
from .config import load_config
from .metrics import evaluate_model
from .preprocessing import FittedPreprocessor
from .training_data import DataPartition, load_prepared_osdr, prepare_training_data


def _same_stats(first, second) -> bool:
    if first is None or second is None:
        return first is None and second is None
    return bool(
        np.array_equal(first.center, second.center)
        and np.array_equal(first.scale, second.scale)
    )


def _assert_preprocessing_matches(
    reconstructed: FittedPreprocessor,
    serialized: FittedPreprocessor,
) -> None:
    if asdict(reconstructed.spec) != asdict(serialized.spec):
        raise ValueError("Reconstructed preprocessing specification changed")
    for name in ("global_stats", "post_harmonization_stats", "final_stats"):
        if not _same_stats(getattr(reconstructed, name), getattr(serialized, name)):
            raise ValueError(
                f"Reconstructed {name} differs from the serialized training fit"
            )
    if set(reconstructed.study_stats) != set(serialized.study_stats):
        raise ValueError("Reconstructed study preprocessing keys changed")
    for study in reconstructed.study_stats:
        if not _same_stats(
            reconstructed.study_stats[study], serialized.study_stats[study]
        ):
            raise ValueError(
                f"Reconstructed preprocessing differs for study {study!r}"
            )
    if reconstructed.harmonizer is not None or serialized.harmonizer is not None:
        raise ValueError(
            "Dedicated harmonizer state cannot be proven equivalent without a "
            "saved prepared matrix"
        )


def _profile_keys(frame: pd.DataFrame) -> list[tuple[str, str]]:
    if not {"profile_id", "accession"}.issubset(frame.columns):
        raise ValueError("Saved observations require profile_id and accession")
    keys = list(
        zip(
            frame["profile_id"].astype(str),
            frame["accession"].astype(str),
        )
    )
    if len(set(keys)) != len(keys):
        raise ValueError("Profile/accession keys are not unique within a partition")
    return keys


def _restore_saved_conditioning(
    partitions: dict[str, DataPartition],
    run_dir: Path,
    adapter: ModelAdapter,
) -> dict[str, DataPartition]:
    encoder = CategoryEncoder.load(run_dir / "categorical_encoder.json")
    if encoder.covariates != adapter.covariates:
        raise ValueError("Saved encoder covariates differ from the trained model")
    if encoder.cardinalities != adapter.cardinalities:
        raise ValueError("Saved encoder cardinalities differ from the trained model")
    restored: dict[str, DataPartition] = {}
    for name, partition in partitions.items():
        obs_path = run_dir / f"{name}_obs.tsv.gz"
        if not obs_path.exists():
            raise FileNotFoundError(f"Saved observations not found: {obs_path}")
        saved_obs = pd.read_csv(obs_path, sep="\t")
        current_keys = _profile_keys(partition.obs)
        saved_keys = _profile_keys(saved_obs)
        positions = {key: index for index, key in enumerate(current_keys)}
        if set(positions) != set(saved_keys):
            raise ValueError(f"Reconstructed {name} profiles differ from training")
        order = np.asarray([positions[key] for key in saved_keys], dtype=int)
        restored[name] = DataPartition(
            name=partition.name,
            matrix=partition.matrix[order],
            obs=saved_obs.reset_index(drop=True),
            categories=encoder.transform(saved_obs),
            weights=partition.weights[order],
        )
    return restored


def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir)
    config = load_config(run_dir / "resolved_config.yaml")
    if args.split == "test" and config.validation.final_test_locked and not args.unlock_test:
        raise SystemExit(
            "The final test split is locked. Re-run with --unlock-test only after "
            "model and preprocessing choices are fixed."
        )
    reconstructed_preprocessor = None
    reconstructed = not (
        (run_dir / "prepared_data.h5").exists()
        or (run_dir / "prepared_osdr.h5").exists()
    )
    if not reconstructed:
        genes, partitions = load_prepared_osdr(run_dir)
    else:
        import yaml

        resolved = yaml.safe_load(
            (run_dir / "resolved_config.yaml").read_text(encoding="utf-8")
        ) or {}
        tissue = str(resolved.get("run", {}).get("tissue_override", "")) or None
        prepared = prepare_training_data(config, tissue=tissue)
        genes, partitions = prepared.genes, prepared.partitions
        reconstructed_preprocessor = prepared.preprocessor
    preprocessor = FittedPreprocessor.load(run_dir)
    if reconstructed_preprocessor is not None:
        _assert_preprocessing_matches(reconstructed_preprocessor, preprocessor)
    adapter = load_adapter(
        run_dir, device_spec=args.device or config.execution.device
    )
    if reconstructed:
        partitions = _restore_saved_conditioning(partitions, run_dir, adapter)
    if genes != adapter.genes:
        raise ValueError("Prepared-data genes differ from the trained model")
    return evaluate_model(
        adapter,
        partitions,
        preprocessor,
        split=args.split,
        output_dir=run_dir / "evaluation",
        seed=config.training.seed,
        max_samples=args.max_samples or config.validation.max_metric_samples,
        save_generated_matrix=(
            args.save_generated or config.execution.save_generated_matrix
        ),
        samples_per_covariate_profile=(
            config.generation.samples_per_covariate_profile
        ),
        synthetic_to_real_ratios=config.generation.synthetic_to_real_ratios,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--split", choices=["validation", "test"], default="validation")
    parser.add_argument("--unlock-test", action="store_true")
    parser.add_argument("--max-samples", type=int, default=0)
    parser.add_argument("--save-generated", action="store_true")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="")
    return parser.parse_args()


def main() -> None:
    path = run(parse_args())
    print(path)


if __name__ == "__main__":
    main()
