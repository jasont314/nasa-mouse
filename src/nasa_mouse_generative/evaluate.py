"""Evaluate a trained run on an accession-held-out validation or test split."""

from __future__ import annotations

import argparse
from pathlib import Path

from .adapters import load_adapter
from .config import load_config
from .metrics import evaluate_model
from .preprocessing import FittedPreprocessor
from .training_data import load_prepared_osdr


def run(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir)
    config = load_config(run_dir / "resolved_config.yaml")
    if args.split == "test" and config.validation.final_test_locked and not args.unlock_test:
        raise SystemExit(
            "The final test split is locked. Re-run with --unlock-test only after "
            "model and preprocessing choices are fixed."
        )
    genes, partitions = load_prepared_osdr(run_dir)
    preprocessor = FittedPreprocessor.load(run_dir)
    adapter = load_adapter(
        run_dir, device_spec=args.device or config.execution.device
    )
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
