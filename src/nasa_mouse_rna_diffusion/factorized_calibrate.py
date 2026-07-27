"""Train-only covariance calibration for factorized DDIM samples."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.metrics import (
    _condition_effect,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)

from .factorized_adapter import load_factorized_role
from .factorized_config import load_factorized_config
from .factorized_evaluate import _class_probe, _per_tissue_effects, _plot_pca


class CovarianceCalibrator:
    """Affine CORAL map fitted without expression pairing or condition labels."""

    def __init__(self, ridge_fraction: float) -> None:
        if float(ridge_fraction) <= 0:
            raise ValueError("ridge_fraction must be positive")
        self.ridge_fraction = float(ridge_fraction)
        self.synthetic_mean: np.ndarray | None = None
        self.real_mean: np.ndarray | None = None
        self.transform: np.ndarray | None = None

    @staticmethod
    def _matrix_power(
        covariance: np.ndarray, power: float, floor: float
    ) -> np.ndarray:
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        powered = np.power(np.maximum(eigenvalues, floor), power)
        return (eigenvectors * powered.reshape(1, -1)) @ eigenvectors.T

    def fit(
        self, real: np.ndarray, synthetic: np.ndarray
    ) -> "CovarianceCalibrator":
        real = np.asarray(real, dtype=np.float64)
        synthetic = np.asarray(synthetic, dtype=np.float64)
        if real.ndim != 2 or synthetic.ndim != 2 or real.shape[1] != synthetic.shape[1]:
            raise ValueError("Calibration matrices must have the same feature width")
        self.real_mean = real.mean(axis=0)
        self.synthetic_mean = synthetic.mean(axis=0)
        real_covariance = np.cov(real, rowvar=False)
        synthetic_covariance = np.cov(synthetic, rowvar=False)
        average_variance = max(
            float(np.mean(np.diag(real_covariance))),
            float(np.mean(np.diag(synthetic_covariance))),
            1e-12,
        )
        floor = self.ridge_fraction * average_variance
        identity = np.eye(real.shape[1], dtype=np.float64)
        real_regularized = real_covariance + floor * identity
        synthetic_regularized = synthetic_covariance + floor * identity
        inverse_synthetic = self._matrix_power(
            synthetic_regularized, -0.5, floor
        )
        real_square_root = self._matrix_power(real_regularized, 0.5, floor)
        self.transform = inverse_synthetic @ real_square_root
        return self

    def apply(self, synthetic: np.ndarray) -> np.ndarray:
        if self.transform is None or self.synthetic_mean is None or self.real_mean is None:
            raise RuntimeError("CovarianceCalibrator must be fit before apply")
        result = (
            (np.asarray(synthetic, dtype=np.float64) - self.synthetic_mean)
            @ self.transform
            + self.real_mean
        )
        return result.astype(np.float32)


def _aligned_expression(
    role: dict[str, object], source_rows: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame]:
    lookup = {
        int(source): index for index, source in enumerate(role["source_row"])
    }
    missing = sorted(set(map(int, source_rows)) - set(lookup))
    if missing:
        raise ValueError(f"Synthetic rows are absent from prepared role: {missing[:5]}")
    indices = np.asarray([lookup[int(source)] for source in source_rows], dtype=int)
    return (
        role["expression"][indices],
        role["samples"].iloc[indices].reset_index(drop=True),
    )


def calibrate_factorized(
    config_path: str | Path,
    *,
    guidance_scale: float = 1.0,
    ridge_fractions: Iterable[float] = (0.001, 0.01, 0.1),
) -> Path:
    config = load_factorized_config(config_path)
    data = config["data"]
    train = load_factorized_role(data["prepared_h5"], data["samples_tsv"], "train")
    validation = load_factorized_role(
        data["prepared_h5"], data["samples_tsv"], "validation"
    )
    evaluation = (
        Path(config["run"]["output_dir"])
        / "evaluation"
        / f"validation_guidance_{float(guidance_scale):g}"
    )
    train_npz = np.load(evaluation / "synthetic_train_expression.npz")
    validation_npz = np.load(evaluation / "synthetic_validation_expression.npz")
    synthetic_train = np.asarray(train_npz["scaled_expression"], dtype=np.float32)
    synthetic_validation = np.asarray(
        validation_npz["scaled_expression"], dtype=np.float32
    )
    real_train, train_samples = _aligned_expression(
        train, np.asarray(train_npz["source_row"], dtype=np.int64)
    )
    real_validation, validation_samples = _aligned_expression(
        validation, np.asarray(validation_npz["source_row"], dtype=np.int64)
    )
    output = evaluation / "covariance_calibration"
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for ridge in map(float, ridge_fractions):
        calibrator = CovarianceCalibrator(ridge).fit(real_train, synthetic_train)
        calibrated_train = calibrator.apply(synthetic_train)
        calibrated_validation = calibrator.apply(synthetic_validation)
        fidelity = generated_quality(
            real_validation,
            calibrated_validation,
            max_pr_samples=len(real_validation),
        )
        memorization = memorization_metrics(
            train["expression"],
            calibrated_validation,
            max_samples=max(len(real_validation), 50),
            seed=int(config["run"]["seed"]),
        )
        selection = fidelity_selection(fidelity, memorization)
        train_fidelity = generated_quality(
            real_train, calibrated_train, max_pr_samples=len(real_train)
        )
        train_memorization = memorization_metrics(
            real_train,
            calibrated_train,
            max_samples=max(len(real_train), 50),
            seed=int(config["run"]["seed"]) + 1,
        )
        train_selection = fidelity_selection(train_fidelity, train_memorization)
        effect = _condition_effect(
            real_validation,
            calibrated_validation,
            validation_samples["condition"].astype(str).to_numpy(),
        )
        effect_gate = conditional_effect_selection(effect)
        per_tissue, tissue_validation = _per_tissue_effects(
            real_validation,
            calibrated_validation,
            validation_samples,
            validation["genes"],
        )
        directory = output / f"ridge_{ridge:g}"
        directory.mkdir(parents=True, exist_ok=True)
        per_tissue.to_csv(
            directory / "per_tissue_condition_recovery.tsv", sep="\t", index=False
        )
        _plot_pca(real_validation, calibrated_validation, validation_samples, directory)
        np.savez_compressed(
            directory / "calibrated_expression.npz",
            validation_expression=calibrated_validation,
            train_expression=calibrated_train,
            validation_source_row=validation_npz["source_row"],
            train_source_row=train_npz["source_row"],
            genes=np.asarray(validation["genes"]),
        )
        summary = {
            "status": "complete",
            "method": "unpaired_train_only_CORAL",
            "ridge_fraction": ridge,
            "guidance_scale": float(guidance_scale),
            "locked_test_opened": False,
            "fit_profiles": int(len(real_train)),
            "validation_profiles": int(len(real_validation)),
            "fidelity": fidelity,
            "memorization": memorization,
            "model_selection": selection,
            "train_fidelity": train_fidelity,
            "train_memorization": train_memorization,
            "train_model_selection": train_selection,
            "pooled_condition_effect": effect,
            "pooled_condition_gate": effect_gate,
            "per_tissue_condition": {
                row["tissue"]: {
                    key: value for key, value in row.items() if key != "tissue"
                }
                for row in per_tissue.to_dict(orient="records")
            },
            "condition_probe": _class_probe(
                real_train,
                train_samples["condition"],
                calibrated_validation,
                validation_samples["condition"],
            ),
            "tissue_probe": _class_probe(
                real_train,
                train_samples["tissue"],
                calibrated_validation,
                validation_samples["tissue"],
            ),
            "limitations": (
                "Calibration matches global training covariance without sample pairing "
                "or condition labels; it is not part of the Lacan et al. architecture."
            ),
        }
        (directory / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        rows.append(
            {
                "ridge_fraction": ridge,
                "correlation": fidelity["correlation_matrix_agreement"],
                "precision": fidelity["precision"],
                "recall": fidelity["recall"],
                "f1": fidelity["f1"],
                "adversarial_accuracy": fidelity["adversarial_accuracy"],
                "frechet_ratio": fidelity["frechet_ratio_to_real_split_p95"],
                "condition_delta_correlation": effect["delta_correlation"],
                "condition_direction_agreement": effect["direction_agreement"],
                "validation_fidelity_gate": selection[
                    "eligible_for_model_selection"
                ],
                "train_fidelity_gate": train_selection[
                    "eligible_for_model_selection"
                ],
                "condition_gate": effect_gate["passed"],
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(output / "calibration_screen.tsv", sep="\t", index=False)
    return output / "calibration_screen.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--guidance-scale", type=float, default=1.0)
    parser.add_argument(
        "--ridge-fractions", nargs="+", type=float, default=[0.001, 0.01, 0.1]
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    calibrate_factorized(
        args.config,
        guidance_scale=args.guidance_scale,
        ridge_fractions=args.ridge_fractions,
    )


if __name__ == "__main__":
    main()
