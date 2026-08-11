"""Condition-blind hierarchical mean calibration for factorized DDIM samples."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.metrics import (
    _condition_effect,
    accession_effect_selection,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)

from .factorized_adapter import load_factorized_role
from .factorized_calibrate import _aligned_expression
from .factorized_config import load_factorized_config
from .factorized_evaluate import (
    _class_probe,
    _per_tissue_effects,
    _plot_pca,
    _variant_label,
)


class HierarchicalMeanCalibrator:
    """Global diagonal alignment plus shrunk condition-blind group shifts."""

    def __init__(
        self,
        group_columns: Iterable[str],
        prior_strength: float,
        *,
        epsilon: float = 1e-5,
    ) -> None:
        self.group_columns = tuple(map(str, group_columns))
        if not self.group_columns:
            raise ValueError("At least one calibration group column is required")
        if "condition" in self.group_columns:
            raise ValueError("Mean calibration cannot use the FLT/GC condition label")
        if float(prior_strength) <= 0:
            raise ValueError("prior_strength must be positive")
        self.prior_strength = float(prior_strength)
        self.epsilon = float(epsilon)
        self.real_mean: np.ndarray | None = None
        self.synthetic_mean: np.ndarray | None = None
        self.global_scale: np.ndarray | None = None
        self.group_shifts: dict[str, np.ndarray] = {}
        self.group_weights: dict[str, float] = {}
        self.group_counts: dict[str, int] = {}

    def _keys(self, metadata: pd.DataFrame) -> pd.Series:
        missing = [column for column in self.group_columns if column not in metadata]
        if missing:
            raise ValueError(f"Calibration metadata lacks columns: {missing}")
        return (
            metadata.loc[:, self.group_columns]
            .fillna("__missing__")
            .astype(str)
            .agg("||".join, axis=1)
        )

    def _global_apply(self, synthetic: np.ndarray) -> np.ndarray:
        if (
            self.real_mean is None
            or self.synthetic_mean is None
            or self.global_scale is None
        ):
            raise RuntimeError("HierarchicalMeanCalibrator must be fit before apply")
        return (
            (np.asarray(synthetic, dtype=np.float64) - self.synthetic_mean)
            * self.global_scale
            + self.real_mean
        )

    def fit(
        self,
        real: np.ndarray,
        synthetic: np.ndarray,
        metadata: pd.DataFrame,
    ) -> "HierarchicalMeanCalibrator":
        real = np.asarray(real, dtype=np.float64)
        synthetic = np.asarray(synthetic, dtype=np.float64)
        if real.shape != synthetic.shape or real.ndim != 2:
            raise ValueError("Real and synthetic calibration matrices must align")
        if len(metadata) != len(real):
            raise ValueError("Calibration metadata must align to expression rows")
        self.real_mean = real.mean(axis=0)
        self.synthetic_mean = synthetic.mean(axis=0)
        self.global_scale = (real.std(axis=0) + self.epsilon) / (
            synthetic.std(axis=0) + self.epsilon
        )
        globally_calibrated = self._global_apply(synthetic)
        keys = self._keys(metadata)
        self.group_shifts = {}
        self.group_weights = {}
        self.group_counts = {}
        for key in sorted(keys.unique()):
            mask = keys.eq(key).to_numpy()
            count = int(mask.sum())
            weight = count / (count + self.prior_strength)
            self.group_shifts[key] = (
                real[mask].mean(axis=0) - globally_calibrated[mask].mean(axis=0)
            )
            self.group_weights[key] = float(weight)
            self.group_counts[key] = count
        return self

    def apply(self, synthetic: np.ndarray, metadata: pd.DataFrame) -> np.ndarray:
        synthetic = np.asarray(synthetic)
        if len(metadata) != len(synthetic):
            raise ValueError("Calibration metadata must align to expression rows")
        calibrated = self._global_apply(synthetic)
        keys = self._keys(metadata)
        for key in sorted(set(keys) & set(self.group_shifts)):
            calibrated[keys.eq(key).to_numpy()] += (
                self.group_weights[key] * self.group_shifts[key]
            )
        return calibrated.astype(np.float32)

    def save(self, directory: str | Path) -> Path:
        if self.global_scale is None:
            raise RuntimeError("Cannot save an unfitted calibrator")
        output = Path(directory)
        output.mkdir(parents=True, exist_ok=True)
        keys = sorted(self.group_shifts)
        np.savez_compressed(
            output / "mean_calibrator.npz",
            real_mean=self.real_mean,
            synthetic_mean=self.synthetic_mean,
            global_scale=self.global_scale,
            group_keys=np.asarray(keys, dtype=str),
            group_shifts=np.asarray([self.group_shifts[key] for key in keys]),
            group_weights=np.asarray([self.group_weights[key] for key in keys]),
            group_counts=np.asarray([self.group_counts[key] for key in keys]),
        )
        manifest = {
            "method": "global_diagonal_plus_shrunk_group_mean",
            "group_columns": list(self.group_columns),
            "prior_strength": self.prior_strength,
            "epsilon": self.epsilon,
            "groups": len(keys),
            "condition_blind": True,
        }
        path = output / "mean_calibrator.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return path

    @classmethod
    def load(cls, directory: str | Path) -> "HierarchicalMeanCalibrator":
        root = Path(directory)
        manifest = json.loads((root / "mean_calibrator.json").read_text())
        arrays = np.load(root / "mean_calibrator.npz")
        calibrator = cls(
            manifest["group_columns"],
            manifest["prior_strength"],
            epsilon=manifest["epsilon"],
        )
        calibrator.real_mean = arrays["real_mean"]
        calibrator.synthetic_mean = arrays["synthetic_mean"]
        calibrator.global_scale = arrays["global_scale"]
        keys = arrays["group_keys"].astype(str).tolist()
        calibrator.group_shifts = dict(zip(keys, arrays["group_shifts"]))
        calibrator.group_weights = dict(zip(keys, arrays["group_weights"].tolist()))
        calibrator.group_counts = dict(zip(keys, arrays["group_counts"].tolist()))
        return calibrator


def calibrate_factorized_means(
    config_path: str | Path,
    *,
    guidance_scale: float = 1.0,
    group_columns: Iterable[str] = ("accession", "tissue"),
    prior_strengths: Iterable[float] = (2.0, 5.0, 10.0),
    evaluation_variant: str = "",
) -> Path:
    config = load_factorized_config(config_path)
    data = config["data"]
    train = load_factorized_role(data["prepared_h5"], data["samples_tsv"], "train")
    validation = load_factorized_role(
        data["prepared_h5"], data["samples_tsv"], "validation"
    )
    variant = _variant_label(evaluation_variant)
    evaluation_name = f"validation_guidance_{float(guidance_scale):g}"
    if variant:
        evaluation_name = f"{evaluation_name}_{variant}"
    evaluation = Path(config["run"]["output_dir"]) / "evaluation" / evaluation_name
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
    columns = tuple(map(str, group_columns))
    output = evaluation / "hierarchical_mean_calibration"
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for prior in map(float, prior_strengths):
        calibrator = HierarchicalMeanCalibrator(columns, prior).fit(
            real_train, synthetic_train, train_samples
        )
        calibrated_train = calibrator.apply(synthetic_train, train_samples)
        calibrated_validation = calibrator.apply(
            synthetic_validation, validation_samples
        )
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
        muscle = tissue_validation.get("skeletal_muscle", {})
        muscle_gate = muscle.get(
            "gate", accession_effect_selection({"accessions": 0})
        )
        gates = {
            "validation_fidelity": bool(selection["eligible_for_model_selection"]),
            "train_fidelity": bool(train_selection["eligible_for_model_selection"]),
            "pooled_condition_effect": bool(effect_gate["passed"]),
            "muscle_accession_effect": bool(muscle_gate["passed"]),
        }
        directory = output / f"prior_{prior:g}"
        directory.mkdir(parents=True, exist_ok=True)
        calibrator.save(directory)
        per_tissue.to_csv(
            directory / "per_tissue_condition_recovery.tsv", sep="\t", index=False
        )
        for tissue, result in tissue_validation.items():
            comparison = result.pop("comparison", None)
            if comparison is not None:
                comparison.to_csv(
                    directory / f"{tissue}_accession_effect_recovery.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
        np.savez_compressed(
            directory / "calibrated_expression.npz",
            validation_expression=calibrated_validation,
            train_expression=calibrated_train,
            validation_source_row=validation_npz["source_row"],
            train_source_row=train_npz["source_row"],
            genes=np.asarray(validation["genes"]),
        )
        _plot_pca(real_validation, calibrated_validation, validation_samples, directory)
        summary = {
            "status": "complete",
            "method": "condition_blind_hierarchical_mean_calibration",
            "fit_split": "train",
            "evaluation_split": "validation",
            "locked_test_opened": False,
            "guidance_scale": float(guidance_scale),
            "evaluation_variant": variant,
            "group_columns": list(columns),
            "prior_strength": prior,
            "fidelity": fidelity,
            "memorization": memorization,
            "model_selection": selection,
            "train_fidelity": train_fidelity,
            "train_memorization": train_memorization,
            "train_model_selection": train_selection,
            "pooled_condition_effect": effect,
            "pooled_condition_gate": effect_gate,
            "skeletal_muscle_accession_validation": muscle,
            "independent_acceptance_gates": gates,
            "all_acceptance_gates_pass": bool(all(gates.values())),
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
            "limitation": (
                "Accession conditioning and calibration support interpolation only for "
                "studies represented during training; condition labels are not used by "
                "the calibrator."
            ),
        }
        (directory / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        rows.append(
            {
                "prior_strength": prior,
                "correlation": fidelity["correlation_matrix_agreement"],
                "correlation_minimum": selection["fidelity_gate"]["requirements"][
                    "correlation_matrix_agreement"
                ]["minimum"],
                "precision": fidelity["precision"],
                "recall": fidelity["recall"],
                "f1": fidelity["f1"],
                "adversarial_accuracy": fidelity["adversarial_accuracy"],
                "frechet_ratio": fidelity["frechet_ratio_to_real_split_p95"],
                "condition_delta_correlation": effect["delta_correlation"],
                "condition_direction_agreement": effect["direction_agreement"],
                "muscle_accession_correlation": muscle.get("summary", {}).get(
                    "meta_effect_correlation", float("nan")
                ),
                "muscle_accession_direction": muscle.get("summary", {}).get(
                    "meta_direction_agreement", float("nan")
                ),
                **{f"gate_{key}": value for key, value in gates.items()},
                "all_gates_pass": bool(all(gates.values())),
            }
        )
    table = pd.DataFrame(rows)
    path = output / "mean_calibration_screen.tsv"
    table.to_csv(path, sep="\t", index=False)
    return path
