"""Repeated-seed train-only distribution calibration for factorized DDIM samples."""

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
from .factorized_evaluate import _per_tissue_effects, _plot_pca, _variant_label
from .factorized_mean_calibrate import HierarchicalMeanCalibrator


def _normalized_variant(value: str) -> str:
    value = str(value).strip()
    if value in {"", "base"}:
        return ""
    return _variant_label(value)


def _evaluation_directory(
    run_output: str | Path, guidance_scale: float, variant: str
) -> Path:
    name = f"validation_guidance_{float(guidance_scale):g}"
    normalized = _normalized_variant(variant)
    if normalized:
        name = f"{name}_{normalized}"
    return Path(run_output) / "evaluation" / name


class PositiveResidualCalibrator:
    """Mean alignment plus Gaussian covariance absent from synthetic training draws."""

    def __init__(
        self,
        group_columns: Iterable[str],
        prior_strength: float,
        residual_scale: float,
        *,
        clip_nonnegative: bool = True,
        noise_group_columns: Iterable[str] = (),
    ) -> None:
        if float(residual_scale) < 0:
            raise ValueError("residual_scale cannot be negative")
        self.mean_calibrator = HierarchicalMeanCalibrator(
            group_columns, prior_strength
        )
        self.residual_scale = float(residual_scale)
        self.clip_nonnegative = bool(clip_nonnegative)
        self.noise_group_columns = tuple(map(str, noise_group_columns))
        self.residual_root: np.ndarray | None = None
        self.positive_trace = 0.0
        self.negative_trace = 0.0

    def fit(
        self,
        real: np.ndarray,
        synthetic: np.ndarray,
        metadata: pd.DataFrame,
    ) -> "PositiveResidualCalibrator":
        real = np.asarray(real, dtype=np.float64)
        synthetic = np.asarray(synthetic, dtype=np.float64)
        self.mean_calibrator.fit(real, synthetic, metadata)
        calibrated = self.mean_calibrator.apply(synthetic, metadata).astype(
            np.float64
        )
        difference = np.cov(real, rowvar=False) - np.cov(
            calibrated, rowvar=False
        )
        difference = (difference + difference.T) / 2.0
        eigenvalues, eigenvectors = np.linalg.eigh(difference)
        tolerance = max(float(eigenvalues.max(initial=0.0)) * 1e-10, 1e-12)
        positive = eigenvalues > tolerance
        self.residual_root = (
            np.sqrt(eigenvalues[positive])[:, None]
            * eigenvectors[:, positive].T
        )
        self.positive_trace = float(eigenvalues[positive].sum())
        self.negative_trace = float(-eigenvalues[eigenvalues < -tolerance].sum())
        return self

    def apply(
        self,
        synthetic: np.ndarray,
        metadata: pd.DataFrame,
        *,
        seed: int,
    ) -> np.ndarray:
        if self.residual_root is None:
            raise RuntimeError("PositiveResidualCalibrator must be fit before apply")
        result = self.mean_calibrator.apply(synthetic, metadata).astype(np.float64)
        if self.residual_scale > 0 and len(self.residual_root):
            rng = np.random.default_rng(int(seed))
            if self.noise_group_columns:
                missing = [
                    column
                    for column in self.noise_group_columns
                    if column not in metadata
                ]
                if missing:
                    raise ValueError(
                        f"Residual-noise metadata lacks columns: {missing}"
                    )
                keys = (
                    metadata.loc[:, self.noise_group_columns]
                    .fillna("__missing__")
                    .astype(str)
                    .agg("||".join, axis=1)
                )
                random = np.empty(
                    (len(result), len(self.residual_root)), dtype=np.float64
                )
                for key in sorted(keys.unique()):
                    indices = np.flatnonzero(keys.eq(key).to_numpy())
                    pairs = len(indices) // 2
                    positive = rng.standard_normal(
                        (pairs, len(self.residual_root))
                    )
                    values = [positive, -positive]
                    if len(indices) % 2:
                        values.append(
                            rng.standard_normal((1, len(self.residual_root)))
                        )
                    random[rng.permutation(indices)] = np.concatenate(values)
            else:
                random = rng.standard_normal(
                    (len(result), len(self.residual_root))
                )
            result += (
                np.sqrt(self.residual_scale) * random @ self.residual_root
            )
        if self.clip_nonnegative:
            np.maximum(result, 0.0, out=result)
        return result.astype(np.float32)

    def save(self, directory: str | Path) -> Path:
        if self.residual_root is None:
            raise RuntimeError("Cannot save an unfitted calibrator")
        output = Path(directory)
        output.mkdir(parents=True, exist_ok=True)
        self.mean_calibrator.save(output)
        np.savez_compressed(
            output / "positive_residual_calibrator.npz",
            residual_root=self.residual_root,
        )
        manifest = {
            "method": "positive_missing_covariance_gaussian",
            "residual_scale": self.residual_scale,
            "clip_nonnegative": self.clip_nonnegative,
            "noise_group_columns": list(self.noise_group_columns),
            "residual_rank": int(len(self.residual_root)),
            "positive_trace": self.positive_trace,
            "negative_trace_discarded": self.negative_trace,
            "fit_condition_blind": True,
            "noise_balancing_uses_condition": (
                "condition" in self.noise_group_columns
            ),
        }
        path = output / "positive_residual_calibrator.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return path

    @classmethod
    def load(cls, directory: str | Path) -> "PositiveResidualCalibrator":
        root = Path(directory)
        manifest = json.loads(
            (root / "positive_residual_calibrator.json").read_text()
        )
        mean = HierarchicalMeanCalibrator.load(root)
        calibrator = cls(
            mean.group_columns,
            mean.prior_strength,
            manifest["residual_scale"],
            clip_nonnegative=manifest["clip_nonnegative"],
            noise_group_columns=manifest.get("noise_group_columns", []),
        )
        calibrator.mean_calibrator = mean
        calibrator.residual_root = np.load(
            root / "positive_residual_calibrator.npz"
        )["residual_root"]
        calibrator.positive_trace = float(manifest["positive_trace"])
        calibrator.negative_trace = float(manifest["negative_trace_discarded"])
        return calibrator


def _generated_role(
    role: dict[str, object], directory: Path, role_name: str
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    arrays = np.load(directory / f"synthetic_{role_name}_expression.npz")
    synthetic = np.asarray(arrays["scaled_expression"], dtype=np.float32)
    real, metadata = _aligned_expression(
        role, np.asarray(arrays["source_row"], dtype=np.int64)
    )
    return real, synthetic, metadata


def _metric_repeat_summary(rows: pd.DataFrame) -> dict[str, object]:
    requirements = {
        "correlation": ("minimum", "correlation_minimum"),
        "precision": ("minimum", 0.95),
        "recall": ("minimum", 0.85),
        "f1": ("minimum", 0.90),
        "adversarial_accuracy": ("range", (0.40, 0.60)),
        "frechet_ratio": ("maximum", 1.0),
    }
    result: dict[str, object] = {}
    for metric, (policy, requirement) in requirements.items():
        values = rows[metric].to_numpy(dtype=float)
        resolved = (
            float(rows[str(requirement)].iloc[0])
            if isinstance(requirement, str)
            else requirement
        )
        if policy == "minimum":
            passed = values >= float(resolved)
        elif policy == "maximum":
            passed = values <= float(resolved)
        else:
            lower, upper = resolved
            passed = (values >= lower) & (values <= upper)
        result[metric] = {
            "mean": float(np.mean(values)),
            "standard_deviation": float(np.std(values, ddof=1)),
            "minimum_observed": float(np.min(values)),
            "maximum_observed": float(np.max(values)),
            "replicates_passed": int(passed.sum()),
            "replicates": int(len(values)),
            "pass_fraction": float(passed.mean()),
            "requirement": resolved,
        }
    return result


def calibrate_factorized_distribution(
    config_path: str | Path,
    *,
    guidance_scale: float = 1.0,
    fit_variants: Iterable[str] = ("base", "seed3022", "seed3023"),
    evaluation_variants: Iterable[str] = (
        "base",
        "seed3021",
        "seed3022",
        "seed3023",
    ),
    group_columns: Iterable[str] = ("accession", "tissue"),
    prior_strength: float = 5.0,
    residual_scale: float = 0.5,
    residual_seed: int = 9100,
    noise_group_columns: Iterable[str] = ("accession", "tissue", "condition"),
    minimum_repeat_pass_fraction: float = 0.75,
) -> Path:
    if not 0 < float(minimum_repeat_pass_fraction) <= 1:
        raise ValueError("minimum_repeat_pass_fraction must be in (0, 1]")
    config = load_factorized_config(config_path)
    data = config["data"]
    run_output = Path(config["run"]["output_dir"])
    train = load_factorized_role(data["prepared_h5"], data["samples_tsv"], "train")
    validation = load_factorized_role(
        data["prepared_h5"], data["samples_tsv"], "validation"
    )
    fit_names = [_normalized_variant(value) for value in fit_variants]
    evaluation_names = [_normalized_variant(value) for value in evaluation_variants]
    if len(set(fit_names)) != len(fit_names) or len(set(evaluation_names)) != len(
        evaluation_names
    ):
        raise ValueError("Calibration variants must be unique")

    real_fit_parts: list[np.ndarray] = []
    synthetic_fit_parts: list[np.ndarray] = []
    metadata_fit_parts: list[pd.DataFrame] = []
    for variant in fit_names:
        directory = _evaluation_directory(run_output, guidance_scale, variant)
        real, synthetic, metadata = _generated_role(train, directory, "train")
        real_fit_parts.append(real)
        synthetic_fit_parts.append(synthetic)
        metadata_fit_parts.append(metadata)
    real_fit = np.concatenate(real_fit_parts)
    synthetic_fit = np.concatenate(synthetic_fit_parts)
    metadata_fit = pd.concat(metadata_fit_parts, ignore_index=True)

    calibrator = PositiveResidualCalibrator(
        group_columns,
        prior_strength,
        residual_scale,
        clip_nonnegative=True,
        noise_group_columns=noise_group_columns,
    ).fit(real_fit, synthetic_fit, metadata_fit)
    output = (
        run_output
        / "evaluation"
        / "repeated_distribution_calibration"
        / f"prior_{float(prior_strength):g}_residual_{float(residual_scale):g}"
    )
    output.mkdir(parents=True, exist_ok=True)
    calibrator.save(output)

    train_rows: list[dict[str, object]] = []
    train_summaries: dict[str, object] = {}
    for index, variant in enumerate(fit_names):
        source = _evaluation_directory(run_output, guidance_scale, variant)
        real, synthetic, metadata = _generated_role(train, source, "train")
        calibrated = calibrator.apply(
            synthetic, metadata, seed=int(residual_seed) + 10_000 + index
        )
        fidelity = generated_quality(real, calibrated, max_pr_samples=len(real))
        memorization = memorization_metrics(
            real,
            calibrated,
            max_samples=len(real),
            seed=int(config["run"]["seed"]) + 1,
        )
        selection = fidelity_selection(fidelity, memorization)
        label = variant or "base"
        train_summaries[label] = {
            "fidelity": fidelity,
            "memorization": memorization,
            "model_selection": selection,
        }
        train_rows.append(
            {
                "variant": label,
                "correlation": fidelity["correlation_matrix_agreement"],
                "correlation_minimum": selection["fidelity_gate"]["requirements"][
                    "correlation_matrix_agreement"
                ]["minimum"],
                "precision": fidelity["precision"],
                "recall": fidelity["recall"],
                "f1": fidelity["f1"],
                "adversarial_accuracy": fidelity["adversarial_accuracy"],
                "frechet_ratio": fidelity["frechet_ratio_to_real_split_p95"],
                "diversity_pass": selection["diversity_gate"]["passed"],
                "memorization_pass": selection["memorization_gate"]["passed"],
                "fidelity_pass": selection["eligible_for_model_selection"],
            }
        )
    train_table = pd.DataFrame(train_rows)
    train_table.to_csv(output / "train_repeat_metrics.tsv", sep="\t", index=False)
    train_metric_stability = _metric_repeat_summary(train_table)
    train_stable = all(
        float(value["pass_fraction"]) >= float(minimum_repeat_pass_fraction)
        for value in train_metric_stability.values()
    ) and bool(
        train_table["diversity_pass"].all()
        and train_table["memorization_pass"].all()
    )

    rows: list[dict[str, object]] = []
    repeat_summaries: dict[str, object] = {}
    for index, variant in enumerate(evaluation_names):
        source = _evaluation_directory(run_output, guidance_scale, variant)
        real, synthetic, metadata = _generated_role(
            validation, source, "validation"
        )
        calibrated = calibrator.apply(
            synthetic, metadata, seed=int(residual_seed) + index
        )
        fidelity = generated_quality(
            real, calibrated, max_pr_samples=len(real)
        )
        memorization = memorization_metrics(
            train["expression"],
            calibrated,
            max_samples=max(len(real), 50),
            seed=int(config["run"]["seed"]),
        )
        selection = fidelity_selection(fidelity, memorization)
        effect = _condition_effect(
            real, calibrated, metadata["condition"].astype(str).to_numpy()
        )
        effect_gate = conditional_effect_selection(effect)
        per_tissue, tissue_validation = _per_tissue_effects(
            real, calibrated, metadata, validation["genes"]
        )
        muscle = tissue_validation.get("skeletal_muscle", {})
        muscle_gate = muscle.get(
            "gate", accession_effect_selection({"accessions": 0})
        )
        label = variant or "base"
        directory = output / label
        directory.mkdir(parents=True, exist_ok=True)
        per_tissue.to_csv(
            directory / "per_tissue_condition_recovery.tsv", sep="\t", index=False
        )
        comparison = muscle.pop("comparison", None)
        if comparison is not None:
            comparison.to_csv(
                directory / "skeletal_muscle_accession_effect_recovery.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
        np.savez_compressed(
            directory / "calibrated_validation_expression.npz",
            scaled_expression=calibrated,
            source_row=np.load(
                source / "synthetic_validation_expression.npz"
            )["source_row"],
            genes=np.asarray(validation["genes"]),
        )
        _plot_pca(real, calibrated, metadata, directory)
        gates = {
            "validation_fidelity": bool(selection["eligible_for_model_selection"]),
            "train_fidelity": bool(train_stable),
            "pooled_condition_effect": bool(effect_gate["passed"]),
            "muscle_accession_effect": bool(muscle_gate["passed"]),
        }
        summary = {
            "status": "complete",
            "split": "validation",
            "locked_test_opened": False,
            "evaluation_variant": label,
            "fit_variants": [value or "base" for value in fit_names],
            "residual_seed": int(residual_seed) + index,
            "calibration_fit_condition_blind": True,
            "residual_noise_balancing_columns": list(
                calibrator.noise_group_columns
            ),
            "fidelity": fidelity,
            "memorization": memorization,
            "model_selection": selection,
            "pooled_condition_effect": effect,
            "pooled_condition_gate": effect_gate,
            "skeletal_muscle_accession_validation": muscle,
            "independent_acceptance_gates": gates,
            "all_acceptance_gates_pass": bool(all(gates.values())),
        }
        (directory / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        repeat_summaries[label] = summary
        rows.append(
            {
                "variant": label,
                "correlation": fidelity["correlation_matrix_agreement"],
                "correlation_minimum": selection["fidelity_gate"]["requirements"][
                    "correlation_matrix_agreement"
                ]["minimum"],
                "precision": fidelity["precision"],
                "recall": fidelity["recall"],
                "f1": fidelity["f1"],
                "adversarial_accuracy": fidelity["adversarial_accuracy"],
                "frechet_ratio": fidelity["frechet_ratio_to_real_split_p95"],
                "negative_fraction": float(np.mean(calibrated < 0)),
                "diversity_pass": selection["diversity_gate"]["passed"],
                "memorization_pass": selection["memorization_gate"]["passed"],
                "validation_fidelity_pass": selection[
                    "eligible_for_model_selection"
                ],
                "condition_delta_correlation": effect["delta_correlation"],
                "condition_direction_agreement": effect["direction_agreement"],
                "condition_effect_pass": effect_gate["passed"],
                "muscle_accession_correlation": muscle.get("summary", {}).get(
                    "meta_effect_correlation", float("nan")
                ),
                "muscle_accession_direction": muscle.get("summary", {}).get(
                    "meta_direction_agreement", float("nan")
                ),
                "muscle_accession_pass": muscle_gate["passed"],
                "all_gates_pass": bool(all(gates.values())),
            }
        )

    table = pd.DataFrame(rows)
    table_path = output / "repeat_metrics.tsv"
    table.to_csv(table_path, sep="\t", index=False)
    metric_stability = _metric_repeat_summary(table)
    required_fraction = float(minimum_repeat_pass_fraction)
    fidelity_stable = all(
        float(value["pass_fraction"]) >= required_fraction
        for value in metric_stability.values()
    ) and bool(table["diversity_pass"].all() and table["memorization_pass"].all())
    effect_stability = {
        "pooled_condition_effect": {
            "pass_fraction": float(table["condition_effect_pass"].mean()),
            "passed": bool(table["condition_effect_pass"].mean() >= required_fraction),
        },
        "muscle_accession_effect": {
            "pass_fraction": float(table["muscle_accession_pass"].mean()),
            "passed": bool(table["muscle_accession_pass"].mean() >= required_fraction),
        },
    }
    stable_gates = {
        "train_fidelity": bool(train_stable),
        "repeated_validation_fidelity": bool(fidelity_stable),
        "repeated_pooled_condition_effect": effect_stability[
            "pooled_condition_effect"
        ]["passed"],
        "repeated_muscle_accession_effect": effect_stability[
            "muscle_accession_effect"
        ]["passed"],
    }
    aggregate = {
        "status": "complete",
        "method": "condition_blind_multi_draw_positive_residual_calibration",
        "fit_split": "train",
        "evaluation_split": "validation",
        "locked_test_opened": False,
        "fit_variants": [value or "base" for value in fit_names],
        "evaluation_variants": [value or "base" for value in evaluation_names],
        "minimum_repeat_pass_fraction": required_fraction,
        "calibration_fit_condition_blind": True,
        "residual_noise_balancing_columns": list(calibrator.noise_group_columns),
        "metric_repeat_stability": metric_stability,
        "train_metric_repeat_stability": train_metric_stability,
        "effect_repeat_stability": effect_stability,
        "train_repeats": train_summaries,
        "independent_stability_gates": stable_gates,
        "all_stability_gates_pass": bool(all(stable_gates.values())),
        "acceptance_rule": (
            "Every metric is gated separately; each must pass in at least the "
            "declared fraction of generation repeats. No composite score is used."
        ),
        "limitation": (
            "Study conditioning and group calibration support interpolation only for "
            "studies represented during training. FLT/GC labels are not used to fit "
            "means or covariance; they only define antithetic noise-balancing strata."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(aggregate, indent=2) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# Repeated distribution calibration\n\n"
        "This train-only artifact aligns global and accession/tissue "
        "means, adds the positive covariance missing from generated training draws, "
        "balances residual noise within covariate/condition strata, and clips scaled "
        "expression to nonnegative values. FLT/GC labels do not fit any mean or "
        "covariance parameter. `repeat_metrics.tsv` "
        "reports every generation seed independently; `summary.json` applies no "
        "composite score. The locked test split was not opened.\n",
        encoding="utf-8",
    )
    return table_path
