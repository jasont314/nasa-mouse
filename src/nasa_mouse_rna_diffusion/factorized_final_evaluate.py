"""One-time locked-test evaluation for a fixed factorized DDIM finalist."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.metrics import (
    _condition_effect,
    accession_effect_selection,
    classifier_utility,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)

from .factorized_adapter import encode_factorized_labels, load_factorized_role
from .factorized_calibrate import _aligned_expression
from .factorized_config import load_factorized_config
from .factorized_distribution_calibrate import (
    PositiveResidualCalibrator,
    _metric_repeat_summary,
)
from .factorized_evaluate import (
    _load_adapter_model,
    _per_tissue_effects,
    _plot_pca,
    _sample,
)
from .upstream import quadratic_beta_schedule


def _require_test_unlock(unlock_test: bool) -> None:
    if not unlock_test:
        raise PermissionError(
            "Locked test evaluation requires the explicit --unlock-test flag"
        )


def evaluate_factorized_finalist_test(
    config_path: str | Path,
    calibrator_dir: str | Path,
    *,
    unlock_test: bool = False,
    sampling_seeds: Iterable[int] = (5020, 5021, 5022, 5023),
    residual_seed: int = 15_020,
    minimum_repeat_pass_fraction: float = 0.75,
) -> Path:
    _require_test_unlock(unlock_test)
    seeds = tuple(map(int, sampling_seeds))
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("Final test sampling seeds must be nonempty and unique")
    if not 0 < float(minimum_repeat_pass_fraction) <= 1:
        raise ValueError("minimum_repeat_pass_fraction must be in (0, 1]")
    config = load_factorized_config(config_path)
    options = config["evaluation"]
    data = config["data"]
    run_output = Path(config["run"]["output_dir"])
    output = run_output / "evaluation" / "final_locked_test"
    if (output / "summary.json").exists():
        raise FileExistsError(
            "Final locked-test summary already exists; refusing to overwrite it"
        )
    output.mkdir(parents=True, exist_ok=False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Finalist ModelDDIM evaluation requires CUDA")
    model, schema, payload = _load_adapter_model(config, device)
    train = load_factorized_role(
        data["prepared_h5"], data["samples_tsv"], "train"
    )
    test = load_factorized_role(data["prepared_h5"], data["samples_tsv"], "test")
    test_labels = encode_factorized_labels(test["samples"], schema)
    calibrator = PositiveResidualCalibrator.load(calibrator_dir)
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)

    rows: list[dict[str, object]] = []
    repeat_summaries: dict[str, object] = {}
    for index, seed in enumerate(seeds):
        synthetic = _sample(
            model,
            schema,
            test_labels,
            genes=len(test["genes"]),
            betas=betas,
            diffusion_timesteps=int(config["model"]["diffusion_timesteps"]),
            sampling_steps=int(options.get("sampling_steps", 1000)),
            batch_size=int(options.get("batch_size", 128)),
            guidance_scale=1.0,
            seed=seed,
            device=device,
        )
        calibrated = calibrator.apply(
            synthetic,
            test["samples"],
            seed=int(residual_seed) + index,
        )
        fidelity = generated_quality(
            test["expression"], calibrated, max_pr_samples=len(calibrated)
        )
        memorization = memorization_metrics(
            train["expression"],
            calibrated,
            max_samples=max(len(calibrated), 50),
            seed=int(config["run"]["seed"]),
        )
        selection = fidelity_selection(fidelity, memorization)
        effect = _condition_effect(
            test["expression"],
            calibrated,
            test["samples"]["condition"].astype(str).to_numpy(),
        )
        effect_gate = conditional_effect_selection(effect)
        per_tissue, tissue_validation = _per_tissue_effects(
            test["expression"], calibrated, test["samples"], test["genes"]
        )
        muscle = tissue_validation.get("skeletal_muscle", {})
        muscle_gate = muscle.get(
            "gate", accession_effect_selection({"accessions": 0})
        )
        label = f"seed{seed}"
        directory = output / label
        directory.mkdir(parents=True)
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
            directory / "calibrated_test_expression.npz",
            scaled_expression=calibrated,
            source_row=test["source_row"],
            genes=np.asarray(test["genes"]),
            sampling_seed=seed,
            residual_seed=int(residual_seed) + index,
        )
        _plot_pca(test["expression"], calibrated, test["samples"], directory)
        summary = {
            "status": "complete",
            "split": "test",
            "locked_test_opened": True,
            "sampling_seed": seed,
            "residual_seed": int(residual_seed) + index,
            "fidelity": fidelity,
            "memorization": memorization,
            "model_selection": selection,
            "pooled_condition_effect": effect,
            "pooled_condition_gate": effect_gate,
            "skeletal_muscle_accession_validation": muscle,
            "skeletal_muscle_accession_gate": muscle_gate,
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
                "fidelity_pass": selection["eligible_for_model_selection"],
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
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(output / "repeat_metrics.tsv", sep="\t", index=False)
    metric_stability = _metric_repeat_summary(table)
    required_fraction = float(minimum_repeat_pass_fraction)
    fidelity_stable = all(
        float(value["pass_fraction"]) >= required_fraction
        for value in metric_stability.values()
    ) and bool(table["diversity_pass"].all() and table["memorization_pass"].all())
    pooled_stable = bool(table["condition_effect_pass"].mean() >= required_fraction)
    muscle_stable = bool(table["muscle_accession_pass"].mean() >= required_fraction)

    base_train_directory = run_output / "evaluation" / "validation_guidance_1"
    raw_train = np.load(base_train_directory / "synthetic_train_expression.npz")
    real_train, train_samples = _aligned_expression(
        train, np.asarray(raw_train["source_row"], dtype=np.int64)
    )
    calibrated_train = calibrator.apply(
        raw_train["scaled_expression"], train_samples, seed=int(residual_seed) + 100
    )
    utility = classifier_utility(
        real_train,
        train_samples["condition"].astype(str).to_numpy(),
        test["expression"],
        test["samples"]["condition"].astype(str).to_numpy(),
        synthetic_train=calibrated_train,
        synthetic_labels=train_samples["condition"].astype(str).to_numpy(),
        allow_augmentation=bool(fidelity_stable and pooled_stable),
    )
    stable_gates = {
        "repeated_test_fidelity": fidelity_stable,
        "repeated_pooled_condition_effect": pooled_stable,
    }
    aggregate = {
        "status": "complete",
        "split": "test",
        "locked_test_opened": True,
        "finalist_scope": "pooled_all_tissue_conditional_generation",
        "sampling_seeds": list(seeds),
        "minimum_repeat_pass_fraction": required_fraction,
        "metric_repeat_stability": metric_stability,
        "pooled_condition_effect_pass_fraction": float(
            table["condition_effect_pass"].mean()
        ),
        "skeletal_muscle_accession_pass_fraction": float(
            table["muscle_accession_pass"].mean()
        ),
        "skeletal_muscle_stable": muscle_stable,
        "classifier_utility": utility,
        "independent_finalist_gates": stable_gates,
        "broad_finalist_pass": bool(all(stable_gates.values())),
        "repeat_summaries": repeat_summaries,
        "model": str(run_output / "model.pt"),
        "model_training_config_sha256": payload["metadata"].get(
            "config_sha256", ""
        ),
        "calibrator": str(Path(calibrator_dir).resolve()),
        "device": torch.cuda.get_device_name(device),
        "acceptance_rule": (
            "Each fidelity metric and pooled FLT/GC recovery must independently pass "
            "in at least the declared fraction of four predeclared test generations."
        ),
        "limitations": [
            "The split measures within-study interpolation, not unseen-study transfer.",
            "Skeletal-muscle accession recovery is reported but is not part of the "
            "broad pooled-generator acceptance rule.",
            "The locked test was opened once after model and calibration selection.",
        ],
    }
    path = output / "summary.json"
    path.write_text(json.dumps(aggregate, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# Final locked-test evaluation\n\n"
        "This directory is the one-time test evaluation of the fixed pooled OSDR "
        "factorized DDIM and its fixed train-only calibrator. Four sampling seeds were "
        "declared before opening test. Metrics are gated independently; no composite "
        "score is used. Results measure within-study interpolation.\n",
        encoding="utf-8",
    )
    return path
