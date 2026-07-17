"""Validation-only evaluation of factorized ModelDDIM adapters."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import torch

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.effect_validation import compare_real_synthetic_effects
from nasa_mouse_generative.metrics import (
    _condition_effect,
    accession_effect_selection,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)

from .factorized_adapter import (
    FactorizedAdapterDDIM,
    FactorizedSchema,
    encode_factorized_labels,
    load_factorized_role,
    neutralize_group,
)
from .factorized_config import load_factorized_config
from .factorized_train import FORMAT, _base_model
from .upstream import ddim_trajectory, quadratic_beta_schedule


def _balanced_indices(
    samples: pd.DataFrame, limit: int, seed: int
) -> np.ndarray:
    if limit <= 0 or len(samples) <= limit:
        return np.arange(len(samples), dtype=int)
    labels = (
        samples[["tissue", "condition"]]
        .fillna("unknown")
        .astype(str)
        .agg("||".join, axis=1)
        .to_numpy()
    )
    rng = np.random.default_rng(seed)
    groups = {
        label: rng.permutation(np.flatnonzero(labels == label))
        for label in sorted(set(labels))
    }
    selected: list[int] = []
    offset = 0
    while len(selected) < limit:
        added = False
        for label in groups:
            if offset < len(groups[label]):
                selected.append(int(groups[label][offset]))
                added = True
                if len(selected) == limit:
                    break
        if not added:
            break
        offset += 1
    return np.asarray(selected, dtype=int)


def _evaluation_sampling_seeds(
    run_seed: int, options: dict[str, object]
) -> tuple[int, int]:
    """Resolve fixed seeds shared by every guidance scale in one screen."""

    return (
        int(options.get("validation_sampling_seed", int(run_seed) + 1000)),
        int(options.get("train_sampling_seed", int(run_seed) + 2000)),
    )


def _variant_label(value: str) -> str:
    value = str(value).strip()
    if not value:
        return ""
    if Path(value).name != value or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
        for character in value
    ):
        raise ValueError("evaluation_variant must contain only letters, digits, _ or -")
    return value


class _GuidedModel(torch.nn.Module):
    def __init__(
        self,
        model: FactorizedAdapterDDIM,
        schema: FactorizedSchema,
        scale: float,
    ) -> None:
        super().__init__()
        self.model = model
        self.schema = schema
        self.scale = float(scale)

    def forward(
        self, expression: torch.Tensor, timesteps: torch.Tensor, labels: torch.Tensor
    ) -> torch.Tensor:
        conditional = self.model(expression, timesteps, labels)
        if self.scale == 1.0:
            return conditional
        neutral = neutralize_group(labels, self.schema, "condition")
        unconditioned = self.model(expression, timesteps, neutral)
        return unconditioned + self.scale * (conditional - unconditioned)


def _sampling_sequence(diffusion_timesteps: int, sampling_steps: int) -> list[int]:
    if sampling_steps <= 0 or sampling_steps > diffusion_timesteps:
        raise ValueError("sampling_steps must be in [1, diffusion_timesteps]")
    if sampling_steps == diffusion_timesteps:
        return list(range(diffusion_timesteps))
    values = np.linspace(0, diffusion_timesteps - 1, sampling_steps)
    return sorted(set(np.rint(values).astype(int).tolist()))


def _sample(
    model: FactorizedAdapterDDIM,
    schema: FactorizedSchema,
    labels: np.ndarray,
    *,
    genes: int,
    betas: torch.Tensor,
    diffusion_timesteps: int,
    sampling_steps: int,
    batch_size: int,
    guidance_scale: float,
    seed: int,
    device: torch.device,
) -> np.ndarray:
    guided = _GuidedModel(model, schema, guidance_scale)
    sequence = _sampling_sequence(diffusion_timesteps, sampling_steps)
    generator = torch.Generator(device=device).manual_seed(int(seed))
    collected: list[np.ndarray] = []
    for start in range(0, len(labels), int(batch_size)):
        end = min(start + int(batch_size), len(labels))
        condition = torch.as_tensor(labels[start:end], device=device)
        noise = torch.randn(
            (end - start, genes), generator=generator, device=device
        )
        generated = ddim_trajectory(
            noise,
            condition,
            guided,
            betas,
            sequence=sequence,
            snapshot_timesteps=(0,),
            eta=0.0,
            generator=generator,
        )[0]
        collected.append(generated.numpy().astype(np.float32))
        print(
            f"[factorized-ddim:evaluate] guidance={guidance_scale:g} "
            f"sampled={end}/{len(labels)}",
            flush=True,
        )
    return np.concatenate(collected)


def _class_probe(
    train: np.ndarray,
    train_labels: Iterable[object],
    synthetic: np.ndarray,
    expected: Iterable[object],
) -> dict[str, float]:
    train_labels = np.asarray(list(map(str, train_labels)))
    expected = np.asarray(list(map(str, expected)))
    represented = np.isin(train_labels, np.unique(expected))
    if len(np.unique(train_labels[represented])) < 2:
        return {"balanced_accuracy": float("nan")}
    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=3000, class_weight="balanced", random_state=0),
    )
    classifier.fit(train[represented], train_labels[represented])
    predicted = classifier.predict(synthetic)
    return {
        "balanced_accuracy": float(balanced_accuracy_score(expected, predicted))
    }


def _per_tissue_effects(
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
    genes: list[str],
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    validation: dict[str, dict[str, object]] = {}
    for tissue in sorted(samples["tissue"].astype(str).unique()):
        mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
        conditions = samples.loc[mask, "condition"].astype(str).to_numpy()
        effect = _condition_effect(real[mask], synthetic[mask], conditions)
        tables, accession_summary = compare_real_synthetic_effects(
            real[mask],
            synthetic[mask],
            samples.loc[mask].reset_index(drop=True),
            genes,
        )
        gate = accession_effect_selection(accession_summary)
        validation[tissue] = {"summary": accession_summary, "gate": gate}
        rows.append(
            {
                "tissue": tissue,
                "profiles": int(mask.sum()),
                **effect,
                "condition_gate_passed": conditional_effect_selection(effect)[
                    "passed"
                ],
                "accessions": accession_summary.get("accessions", 0),
                "accession_meta_effect_correlation": accession_summary.get(
                    "meta_effect_correlation", float("nan")
                ),
                "accession_meta_direction_agreement": accession_summary.get(
                    "meta_direction_agreement", float("nan")
                ),
                "accession_gate_passed": gate["passed"],
            }
        )
        comparison = tables.get("comparison", pd.DataFrame())
        if not comparison.empty:
            validation[tissue]["comparison"] = comparison
    return pd.DataFrame(rows), validation


def _plot_pca(
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
    output: Path,
) -> None:
    coordinates = PCA(n_components=2, random_state=0).fit_transform(
        np.concatenate([real, synthetic])
    )
    real_coordinates = coordinates[: len(real)]
    synthetic_coordinates = coordinates[len(real) :]
    tissues = sorted(samples["tissue"].astype(str).unique())
    colors = plt.get_cmap("tab20", max(len(tissues), 1))
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.4))
    for index, tissue in enumerate(tissues):
        mask = samples["tissue"].astype(str).eq(tissue).to_numpy()
        axes[0].scatter(
            real_coordinates[mask, 0],
            real_coordinates[mask, 1],
            s=20,
            alpha=0.45,
            color=colors(index),
            edgecolors="none",
        )
        axes[0].scatter(
            synthetic_coordinates[mask, 0],
            synthetic_coordinates[mask, 1],
            s=24,
            marker="x",
            alpha=0.78,
            color=colors(index),
            label=tissue,
        )
    condition_colors = {"flight": "#C14924", "ground_control": "#176B87"}
    for condition, color in condition_colors.items():
        mask = samples["condition"].astype(str).eq(condition).to_numpy()
        axes[1].scatter(
            real_coordinates[mask, 0],
            real_coordinates[mask, 1],
            s=20,
            alpha=0.40,
            color=color,
            edgecolors="none",
        )
        axes[1].scatter(
            synthetic_coordinates[mask, 0],
            synthetic_coordinates[mask, 1],
            s=24,
            marker="x",
            alpha=0.78,
            color=color,
            label=condition,
        )
    axes[0].set_title("Tissue: real circles, synthetic crosses")
    axes[1].set_title("Condition: real circles, synthetic crosses")
    axes[0].legend(frameon=False, fontsize=6, ncol=2, bbox_to_anchor=(1.0, 1.0))
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel("PCA 1")
        axis.set_ylabel("PCA 2")
        axis.grid(alpha=0.15)
    figure.tight_layout()
    figure.savefig(output / "real_vs_synthetic_pca.png", dpi=220, bbox_inches="tight")
    figure.savefig(output / "real_vs_synthetic_pca.pdf", bbox_inches="tight")
    plt.close(figure)


def _load_adapter_model(
    config: dict[str, Any], device: torch.device, artifact_name: str = "model.pt"
) -> tuple[FactorizedAdapterDDIM, FactorizedSchema, dict[str, Any]]:
    output = Path(config["run"]["output_dir"])
    payload = torch.load(
        output / artifact_name, map_location="cpu", weights_only=False
    )
    if payload.get("format") != FORMAT:
        raise ValueError("Incompatible factorized adapter model")
    schema = FactorizedSchema.from_dict(payload["metadata"]["schema"])
    base, _, classes = _base_model(
        config["data"]["pretrained_model"],
        config["model"],
        len(payload["metadata"]["genes"]),
    )
    if tuple(classes) != schema.base_classes:
        raise ValueError("Pretrained classes differ from the adapter schema")
    adapter_options = config.get("adapter", {})
    model = FactorizedAdapterDDIM(
        base,
        schema,
        domain_lora_rank=int(adapter_options.get("domain_lora_rank", 0)),
        domain_lora_alpha=float(adapter_options.get("domain_lora_alpha", 1.0)),
    )
    model.load_adapter_state_dict(payload["adapter_state_dict"])
    model.to(device).eval()
    return model, schema, payload


def evaluate_factorized(
    config_path: str | Path,
    *,
    guidance_scales: Iterable[float] | None = None,
    model_artifact: str = "model.pt",
    validation_sampling_seed: int | None = None,
    train_sampling_seed: int | None = None,
    evaluation_variant: str = "",
) -> Path:
    config = load_factorized_config(config_path)
    options = config["evaluation"]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Factorized paper ModelDDIM evaluation requires CUDA")
    if Path(model_artifact).name != model_artifact:
        raise ValueError("model_artifact must be a filename under the run directory")
    variant = _variant_label(evaluation_variant)
    model, schema, payload = _load_adapter_model(
        config, device, artifact_name=model_artifact
    )
    data_options = config["data"]
    train = load_factorized_role(
        data_options["prepared_h5"], data_options["samples_tsv"], "train"
    )
    validation = load_factorized_role(
        data_options["prepared_h5"], data_options["samples_tsv"], "validation"
    )
    train_labels = encode_factorized_labels(train["samples"], schema)
    validation_labels = encode_factorized_labels(validation["samples"], schema)
    validation_indices = _balanced_indices(
        validation["samples"],
        int(options.get("metric_samples", len(validation["expression"]))),
        int(config["run"]["seed"]) + 1,
    )
    train_indices = _balanced_indices(
        train["samples"],
        int(options.get("train_metric_samples", 240)),
        int(config["run"]["seed"]) + 2,
    )
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)
    scales = list(
        map(
            float,
            guidance_scales
            if guidance_scales is not None
            else options.get("guidance_scales", [1.0]),
        )
    )
    run_output = Path(config["run"]["output_dir"])
    summaries: list[dict[str, object]] = []
    summary_paths: list[Path] = []
    default_validation_seed, default_train_seed = _evaluation_sampling_seeds(
        int(config["run"]["seed"]), options
    )
    validation_sampling_seed = (
        default_validation_seed
        if validation_sampling_seed is None
        else int(validation_sampling_seed)
    )
    train_sampling_seed = (
        default_train_seed if train_sampling_seed is None else int(train_sampling_seed)
    )
    for scale in scales:
        artifact_label = Path(model_artifact).stem
        output_label = (
            f"validation_guidance_{scale:g}"
            if model_artifact == "model.pt"
            else f"validation_{artifact_label}_guidance_{scale:g}"
        )
        if variant:
            output_label = f"{output_label}_{variant}"
        output = run_output / "evaluation" / output_label
        output.mkdir(parents=True, exist_ok=True)
        started = time.time()
        synthetic_all = _sample(
            model,
            schema,
            validation_labels,
            genes=len(validation["genes"]),
            betas=betas,
            diffusion_timesteps=int(config["model"]["diffusion_timesteps"]),
            sampling_steps=int(options.get("sampling_steps", 1000)),
            batch_size=int(options.get("batch_size", 128)),
            guidance_scale=scale,
            seed=validation_sampling_seed,
            device=device,
        )
        synthetic_train = _sample(
            model,
            schema,
            train_labels[train_indices],
            genes=len(train["genes"]),
            betas=betas,
            diffusion_timesteps=int(config["model"]["diffusion_timesteps"]),
            sampling_steps=int(options.get("sampling_steps", 1000)),
            batch_size=int(options.get("batch_size", 128)),
            guidance_scale=scale,
            seed=train_sampling_seed,
            device=device,
        )
        real = validation["expression"][validation_indices]
        synthetic = synthetic_all[validation_indices]
        samples = validation["samples"].iloc[validation_indices].reset_index(drop=True)
        fidelity = generated_quality(real, synthetic, max_pr_samples=len(real))
        memorization = memorization_metrics(
            train["expression"],
            synthetic,
            max_samples=max(len(real), 50),
            seed=int(config["run"]["seed"]),
        )
        selection = fidelity_selection(fidelity, memorization)
        train_fidelity = generated_quality(
            train["expression"][train_indices],
            synthetic_train,
            max_pr_samples=len(train_indices),
        )
        train_memorization = memorization_metrics(
            train["expression"],
            synthetic_train,
            max_samples=max(len(train_indices), 50),
            seed=int(config["run"]["seed"]) + 1,
        )
        train_selection = fidelity_selection(train_fidelity, train_memorization)
        effect = _condition_effect(
            real,
            synthetic,
            samples["condition"].astype(str).to_numpy(),
        )
        effect_gate = conditional_effect_selection(effect)
        per_tissue, tissue_validation = _per_tissue_effects(
            real, synthetic, samples, validation["genes"]
        )
        per_tissue.to_csv(
            output / "per_tissue_condition_recovery.tsv", sep="\t", index=False
        )
        for tissue, result in tissue_validation.items():
            comparison = result.pop("comparison", None)
            if comparison is not None:
                comparison.to_csv(
                    output / f"{tissue}_accession_effect_recovery.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
        muscle = tissue_validation.get("skeletal_muscle", {})
        muscle_gate = muscle.get(
            "gate",
            accession_effect_selection({"accessions": 0}),
        )
        condition_probe = _class_probe(
            train["expression"],
            train["samples"]["condition"],
            synthetic,
            samples["condition"],
        )
        tissue_probe = _class_probe(
            train["expression"],
            train["samples"]["tissue"],
            synthetic,
            samples["tissue"],
        )
        all_gates = {
            "validation_fidelity": bool(selection["eligible_for_model_selection"]),
            "train_fidelity": bool(train_selection["eligible_for_model_selection"]),
            "pooled_condition_effect": bool(effect_gate["passed"]),
            "muscle_accession_effect": bool(muscle_gate["passed"]),
        }
        _plot_pca(real, synthetic, samples, output)
        np.savez_compressed(
            output / "synthetic_validation_expression.npz",
            scaled_expression=synthetic_all,
            source_row=validation["source_row"],
            genes=np.asarray(validation["genes"]),
            guidance_scale=scale,
        )
        np.savez_compressed(
            output / "synthetic_train_expression.npz",
            scaled_expression=synthetic_train,
            source_row=train["source_row"][train_indices],
            genes=np.asarray(train["genes"]),
            guidance_scale=scale,
        )
        summary: dict[str, object] = {
            "status": "complete",
            "split": "validation",
            "locked_test_opened": False,
            "guidance_scale": scale,
            "sampling_steps": int(options.get("sampling_steps", 1000)),
            "validation_sampling_seed": validation_sampling_seed,
            "train_sampling_seed": train_sampling_seed,
            "evaluation_variant": variant,
            "sampling_seconds": float(time.time() - started),
            "profiles": int(len(real)),
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
            "skeletal_muscle_accession_validation": muscle,
            "condition_probe": condition_probe,
            "tissue_probe": tissue_probe,
            "independent_acceptance_gates": all_gates,
            "all_acceptance_gates_pass": bool(all(all_gates.values())),
            "acceptance_rule": (
                "validation fidelity, paper-train fidelity, pooled condition effect, "
                "and skeletal-muscle accession effect must pass independently"
            ),
            "device": torch.cuda.get_device_name(device),
            "model": str(run_output / model_artifact),
            "pretrained_model": payload["metadata"]["pretrained_model"],
        }
        summary_path = output / "summary.json"
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        summary_paths.append(summary_path)
        summaries.append(
            {
                "guidance_scale": scale,
                "correlation": fidelity["correlation_matrix_agreement"],
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
                "condition_probe_balanced_accuracy": condition_probe[
                    "balanced_accuracy"
                ],
                "tissue_probe_balanced_accuracy": tissue_probe["balanced_accuracy"],
                **{f"gate_{key}": value for key, value in all_gates.items()},
                "all_gates_pass": bool(all(all_gates.values())),
            }
        )
    table = pd.DataFrame(summaries)
    screen_name = (
        "guidance_screen.tsv"
        if model_artifact == "model.pt"
        else f"{Path(model_artifact).stem}_guidance_screen.tsv"
    )
    if variant:
        screen_name = f"{Path(screen_name).stem}_{variant}.tsv"
    table.to_csv(run_output / "evaluation" / screen_name, sep="\t", index=False)
    return summary_paths[0]
