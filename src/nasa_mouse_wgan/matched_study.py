"""Matched study-conditioned WGAN-GP for the OSDR diffusion finalist split."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import shutil
import time
from typing import Any, Iterable

import h5py
import numpy as np
import pandas as pd
import torch
import yaml

from nasa_mouse_diffusion.evaluate import generated_quality
from nasa_mouse_generative.adapters.wgan import WGANAdapter
from nasa_mouse_generative.adapters.base import seed_everything
from nasa_mouse_generative.conditioning import CategoryEncoder
from nasa_mouse_generative.metrics import (
    _condition_effect,
    classifier_utility,
    conditional_effect_selection,
    fidelity_selection,
    memorization_metrics,
)
from nasa_mouse_generative.training_data import DataPartition
from nasa_mouse_rna_diffusion.factorized_adapter import load_factorized_role
from nasa_mouse_rna_diffusion.factorized_distribution_calibrate import (
    PositiveResidualCalibrator,
    _metric_repeat_summary,
)
from nasa_mouse_rna_diffusion.factorized_evaluate import (
    _per_tissue_effects,
    _plot_pca,
)


FORMAT = "vinas_wgan_gp_matched_study_v1"
REFERENCE_VALUE = "archs4_reference"
TISSUE_FALLBACKS = {
    "cecum": "colon",
    "cells": "cultured_cells",
    "eye": "retina",
    "optic_nerve": "retina",
}


def _decode(values: Iterable[object]) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    config = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    if config.get("contract") != FORMAT:
        raise ValueError(f"Expected contract={FORMAT!r}")
    required = {
        "run": ("output_dir", "seed", "source_root"),
        "data": ("archs4_h5", "osdr_h5", "osdr_samples"),
        "conditioning": ("covariates",),
        "model": (
            "noise_dim",
            "hidden_dims",
            "batch_size",
            "reference_epochs",
            "finetune_epochs",
        ),
        "training": ("checkpoint_every_epochs", "device"),
        "evaluation": (
            "train_generation_seeds",
            "validation_generation_seeds",
            "test_generation_seeds",
            "minimum_repeat_pass_fraction",
        ),
    }
    for section, names in required.items():
        if section not in config:
            raise ValueError(f"Missing config section {section!r}")
        missing = [name for name in names if name not in config[section]]
        if missing:
            raise ValueError(f"Config section {section!r} lacks {missing}")
    covariates = tuple(map(str, config["conditioning"]["covariates"]))
    expected = {
        "tissue",
        "condition",
        "study",
        "material_type",
        "sex",
        "muscle_group",
    }
    if set(covariates) != expected:
        raise ValueError(
            "Matched conditioning must contain tissue, condition, study, "
            "material_type, sex, and muscle_group exactly once"
        )
    fraction = float(config["evaluation"]["minimum_repeat_pass_fraction"])
    if not 0 < fraction <= 1:
        raise ValueError("minimum_repeat_pass_fraction must be in (0, 1]")
    config["_config_path"] = str(source.resolve())
    config["_config_sha256"] = _sha256(source)
    return config


@dataclass
class NativeLogZScaler:
    """WGAN paper transform: log1p expression and train-gene z-score."""

    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, tpm: np.ndarray) -> "NativeLogZScaler":
        logged = np.log1p(np.maximum(np.asarray(tpm, dtype=np.float64), 0.0))
        mean = logged.mean(axis=0)
        scale = logged.std(axis=0)
        scale[~np.isfinite(scale) | (scale < 1e-6)] = 1.0
        return cls(mean=mean.astype(np.float32), scale=scale.astype(np.float32))

    def transform(self, tpm: np.ndarray) -> np.ndarray:
        logged = np.log1p(np.maximum(np.asarray(tpm, dtype=np.float32), 0.0))
        return ((logged - self.mean) / self.scale).astype(np.float32)

    def inverse_to_scaled(
        self, transformed: np.ndarray, maxabs_scale: np.ndarray
    ) -> np.ndarray:
        logged = (
            np.asarray(transformed, dtype=np.float64) * self.scale + self.mean
        )
        # The source workflow clips exponent inputs before reversing log1p.
        tpm = np.expm1(np.clip(logged, -30.0, 30.0))
        np.maximum(tpm, 0.0, out=tpm)
        denominator = np.maximum(np.asarray(maxabs_scale, dtype=np.float64), 1e-8)
        return (tpm / denominator).astype(np.float32)

    def save(self, path: str | Path) -> Path:
        output = Path(path)
        np.savez_compressed(output, mean=self.mean, scale=self.scale)
        return output

    @classmethod
    def load(cls, path: str | Path) -> "NativeLogZScaler":
        arrays = np.load(path)
        return cls(
            mean=np.asarray(arrays["mean"], dtype=np.float32),
            scale=np.asarray(arrays["scale"], dtype=np.float32),
        )


@dataclass
class MatchedData:
    genes: list[str]
    maxabs_scale: np.ndarray
    scaler: NativeLogZScaler
    encoder: CategoryEncoder
    reference: dict[str, DataPartition]
    osdr: dict[str, DataPartition]
    osdr_scaled: dict[str, np.ndarray]
    source_manifest: dict[str, object]


def _reference_role(
    handle: h5py.File, role: str, classes: list[str], maxabs_scale: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame]:
    group = handle[role]
    expression = np.asarray(group["expression"][:], dtype=np.float32)
    tpm = expression * maxabs_scale
    class_index = np.asarray(group["class_index"][:], dtype=np.int64)
    source_row = np.asarray(group["source_row"][:], dtype=np.int64)
    tissues = np.asarray(classes, dtype=object)[class_index]
    obs = pd.DataFrame(
        {
            "profile_id": [f"archs4:{value}" for value in source_row],
            "accession": REFERENCE_VALUE,
            "condition": REFERENCE_VALUE,
            "tissue": tissues,
            "material_type": REFERENCE_VALUE,
            "study": REFERENCE_VALUE,
            "sex": "unknown_sex",
            "muscle_group": "not_applicable",
            "source": "archs4",
            "source_row": source_row,
        }
    )
    return tpm, obs


def _partition(
    name: str,
    tpm: np.ndarray,
    obs: pd.DataFrame,
    scaler: NativeLogZScaler,
    encoder: CategoryEncoder,
) -> DataPartition:
    weights = np.ones(len(obs), dtype=np.float32)
    if len(weights):
        weights /= weights.sum()
    return DataPartition(
        name=name,
        matrix=scaler.transform(tpm),
        obs=obs.reset_index(drop=True),
        categories=encoder.transform(obs),
        weights=weights,
    )


def prepare_matched_data(
    config: dict[str, Any], *, include_test: bool = False
) -> MatchedData:
    data = config["data"]
    archs4_path = Path(data["archs4_h5"])
    osdr_path = Path(data["osdr_h5"])
    sample_path = Path(data["osdr_samples"])
    for path in (archs4_path, osdr_path, sample_path):
        if not path.exists():
            raise FileNotFoundError(path)

    with h5py.File(archs4_path, "r") as handle:
        genes = _decode(handle["genes"][:])
        classes = _decode(handle["classes"][:])
        maxabs_scale = np.asarray(handle["maxabs_scale"][:], dtype=np.float32)
        reference_raw = {
            role: _reference_role(handle, role, classes, maxabs_scale)
            for role in ("train", "validation")
        }

    roles = ["train", "validation"] + (["test"] if include_test else [])
    osdr_raw: dict[str, dict[str, object]] = {}
    for role in roles:
        loaded = load_factorized_role(osdr_path, sample_path, role)
        if list(loaded["genes"]) != genes:
            raise ValueError(f"OSDR {role} genes differ from ARCHS4 reference")
        if not np.allclose(loaded["maxabs_scale"], maxabs_scale, rtol=1e-6, atol=1e-8):
            raise ValueError("OSDR and ARCHS4 MaxAbs scales differ")
        osdr_raw[role] = loaded

    scaler = NativeLogZScaler.fit(reference_raw["train"][0])
    covariates = tuple(map(str, config["conditioning"]["covariates"]))
    encoder = CategoryEncoder.fit(
        [reference_raw["train"][1], osdr_raw["train"]["samples"]], covariates
    )
    reference = {
        role: _partition(role, tpm, obs, scaler, encoder)
        for role, (tpm, obs) in reference_raw.items()
    }
    osdr: dict[str, DataPartition] = {}
    osdr_scaled: dict[str, np.ndarray] = {}
    for role, loaded in osdr_raw.items():
        tpm = np.asarray(loaded["analysis_expression"], dtype=np.float32)
        obs = loaded["samples"].copy()
        osdr[role] = _partition(role, tpm, obs, scaler, encoder)
        osdr_scaled[role] = np.asarray(loaded["expression"], dtype=np.float32)
        for index, covariate in enumerate(encoder.covariates):
            unknown = encoder.vocabularies[covariate].index("__unknown__")
            if np.any(osdr[role].categories[:, index] == unknown):
                raise ValueError(
                    f"OSDR {role} contains unseen {covariate!r} categories"
                )

    return MatchedData(
        genes=genes,
        maxabs_scale=maxabs_scale,
        scaler=scaler,
        encoder=encoder,
        reference=reference,
        osdr=osdr,
        osdr_scaled=osdr_scaled,
        source_manifest={
            "raw_integrated_osdr_h5_used": False,
            "archs4_h5": str(archs4_path.resolve()),
            "archs4_h5_sha256": _sha256(archs4_path),
            "osdr_h5": str(osdr_path.resolve()),
            "osdr_h5_sha256": _sha256(osdr_path),
            "osdr_samples": str(sample_path.resolve()),
            "osdr_samples_sha256": _sha256(sample_path),
            "normalization": "full_transcriptome_tpm_then_log1p_archs4_train_zscore",
            "split": "within_accession_tissue_condition",
        },
    )


def _copy_embedding_rows(
    embedding: torch.nn.Embedding,
    vocabulary: list[str],
    source: str,
    targets: Iterable[str],
) -> int:
    lookup = {value: index for index, value in enumerate(vocabulary)}
    if source not in lookup:
        raise ValueError(f"Embedding source category {source!r} is missing")
    copied = 0
    with torch.no_grad():
        for target in targets:
            if target in lookup and target != source:
                embedding.weight[lookup[target]].copy_(embedding.weight[lookup[source]])
                copied += 1
    return copied


def initialize_query_embeddings(
    adapter: WGANAdapter, encoder: CategoryEncoder
) -> dict[str, int]:
    """Make every new OSDR covariate profile equal to a pretrained profile."""

    modules = {
        "generator": adapter.model.generator.covariates.embeddings,
        "critic": adapter.model.critic.covariates.embeddings,
    }
    counts: dict[str, int] = {}
    for module_name, embeddings in modules.items():
        for index, covariate in enumerate(encoder.covariates):
            vocabulary = encoder.vocabularies[covariate]
            if covariate == "tissue":
                copied = 0
                for target, source in TISSUE_FALLBACKS.items():
                    copied += _copy_embedding_rows(
                        embeddings[index], vocabulary, source, (target,)
                    )
            else:
                source = {
                    "condition": REFERENCE_VALUE,
                    "study": REFERENCE_VALUE,
                    "material_type": REFERENCE_VALUE,
                    "sex": "unknown_sex",
                    "muscle_group": "not_applicable",
                }[covariate]
                copied = _copy_embedding_rows(
                    embeddings[index], vocabulary, source, vocabulary
                )
            counts[f"{module_name}.{covariate}"] = copied
    return counts


def _write_prepared_artifacts(
    output: Path, config: dict[str, Any], data: MatchedData
) -> None:
    data.encoder.save(output / "categorical_encoder.json")
    data.scaler.save(output / "native_log_z_scaler.npz")
    (output / "genes.tsv").write_text(
        "gene_id\n" + "\n".join(data.genes) + "\n", encoding="utf-8"
    )
    payload = {
        "contract": FORMAT,
        "config": config["_config_path"],
        "config_sha256": config["_config_sha256"],
        "genes": len(data.genes),
        "covariates": list(data.encoder.covariates),
        "cardinalities": data.encoder.cardinalities,
        "reference_profiles": {
            key: len(value) for key, value in data.reference.items()
        },
        "osdr_profiles_loaded": {key: len(value) for key, value in data.osdr.items()},
        "test_loaded": "test" in data.osdr,
        "sources": data.source_manifest,
    }
    (output / "data_manifest.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def train(config_path: str | Path) -> Path:
    config = load_config(config_path)
    output = Path(config["run"]["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "training_summary.json"
    if summary_path.exists() and (output / "model.pt").exists():
        return summary_path
    data = prepare_matched_data(config, include_test=False)
    _write_prepared_artifacts(output, config, data)
    shutil.copy2(config_path, output / "resolved_config.yaml")

    parameters = dict(config["model"])
    parameters["_paper_native"] = True
    adapter = WGANAdapter(
        genes=data.genes,
        cardinalities=data.encoder.cardinalities,
        covariates=data.encoder.covariates,
        parameters=parameters,
        device_spec=str(config["training"]["device"]),
        output_dir=output,
        checkpoint_every=int(config["training"]["checkpoint_every_epochs"]),
        resume=True,
        seed=int(config["run"]["seed"]),
        num_workers=int(config["training"].get("num_workers", 0)),
        source_path=str(config["run"]["source_root"]),
        validation_partition=data.reference["validation"],
    )
    reference_initializer = str(
        config["training"].get("reference_initializer", "")
    ).strip()
    initialization_lineage: dict[str, object] = {}
    if reference_initializer and not adapter.checkpoint_path.exists():
        initializer_path = Path(reference_initializer)
        if not initializer_path.exists():
            raise FileNotFoundError(initializer_path)
        payload = torch.load(
            initializer_path, map_location=adapter.device, weights_only=False
        )
        adapter._validate_payload(payload)
        adapter.model.load_state_dict(payload["model_state_dict"])
        adapter._restore_common(payload)
        adapter.early_stopped_stages = set(
            map(str, payload.get("early_stopped_stages", []))
        )
        if "reference" not in adapter.early_stopped_stages:
            raise ValueError(
                "reference_initializer must contain a completed reference stage"
            )
        adapter._resume_payload = None
        seed_everything(adapter.seed)
        shutil.copy2(initializer_path, output / "reference_model.pt")
        initialization_lineage = {
            "reference_initializer": str(initializer_path.resolve()),
            "reference_initializer_sha256": _sha256(initializer_path),
            "imported_completed_epochs": dict(adapter.state.completed_epochs),
            "imported_early_stopped_stages": sorted(adapter.early_stopped_stages),
        }
        (output / "reference_initialization.json").write_text(
            json.dumps(initialization_lineage, indent=2) + "\n",
            encoding="utf-8",
        )
    elif (output / "reference_initialization.json").exists():
        initialization_lineage = json.loads(
            (output / "reference_initialization.json").read_text(
                encoding="utf-8"
            )
        )
    if adapter.device.type != "cuda":
        raise RuntimeError("Matched WGAN training requires CUDA")
    torch.cuda.reset_peak_memory_stats(adapter.device)
    started = time.time()
    runtimes: dict[str, float] = {}

    stage_started = time.time()
    adapter.fit_stage(
        data.reference["train"],
        stage="reference",
        epochs=int(parameters["reference_epochs"]),
        learning_rate=float(parameters["learning_rate"]),
    )
    runtimes["reference"] = float(time.time() - stage_started)

    initialization: dict[str, int] = {}
    if int(adapter.state.completed_epochs.get("osdr_finetune", 0)) == 0:
        reference_path = output / "reference_model.pt"
        if not reference_path.exists():
            adapter.save_final()
            shutil.copy2(output / "model.pt", reference_path)
        initialization = initialize_query_embeddings(adapter, data.encoder)
        (output / "query_embedding_initialization.json").write_text(
            json.dumps(
                {
                    "method": "copy_reference_rows_with_tissue_fallbacks",
                    "copied_rows": initialization,
                    "tissue_fallbacks": TISSUE_FALLBACKS,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    adapter.validation_partition = data.osdr["validation"]
    stage_started = time.time()
    adapter.fit_stage(
        data.osdr["train"],
        stage="osdr_finetune",
        epochs=int(parameters["finetune_epochs"]),
        learning_rate=float(parameters["finetune_learning_rate"]),
    )
    runtimes["osdr_finetune"] = float(time.time() - stage_started)
    model_path = adapter.save_final()
    summary = {
        "status": "complete",
        "contract": FORMAT,
        "model": str(model_path),
        "device": adapter.device_summary(),
        "genes": len(data.genes),
        "conditioning_covariates": list(data.encoder.covariates),
        "completed_epochs": adapter.state.completed_epochs,
        "early_stopped_stages": sorted(adapter.early_stopped_stages),
        "early_stopping_state": adapter.early_stopping_state,
        "stage_training_seconds": runtimes,
        "training_seconds": float(time.time() - started),
        "cuda_peak_memory_gb": float(
            torch.cuda.max_memory_allocated(adapter.device) / 1024**3
        ),
        "query_embedding_initialization": initialization,
        "initialization_lineage": initialization_lineage,
        "test_loaded": False,
        "data": data.source_manifest,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path


def _load_model_and_data(
    config: dict[str, Any], *, include_test: bool
) -> tuple[WGANAdapter, MatchedData]:
    output = Path(config["run"]["output_dir"])
    if not (output / "training_summary.json").exists():
        raise FileNotFoundError("Train the matched WGAN before evaluation")
    data = prepare_matched_data(config, include_test=include_test)
    adapter = WGANAdapter.load(
        output, device_spec=str(config["training"]["device"])
    )
    if adapter.covariates != data.encoder.covariates:
        raise ValueError("Saved WGAN covariates differ from the matched data")
    return adapter, data


def _generate_scaled(
    adapter: WGANAdapter,
    partition: DataPartition,
    scaler: NativeLogZScaler,
    maxabs_scale: np.ndarray,
    *,
    seed: int,
    batch_size: int,
) -> np.ndarray:
    transformed = adapter.generate(
        partition.categories, seed=int(seed), batch_size=int(batch_size)
    )
    return scaler.inverse_to_scaled(transformed, maxabs_scale)


def _repeat_result(
    *,
    label: str,
    split: str,
    seed: int,
    real: np.ndarray,
    synthetic: np.ndarray,
    samples: pd.DataFrame,
    genes: list[str],
    train_real: np.ndarray,
    run_seed: int,
    output: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    output.mkdir(parents=True, exist_ok=True)
    fidelity = generated_quality(real, synthetic, max_pr_samples=len(real))
    memorization = memorization_metrics(
        train_real,
        synthetic,
        max_samples=max(len(synthetic), 50),
        seed=int(run_seed),
    )
    selection = fidelity_selection(fidelity, memorization)
    effect = _condition_effect(
        real, synthetic, samples["condition"].astype(str).to_numpy()
    )
    effect_gate = conditional_effect_selection(effect)
    per_tissue, tissue_validation = _per_tissue_effects(
        real, synthetic, samples, genes
    )
    per_tissue.to_csv(
        output / "per_tissue_condition_recovery.tsv", sep="\t", index=False
    )
    muscle = dict(tissue_validation.get("skeletal_muscle", {}))
    comparison = muscle.pop("comparison", None)
    if comparison is not None:
        comparison.to_csv(
            output / "skeletal_muscle_accession_effect_recovery.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    muscle_gate = muscle.get("gate", {"passed": False})
    np.savez_compressed(
        output / "synthetic_scaled_expression.npz",
        scaled_expression=synthetic,
        sampling_seed=int(seed),
    )
    _plot_pca(real, synthetic, samples, output)
    summary = {
        "status": "complete",
        "split": split,
        "variant": label,
        "sampling_seed": int(seed),
        "fidelity": fidelity,
        "memorization": memorization,
        "model_selection": selection,
        "pooled_condition_effect": effect,
        "pooled_condition_gate": effect_gate,
        "skeletal_muscle_accession_validation": muscle,
        "skeletal_muscle_accession_gate": muscle_gate,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    row = {
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
        "negative_fraction": float(np.mean(synthetic < 0)),
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
        "muscle_accession_pass": muscle_gate.get("passed", False),
    }
    return summary, row


def _stability(
    table: pd.DataFrame, minimum_fraction: float
) -> tuple[dict[str, object], bool, bool, bool]:
    metrics = _metric_repeat_summary(table)
    fidelity = all(
        float(value["pass_fraction"]) >= minimum_fraction
        for value in metrics.values()
    ) and bool(table["diversity_pass"].all() and table["memorization_pass"].all())
    condition = bool(table["condition_effect_pass"].mean() >= minimum_fraction)
    muscle = bool(table["muscle_accession_pass"].mean() >= minimum_fraction)
    return metrics, fidelity, condition, muscle


def evaluate_validation(config_path: str | Path) -> Path:
    config = load_config(config_path)
    adapter, data = _load_model_and_data(config, include_test=False)
    options = config["evaluation"]
    output = Path(config["run"]["output_dir"]) / "evaluation" / "matched_validation"
    summary_path = output / "summary.json"
    if summary_path.exists():
        raise FileExistsError(f"Validation result already exists: {summary_path}")
    output.mkdir(parents=True, exist_ok=False)
    batch_size = int(options["batch_size"])

    real_fit_parts: list[np.ndarray] = []
    synthetic_fit_parts: list[np.ndarray] = []
    metadata_fit_parts: list[pd.DataFrame] = []
    for seed in map(int, options["train_generation_seeds"]):
        generated = _generate_scaled(
            adapter,
            data.osdr["train"],
            data.scaler,
            data.maxabs_scale,
            seed=seed,
            batch_size=batch_size,
        )
        real_fit_parts.append(data.osdr_scaled["train"])
        synthetic_fit_parts.append(generated)
        metadata_fit_parts.append(data.osdr["train"].obs.copy())

    calibrator = PositiveResidualCalibrator(
        ("accession", "tissue"),
        float(options["calibration_prior_strength"]),
        float(options["calibration_residual_scale"]),
        clip_nonnegative=True,
        noise_group_columns=("accession", "tissue", "condition"),
    ).fit(
        np.concatenate(real_fit_parts),
        np.concatenate(synthetic_fit_parts),
        pd.concat(metadata_fit_parts, ignore_index=True),
    )
    calibrator_dir = output / "train_only_calibrator"
    calibrator.save(calibrator_dir)

    raw_rows: list[dict[str, object]] = []
    calibrated_rows: list[dict[str, object]] = []
    raw_summaries: dict[str, object] = {}
    calibrated_summaries: dict[str, object] = {}
    residual_seed = int(options["calibration_residual_seed"])
    for index, seed in enumerate(map(int, options["validation_generation_seeds"])):
        label = f"seed{seed}"
        raw = _generate_scaled(
            adapter,
            data.osdr["validation"],
            data.scaler,
            data.maxabs_scale,
            seed=seed,
            batch_size=batch_size,
        )
        raw_summary, raw_row = _repeat_result(
            label=label,
            split="validation",
            seed=seed,
            real=data.osdr_scaled["validation"],
            synthetic=raw,
            samples=data.osdr["validation"].obs,
            genes=data.genes,
            train_real=data.osdr_scaled["train"],
            run_seed=int(config["run"]["seed"]),
            output=output / "raw" / label,
        )
        calibrated = calibrator.apply(
            raw, data.osdr["validation"].obs, seed=residual_seed + index
        )
        calibrated_summary, calibrated_row = _repeat_result(
            label=label,
            split="validation",
            seed=seed,
            real=data.osdr_scaled["validation"],
            synthetic=calibrated,
            samples=data.osdr["validation"].obs,
            genes=data.genes,
            train_real=data.osdr_scaled["train"],
            run_seed=int(config["run"]["seed"]),
            output=output / "calibrated" / label,
        )
        raw_rows.append(raw_row)
        calibrated_rows.append(calibrated_row)
        raw_summaries[label] = raw_summary
        calibrated_summaries[label] = calibrated_summary

    raw_table = pd.DataFrame(raw_rows)
    calibrated_table = pd.DataFrame(calibrated_rows)
    raw_table.to_csv(output / "raw_repeat_metrics.tsv", sep="\t", index=False)
    calibrated_table.to_csv(
        output / "calibrated_repeat_metrics.tsv", sep="\t", index=False
    )
    minimum = float(options["minimum_repeat_pass_fraction"])
    raw_metrics, raw_fidelity, raw_condition, raw_muscle = _stability(
        raw_table, minimum
    )
    metrics, fidelity, condition, muscle = _stability(calibrated_table, minimum)
    summary = {
        "status": "complete",
        "split": "validation",
        "locked_test_opened": False,
        "model": str(Path(config["run"]["output_dir"]) / "model.pt"),
        "sampling_seeds": list(map(int, options["validation_generation_seeds"])),
        "minimum_repeat_pass_fraction": minimum,
        "raw": {
            "metric_repeat_stability": raw_metrics,
            "repeated_fidelity": raw_fidelity,
            "pooled_condition_effect": raw_condition,
            "skeletal_muscle_diagnostic": raw_muscle,
            "repeat_summaries": raw_summaries,
        },
        "calibrated": {
            "metric_repeat_stability": metrics,
            "repeated_fidelity": fidelity,
            "pooled_condition_effect": condition,
            "skeletal_muscle_diagnostic": muscle,
            "repeat_summaries": calibrated_summaries,
        },
        "validation_candidate_pass": bool(fidelity and condition),
        "acceptance_scope": "pooled within-study conditional generation",
        "test_loaded": False,
        "calibrator": str(calibrator_dir),
        "selection_rule": (
            "Every fidelity metric and pooled FLT/GC effect must independently pass "
            "in at least the declared fraction of generation repeats; no composite."
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# Matched WGAN validation\n\n"
        "Repeated validation of the study-conditioned WGAN on the diffusion "
        "finalist split. Calibration is fit on OSDR training draws only. The locked "
        "test is not loaded.\n",
        encoding="utf-8",
    )
    return summary_path


def screen_calibration(config_path: str | Path) -> Path:
    """Select a train-only calibrator using repeated validation, never test."""

    config = load_config(config_path)
    run_output = Path(config["run"]["output_dir"])
    validation = run_output / "evaluation" / "matched_validation"
    base_summary = json.loads((validation / "summary.json").read_text())
    if base_summary.get("locked_test_opened", True):
        raise RuntimeError("Calibration selection requires validation-only inputs")
    output = validation / "calibration_screen"
    summary_path = output / "summary.json"
    if summary_path.exists() or output.exists():
        raise FileExistsError(f"Calibration screen already exists: {output}")
    output.mkdir(parents=True, exist_ok=False)
    adapter, data = _load_model_and_data(config, include_test=False)
    options = config["evaluation"]
    batch_size = int(options["batch_size"])
    train_seeds = tuple(map(int, options["train_generation_seeds"]))
    validation_seeds = tuple(map(int, options["validation_generation_seeds"]))
    priors = tuple(map(float, options["calibration_screen_prior_strengths"]))
    residual_scales = tuple(
        map(float, options["calibration_screen_residual_scales"])
    )
    if not priors or not residual_scales:
        raise ValueError("Calibration screen grid cannot be empty")

    real_fit_parts: list[np.ndarray] = []
    synthetic_fit_parts: list[np.ndarray] = []
    metadata_fit_parts: list[pd.DataFrame] = []
    train_draws = output / "train_draws"
    train_draws.mkdir()
    for seed in train_seeds:
        generated = _generate_scaled(
            adapter,
            data.osdr["train"],
            data.scaler,
            data.maxabs_scale,
            seed=seed,
            batch_size=batch_size,
        )
        np.savez_compressed(
            train_draws / f"seed{seed}.npz",
            scaled_expression=generated,
            sampling_seed=seed,
        )
        real_fit_parts.append(data.osdr_scaled["train"])
        synthetic_fit_parts.append(generated)
        metadata_fit_parts.append(data.osdr["train"].obs.copy())
    real_fit = np.concatenate(real_fit_parts)
    synthetic_fit = np.concatenate(synthetic_fit_parts)
    metadata_fit = pd.concat(metadata_fit_parts, ignore_index=True)

    validation_draws = {
        seed: np.asarray(
            np.load(validation / "raw" / f"seed{seed}" / "synthetic_scaled_expression.npz")[
                "scaled_expression"
            ],
            dtype=np.float32,
        )
        for seed in validation_seeds
    }
    residual_seed = int(options["calibration_residual_seed"])
    all_rows: list[dict[str, object]] = []
    variant_rows: list[dict[str, object]] = []
    calibrators: dict[float, PositiveResidualCalibrator] = {}
    minimum = float(options["minimum_repeat_pass_fraction"])
    for prior in priors:
        calibrator = PositiveResidualCalibrator(
            ("accession", "tissue"),
            prior,
            0.0,
            clip_nonnegative=True,
            noise_group_columns=("accession", "tissue", "condition"),
        ).fit(real_fit, synthetic_fit, metadata_fit)
        calibrators[prior] = calibrator
        for residual_scale in residual_scales:
            calibrator.residual_scale = residual_scale
            variant = f"prior_{prior:g}_residual_{residual_scale:g}"
            rows: list[dict[str, object]] = []
            for index, seed in enumerate(validation_seeds):
                synthetic = calibrator.apply(
                    validation_draws[seed],
                    data.osdr["validation"].obs,
                    seed=residual_seed + index,
                )
                fidelity = generated_quality(
                    data.osdr_scaled["validation"],
                    synthetic,
                    max_pr_samples=len(synthetic),
                )
                memorization = memorization_metrics(
                    data.osdr_scaled["train"],
                    synthetic,
                    max_samples=max(len(synthetic), 50),
                    seed=int(config["run"]["seed"]),
                )
                selection = fidelity_selection(fidelity, memorization)
                effect = _condition_effect(
                    data.osdr_scaled["validation"],
                    synthetic,
                    data.osdr["validation"].obs["condition"].astype(str).to_numpy(),
                )
                effect_gate = conditional_effect_selection(effect)
                row = {
                    "calibration_variant": variant,
                    "sampling_seed": seed,
                    "correlation": fidelity["correlation_matrix_agreement"],
                    "correlation_minimum": selection["fidelity_gate"][
                        "requirements"
                    ]["correlation_matrix_agreement"]["minimum"],
                    "precision": fidelity["precision"],
                    "recall": fidelity["recall"],
                    "f1": fidelity["f1"],
                    "adversarial_accuracy": fidelity["adversarial_accuracy"],
                    "frechet_ratio": fidelity["frechet_ratio_to_real_split_p95"],
                    "negative_fraction": float(np.mean(synthetic < 0)),
                    "diversity_pass": selection["diversity_gate"]["passed"],
                    "memorization_pass": selection["memorization_gate"]["passed"],
                    "fidelity_pass": selection["eligible_for_model_selection"],
                    "condition_delta_correlation": effect["delta_correlation"],
                    "condition_direction_agreement": effect["direction_agreement"],
                    "condition_effect_pass": effect_gate["passed"],
                }
                rows.append(row)
                all_rows.append(row)
            table = pd.DataFrame(rows)
            metrics = _metric_repeat_summary(table)
            fidelity_stable = all(
                float(value["pass_fraction"]) >= minimum
                for value in metrics.values()
            ) and bool(
                table["diversity_pass"].all()
                and table["memorization_pass"].all()
            )
            condition_stable = bool(
                table["condition_effect_pass"].mean() >= minimum
            )
            variant_rows.append(
                {
                    "calibration_variant": variant,
                    "prior_strength": prior,
                    "residual_scale": residual_scale,
                    "mean_correlation": table["correlation"].mean(),
                    "correlation_pass_fraction": metrics["correlation"][
                        "pass_fraction"
                    ],
                    "mean_precision": table["precision"].mean(),
                    "mean_recall": table["recall"].mean(),
                    "mean_f1": table["f1"].mean(),
                    "mean_adversarial_accuracy": table[
                        "adversarial_accuracy"
                    ].mean(),
                    "adversarial_accuracy_pass_fraction": metrics[
                        "adversarial_accuracy"
                    ]["pass_fraction"],
                    "mean_frechet_ratio": table["frechet_ratio"].mean(),
                    "fidelity_stable": fidelity_stable,
                    "condition_effect_pass_fraction": table[
                        "condition_effect_pass"
                    ].mean(),
                    "condition_stable": condition_stable,
                    "candidate_pass": bool(fidelity_stable and condition_stable),
                }
            )

    repeats = pd.DataFrame(all_rows)
    variants = pd.DataFrame(variant_rows)
    repeats.to_csv(output / "repeat_metrics.tsv", sep="\t", index=False)
    variants.to_csv(output / "variant_summary.tsv", sep="\t", index=False)
    eligible = variants.loc[variants["candidate_pass"]].copy()
    selected: dict[str, object] | None = None
    if not eligible.empty:
        eligible["aa_distance_from_chance"] = (
            eligible["mean_adversarial_accuracy"] - 0.5
        ).abs()
        eligible = eligible.sort_values(
            ["aa_distance_from_chance", "mean_correlation"],
            ascending=[True, False],
            kind="stable",
        )
        selected = eligible.iloc[0].to_dict()
        selected_prior = float(selected["prior_strength"])
        selected_calibrator = calibrators[selected_prior]
        selected_calibrator.residual_scale = float(selected["residual_scale"])
        selected_directory = output / "selected_calibrator"
        selected_calibrator.save(selected_directory)
        selected["calibrator"] = str(selected_directory)

    summary = {
        "status": "complete",
        "split": "validation",
        "locked_test_opened": False,
        "grid": {
            "prior_strengths": list(priors),
            "residual_scales": list(residual_scales),
            "generation_seeds": list(validation_seeds),
        },
        "minimum_repeat_pass_fraction": minimum,
        "eligible_variants": int(variants["candidate_pass"].sum()),
        "validation_candidate_pass": selected is not None,
        "selected": selected,
        "selection_rule": (
            "First require every fidelity metric and pooled condition recovery to "
            "pass independently in the declared repeat fraction. Among eligible "
            "variants only, minimize mean |AA-0.5| and then maximize mean Corr."
        ),
        "test_loaded": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    finalist_path = validation / "finalist_selection.json"
    finalist_path.write_text(
        json.dumps(
            {
                "validation_candidate_pass": selected is not None,
                "selected": selected,
                "selection_summary": str(summary_path),
                "locked_test_opened": False,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return summary_path


def _require_test_unlock(unlock_test: bool) -> None:
    if not unlock_test:
        raise PermissionError("Matched WGAN test evaluation requires --unlock-test")


def evaluate_test(config_path: str | Path, *, unlock_test: bool = False) -> Path:
    _require_test_unlock(unlock_test)
    config = load_config(config_path)
    run_output = Path(config["run"]["output_dir"])
    validation_path = (
        run_output
        / "evaluation"
        / "matched_validation"
        / "finalist_selection.json"
    )
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if not validation.get("validation_candidate_pass", False):
        raise RuntimeError("Validation candidate did not pass; locked test stays closed")
    output = run_output / "evaluation" / "final_locked_test"
    summary_path = output / "summary.json"
    if summary_path.exists() or output.exists():
        raise FileExistsError(f"Refusing to overwrite locked test output: {output}")
    adapter, data = _load_model_and_data(config, include_test=True)
    output.mkdir(parents=True, exist_ok=False)
    options = config["evaluation"]
    calibrator = PositiveResidualCalibrator.load(
        validation["selected"]["calibrator"]
    )
    batch_size = int(options["batch_size"])
    residual_seed = int(options["calibration_residual_seed"]) + 10_000
    rows: list[dict[str, object]] = []
    summaries: dict[str, object] = {}
    for index, seed in enumerate(map(int, options["test_generation_seeds"])):
        label = f"seed{seed}"
        raw = _generate_scaled(
            adapter,
            data.osdr["test"],
            data.scaler,
            data.maxabs_scale,
            seed=seed,
            batch_size=batch_size,
        )
        calibrated = calibrator.apply(
            raw, data.osdr["test"].obs, seed=residual_seed + index
        )
        summary, row = _repeat_result(
            label=label,
            split="test",
            seed=seed,
            real=data.osdr_scaled["test"],
            synthetic=calibrated,
            samples=data.osdr["test"].obs,
            genes=data.genes,
            train_real=data.osdr_scaled["train"],
            run_seed=int(config["run"]["seed"]),
            output=output / label,
        )
        rows.append(row)
        summaries[label] = summary

    table = pd.DataFrame(rows)
    table.to_csv(output / "repeat_metrics.tsv", sep="\t", index=False)
    minimum = float(options["minimum_repeat_pass_fraction"])
    metrics, fidelity, condition, muscle = _stability(table, minimum)

    train_seed = int(options["train_generation_seeds"][0])
    synthetic_train = _generate_scaled(
        adapter,
        data.osdr["train"],
        data.scaler,
        data.maxabs_scale,
        seed=train_seed,
        batch_size=batch_size,
    )
    synthetic_train = calibrator.apply(
        synthetic_train, data.osdr["train"].obs, seed=residual_seed + 100
    )
    utility = classifier_utility(
        data.osdr_scaled["train"],
        data.osdr["train"].obs["condition"].astype(str).to_numpy(),
        data.osdr_scaled["test"],
        data.osdr["test"].obs["condition"].astype(str).to_numpy(),
        synthetic_train=synthetic_train,
        synthetic_labels=data.osdr["train"].obs["condition"].astype(str).to_numpy(),
        allow_augmentation=bool(fidelity and condition),
    )
    aggregate = {
        "status": "complete",
        "split": "test",
        "locked_test_opened": True,
        "sampling_seeds": list(map(int, options["test_generation_seeds"])),
        "minimum_repeat_pass_fraction": minimum,
        "metric_repeat_stability": metrics,
        "repeated_test_fidelity": fidelity,
        "repeated_pooled_condition_effect": condition,
        "skeletal_muscle_diagnostic": muscle,
        "broad_finalist_pass": bool(fidelity and condition),
        "classifier_utility": utility,
        "repeat_summaries": summaries,
        "device": adapter.device_summary(),
        "limitations": [
            "The split measures within-study interpolation, not unseen-study transfer.",
            "The test was opened once after WGAN validation selection.",
        ],
    }
    summary_path.write_text(
        json.dumps(aggregate, indent=2) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# Final matched WGAN locked test\n\n"
        "One-time test evaluation of the fixed study-conditioned WGAN. Metrics are "
        "gated independently and no composite score is used.\n",
        encoding="utf-8",
    )
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=(
            "train",
            "evaluate-validation",
            "screen-calibration",
            "evaluate-test",
        ),
    )
    parser.add_argument(
        "--config",
        default="configs/generative/wgan_matched_study_conditioned.yaml",
    )
    parser.add_argument("--unlock-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "train":
        path = train(args.config)
    elif args.command == "evaluate-validation":
        path = evaluate_validation(args.config)
    elif args.command == "screen-calibration":
        path = screen_calibration(args.config)
    else:
        path = evaluate_test(args.config, unlock_test=args.unlock_test)
    print(path)


if __name__ == "__main__":
    main()
