"""Train staged residual adapters on top of a frozen ARCHS4 ModelDDIM."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import random
import shutil
import time
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler
import yaml

from .factorized_adapter import (
    FactorizedAdapterDDIM,
    FactorizedSchema,
    build_factorized_schema,
    encode_factorized_labels,
    load_factorized_role,
)
from .factorized_config import load_factorized_config
from .train import _atomic_torch_save, _seed_everything
from .upstream import (
    EMA,
    antithetic_timesteps,
    model_config,
    noise_estimation_loss,
    quadratic_beta_schedule,
    upstream_model_class,
    verify_source,
)


FORMAT = "nasa_mouse_lacan_factorized_adapter_v1"
CHECKPOINT_FORMAT = "nasa_mouse_lacan_factorized_adapter_checkpoint_v1"


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _base_model(
    pretrained_path: str | Path,
    model_options: dict[str, Any],
    expression_dim: int,
) -> tuple[torch.nn.Module, dict[str, Any], list[str]]:
    payload = torch.load(pretrained_path, map_location="cpu", weights_only=False)
    metadata = payload.get("metadata", {})
    classes = list(map(str, metadata.get("classes", [])))
    if not classes:
        raise ValueError("Pretrained ARCHS4 checkpoint has no tissue classes")
    base = upstream_model_class()(
        model_config(
            expression_dim=expression_dim,
            num_classes=len(classes),
            model=model_options,
        )
    )
    source_state = payload.get("ema_state_dict", payload.get("model_state_dict"))
    if source_state is None:
        raise ValueError("Pretrained checkpoint has no model or EMA state")
    base.load_state_dict(source_state)
    return base, payload, classes


def _sampling_weights(samples: pd.DataFrame, columns: list[str]) -> torch.Tensor:
    missing = [column for column in columns if column not in samples]
    if missing:
        raise ValueError(f"Balanced-sampling columns are absent: {missing}")
    keys = samples.loc[:, columns].fillna("unknown").astype(str).agg("||".join, axis=1)
    counts = keys.map(keys.value_counts()).to_numpy(dtype=np.float64)
    weights = 1.0 / np.maximum(counts, 1.0)
    weights /= weights.mean()
    return torch.as_tensor(weights, dtype=torch.double)


def _loader(
    expression: np.ndarray,
    labels: np.ndarray,
    samples: pd.DataFrame,
    *,
    batch_size: int,
    steps: int,
    strata: list[str],
    workers: int,
    seed: int,
) -> DataLoader:
    dataset = TensorDataset(
        torch.from_numpy(np.asarray(expression, dtype=np.float32)),
        torch.from_numpy(np.asarray(labels, dtype=np.int64)),
    )
    generator = torch.Generator().manual_seed(int(seed))
    sampler = WeightedRandomSampler(
        _sampling_weights(samples, strata),
        num_samples=int(batch_size) * int(steps),
        replacement=True,
        generator=generator,
    )
    options: dict[str, Any] = {
        "batch_size": int(batch_size),
        "sampler": sampler,
        "num_workers": int(workers),
        "pin_memory": True,
        "drop_last": True,
    }
    if workers:
        options.update(prefetch_factor=2, persistent_workers=True)
    return DataLoader(dataset, **options)


def _drop_condition_rows(
    labels: torch.Tensor,
    schema: FactorizedSchema,
    probability: float,
) -> torch.Tensor:
    if probability <= 0:
        return labels
    result = labels.clone()
    dropped = torch.rand(len(result), device=result.device) < float(probability)
    bounds = schema.group_slices()["condition"]
    start = schema.base_width + int(bounds.start)
    stop = schema.base_width + int(bounds.stop)
    result[dropped, start:stop] = 0
    return result


def _correlation_structure_loss(
    clean: torch.Tensor, reconstructed: torch.Tensor, gene_indices: torch.Tensor
) -> torch.Tensor:
    """Compare gene-correlation structure across profiles in one batch."""

    if len(clean) < 3:
        return reconstructed.sum() * 0.0
    expected = clean.index_select(1, gene_indices).float()
    observed = reconstructed.index_select(1, gene_indices).float()

    def correlation(values: torch.Tensor) -> torch.Tensor:
        centered = values - values.mean(dim=0, keepdim=True)
        norms = torch.linalg.vector_norm(centered, dim=0).clamp_min(1e-6)
        return centered.T.mm(centered) / torch.outer(norms, norms)

    expected_correlation = correlation(expected).detach()
    observed_correlation = correlation(observed)
    upper = torch.triu_indices(
        observed_correlation.shape[0],
        observed_correlation.shape[1],
        offset=1,
        device=observed.device,
    )
    return torch.mean(
        (
            observed_correlation[upper[0], upper[1]]
            - expected_correlation[upper[0], upper[1]]
        ).square()
    )


def _regularized_noise_loss(
    model: FactorizedAdapterDDIM,
    clean: torch.Tensor,
    timesteps: torch.Tensor,
    noise: torch.Tensor,
    betas: torch.Tensor,
    labels: torch.Tensor,
    options: dict[str, object] | None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if not options:
        loss, error = noise_estimation_loss(
            model, clean, timesteps, noise, betas, labels
        )
        return loss, error, loss.detach() * 0.0
    alpha = (1 - betas).cumprod(dim=0).index_select(0, timesteps).view(-1, 1)
    noisy = clean * alpha.sqrt() + noise * (1.0 - alpha).sqrt()
    prediction = model(noisy, timesteps, labels)
    residual = noise - prediction
    noise_loss = residual.square().sum(dim=1).mean(dim=0)
    error = residual.abs().mean(dim=1).mean(dim=0)
    eligible = timesteps <= int(options["max_timestep"])
    if int(eligible.sum().item()) < 3:
        return noise_loss, error, noise_loss.detach() * 0.0
    reconstructed = (
        noisy[eligible]
        - prediction[eligible] * (1.0 - alpha[eligible]).sqrt()
    ) / alpha[eligible].sqrt().clamp_min(1e-6)
    gene_count = min(int(options["genes"]), clean.shape[1])
    gene_indices = torch.randperm(clean.shape[1], device=clean.device)[:gene_count]
    correlation_loss = _correlation_structure_loss(
        clean[eligible], reconstructed, gene_indices
    )
    return (
        noise_loss + float(options["weight"]) * correlation_loss,
        error,
        correlation_loss,
    )


@torch.no_grad()
def _validation_loss(
    model: FactorizedAdapterDDIM,
    expression: np.ndarray,
    labels: np.ndarray,
    betas: torch.Tensor,
    *,
    batch_size: int,
    seed: int,
    device: torch.device,
) -> dict[str, float]:
    model.eval()
    generator = torch.Generator(device=device).manual_seed(int(seed))
    total_loss = 0.0
    total_error = 0.0
    profiles = 0
    for start in range(0, len(expression), int(batch_size)):
        end = min(start + int(batch_size), len(expression))
        clean = torch.as_tensor(expression[start:end], device=device)
        condition = torch.as_tensor(labels[start:end], device=device)
        noise = torch.randn(clean.shape, generator=generator, device=device)
        timesteps = torch.randint(
            0,
            len(betas),
            (len(clean),),
            generator=generator,
            device=device,
        )
        with torch.autocast("cuda", dtype=torch.float16, enabled=True):
            loss, error = noise_estimation_loss(
                model, clean, timesteps, noise, betas, condition
            )
        total_loss += float(loss.cpu()) * len(clean)
        total_error += float(error.cpu()) * len(clean)
        profiles += len(clean)
    model.train()
    return {
        "validation_loss": total_loss / profiles,
        "validation_noise_absolute_error": total_error / profiles,
    }


def _checkpoint_payload(
    model: FactorizedAdapterDDIM,
    ema: EMA,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LRScheduler,
    scaler: torch.amp.GradScaler,
    *,
    stage: str,
    step: int,
    history: list[dict[str, object]],
    metadata: dict[str, object],
) -> dict[str, object]:
    return {
        "format": CHECKPOINT_FORMAT,
        "adapter_state_dict": model.adapter_state_dict(),
        "ema_state_dict": {
            key: value.detach().cpu() for key, value in ema.state_dict().items()
        },
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "scaler_state_dict": scaler.state_dict(),
        "stage": stage,
        "step": int(step),
        "history": history,
        "metadata": metadata,
        "torch_rng_state": torch.get_rng_state(),
        "cuda_rng_state": torch.cuda.get_rng_state_all(),
        "numpy_rng_state": np.random.get_state(),
        "python_rng_state": random.getstate(),
    }


def _restore_rng(payload: dict[str, Any]) -> None:
    torch.set_rng_state(payload["torch_rng_state"].cpu())
    torch.cuda.set_rng_state_all(
        [value.cpu() for value in payload.get("cuda_rng_state", [])]
    )
    np.random.set_state(payload["numpy_rng_state"])
    random.setstate(payload["python_rng_state"])


def _train_stage(
    model: FactorizedAdapterDDIM,
    schema: FactorizedSchema,
    train: dict[str, object],
    validation: dict[str, object],
    labels: np.ndarray,
    validation_labels: np.ndarray,
    betas: torch.Tensor,
    stage: str,
    options: dict[str, Any],
    common: dict[str, Any],
    output: Path,
    metadata: dict[str, object],
    *,
    seed: int,
    restart: bool,
    device: torch.device,
) -> list[dict[str, object]]:
    artifact = output / f"{stage}_adapter.pt"
    checkpoint = output / "checkpoints" / f"{stage}_latest.pt"
    model.set_trainable_groups([stage])
    if artifact.exists() and not restart:
        payload = torch.load(artifact, map_location="cpu", weights_only=False)
        model.load_adapter_state_dict(payload["adapter_state_dict"])
        print(f"[factorized-ddim] loaded completed stage={stage}", flush=True)
        return list(payload.get("history", []))

    steps = int(options["steps"])
    batch_size = int(common["batch_size"])
    loader = _loader(
        train["expression"],
        labels,
        train["samples"],
        batch_size=batch_size,
        steps=steps,
        strata=list(map(str, options["sampling_strata"])),
        workers=int(common.get("num_workers", 0)),
        seed=seed,
    )
    parameters = [value for value in model.parameters() if value.requires_grad]
    optimizer = torch.optim.Adam(
        parameters,
        lr=float(options["learning_rate"]),
        betas=(float(common.get("beta1", 0.9)), float(common.get("beta2", 0.999))),
        eps=float(common.get("epsilon", 1e-8)),
    )
    warmup_steps = min(int(options.get("warmup_steps", 100)), steps - 1)
    scheduler = torch.optim.lr_scheduler.OneCycleLR(
        optimizer,
        max_lr=float(options["learning_rate"]),
        total_steps=steps,
        pct_start=max(warmup_steps / steps, 1.0 / steps),
    )
    scaler = torch.amp.GradScaler(
        "cuda", enabled=bool(common.get("amp", True))
    )
    ema = EMA(model, float(common.get("ema_decay", 0.999)))
    start_step = 0
    history: list[dict[str, object]] = []
    if checkpoint.exists() and not restart:
        payload = torch.load(checkpoint, map_location=device, weights_only=False)
        if payload.get("format") != CHECKPOINT_FORMAT or payload.get("stage") != stage:
            raise ValueError(f"Incompatible factorized checkpoint: {checkpoint}")
        model.load_adapter_state_dict(payload["adapter_state_dict"])
        optimizer.load_state_dict(payload["optimizer_state_dict"])
        scheduler.load_state_dict(payload["scheduler_state_dict"])
        scaler.load_state_dict(payload["scaler_state_dict"])
        ema.load_state_dict(payload["ema_state_dict"], device)
        start_step = int(payload["step"])
        history = list(payload.get("history", []))
        _restore_rng(payload)

    condition_dropout = (
        float(common.get("condition_dropout", 0.0)) if stage == "condition" else 0.0
    )
    print(
        f"[factorized-ddim] stage={stage} device={torch.cuda.get_device_name(device)} "
        f"trainable={model.trainable_parameter_count():,} profiles={len(train['expression'])} "
        f"batch={batch_size} steps={steps} lr={options['learning_rate']}",
        flush=True,
    )
    model.train()
    torch.cuda.reset_peak_memory_stats(device)
    interval_loss = 0.0
    interval_error = 0.0
    interval_correlation_loss = 0.0
    interval_profiles = 0
    interval_started = time.time()
    log_every = int(common.get("log_every_steps", 100))
    validate_every = int(common.get("validate_every_steps", 500))
    checkpoint_every = int(common.get("checkpoint_every_steps", 500))
    for step_index, (clean, condition) in enumerate(loader, start=1):
        if step_index <= start_step:
            continue
        clean = clean.to(device, non_blocking=True)
        condition = condition.to(device, non_blocking=True)
        condition = _drop_condition_rows(
            condition, schema, condition_dropout
        )
        noise = torch.randn_like(clean)
        timesteps = antithetic_timesteps(len(clean), len(betas), device)
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast(
            "cuda", dtype=torch.float16, enabled=bool(common.get("amp", True))
        ):
            loss, error, correlation_loss = _regularized_noise_loss(
                model,
                clean,
                timesteps,
                noise,
                betas,
                condition,
                options.get("correlation_regularization"),
            )
        if not torch.isfinite(loss):
            raise FloatingPointError(f"Non-finite {stage} loss at step {step_index}")
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        scheduler.step()
        ema.update(model)
        interval_loss += float(loss.detach().cpu()) * len(clean)
        interval_error += float(error.detach().cpu()) * len(clean)
        interval_correlation_loss += float(correlation_loss.detach().cpu()) * len(
            clean
        )
        interval_profiles += len(clean)

        should_log = step_index == 1 or step_index % log_every == 0 or step_index == steps
        if should_log:
            row: dict[str, object] = {
                "stage": stage,
                "step": step_index,
                "loss": interval_loss / interval_profiles,
                "noise_absolute_error": interval_error / interval_profiles,
                "correlation_structure_loss": (
                    interval_correlation_loss / interval_profiles
                ),
                "learning_rate": float(optimizer.param_groups[0]["lr"]),
                "interval_seconds": float(time.time() - interval_started),
                "cuda_peak_memory_gb": float(
                    torch.cuda.max_memory_allocated(device) / 1024**3
                ),
            }
            if step_index % validate_every == 0 or step_index == steps:
                row.update(
                    _validation_loss(
                        model,
                        validation["expression"],
                        validation_labels,
                        betas,
                        batch_size=batch_size,
                        seed=seed + step_index,
                        device=device,
                    )
                )
            history.append(row)
            print(
                f"[factorized-ddim] stage={stage} step={step_index}/{steps} "
                f"loss={row['loss']:.5f} mae={row['noise_absolute_error']:.5f} "
                f"lr={row['learning_rate']:.3g}",
                flush=True,
            )
            pd.DataFrame(history).to_csv(
                output / f"{stage}_training_history.tsv", sep="\t", index=False
            )
            interval_loss = 0.0
            interval_error = 0.0
            interval_correlation_loss = 0.0
            interval_profiles = 0
            interval_started = time.time()
        if step_index % checkpoint_every == 0 and step_index < steps:
            _atomic_torch_save(
                _checkpoint_payload(
                    model,
                    ema,
                    optimizer,
                    scheduler,
                    scaler,
                    stage=stage,
                    step=step_index,
                    history=history,
                    metadata=metadata,
                ),
                checkpoint,
            )

    ema.copy_to(model)
    _atomic_torch_save(
        {
            "format": FORMAT,
            "stage": stage,
            "adapter_state_dict": model.adapter_state_dict(),
            "history": history,
            "metadata": metadata,
        },
        artifact,
    )
    checkpoint.unlink(missing_ok=True)
    return history


def train_factorized(config_path: str | Path, *, restart: bool = False) -> Path:
    config = load_factorized_config(config_path)
    run = config["run"]
    data_options = config["data"]
    common = config["training"]
    output = Path(run["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    if restart:
        for path in (
            output / "domain_adapter.pt",
            output / "condition_adapter.pt",
            output / "model.pt",
        ):
            path.unlink(missing_ok=True)
        shutil.rmtree(output / "checkpoints", ignore_errors=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Factorized paper ModelDDIM adapters require CUDA")
    seed = int(run["seed"])
    _seed_everything(seed)
    torch.backends.cudnn.benchmark = True
    source_manifest = verify_source(run.get("source_root"))
    train = load_factorized_role(
        data_options["prepared_h5"], data_options["samples_tsv"], "train"
    )
    validation = load_factorized_role(
        data_options["prepared_h5"], data_options["samples_tsv"], "validation"
    )
    base, pretrained_payload, classes = _base_model(
        data_options["pretrained_model"], config["model"], len(train["genes"])
    )
    conditioning = config.get("conditioning", {})
    schema = build_factorized_schema(
        train["samples"],
        classes,
        include_study=bool(conditioning.get("study", False)),
        include_material_type=bool(conditioning.get("material_type", False)),
    )
    train_labels = encode_factorized_labels(train["samples"], schema)
    validation_labels = encode_factorized_labels(validation["samples"], schema)
    adapter_options = config.get("adapter", {})
    model = FactorizedAdapterDDIM(
        base,
        schema,
        domain_lora_rank=int(adapter_options.get("domain_lora_rank", 0)),
        domain_lora_alpha=float(adapter_options.get("domain_lora_alpha", 1.0)),
    ).to(device)
    initial_model = str(adapter_options.get("initial_model", ""))
    if initial_model:
        initial_payload = torch.load(
            initial_model, map_location="cpu", weights_only=False
        )
        if initial_payload.get("format") != FORMAT:
            raise ValueError("adapter.initial_model has an incompatible format")
        if initial_payload.get("metadata", {}).get("schema") != schema.as_dict():
            raise ValueError("adapter.initial_model uses a different factorized schema")
        model.load_adapter_state_dict(initial_payload["adapter_state_dict"])
    betas = quadratic_beta_schedule(
        beta_start=float(config["model"]["beta_start"]),
        beta_end=float(config["model"]["beta_end"]),
        timesteps=int(config["model"]["diffusion_timesteps"]),
    ).to(device)
    metadata: dict[str, object] = {
        "config": str(Path(config_path).resolve()),
        "config_sha256": _sha256(config_path),
        "source": source_manifest,
        "prepared_h5": str(Path(data_options["prepared_h5"]).resolve()),
        "samples_tsv": str(Path(data_options["samples_tsv"]).resolve()),
        "pretrained_model": str(Path(data_options["pretrained_model"]).resolve()),
        "pretrained_model_sha256": _sha256(data_options["pretrained_model"]),
        "pretrained_epoch": int(pretrained_payload.get("epoch", 0)),
        "pretrained_global_step": int(pretrained_payload.get("global_step", 0)),
        "pretrained_state": "ema_state_dict",
        "initial_adapter_model": (
            str(Path(initial_model).resolve()) if initial_model else ""
        ),
        "initial_adapter_model_sha256": (
            _sha256(initial_model) if initial_model else ""
        ),
        "genes": train["genes"],
        "schema": schema.as_dict(),
        "train_profiles": int(len(train["expression"])),
        "validation_profiles": int(len(validation["expression"])),
        "locked_test_opened": False,
        "device": torch.cuda.get_device_name(device),
        "torch_version": str(torch.__version__),
        "implementation_contract": {
            "pinned_base_model": True,
            "base_parameters_frozen": True,
            "zero_initialization_preserves_base_function": True,
            "factorized_tissue_condition_covariates": True,
            "study_conditioning": bool(conditioning.get("study", False)),
            "material_type_conditioning": bool(
                conditioning.get("material_type", False)
            ),
            "domain_lora_rank": int(
                adapter_options.get("domain_lora_rank", 0)
            ),
            "domain_lora_alpha": float(
                adapter_options.get("domain_lora_alpha", 1.0)
            ),
            "domain_then_condition_staging": True,
            "condition_dropout_for_guidance": float(
                common.get("condition_dropout", 0.0)
            ),
            "test_expression_loaded": False,
        },
    }
    resolved = dict(config)
    resolved.pop("_config_path", None)
    resolved["resolved"] = metadata
    (output / "resolved_config.yaml").write_text(
        yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8"
    )

    histories: dict[str, list[dict[str, object]]] = {}
    for offset, stage in enumerate(("domain", "condition")):
        histories[stage] = _train_stage(
            model,
            schema,
            train,
            validation,
            train_labels,
            validation_labels,
            betas,
            stage,
            common["stages"][stage],
            common,
            output,
            metadata,
            seed=seed + offset * 100_000,
            restart=restart,
            device=device,
        )
    final_path = output / "model.pt"
    _atomic_torch_save(
        {
            "format": FORMAT,
            "adapter_state_dict": model.adapter_state_dict(),
            "metadata": metadata,
            "history": histories,
        },
        final_path,
    )
    summary = {
        "status": "complete",
        "model": str(final_path),
        "device": metadata["device"],
        "base_parameter_count": int(
            sum(parameter.numel() for parameter in model.base_model.parameters())
        ),
        "adapter_parameter_count": int(
            sum(
                parameter.numel()
                for name, parameter in model.named_parameters()
                if not name.startswith("base_model.")
            )
        ),
        "profiles": {
            "train": int(len(train["expression"])),
            "validation": int(len(validation["expression"])),
        },
        "stages": {
            stage: {
                "steps": int(common["stages"][stage]["steps"]),
                "learning_rate": float(
                    common["stages"][stage]["learning_rate"]
                ),
                "final_loss": float(histories[stage][-1]["loss"]),
            }
            for stage in ("domain", "condition")
        },
        "locked_test_opened": False,
    }
    (output / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return final_path
