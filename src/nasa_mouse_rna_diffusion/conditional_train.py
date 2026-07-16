"""Train the unmodified upstream ModelDDIM on API-derived OSDR conditions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
import yaml

from .conditional_config import load_conditional_config
from .conditional_data import load_conditional_prepared, prepare_conditional
from .train import (
    _atomic_torch_save,
    _checkpoint_payload,
    _restore_rng,
    _seed_everything,
)
from .upstream import (
    EMA,
    antithetic_timesteps,
    model_config,
    noise_estimation_loss,
    quadratic_beta_schedule,
    upstream_model_class,
    verify_source,
)


CHECKPOINT_FORMAT = "nasa_mouse_lacan_conditional_osdr_checkpoint_v1"
MODEL_FORMAT = "nasa_mouse_lacan_conditional_osdr_model_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checkpoint(
    *,
    model,
    optimizer,
    scheduler,
    scaler,
    ema,
    epoch: int,
    step: int,
    history: list[dict[str, Any]],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    return _checkpoint_payload(
        model=model,
        optimizer=optimizer,
        scheduler=scheduler,
        scaler=scaler,
        ema=ema,
        epoch=epoch,
        step=step,
        history=history,
        metadata=metadata,
        format_name=CHECKPOINT_FORMAT,
    )


def _expanded_condition_state(
    template: dict[str, torch.Tensor],
    source: dict[str, torch.Tensor],
    *,
    old_classes: list[str],
    new_classes: list[str],
    embedding_dim: int,
) -> tuple[dict[str, torch.Tensor], dict[str, Any]]:
    """Map tissue-class weights into tissue/reference columns of a larger model."""

    expanded = {
        key: value.detach().cpu().clone() for key, value in template.items()
    }
    for key, value in source.items():
        if (
            key != "y_emb.weight"
            and key in expanded
            and tuple(value.shape) == tuple(expanded[key].shape)
        ):
            expanded[key] = value.detach().cpu().clone()

    new_map = {label: index for index, label in enumerate(new_classes)}
    class_mapping = {
        old_index: new_map[f"tissue={tissue}||condition=reference"]
        for old_index, tissue in enumerate(old_classes)
    }
    if "y_emb.weight" in source:
        for old_index, new_index in class_mapping.items():
            expanded["y_emb.weight"][new_index] = source["y_emb.weight"][
                old_index
            ].cpu()
    condition_keys: list[str] = []
    old_width = len(old_classes) * int(embedding_dim)
    new_width = len(new_classes) * int(embedding_dim)
    for key, old_value in source.items():
        if key == "y_emb.weight" or key not in expanded or old_value.ndim != 2:
            continue
        new_value = expanded[key]
        if old_value.shape[1] - old_width != new_value.shape[1] - new_width:
            continue
        base = int(old_value.shape[1] - old_width)
        new_value[:, :base] = old_value[:, :base].cpu()
        for old_index, new_index in class_mapping.items():
            old_start = base + old_index * int(embedding_dim)
            new_start = base + new_index * int(embedding_dim)
            new_value[
                :, new_start : new_start + int(embedding_dim)
            ] = old_value[
                :, old_start : old_start + int(embedding_dim)
            ].cpu()
        condition_keys.append(key)
    return expanded, {
        "old_classes": old_classes,
        "new_classes": new_classes,
        "mapped_classes": len(class_mapping),
        "mapped_condition_weight_keys": sorted(condition_keys),
        "new_flt_gc_condition_columns_randomly_initialized": True,
    }


def train_conditional(config_path: str | Path, *, restart: bool = False) -> Path:
    config = load_conditional_config(config_path)
    run = config["run"]
    training = config["training"]
    model_options = config["model"]
    regime = str(training["regime"])
    output = Path(run["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output / "checkpoints/latest.pt"
    if restart and checkpoint_path.exists():
        checkpoint_path.unlink()

    source_manifest = verify_source(run["source_root"])
    prepared_path = prepare_conditional(config_path)
    prepared = load_conditional_prepared(prepared_path)
    train_expression = torch.from_numpy(prepared["train"]["expression"])
    train_class = torch.from_numpy(prepared["train"]["class_index"])
    train_labels = torch.nn.functional.one_hot(
        train_class, num_classes=len(prepared["classes"])
    ).long()
    dataset = TensorDataset(train_expression, train_labels)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("The 228M-parameter upstream ModelDDIM requires CUDA")

    seed = int(run["seed"])
    _seed_everything(seed)
    torch.backends.cudnn.benchmark = True
    workers = int(training["num_workers"])
    loader_options: dict[str, Any] = {
        "batch_size": int(training["batch_size"]),
        "shuffle": True,
        "num_workers": workers,
        "pin_memory": True,
        "drop_last": False,
    }
    if workers:
        loader_options.update(prefetch_factor=2, persistent_workers=True)
    loader = DataLoader(dataset, **loader_options)
    namespace = model_config(
        expression_dim=len(prepared["genes"]),
        num_classes=len(prepared["classes"]),
        model=model_options,
    )
    model = upstream_model_class()(namespace).to(device)
    parameter_count = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    pretraining_audit: dict[str, Any] = {}
    pretrained_payload: dict[str, Any] | None = None
    if regime == "archs4_pretrain_osdr_finetune":
        pretrained_path = Path(training["pretrained_model"])
        pretrained_payload = torch.load(
            pretrained_path, map_location="cpu", weights_only=False
        )
        old_classes = list(map(str, pretrained_payload["metadata"]["classes"]))
        expanded, pretraining_audit = _expanded_condition_state(
            model.state_dict(),
            pretrained_payload["model_state_dict"],
            old_classes=old_classes,
            new_classes=prepared["classes"],
            embedding_dim=int(model_options["tissue_embedding_dim"]),
        )
        model.load_state_dict(expanded)
        pretraining_audit.update(
            {
                "pretrained_model": str(pretrained_path.resolve()),
                "pretrained_model_sha256": _sha256(pretrained_path),
                "pretrained_epoch": int(pretrained_payload.get("epoch", 0)),
                "pretrained_global_step": int(
                    pretrained_payload.get("global_step", 0)
                ),
            }
        )
        stage_epochs = int(training["finetune_epochs"])
        stage_learning_rate = float(training["finetune_learning_rate"])
        peak_step = int(training["finetune_one_cycle_peak_step"])
    else:
        stage_epochs = int(training["epochs"])
        stage_learning_rate = float(training["learning_rate"])
        peak_step = int(training["one_cycle_peak_step"])
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=stage_learning_rate,
        betas=(float(training["beta1"]), float(training["beta2"])),
        eps=float(training["epsilon"]),
        weight_decay=float(training["weight_decay"]),
    )
    total_steps = len(loader) * stage_epochs
    if peak_step >= total_steps:
        raise ValueError("OneCycle peak step must occur before the final optimizer step")
    scheduler = torch.optim.lr_scheduler.OneCycleLR(
        optimizer,
        max_lr=stage_learning_rate,
        steps_per_epoch=len(loader),
        epochs=stage_epochs,
        pct_start=peak_step / total_steps,
    )
    amp = bool(training["amp"])
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    ema = EMA(model, float(model_options["ema_decay"]))
    if pretrained_payload is not None:
        expanded_ema, _ = _expanded_condition_state(
            model.state_dict(),
            pretrained_payload["ema_state_dict"],
            old_classes=list(map(str, pretrained_payload["metadata"]["classes"])),
            new_classes=prepared["classes"],
            embedding_dim=int(model_options["tissue_embedding_dim"]),
        )
        ema.load_state_dict(expanded_ema, device)
        del expanded_ema, expanded, pretrained_payload
    betas = quadratic_beta_schedule(
        beta_start=float(model_options["beta_start"]),
        beta_end=float(model_options["beta_end"]),
        timesteps=int(model_options["diffusion_timesteps"]),
    ).to(device)
    metadata = {
        "config": str(Path(config_path).resolve()),
        "source": source_manifest,
        "prepared_data": str(Path(prepared_path).resolve()),
        "prepared_data_sha256": json.loads(
            Path(prepared_path).with_suffix(".manifest.json").read_text(encoding="utf-8")
        )["prepared_h5_sha256"],
        "genes": prepared["genes"],
        "classes": prepared["classes"],
        "conditioning_covariates": prepared["conditioning_covariates"],
        "regime": regime,
        "pretraining": pretraining_audit,
        "parameter_count": int(parameter_count),
        "device": torch.cuda.get_device_name(device),
        "torch_version": str(torch.__version__),
        "batches_per_epoch": len(loader),
        "total_optimizer_steps": total_steps,
        "one_cycle_peak_step": peak_step,
        "stage_epochs": stage_epochs,
        "stage_learning_rate": stage_learning_rate,
        "regime": regime,
        "pretraining": pretraining_audit,
        "extension_label": "unmodified_ModelDDIM_NASA_OSDR_condition_extension",
        "implementation_contract": {
            "upstream_model_class": True,
            "upstream_noise_objective": True,
            "antithetic_timesteps": True,
            "parameter_only_ema": True,
            "weighted_sampling": False,
            "gradient_clipping": False,
            "full_transcriptome_tpm_before_landmarks": True,
            "accession_grouped_split": True,
        },
    }
    resolved = dict(config)
    resolved.pop("_config_path", None)
    resolved["resolved"] = metadata
    (output / "resolved_config.yaml").write_text(
        yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8"
    )

    start_epoch = 0
    global_step = 0
    history: list[dict[str, Any]] = []
    if checkpoint_path.exists():
        payload = torch.load(
            checkpoint_path, map_location=device, weights_only=False
        )
        if payload.get("format") != CHECKPOINT_FORMAT:
            raise ValueError(f"Incompatible checkpoint: {checkpoint_path}")
        observed = payload.get("metadata", {})
        for key in ("genes", "classes", "parameter_count", "prepared_data_sha256"):
            if observed.get(key) != metadata.get(key):
                raise ValueError(f"Checkpoint {key} differs from this run")
        model.load_state_dict(payload["model_state_dict"])
        optimizer.load_state_dict(payload["optimizer_state_dict"])
        scheduler.load_state_dict(payload["scheduler_state_dict"])
        scaler.load_state_dict(payload["scaler_state_dict"])
        ema.load_state_dict(payload["ema_state_dict"], device)
        start_epoch = int(payload["epoch"])
        global_step = int(payload["global_step"])
        history = list(payload.get("history", []))
        _restore_rng(payload)
        print(
            f"[conditional-ddim:train] resumed epoch={start_epoch} step={global_step}",
            flush=True,
        )

    print(
        f"[conditional-ddim:train] device={metadata['device']} params={parameter_count:,} "
        f"profiles={len(dataset):,} classes={len(prepared['classes'])} "
        f"batches={len(loader)} stage={regime} epochs={stage_epochs}",
        flush=True,
    )
    torch.cuda.reset_peak_memory_stats(device)
    started = time.time()
    completed_epoch = start_epoch
    try:
        for epoch_index in range(start_epoch, stage_epochs):
            model.train()
            epoch_loss = 0.0
            epoch_error = 0.0
            epoch_profiles = 0
            epoch_started = time.time()
            for clean, labels in loader:
                clean = clean.to(device, non_blocking=True)
                labels = labels.to(device, non_blocking=True)
                batch_size = len(clean)
                noise = torch.randn_like(clean)
                timesteps = antithetic_timesteps(
                    batch_size,
                    int(model_options["diffusion_timesteps"]),
                    device,
                )
                optimizer.zero_grad(set_to_none=False)
                with torch.autocast("cuda", dtype=torch.float16, enabled=amp):
                    loss, absolute_error = noise_estimation_loss(
                        model, clean, timesteps, noise, betas, labels
                    )
                if not torch.isfinite(loss):
                    raise FloatingPointError(
                        f"Non-finite loss at epoch={epoch_index + 1} step={global_step + 1}"
                    )
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
                ema.update(model)
                scheduler.step()
                global_step += 1
                epoch_profiles += batch_size
                epoch_loss += float(loss.detach().cpu()) * batch_size
                epoch_error += float(absolute_error.detach().cpu()) * batch_size
            completed_epoch = epoch_index + 1
            row = {
                "epoch": completed_epoch,
                "global_step": global_step,
                "loss": epoch_loss / epoch_profiles,
                "noise_absolute_error": epoch_error / epoch_profiles,
                "learning_rate": float(optimizer.param_groups[0]["lr"]),
                "epoch_seconds": float(time.time() - epoch_started),
                "cuda_peak_memory_gb": float(
                    torch.cuda.max_memory_allocated(device) / 1024**3
                ),
            }
            history.append(row)
            log_every = int(training["log_every_epochs"])
            if completed_epoch == 1 or completed_epoch % log_every == 0:
                elapsed = time.time() - started
                rate = (completed_epoch - start_epoch) / max(elapsed, 1e-8)
                remaining = (stage_epochs - completed_epoch) / max(
                    rate, 1e-8
                )
                print(
                    f"[conditional-ddim:train] epoch={completed_epoch}/"
                    f"{stage_epochs} loss={row['loss']:.6f} "
                    f"mae={row['noise_absolute_error']:.6f} "
                    f"lr={row['learning_rate']:.8g} eta_hours={remaining / 3600:.2f}",
                    flush=True,
                )
            if completed_epoch % int(training["checkpoint_every_epochs"]) == 0:
                _atomic_torch_save(
                    _checkpoint(
                        model=model,
                        optimizer=optimizer,
                        scheduler=scheduler,
                        scaler=scaler,
                        ema=ema,
                        epoch=completed_epoch,
                        step=global_step,
                        history=history,
                        metadata=metadata,
                    ),
                    checkpoint_path,
                )
                pd.DataFrame(history).to_csv(
                    output / "training_history.tsv", sep="\t", index=False
                )
    except BaseException:
        if completed_epoch > start_epoch:
            _atomic_torch_save(
                _checkpoint(
                    model=model,
                    optimizer=optimizer,
                    scheduler=scheduler,
                    scaler=scaler,
                    ema=ema,
                    epoch=completed_epoch,
                    step=global_step,
                    history=history,
                    metadata=metadata,
                ),
                checkpoint_path,
            )
            pd.DataFrame(history).to_csv(
                output / "training_history.tsv", sep="\t", index=False
            )
        raise

    final_path = output / "model.pt"
    _atomic_torch_save(
        {
            "format": MODEL_FORMAT,
            "model_state_dict": model.state_dict(),
            "ema_state_dict": ema.state_dict(),
            "metadata": metadata,
            "epoch": completed_epoch,
            "global_step": global_step,
        },
        final_path,
    )
    pd.DataFrame(history).to_csv(
        output / "training_history.tsv", sep="\t", index=False
    )
    if not bool(training.get("retain_training_checkpoint", False)):
        checkpoint_path.unlink(missing_ok=True)
        checkpoint_path.parent.rmdir()
    summary = {
        "status": "complete",
        "model": "Lacan et al. upstream ModelDDIM OSDR condition extension",
        "source": source_manifest,
        "run_dir": str(output),
        "model_path": str(final_path),
        "prepared_data": str(prepared_path),
        "device": metadata["device"],
        "parameter_count": parameter_count,
        "classes": prepared["classes"],
        "conditioning_covariates": prepared["conditioning_covariates"],
        "profiles": {
            role: len(prepared[role]["expression"])
            for role in ("train", "validation", "test")
        },
        "epochs": completed_epoch,
        "global_steps": global_step,
        "training_seconds_this_invocation": float(time.time() - started),
        "final_loss": history[-1]["loss"],
        "final_noise_absolute_error": history[-1]["noise_absolute_error"],
        "cuda_peak_memory_gb": float(
            torch.cuda.max_memory_allocated(device) / 1024**3
        ),
        "training_checkpoint_retained": bool(
            training.get("retain_training_checkpoint", False)
        ),
    }
    (output / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# OSDR Conditional Upstream ModelDDIM\n\n"
        "This NASA extension retains the pinned Lacan et al. ModelDDIM architecture, "
        "loss, optimizer, schedule, and 15,000-epoch duration. The substituted data "
        "are API-derived mouse OSDR profiles conditioned jointly on tissue and "
        "flight/ground-control state. It is not a reproduction of the GTEx cohort.\n\n"
        f"Training regime: `{regime}`. ARCHS4 pretraining, when selected, is loaded "
        "from the completed exact tissue checkpoint and transferred with an explicit "
        "class-column mapping before OSDR fine-tuning.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2), flush=True)
    return final_path
