"""Train the pinned Lacan et al. DDIM on paper-matched ARCHS4 mouse data."""

from __future__ import annotations

import json
from pathlib import Path
import random
import time
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
import yaml

from .config import load_config
from .data import load_prepared, prepare
from .upstream import (
    EMA,
    antithetic_timesteps,
    model_config,
    noise_estimation_loss,
    quadratic_beta_schedule,
    upstream_model_class,
    verify_source,
)


def _seed_everything(seed: int) -> None:
    random.seed(int(seed))
    np.random.seed(int(seed))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))


def _atomic_torch_save(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def _checkpoint_payload(
    *,
    model: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LRScheduler,
    scaler: torch.amp.GradScaler,
    ema: EMA,
    epoch: int,
    step: int,
    history: list[dict[str, Any]],
    metadata: dict[str, Any],
    format_name: str = "nasa_mouse_lacan_paper_parity_checkpoint_v1",
) -> dict[str, Any]:
    return {
        "format": str(format_name),
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "scaler_state_dict": scaler.state_dict(),
        "ema_state_dict": ema.state_dict(),
        "epoch": int(epoch),
        "global_step": int(step),
        "history": history,
        "torch_rng_state": torch.get_rng_state(),
        "cuda_rng_state": torch.cuda.get_rng_state_all(),
        "numpy_rng_state": np.random.get_state(),
        "python_rng_state": random.getstate(),
        "metadata": metadata,
    }


def _restore_rng(payload: dict[str, Any]) -> None:
    torch.set_rng_state(payload["torch_rng_state"].cpu())
    if torch.cuda.is_available() and payload.get("cuda_rng_state"):
        torch.cuda.set_rng_state_all(
            [value.cpu() for value in payload["cuda_rng_state"]]
        )
    np.random.set_state(payload["numpy_rng_state"])
    random.setstate(payload["python_rng_state"])


def train(config_path: str | Path, *, restart: bool = False) -> Path:
    config = load_config(config_path)
    run = config["run"]
    training = config["training"]
    model_options = config["model"]
    output = Path(run["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output / "checkpoints/latest.pt"
    if restart and checkpoint_path.exists():
        checkpoint_path.unlink()

    source_manifest = verify_source(run["source_root"])
    prepared_path = prepare(config_path)
    prepared = load_prepared(prepared_path)
    train_expression = torch.from_numpy(prepared["train"]["expression"])
    train_class = torch.from_numpy(prepared["train"]["class_index"])
    train_labels = torch.nn.functional.one_hot(
        train_class, num_classes=len(prepared["classes"])
    ).long()
    dataset = TensorDataset(train_expression, train_labels)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("The paper-parity 228M-parameter DDIM requires a CUDA GPU")

    seed = int(run["seed"])
    _seed_everything(seed)
    torch.backends.cudnn.benchmark = True
    loader = DataLoader(
        dataset,
        batch_size=int(training["batch_size"]),
        shuffle=True,
        num_workers=int(training["num_workers"]),
        pin_memory=True,
        prefetch_factor=2,
        persistent_workers=int(training["num_workers"]) > 0,
        drop_last=False,
    )
    namespace = model_config(
        expression_dim=len(prepared["genes"]),
        num_classes=len(prepared["classes"]),
        model=model_options,
    )
    model = upstream_model_class()(namespace).to(device)
    parameter_count = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(training["learning_rate"]),
        betas=(float(training["beta1"]), float(training["beta2"])),
        eps=float(training["epsilon"]),
        weight_decay=float(training["weight_decay"]),
    )
    total_steps = len(loader) * int(training["epochs"])
    peak_step = int(training["one_cycle_peak_step"])
    scheduler = torch.optim.lr_scheduler.OneCycleLR(
        optimizer,
        max_lr=float(training["learning_rate"]),
        steps_per_epoch=len(loader),
        epochs=int(training["epochs"]),
        pct_start=peak_step / total_steps,
    )
    amp = bool(training["amp"])
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    ema = EMA(model, float(model_options["ema_decay"]))
    betas = quadratic_beta_schedule(
        beta_start=float(model_options["beta_start"]),
        beta_end=float(model_options["beta_end"]),
        timesteps=int(model_options["diffusion_timesteps"]),
    ).to(device)
    metadata = {
        "config": str(Path(config_path).resolve()),
        "source": source_manifest,
        "prepared_data": str(Path(prepared_path).resolve()),
        "genes": prepared["genes"],
        "classes": prepared["classes"],
        "parameter_count": int(parameter_count),
        "device": torch.cuda.get_device_name(device),
        "torch_version": str(torch.__version__),
        "batches_per_epoch": len(loader),
        "total_optimizer_steps": total_steps,
        "one_cycle_peak_step": peak_step,
        "implementation_contract": {
            "upstream_model_class": True,
            "upstream_noise_objective": True,
            "antithetic_timesteps": True,
            "parameter_only_ema": True,
            "weighted_sampling": False,
            "gradient_clipping": False,
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
        if payload.get("format") != "nasa_mouse_lacan_paper_parity_checkpoint_v1":
            raise ValueError(f"Incompatible checkpoint: {checkpoint_path}")
        observed = payload.get("metadata", {})
        for key in ("genes", "classes", "parameter_count"):
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
            f"[rna-diffusion:train] resumed epoch={start_epoch} step={global_step}",
            flush=True,
        )

    print(
        f"[rna-diffusion:train] device={metadata['device']} params={parameter_count:,} "
        f"profiles={len(dataset):,} batches={len(loader)} epochs={training['epochs']}",
        flush=True,
    )
    started = time.time()
    epoch_started = started
    completed_epoch = start_epoch
    try:
        for epoch_index in range(start_epoch, int(training["epochs"])):
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
                remaining = (int(training["epochs"]) - completed_epoch) / max(rate, 1e-8)
                print(
                    f"[rna-diffusion:train] epoch={completed_epoch}/{training['epochs']} "
                    f"step={global_step}/{total_steps} loss={row['loss']:.6f} "
                    f"mae={row['noise_absolute_error']:.6f} "
                    f"lr={row['learning_rate']:.8g} eta_hours={remaining / 3600:.2f}",
                    flush=True,
                )
            checkpoint_every = int(training["checkpoint_every_epochs"])
            if completed_epoch % checkpoint_every == 0:
                _atomic_torch_save(
                    _checkpoint_payload(
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
                print(
                    f"[rna-diffusion:train] checkpoint epoch={completed_epoch}",
                    flush=True,
                )
    except BaseException:
        if completed_epoch > start_epoch:
            _atomic_torch_save(
                _checkpoint_payload(
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

    _atomic_torch_save(
        _checkpoint_payload(
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
    final_path = output / "model.pt"
    _atomic_torch_save(
        {
            "format": "nasa_mouse_lacan_paper_parity_model_v1",
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
    summary = {
        "status": "complete",
        "model": "Lacan et al. upstream ModelDDIM",
        "source": source_manifest,
        "run_dir": str(output),
        "model_path": str(final_path),
        "checkpoint_path": str(checkpoint_path),
        "prepared_data": str(prepared_path),
        "device": metadata["device"],
        "parameter_count": parameter_count,
        "classes": prepared["classes"],
        "profiles": {role: len(prepared[role]["expression"]) for role in ("train", "validation", "test")},
        "epochs": completed_epoch,
        "global_steps": global_step,
        "training_seconds_this_invocation": float(time.time() - started),
        "final_loss": history[-1]["loss"],
        "final_noise_absolute_error": history[-1]["noise_absolute_error"],
        "cuda_peak_memory_gb": history[-1]["cuda_peak_memory_gb"],
    }
    (output / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# ARCHS4 Mouse DDIM Reference\n\n"
        "This run uses the pinned, unmodified Lacan et al. `ModelDDIM` architecture "
        "with the retained GTEx hyperparameters. The substituted data are "
        f"{sum(len(prepared[role]['expression']) for role in ('train', 'validation', 'test')):,} "
        "healthy-preferred mouse ARCHS4 bulk profiles represented by a deterministic "
        "974-gene mouse landmark panel. The exact cohort and split are recorded in "
        "the prepared-data manifest.\n\n"
        "See `run_summary.json`, `resolved_config.yaml`, and `training_history.tsv`.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2), flush=True)
    return final_path
