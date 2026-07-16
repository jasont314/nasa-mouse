"""Configuration helpers for the paper-parity mouse DDIM experiment."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


REQUIRED_SPLIT = {"train": 9796, "validation": 2448, "test": 5000}


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("RNA diffusion configuration must be a mapping")
    data = payload.get("data", {})
    model = payload.get("model", {})
    training = payload.get("training", {})
    split = data.get("split", {})
    if split != REQUIRED_SPLIT:
        raise ValueError(
            f"Paper-parity split must be {REQUIRED_SPLIT}, observed {split}"
        )
    if int(data.get("profiles", 0)) != sum(REQUIRED_SPLIT.values()):
        raise ValueError("Paper-parity data must contain exactly 17,244 profiles")
    expected_model = {
        "hidden_dims": [8192, 8192],
        "dropout": 0.1,
        "time_embedding_dim": 1,
        "tissue_embedding_dim": 2,
        "sinusoidal_time": False,
        "diffusion_timesteps": 1000,
        "beta_schedule": "quad",
        "beta_start": 0.0001,
        "beta_end": 0.02,
        "ema_decay": 0.999,
    }
    for key, expected in expected_model.items():
        if model.get(key) != expected:
            raise ValueError(
                f"Paper-parity model.{key} must be {expected!r}, "
                f"observed {model.get(key)!r}"
            )
    expected_training = {
        "epochs": 15000,
        "batch_size": 2048,
        "optimizer": "Adam",
        "learning_rate": 0.0004783833151836702,
        "weight_decay": 0.0,
        "one_cycle_peak_step": 1000,
        "gradient_clipping": False,
        "weighted_sampling": False,
        "amp": True,
    }
    for key, expected in expected_training.items():
        if training.get(key) != expected:
            raise ValueError(
                f"Paper-parity training.{key} must be {expected!r}, "
                f"observed {training.get(key)!r}"
            )
    payload["_config_path"] = str(config_path.resolve())
    return payload
