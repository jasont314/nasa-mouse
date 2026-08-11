"""Configuration helpers for the paper-parity mouse DDIM experiment."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


REQUIRED_SPLIT = {"train": 9796, "validation": 2448, "test": 5000}
CUSTOM_REFERENCE_CONTRACT = "lacan_modelddim_reference_v1"
DISJOINT_REFERENCE_CONTRACT = "lacan_modelddim_osdr_disjoint_reference_v1"
EXPECTED_MODEL = {
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
EXPECTED_TRAINING = {
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


def validate_paper_model_training(payload: dict[str, Any]) -> None:
    model = payload.get("model", {})
    training = payload.get("training", {})
    for key, expected in EXPECTED_MODEL.items():
        if model.get(key) != expected:
            raise ValueError(
                f"Paper-parity model.{key} must be {expected!r}, "
                f"observed {model.get(key)!r}"
            )
    for key, expected in EXPECTED_TRAINING.items():
        if training.get(key) != expected:
            raise ValueError(
                f"Paper-parity training.{key} must be {expected!r}, "
                f"observed {training.get(key)!r}"
            )


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("RNA diffusion configuration must be a mapping")
    data = payload.get("data", {})
    contract = payload.get("contract")
    if contract in {CUSTOM_REFERENCE_CONTRACT, DISJOINT_REFERENCE_CONTRACT}:
        profiles = int(data.get("profiles", 0))
        if contract == CUSTOM_REFERENCE_CONTRACT:
            quotas = {
                str(tissue): int(count)
                for tissue, count in data.get("profiles_per_tissue", {}).items()
            }
            if profiles <= 0 or sum(quotas.values()) != profiles:
                raise ValueError(
                    "Custom reference profiles_per_tissue must contain positive quotas "
                    "that sum to data.profiles"
                )
            if len(quotas) < 2 or any(count <= 0 for count in quotas.values()):
                raise ValueError(
                    "The released one-hot embedding requires at least two populated "
                    "reference classes"
                )
        elif profiles != sum(REQUIRED_SPLIT.values()):
            raise ValueError(
                "OSDR-disjoint paper-parity data must contain exactly 17,244 profiles"
            )
        elif not (
            data.get("exclude_series_ids")
            or data.get("exclude_series_ids_file")
        ):
            raise ValueError(
                "OSDR-disjoint reference requires excluded GEO series"
            )
        if data.get("split_strategy") != "series_holdout":
            raise ValueError("Reference split_strategy must be series_holdout")
        fractions = data.get("split_fractions", {})
        if set(fractions) != {"train", "validation", "test"}:
            raise ValueError(
                "Reference split_fractions must define train, validation, and test"
            )
        if any(float(value) <= 0 for value in fractions.values()) or abs(
            sum(map(float, fractions.values())) - 1.0
        ) > 1e-8:
            raise ValueError("Reference split_fractions must be positive and sum to 1")
    else:
        split = data.get("split", {})
        if split != REQUIRED_SPLIT:
            raise ValueError(
                f"Paper-parity split must be {REQUIRED_SPLIT}, observed {split}"
            )
        if int(data.get("profiles", 0)) != sum(REQUIRED_SPLIT.values()):
            raise ValueError("Paper-parity data must contain exactly 17,244 profiles")
    validate_paper_model_training(payload)
    payload["_config_path"] = str(config_path.resolve())
    return payload
