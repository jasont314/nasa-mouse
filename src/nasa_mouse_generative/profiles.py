"""Resolve paper/native model profiles and per-run parameter overrides."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import yaml

from .config import BenchmarkConfig


def load_model_parameters(config: BenchmarkConfig) -> dict[str, Any]:
    path = Path(config.execution.model_profiles)
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    profiles = payload.get("profiles", {})
    model_profiles = profiles.get(config.training.model)
    if not isinstance(model_profiles, dict):
        raise ValueError(
            f"No model profiles for {config.training.model!r} in {path}"
        )
    profile = model_profiles.get(config.training.model_profile)
    if not isinstance(profile, dict):
        raise ValueError(
            f"Unknown profile {config.training.model_profile!r} for "
            f"{config.training.model}; choose from {sorted(model_profiles)}"
        )
    resolved = dict(profile)
    resolved.update(config.training.model_parameters)
    return resolved


def resolve_preprocessing_profile(config: BenchmarkConfig) -> BenchmarkConfig:
    name = config.preprocessing.profile
    if name in {"", "custom"}:
        return config
    path = Path(config.execution.preprocessing_profiles)
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if name == "model_native":
        profile = payload.get("paper_native_profiles", {}).get(config.training.model)
    else:
        profile = payload.get("shared_profiles", {}).get(name)
    if not isinstance(profile, dict):
        raise ValueError(f"Unknown preprocessing profile {name!r} in {path}")
    supported = {
        key: profile[key]
        for key in (
            "input_units",
            "library_normalization",
            "transform",
            "scaler",
        )
        if key in profile
    }
    resolved = replace(
        config,
        preprocessing=replace(config.preprocessing, **supported),
    )
    resolved.validate()
    return resolved


def epochs_for_stage(parameters: dict[str, Any], stage: str) -> int:
    default = int(parameters.get("epochs", 100))
    if stage == "reference":
        return int(parameters.get("reference_epochs", default))
    if stage == "osdr_finetune":
        return int(parameters.get("finetune_epochs", default))
    if stage == "osdr":
        return int(parameters.get("osdr_epochs", default))
    raise ValueError(f"Unknown training stage: {stage}")


def learning_rate_for_stage(parameters: dict[str, Any], stage: str) -> float:
    default = float(parameters.get("learning_rate", 1e-4))
    if stage == "osdr_finetune":
        return float(parameters.get("finetune_learning_rate", default))
    return default
