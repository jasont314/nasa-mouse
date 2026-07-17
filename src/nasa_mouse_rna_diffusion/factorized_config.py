"""Configuration validation for factorized ModelDDIM residual adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from .config import EXPECTED_MODEL


CONTRACT = "lacan_modelddim_factorized_adapter_v1"


def load_factorized_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict) or payload.get("contract") != CONTRACT:
        raise ValueError(f"Factorized adapter config must declare {CONTRACT}")
    model = payload.get("model", {})
    for key, expected in EXPECTED_MODEL.items():
        if model.get(key) != expected:
            raise ValueError(
                f"Factorized adapter retains paper model.{key}={expected!r}"
            )
    data = payload.get("data", {})
    for key in ("prepared_h5", "samples_tsv", "pretrained_model"):
        if not data.get(key):
            raise ValueError(f"Factorized adapter requires data.{key}")
    training = payload.get("training", {})
    if set(training.get("stages", {})) != {"domain", "condition"}:
        raise ValueError("training.stages must define domain and condition")
    for name, stage in training["stages"].items():
        if int(stage.get("steps", 0)) <= 0 or float(stage.get("learning_rate", 0)) <= 0:
            raise ValueError(f"Stage {name} requires positive steps and learning_rate")
    if not 0.0 <= float(training.get("condition_dropout", 0.0)) < 1.0:
        raise ValueError("training.condition_dropout must be in [0, 1)")
    conditioning = payload.get("conditioning", {})
    unknown_conditioning = set(conditioning).difference(
        {"study", "material_type"}
    )
    if unknown_conditioning:
        raise ValueError(
            f"Unsupported factorized conditioning keys: {sorted(unknown_conditioning)}"
        )
    adapter = payload.get("adapter", {})
    if int(adapter.get("domain_lora_rank", 0)) < 0:
        raise ValueError("adapter.domain_lora_rank cannot be negative")
    if float(adapter.get("domain_lora_alpha", 1.0)) <= 0:
        raise ValueError("adapter.domain_lora_alpha must be positive")
    if adapter.get("initial_model") and not isinstance(
        adapter.get("initial_model"), str
    ):
        raise ValueError("adapter.initial_model must be a path string")
    for name, stage in training["stages"].items():
        regularization = stage.get("correlation_regularization")
        if regularization is not None:
            if name != "domain":
                raise ValueError(
                    "correlation_regularization is supported only for the domain stage"
                )
            if float(regularization.get("weight", 0.0)) <= 0.0:
                raise ValueError("correlation_regularization.weight must be positive")
            if int(regularization.get("genes", 0)) < 2:
                raise ValueError("correlation_regularization.genes must be at least two")
            if int(regularization.get("max_timestep", -1)) < 0:
                raise ValueError(
                    "correlation_regularization.max_timestep cannot be negative"
                )
        effect_regularization = stage.get("effect_regularization")
        if effect_regularization is not None:
            if name != "condition":
                raise ValueError(
                    "effect_regularization is supported only for the condition stage"
                )
            if float(effect_regularization.get("weight", 0.0)) <= 0.0:
                raise ValueError("effect_regularization.weight must be positive")
            if int(effect_regularization.get("genes", 0)) < 2:
                raise ValueError("effect_regularization.genes must be at least two")
            if int(effect_regularization.get("max_timestep", -1)) < 0:
                raise ValueError("effect_regularization.max_timestep cannot be negative")
    evaluation = payload.get("evaluation", {})
    if evaluation.get("split", "validation") != "validation":
        raise ValueError("Development configs must evaluate validation, not locked test")
    if evaluation.get("sampling_noise", "pseudo_random") not in {
        "pseudo_random",
        "stratified_antithetic",
    }:
        raise ValueError("evaluation.sampling_noise is unsupported")
    payload["_config_path"] = str(config_path.resolve())
    return payload
