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
    evaluation = payload.get("evaluation", {})
    if evaluation.get("split", "validation") != "validation":
        raise ValueError("Development configs must evaluate validation, not locked test")
    payload["_config_path"] = str(config_path.resolve())
    return payload
