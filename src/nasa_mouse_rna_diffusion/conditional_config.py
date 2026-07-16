"""Configuration contract for the OSDR extension of upstream ModelDDIM."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from .config import validate_paper_model_training


def load_conditional_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("Conditional RNA diffusion configuration must be a mapping")
    if payload.get("contract") != "lacan_modelddim_osdr_extension_v1":
        raise ValueError("Conditional DDIM config lacks its declared extension contract")
    validate_paper_model_training(payload)
    data = payload.get("data", {})
    if data.get("expression_source") != "nasa_osdr_api":
        raise ValueError("Conditional DDIM expression_source must be nasa_osdr_api")
    covariates = tuple(map(str, data.get("conditioning_covariates", [])))
    if not covariates or "condition" not in covariates:
        raise ValueError("Conditional DDIM requires condition in conditioning_covariates")
    if "study" in covariates and data.get("unseen_study_policy") != "unknown_class":
        raise ValueError(
            "Study conditioning with accession holdout requires unseen_study_policy=unknown_class"
        )
    if int(data.get("landmark_dimensions", 0)) != 974:
        raise ValueError("The Lacan paper architecture requires 974 landmark genes")
    if payload.get("training", {}).get("regime") not in {
        "osdr_only",
        "archs4_pretrain_osdr_finetune",
    }:
        raise ValueError("Unsupported conditional DDIM training regime")
    payload["_config_path"] = str(config_path.resolve())
    return payload
