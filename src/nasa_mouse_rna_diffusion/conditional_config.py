"""Configuration contract for the OSDR extension of upstream ModelDDIM."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from .config import validate_paper_model_training
from nasa_mouse_generative.config import PreprocessingConfig


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
    split_strategy = str(data.get("split_strategy", "accession_holdout"))
    if split_strategy not in {"accession_holdout", "within_study_stratified"}:
        raise ValueError("Unsupported conditional DDIM data.split_strategy")
    data["split_strategy"] = split_strategy
    validation_accessions = tuple(
        map(str, data.get("validation_accessions", []))
    )
    test_accessions = tuple(map(str, data.get("test_accessions", [])))
    if validation_accessions or test_accessions:
        if split_strategy != "accession_holdout":
            raise ValueError(
                "Explicit validation/test accessions require accession_holdout"
            )
        if not validation_accessions or not test_accessions:
            raise ValueError(
                "Explicit accession splitting requires both validation_accessions "
                "and test_accessions"
            )
        overlap = set(validation_accessions) & set(test_accessions)
        if overlap:
            raise ValueError(
                f"Validation and test accessions overlap: {sorted(overlap)}"
            )
        data["validation_accessions"] = list(validation_accessions)
        data["test_accessions"] = list(test_accessions)
    if (
        "study" in covariates
        and split_strategy == "accession_holdout"
        and data.get("unseen_study_policy") != "unknown_class"
    ):
        raise ValueError(
            "Study conditioning with accession holdout requires unseen_study_policy=unknown_class"
        )
    if int(data.get("landmark_dimensions", 0)) != 974:
        raise ValueError("The Lacan paper architecture requires 974 landmark genes")
    expression_representation = str(
        data.get("expression_representation", "full_transcriptome_tpm")
    )
    if expression_representation not in {
        "full_transcriptome_tpm",
        "deseq2_median_of_ratios_by_study",
        "deseq2_median_of_ratios_pooled",
    }:
        raise ValueError(
            "Unsupported conditional DDIM data.expression_representation"
        )
    data["expression_representation"] = expression_representation
    regime = payload.get("training", {}).get("regime")
    if regime not in {
        "osdr_only",
        "archs4_pretrain_osdr_finetune",
    }:
        raise ValueError("Unsupported conditional DDIM training regime")
    if regime == "archs4_pretrain_osdr_finetune":
        training = payload["training"]
        for key in (
            "pretrained_model",
            "finetune_epochs",
            "finetune_learning_rate",
            "finetune_one_cycle_peak_step",
        ):
            if key not in training:
                raise ValueError(f"Pretrain/fine-tune config requires training.{key}")
        initialization = str(
            training.get("pretrained_condition_initialization", "reference_only")
        )
        if initialization not in {
            "reference_only",
            "function_preserving_tissue",
        }:
            raise ValueError(
                "Unsupported training.pretrained_condition_initialization"
            )
        training["pretrained_condition_initialization"] = initialization
        if data.get("preprocessing"):
            raise ValueError(
                "Custom OSDR preprocessing is not compatible with the existing "
                "MaxAbs ARCHS4 checkpoint; pretrain a matching reference first or "
                "use osdr_only."
            )
    preprocessing = data.get("preprocessing")
    if preprocessing:
        options = dict(preprocessing)
        if "harmonization_covariates" in options:
            options["harmonization_covariates"] = tuple(
                options["harmonization_covariates"]
            )
        PreprocessingConfig(**options)
    elif expression_representation != "full_transcriptome_tpm":
        raise ValueError(
            "DESeq2 median-of-ratios representations require data.preprocessing"
        )
    payload["_config_path"] = str(config_path.resolve())
    return payload
