"""Configuration loading and validation for one benchmark run."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .models import MODEL_REGISTRY, require_generation_capability


ALLOWED_TRAINING_REGIMES = {
    "osdr_only",
    "archs4_only",
    "archs4_pretrain_osdr_finetune",
}
ALLOWED_TISSUE_MODES = {"pooled_conditioned", "per_tissue"}
ALLOWED_INPUT_UNITS = {"raw_counts", "cpm", "tpm"}
ALLOWED_LIBRARY_NORMALIZATIONS = {"none", "cpm", "tpm"}
ALLOWED_TRANSFORMS = {"none", "log1p", "log2p1"}
ALLOWED_SCALERS = {"none", "zscore", "robust", "maxabs"}
ALLOWED_HARMONIZERS = {
    "none",
    "within_study_zscore",
    "within_study_then_global_zscore",
    "combat",
    "combat_seq",
    "mober",
}
ALLOWED_STUDY_POLICIES = {"not_conditioned", "conditioned"}
ALLOWED_TASKS = {"conditional_generation", "representation"}
ALLOWED_UNSEEN_STUDY_POLICIES = {
    "global_train_fallback",
    "transductive_unlabeled",
}
ALLOWED_ARCHS4_COHORTS = {"control_only", "healthy_preferred", "broad"}
ALLOWED_ACCESSION_SCOPES = {"all_eligible", "single", "selected"}
ALLOWED_FEATURE_SPACES = {
    "all_shared",
    "hvg",
    "reactome_shared",
    "l1000_landmarks",
}
ALLOWED_CONDITIONING_COVARIATES = {
    "condition",
    "tissue",
    "material_type",
    "muscle_group",
    "study",
    "sex",
    "assay",
    "platform",
    "data_source",
}
ALLOWED_TECHNICAL_REPLICATE_POLICIES = {"keep", "sum", "mean"}


@dataclass(frozen=True)
class PreprocessingConfig:
    input_units: str = "raw_counts"
    library_normalization: str = "cpm"
    transform: str = "log1p"
    scaler: str = "zscore"
    harmonization: str = "none"
    unseen_study_policy: str = "global_train_fallback"
    gene_lengths: str = ""


@dataclass(frozen=True)
class DataConfig:
    osdr_metadata: str = "data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
    archs4_h5: str = "assets/archs4/mouse_gene_v2.5.h5"
    archs4_cohort: str = "healthy_preferred"
    archs4_max_per_tissue: int = 10000
    archs4_max_per_series: int = 100
    osdr_accession_scope: str = "all_eligible"
    osdr_include_accessions: tuple[str, ...] = ()
    osdr_exclude_accessions: tuple[str, ...] = ()
    osdr_tissues: tuple[str, ...] = ()
    technical_replicate_policy: str = "sum"
    min_confirmatory_total: int = 60
    min_confirmatory_per_condition: int = 20
    min_confirmatory_accessions: int = 5
    min_exploratory_total: int = 30
    min_exploratory_per_condition: int = 10
    min_exploratory_accessions: int = 2


@dataclass(frozen=True)
class TrainingConfig:
    model: str = "vinas_wgan_gp"
    task: str = "conditional_generation"
    regime: str = "archs4_pretrain_osdr_finetune"
    tissue_mode: str = "pooled_conditioned"
    condition_on_flight: bool = True
    study_policy: str = "not_conditioned"
    conditioning_covariates: tuple[str, ...] = (
        "condition",
        "tissue",
        "material_type",
        "muscle_group",
        "sex",
        "assay",
        "platform",
        "data_source",
    )
    seed: int = 2020
    repeats: int = 5


@dataclass(frozen=True)
class FeatureConfig:
    space: str = "all_shared"
    max_genes: int = 0
    hvg_genes: int = 2000
    l1000_map: str = "data/diffusion/l1000_human_to_mouse_ensembl.tsv"
    reactome_gmt: str = "data/pathways/reactome_current_mouse_ensembl.gmt"


@dataclass(frozen=True)
class GenerationConfig:
    samples_per_covariate_profile: int = 100
    synthetic_to_real_ratios: tuple[float, ...] = (0.5, 1.0, 2.0)
    paired_counterfactual: bool = False


@dataclass(frozen=True)
class ValidationConfig:
    split_unit: str = "accession"
    selection_metric: str = "heldout_fidelity_composite"
    downstream_flt_gc_secondary: bool = True
    final_test_locked: bool = True
    allow_transductive_preprocessing: bool = False
    pooled_validation_fraction: float = 0.15
    pooled_test_fraction: float = 0.15


@dataclass(frozen=True)
class BenchmarkConfig:
    version: int = 1
    output_root: str = "outputs/generative_benchmark"
    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    data: DataConfig = field(default_factory=DataConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    generation: GenerationConfig = field(default_factory=GenerationConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)

    def validate(self) -> None:
        p = self.preprocessing
        d = self.data
        f = self.features
        t = self.training
        g = self.generation
        v = self.validation
        if self.version != 1:
            raise ValueError(f"Unsupported configuration version: {self.version}")
        if t.model not in MODEL_REGISTRY:
            raise ValueError(f"Unknown model {t.model!r}; choose from {sorted(MODEL_REGISTRY)}")
        if t.task not in ALLOWED_TASKS:
            raise ValueError(f"Unsupported task: {t.task}")
        if t.task == "conditional_generation":
            require_generation_capability(t.model)
        if t.regime not in ALLOWED_TRAINING_REGIMES:
            raise ValueError(f"Unsupported training regime: {t.regime}")
        if t.tissue_mode not in ALLOWED_TISSUE_MODES:
            raise ValueError(f"Unsupported tissue mode: {t.tissue_mode}")
        if t.study_policy not in ALLOWED_STUDY_POLICIES:
            raise ValueError(f"Unsupported study policy: {t.study_policy}")
        unknown_covariates = set(t.conditioning_covariates).difference(
            ALLOWED_CONDITIONING_COVARIATES
        )
        if unknown_covariates:
            raise ValueError(
                f"Unsupported conditioning covariates: {sorted(unknown_covariates)}"
            )
        if t.condition_on_flight and "condition" not in t.conditioning_covariates:
            raise ValueError(
                "condition_on_flight=true requires condition in conditioning_covariates"
            )
        if t.study_policy == "conditioned" and "study" not in t.conditioning_covariates:
            raise ValueError(
                "study_policy=conditioned requires study in conditioning_covariates"
            )
        if t.study_policy == "not_conditioned" and "study" in t.conditioning_covariates:
            raise ValueError(
                "Remove study from conditioning_covariates when study_policy=not_conditioned"
            )
        if p.input_units not in ALLOWED_INPUT_UNITS:
            raise ValueError(f"Unsupported input units: {p.input_units}")
        if p.library_normalization not in ALLOWED_LIBRARY_NORMALIZATIONS:
            raise ValueError(
                f"Unsupported library normalization: {p.library_normalization}"
            )
        if p.transform not in ALLOWED_TRANSFORMS:
            raise ValueError(f"Unsupported transform: {p.transform}")
        if p.scaler not in ALLOWED_SCALERS:
            raise ValueError(f"Unsupported scaler: {p.scaler}")
        if p.harmonization not in ALLOWED_HARMONIZERS:
            raise ValueError(f"Unsupported harmonization: {p.harmonization}")
        if p.unseen_study_policy not in ALLOWED_UNSEEN_STUDY_POLICIES:
            raise ValueError(
                f"Unsupported unseen-study policy: {p.unseen_study_policy}"
            )
        if p.library_normalization != "none" and p.input_units != "raw_counts":
            raise ValueError(
                "Library normalization can only be fitted to raw-count input; "
                "set library_normalization=none for precomputed CPM or TPM."
            )
        if p.library_normalization == "tpm" and not p.gene_lengths:
            raise ValueError("TPM normalization requires preprocessing.gene_lengths")
        if p.unseen_study_policy == "transductive_unlabeled" and not v.allow_transductive_preprocessing:
            raise ValueError(
                "Transductive unseen-study scaling must be explicitly enabled in validation."
            )
        if t.condition_on_flight is False and t.task == "conditional_generation":
            # This is a supported negative control, not the default scientific run.
            pass
        if d.archs4_cohort not in ALLOWED_ARCHS4_COHORTS:
            raise ValueError(
                f"data.archs4_cohort must be one of {sorted(ALLOWED_ARCHS4_COHORTS)}"
            )
        if d.osdr_accession_scope not in ALLOWED_ACCESSION_SCOPES:
            raise ValueError(
                f"Unsupported OSDR accession scope: {d.osdr_accession_scope}"
            )
        if d.osdr_accession_scope == "single" and len(d.osdr_include_accessions) != 1:
            raise ValueError(
                "osdr_accession_scope=single requires exactly one included accession"
            )
        if d.osdr_accession_scope == "selected" and not d.osdr_include_accessions:
            raise ValueError(
                "osdr_accession_scope=selected requires osdr_include_accessions"
            )
        if d.technical_replicate_policy not in ALLOWED_TECHNICAL_REPLICATE_POLICIES:
            raise ValueError(
                "Unsupported technical replicate policy: "
                f"{d.technical_replicate_policy}"
            )
        if f.space not in ALLOWED_FEATURE_SPACES:
            raise ValueError(f"Unsupported feature space: {f.space}")
        if f.max_genes < 0 or f.hvg_genes < 1:
            raise ValueError("Feature counts must be non-negative, with hvg_genes >= 1")
        if g.samples_per_covariate_profile < 1:
            raise ValueError("samples_per_covariate_profile must be positive")
        if not g.synthetic_to_real_ratios or any(
            ratio <= 0 for ratio in g.synthetic_to_real_ratios
        ):
            raise ValueError("synthetic_to_real_ratios must contain positive values")
        if g.paired_counterfactual:
            model = MODEL_REGISTRY[t.model]
            if not model.supports_counterfactual_pairing:
                raise ValueError(
                    f"{model.display_name} does not support paired counterfactual generation"
                )
            if not t.condition_on_flight:
                raise ValueError(
                    "Paired FLT/GC counterfactual generation requires condition conditioning"
                )
        if not 0 < v.pooled_validation_fraction < 0.5:
            raise ValueError("pooled_validation_fraction must be between 0 and 0.5")
        if not 0 < v.pooled_test_fraction < 0.5:
            raise ValueError("pooled_test_fraction must be between 0 and 0.5")
        if v.pooled_validation_fraction + v.pooled_test_fraction >= 0.8:
            raise ValueError("Pooled validation and test fractions leave too little training data")
        for name, value in {
            "min_confirmatory_total": d.min_confirmatory_total,
            "min_confirmatory_per_condition": d.min_confirmatory_per_condition,
            "min_confirmatory_accessions": d.min_confirmatory_accessions,
            "min_exploratory_total": d.min_exploratory_total,
            "min_exploratory_per_condition": d.min_exploratory_per_condition,
            "min_exploratory_accessions": d.min_exploratory_accessions,
            "repeats": t.repeats,
        }.items():
            if value < 1:
                raise ValueError(f"{name} must be positive")


def _construct(data: dict[str, Any]) -> BenchmarkConfig:
    data_options = dict(data.get("data", {}))
    for key in ("osdr_include_accessions", "osdr_exclude_accessions", "osdr_tissues"):
        if key in data_options:
            data_options[key] = tuple(map(str, data_options[key] or ()))
    training_options = dict(data.get("training", {}))
    if "conditioning_covariates" in training_options:
        training_options["conditioning_covariates"] = tuple(
            map(str, training_options["conditioning_covariates"] or ())
        )
    generation_options = dict(data.get("generation", {}))
    if "synthetic_to_real_ratios" in generation_options:
        generation_options["synthetic_to_real_ratios"] = tuple(
            map(float, generation_options["synthetic_to_real_ratios"] or ())
        )
    return BenchmarkConfig(
        version=int(data.get("version", 1)),
        output_root=str(data.get("output_root", "outputs/generative_benchmark")),
        preprocessing=PreprocessingConfig(**data.get("preprocessing", {})),
        data=DataConfig(**data_options),
        features=FeatureConfig(**data.get("features", {})),
        training=TrainingConfig(**training_options),
        generation=GenerationConfig(**generation_options),
        validation=ValidationConfig(**data.get("validation", {})),
    )


def load_config(path: str | Path) -> BenchmarkConfig:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("Configuration root must be a mapping")
    config = _construct(payload)
    config.validate()
    return config
