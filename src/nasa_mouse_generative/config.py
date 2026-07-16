"""Configuration loading and validation for one benchmark run."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
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
ALLOWED_SCALERS = {
    "none",
    "zscore",
    "global_zscore",
    "nonzero_global_zscore",
    "robust",
    "maxabs",
}
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
    "age",
    "assay",
    "platform",
    "data_source",
}
ALLOWED_TECHNICAL_REPLICATE_POLICIES = {"keep", "sum", "mean"}
ALLOWED_DEVICES = {"auto", "cpu", "cuda"}
ALLOWED_FEATURE_SELECTION_SOURCES = {"auto", "osdr_train", "archs4_reference"}


@dataclass(frozen=True)
class PreprocessingConfig:
    profile: str = "custom"
    input_units: str = "raw_counts"
    library_normalization: str = "cpm"
    transform: str = "log1p"
    scaler: str = "zscore"
    harmonization: str = "none"
    harmonization_covariates: tuple[str, ...] = ("condition", "tissue", "sex")
    harmonization_parameters: dict[str, Any] = field(default_factory=dict)
    unseen_study_policy: str = "global_train_fallback"
    gene_lengths: str = ""


@dataclass(frozen=True)
class DataConfig:
    osdr_metadata: str = "data/osdr_api/osdr_api_mouse_bulk_rnaseq_flt_gc_metadata.tsv"
    osdr_h5ad: str = "outputs/generative_benchmark/data/osdr/osdr_api_raw_counts.h5ad"
    archs4_h5: str = "assets/archs4/mouse_gene_v2.5.h5"
    archs4_catalog_dir: str = "outputs/generative_benchmark/data_audit/archs4"
    split_dir: str = "outputs/generative_benchmark/splits"
    archs4_cohort: str = "healthy_preferred"
    archs4_max_per_tissue: int = 10000
    archs4_max_per_series: int = 100
    archs4_max_corrupt_profiles: int = 100
    archs4_sample_limit: int = 0
    osdr_sample_limit: int = 0
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
    model_profile: str = "practical_screen"
    model_parameters: dict[str, Any] = field(default_factory=dict)
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
    repeats: int = 1


@dataclass(frozen=True)
class FeatureConfig:
    space: str = "all_shared"
    selection_source: str = "auto"
    selection_sample_limit: int = 5000
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
    fold_id: str = ""
    max_metric_samples: int = 2000


@dataclass(frozen=True)
class ExecutionConfig:
    device: str = "auto"
    resume: bool = True
    checkpoint_every_epochs: int = 100
    num_workers: int = 0
    cache_archs4: bool = True
    evaluate_after_training: bool = True
    save_generated_matrix: bool = False
    save_prepared_data: bool = False
    model_profiles: str = "configs/generative/model_profiles.yaml"
    preprocessing_profiles: str = "configs/generative/preprocessing_profiles.yaml"
    genejepa_source: str = "assets/model_sources/GeneJEPA"
    wgan_source: str = "assets/model_sources/adversarial-gene-expression"
    diffusion_source: str = "assets/model_sources/rna-diffusion"
    retain_training_checkpoint: bool = False
    min_free_space_gb: float = 8.0
    max_run_output_gb: float = 12.0


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
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)

    def validate(self) -> None:
        p = self.preprocessing
        d = self.data
        f = self.features
        t = self.training
        g = self.generation
        v = self.validation
        e = self.execution
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
        if t.regime == "archs4_only" and t.condition_on_flight:
            raise ValueError(
                "ARCHS4 has no flight/ground-control labels. Use "
                "condition_on_flight=false for the ARCHS4-only tissue baseline; "
                "FLT/GC generation requires OSDR fine-tuning."
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
        unknown_harmonization_covariates = set(
            p.harmonization_covariates
        ).difference(ALLOWED_CONDITIONING_COVARIATES)
        if unknown_harmonization_covariates:
            raise ValueError(
                "Unsupported harmonization covariates: "
                f"{sorted(unknown_harmonization_covariates)}"
            )
        if p.harmonization in {"combat", "combat_seq"} and not (
            v.allow_transductive_preprocessing
        ):
            raise ValueError(
                f"{p.harmonization} has no frozen transform for unseen batches. "
                "Set validation.allow_transductive_preprocessing=true and treat "
                "the run as a transductive sensitivity analysis."
            )
        if p.harmonization == "combat_seq" and p.input_units != "raw_counts":
            raise ValueError("ComBat-seq requires raw-count input")
        if (
            p.harmonization in {"combat", "combat_seq"}
            and not t.condition_on_flight
            and "condition" in p.harmonization_covariates
        ):
            raise ValueError(
                "An unconditional negative-control run cannot preserve condition "
                "inside ComBat preprocessing; remove condition from "
                "preprocessing.harmonization_covariates."
            )
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
        if f.selection_source not in ALLOWED_FEATURE_SELECTION_SOURCES:
            raise ValueError(
                f"Unsupported feature-selection source: {f.selection_source}"
            )
        if f.selection_sample_limit < 0:
            raise ValueError("Feature-selection sample limit cannot be negative")
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
        if v.max_metric_samples < 5:
            raise ValueError("validation.max_metric_samples must be at least 5")
        if e.device not in ALLOWED_DEVICES:
            raise ValueError(f"Unsupported execution device: {e.device}")
        if e.checkpoint_every_epochs < 1:
            raise ValueError("execution.checkpoint_every_epochs must be positive")
        if e.num_workers < 0:
            raise ValueError("execution.num_workers cannot be negative")
        if e.min_free_space_gb < 0 or e.max_run_output_gb <= 0:
            raise ValueError(
                "Storage guards require min_free_space_gb >= 0 and max_run_output_gb > 0"
            )
        if d.archs4_max_corrupt_profiles < 0:
            raise ValueError("ARCHS4 corrupt-profile cap cannot be negative")
        if d.archs4_sample_limit < 0 or d.osdr_sample_limit < 0:
            raise ValueError("Runtime sample limits cannot be negative")
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

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


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
    if "model_parameters" in training_options:
        training_options["model_parameters"] = dict(
            training_options["model_parameters"] or {}
        )
    preprocessing_options = dict(data.get("preprocessing", {}))
    if "harmonization_covariates" in preprocessing_options:
        preprocessing_options["harmonization_covariates"] = tuple(
            map(str, preprocessing_options["harmonization_covariates"] or ())
        )
    if "harmonization_parameters" in preprocessing_options:
        preprocessing_options["harmonization_parameters"] = dict(
            preprocessing_options["harmonization_parameters"] or {}
        )
    generation_options = dict(data.get("generation", {}))
    if "synthetic_to_real_ratios" in generation_options:
        generation_options["synthetic_to_real_ratios"] = tuple(
            map(float, generation_options["synthetic_to_real_ratios"] or ())
        )
    return BenchmarkConfig(
        version=int(data.get("version", 1)),
        output_root=str(data.get("output_root", "outputs/generative_benchmark")),
        preprocessing=PreprocessingConfig(**preprocessing_options),
        data=DataConfig(**data_options),
        features=FeatureConfig(**data.get("features", {})),
        training=TrainingConfig(**training_options),
        generation=GenerationConfig(**generation_options),
        validation=ValidationConfig(**data.get("validation", {})),
        execution=ExecutionConfig(**data.get("execution", {})),
    )


def load_config(path: str | Path) -> BenchmarkConfig:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("Configuration root must be a mapping")
    config = _construct(payload)
    config.validate()
    return config


def load_config_with_overrides(
    path: str | Path, overrides: list[str] | tuple[str, ...]
) -> BenchmarkConfig:
    """Load YAML and apply dotted ``key=value`` overrides before validation."""

    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        raise ValueError("Configuration root must be a mapping")
    for override in overrides:
        if "=" not in override:
            raise ValueError(f"Expected dotted key=value override, got {override!r}")
        dotted_key, raw_value = override.split("=", 1)
        keys = [key for key in dotted_key.split(".") if key]
        if not keys:
            raise ValueError(f"Invalid override key: {dotted_key!r}")
        cursor = payload
        for key in keys[:-1]:
            value = cursor.setdefault(key, {})
            if not isinstance(value, dict):
                raise ValueError(
                    f"Cannot assign {dotted_key!r}; {key!r} is not a mapping"
                )
            cursor = value
        cursor[keys[-1]] = yaml.safe_load(raw_value)
    config = _construct(payload)
    config.validate()
    return config
