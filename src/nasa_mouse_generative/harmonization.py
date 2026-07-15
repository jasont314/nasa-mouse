"""Harmonization method registry and leakage constraints."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HarmonizationMethod:
    method_id: str
    label: str
    fold_behavior: str
    implementation: str
    notes: str


HARMONIZATION_REGISTRY = {
    "none": HarmonizationMethod(
        "none",
        "No cross-study harmonization",
        "inductive",
        "built_in",
        "Model may still receive study as a condition when configured.",
    ),
    "within_study_zscore": HarmonizationMethod(
        "within_study_zscore",
        "Study-wise log transform and gene z-score",
        "inductive_with_global_fallback_or_explicit_transductive_mode",
        "built_in",
        "Matches the core Ilangovan et al. study-by-study scaling design.",
    ),
    "within_study_then_global_zscore": HarmonizationMethod(
        "within_study_then_global_zscore",
        "Study-wise gene z-score followed by pooled gene z-score",
        "inductive_with_global_fallback_or_explicit_transductive_mode",
        "built_in",
        "Mentor-proposed two-stage standardization; all statistics are fold-scoped.",
    ),
    "combat": HarmonizationMethod(
        "combat",
        "ComBat",
        "transductive_for_unseen_batches",
        "built_in_parametric_empirical_bayes_adapter",
        "Frozen parameters are used for known batches; unseen held-out batches "
        "require an explicitly enabled transductive estimate.",
    ),
    "combat_seq": HarmonizationMethod(
        "combat_seq",
        "ComBat-seq",
        "transductive_for_unseen_batches",
        "R_sva_adapter",
        "Calls Bioconductor sva::ComBat_seq on integer counts; fractional "
        "count-like input and singleton batches require explicit policies.",
    ),
    "mober": HarmonizationMethod(
        "mober",
        "MOBER batch-aware VAE projection",
        "inductive_projection_to_trained_source",
        "official_architecture_compatible_pytorch_adapter",
        "Projects held-out samples through a frozen encoder and target-batch "
        "decoder; validate preservation of FLT/GC by accession.",
    ),
}
