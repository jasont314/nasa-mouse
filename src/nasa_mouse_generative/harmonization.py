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
        "scanpy_or_R_sva_adapter",
        "Frontiers 2023 found library-preparation/ComBat best for its seven-study liver case; not a universal choice.",
    ),
    "combat_seq": HarmonizationMethod(
        "combat_seq",
        "ComBat-seq",
        "transductive_for_unseen_batches",
        "R_sva_adapter",
        "Requires integer counts and an R environment with sva; unavailable in the current environment.",
    ),
    "mober": HarmonizationMethod(
        "mober",
        "MOBER batch-aware VAE projection",
        "inductive_projection_to_trained_source",
        "official_mober_adapter",
        "Must preserve FLT/GC while adversarially reducing source information; validate by held-out accession.",
    ),
}
