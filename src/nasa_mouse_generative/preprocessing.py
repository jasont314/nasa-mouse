"""Fold-aware preprocessing primitives shared by model adapters."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from .config import PreprocessingConfig


EPSILON = 1e-8


@dataclass
class ScaleStats:
    center: np.ndarray
    scale: np.ndarray


@dataclass
class FittedPreprocessor:
    spec: PreprocessingConfig
    global_stats: ScaleStats | None = None
    study_stats: dict[str, ScaleStats] = field(default_factory=dict)
    post_harmonization_stats: ScaleStats | None = None
    final_stats: ScaleStats | None = None

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        *,
        gene_lengths: np.ndarray | None = None,
    ) -> np.ndarray:
        studies_array = np.asarray(list(studies), dtype=str)
        values = self._base_transform(matrix, gene_lengths=gene_lengths)
        if values.shape[0] != studies_array.shape[0]:
            raise ValueError("studies must contain one value per sample")

        if self.spec.harmonization in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            self.global_stats = fit_stats(values, "zscore")
            harmonized = np.empty_like(values, dtype=np.float32)
            for study in sorted(set(studies_array)):
                mask = studies_array == study
                stats = fit_stats(values[mask], "zscore")
                self.study_stats[study] = stats
                harmonized[mask] = apply_stats(values[mask], stats)
            if self.spec.harmonization == "within_study_then_global_zscore":
                self.post_harmonization_stats = fit_stats(harmonized, "zscore")
                harmonized = apply_stats(harmonized, self.post_harmonization_stats)
            return self._fit_final_scaler(harmonized)

        if self.spec.harmonization not in {"none", "combat", "combat_seq", "mober"}:
            raise ValueError(f"Unsupported harmonization: {self.spec.harmonization}")
        if self.spec.harmonization != "none":
            raise NotImplementedError(
                f"{self.spec.harmonization} is a model-based/transductive harmonizer and "
                "must be invoked through its dedicated adapter."
            )
        return self._fit_final_scaler(values)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        *,
        gene_lengths: np.ndarray | None = None,
        allow_transductive: bool = False,
    ) -> np.ndarray:
        studies_array = np.asarray(list(studies), dtype=str)
        values = self._base_transform(matrix, gene_lengths=gene_lengths)
        if values.shape[0] != studies_array.shape[0]:
            raise ValueError("studies must contain one value per sample")
        if self.spec.harmonization in {
            "within_study_zscore",
            "within_study_then_global_zscore",
        }:
            if self.global_stats is None:
                raise RuntimeError("Preprocessor is not fitted")
            harmonized = np.empty_like(values, dtype=np.float32)
            for study in sorted(set(studies_array)):
                mask = studies_array == study
                stats = self.study_stats.get(study)
                if stats is None:
                    if self.spec.unseen_study_policy == "transductive_unlabeled":
                        if not allow_transductive:
                            raise ValueError(
                                "Unseen-study transductive scaling was requested but not allowed."
                            )
                        stats = fit_stats(values[mask], "zscore")
                    else:
                        stats = self.global_stats
                harmonized[mask] = apply_stats(values[mask], stats)
            if self.post_harmonization_stats is not None:
                harmonized = apply_stats(harmonized, self.post_harmonization_stats)
            return self._apply_final_scaler(harmonized)
        return self._apply_final_scaler(values)

    def _base_transform(
        self, matrix: np.ndarray, *, gene_lengths: np.ndarray | None
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        if values.ndim != 2:
            raise ValueError("expression matrix must be samples x genes")
        if not np.isfinite(values).all():
            raise ValueError("expression matrix contains non-finite values")
        if (values < 0).any() and self.spec.input_units == "raw_counts":
            raise ValueError("raw counts cannot be negative")

        normalization = self.spec.library_normalization
        if normalization == "none":
            pass
        elif normalization == "cpm":
            library = values.sum(axis=1, keepdims=True)
            values = values / np.maximum(library, EPSILON) * 1_000_000.0
        elif normalization == "tpm":
            if gene_lengths is None:
                raise ValueError("TPM normalization requires gene lengths")
            lengths = np.asarray(gene_lengths, dtype=np.float64).reshape(1, -1)
            if lengths.shape[1] != values.shape[1] or (lengths <= 0).any():
                raise ValueError("gene lengths must be positive and aligned to matrix columns")
            rates = values / (lengths / 1000.0)
            values = rates / np.maximum(rates.sum(axis=1, keepdims=True), EPSILON) * 1_000_000.0
        else:
            raise ValueError(f"Unsupported library normalization: {normalization}")

        if self.spec.transform == "log1p":
            values = np.log1p(values)
        elif self.spec.transform == "log2p1":
            values = np.log2(values + 1.0)
        elif self.spec.transform != "none":
            raise ValueError(f"Unsupported transform: {self.spec.transform}")
        return values.astype(np.float32)

    def _fit_final_scaler(self, values: np.ndarray) -> np.ndarray:
        if self.spec.scaler == "none":
            return values.astype(np.float32, copy=False)
        stats = fit_stats(values, self.spec.scaler)
        # Keep final scaling separate from the mentor method's second z-score.
        self.final_stats = stats
        return apply_stats(values, stats)

    def _apply_final_scaler(self, values: np.ndarray) -> np.ndarray:
        if self.spec.scaler == "none":
            return values.astype(np.float32, copy=False)
        if self.final_stats is None:
            raise RuntimeError("Preprocessor is not fitted")
        return apply_stats(values, self.final_stats)


def fit_stats(values: np.ndarray, method: str) -> ScaleStats:
    array = np.asarray(values, dtype=np.float64)
    if method == "zscore":
        center = array.mean(axis=0)
        scale = array.std(axis=0)
    elif method == "robust":
        center = np.median(array, axis=0)
        q25, q75 = np.percentile(array, [25, 75], axis=0)
        scale = q75 - q25
    elif method == "maxabs":
        center = np.zeros(array.shape[1], dtype=np.float64)
        scale = np.max(np.abs(array), axis=0)
    else:
        raise ValueError(f"Unsupported scaler: {method}")
    scale = np.where(np.abs(scale) < EPSILON, 1.0, scale)
    return ScaleStats(center=center.astype(np.float32), scale=scale.astype(np.float32))


def apply_stats(values: np.ndarray, stats: ScaleStats) -> np.ndarray:
    return ((np.asarray(values, dtype=np.float32) - stats.center) / stats.scale).astype(
        np.float32
    )
