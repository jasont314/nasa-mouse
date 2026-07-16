"""Adapters for the correction functions distributed by the official MBatch repo."""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .harmonizers import Harmonizer, _batch_values, _metadata_frame, _resolve_batch_key


METHOD_FILES = {
    "mbatch_median_polish": "BEA_CorrectionsMP.R",
    "mbatch_empirical_bayes": "BEA_CorrectionsEB.R",
    "mbatch_anova": "BEA_CorrectionsAN.R",
}


class MBatchHarmonizer(Harmonizer):
    """Run official MBatch R functions with fold-scoped training anchors."""

    is_transductive = True

    def __init__(self, *, method_id: str, **kwargs: Any) -> None:
        if method_id not in METHOD_FILES:
            raise ValueError(f"Unsupported MBatch correction: {method_id}")
        super().__init__(**kwargs)
        self.method_id = method_id
        self.batch_key = ""
        self.batch_levels: list[str] = []
        self.anchor_matrix: np.ndarray | None = None
        self.anchor_metadata: pd.DataFrame | None = None
        self.r_runtime: dict[str, str] = {}
        self.source_identity: dict[str, object] = {}
        self.nonfinite_audit: list[dict[str, object]] = []

    def _source_path(self) -> Path:
        root = Path(
            self.parameters.get(
                "source_root", "assets/model_sources/MBatch/apps/MBatch/R"
            )
        )
        path = (root / METHOD_FILES[self.method_id]).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Official MBatch source file is absent: {path}")
        return path

    def _run(self, matrix: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError("MBatch correction requires finite expression")
        if self.covariates:
            raise ValueError(
                "The matched MBatch adapters are outcome-blind; "
                "harmonization_covariates must be empty"
            )
        rscript = str(self.parameters.get("rscript", "Rscript"))
        executable = shutil.which(rscript) if "/" not in rscript else rscript
        if not executable and rscript == "Rscript":
            environment_rscript = Path(sys.executable).with_name("Rscript")
            executable = (
                str(environment_rscript) if environment_rscript.exists() else None
            )
        if not executable or not Path(executable).exists():
            raise RuntimeError("MBatch correction requires a working Rscript")
        source_path = self._source_path()
        digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
        self.source_identity = {
            "path": str(source_path),
            "sha256": digest,
        }
        batches = _batch_values(frame, self.batch_key)
        metadata = pd.DataFrame({"batch": batches})
        script = r'''args <- commandArgs(trailingOnly=TRUE)
logDebug <- function(...) invisible(NULL)
logInfo <- function(...) invisible(NULL)
logWarn <- function(...) invisible(NULL)
stopifnotWithLogging <- function(message, condition) {
  if (!isTRUE(condition)) stop(message)
}
source(args[4], local=environment())
samples_by_genes <- as.matrix(read.table(
  args[1], header=FALSE, sep="\t", check.names=FALSE
))
meta <- read.delim(args[2], stringsAsFactors=FALSE, check.names=FALSE)
sample_names <- sprintf("sample_%06d", seq_len(nrow(samples_by_genes)))
dat <- t(samples_by_genes)
colnames(dat) <- sample_names
si <- data.frame(Batch=factor(meta$batch), row.names=sample_names)
method <- args[5]
adjusted <- switch(
  method,
  mbatch_median_polish=MP(dat, si, by="Batch", overall=FALSE),
  mbatch_empirical_bayes=EB(
    dat, si, par.prior=TRUE, by="Batch", covariates=NULL,
    prior.plots=FALSE, theNumberOfThreads=1
  ),
  mbatch_anova=AN(dat, si, by="Batch", var.adj=TRUE),
  stop(paste("unsupported method", method))
)
if (is.null(adjusted)) stop("MBatch returned NULL")
nonfinite_genes <- apply(!is.finite(adjusted), 1, any)
nonfinite_gene_count <- sum(nonfinite_genes)
nonfinite_policy <- args[7]
if (nonfinite_gene_count > 0) {
  if (nonfinite_policy != "identity_gene") {
    stop(paste(
      "MBatch produced non-finite values for", nonfinite_gene_count,
      "genes; set nonfinite_policy=identity_gene for an audited fallback"
    ))
  }
  adjusted[nonfinite_genes, ] <- dat[nonfinite_genes, , drop=FALSE]
}
write.table(
  t(adjusted), file=args[3], sep="\t", row.names=FALSE,
  col.names=FALSE, quote=FALSE
)
writeLines(
  c(R.version.string, as.character(nonfinite_gene_count), nonfinite_policy),
  con=args[6]
)
'''
        with tempfile.TemporaryDirectory(prefix="nasa_mouse_mbatch_") as tmp:
            root = Path(tmp)
            matrix_path = root / "matrix.tsv"
            metadata_path = root / "metadata.tsv"
            output_path = root / "adjusted.tsv"
            script_path = root / "mbatch_adapter.R"
            runtime_path = root / "runtime.txt"
            np.savetxt(matrix_path, values, delimiter="\t", fmt="%.10g")
            metadata.to_csv(metadata_path, sep="\t", index=False)
            script_path.write_text(script, encoding="utf-8")
            completed = subprocess.run(
                [
                    str(executable),
                    str(script_path),
                    str(matrix_path),
                    str(metadata_path),
                    str(output_path),
                    str(source_path),
                    self.method_id,
                    str(runtime_path),
                    str(self.parameters.get("nonfinite_policy", "error")),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=int(self.parameters.get("timeout_seconds", 7200)),
            )
            if completed.returncode != 0:
                detail = (completed.stderr or completed.stdout)[-4000:]
                raise RuntimeError(f"{self.method_id} failed:\n{detail}")
            adjusted = np.loadtxt(output_path, delimiter="\t", ndmin=2)
            runtime_lines = runtime_path.read_text(encoding="utf-8").splitlines()
            self.r_runtime = {
                "r": runtime_lines[0] if runtime_lines else "unknown"
            }
            self.nonfinite_audit.append(
                {
                    "genes_restored_to_input": (
                        int(runtime_lines[1]) if len(runtime_lines) > 1 else 0
                    ),
                    "policy": (
                        runtime_lines[2] if len(runtime_lines) > 2 else "unknown"
                    ),
                    "samples": int(len(values)),
                    "genes": int(values.shape[1]),
                }
            )
        if adjusted.shape != values.shape:
            raise RuntimeError(
                f"{self.method_id} returned {adjusted.shape}, expected {values.shape}"
            )
        if not np.isfinite(adjusted).all():
            raise RuntimeError(f"{self.method_id} produced non-finite expression")
        return adjusted.astype(np.float32)

    def _select_anchors(
        self, matrix: np.ndarray, frame: pd.DataFrame
    ) -> tuple[np.ndarray, pd.DataFrame]:
        limit = int(self.parameters.get("anchor_samples", 64))
        if limit <= 0 or len(frame) <= limit:
            return matrix.copy(), frame.copy()
        rng = np.random.default_rng(self.seed)
        batches = _batch_values(frame, self.batch_key)
        groups = {
            level: np.flatnonzero(batches == level)
            for level in sorted(set(batches))
        }
        selected: list[int] = []
        while len(selected) < limit and any(len(items) for items in groups.values()):
            for level, items in groups.items():
                if not len(items) or len(selected) >= limit:
                    continue
                choice = int(rng.choice(items))
                selected.append(choice)
                groups[level] = items[items != choice]
        selected = sorted(selected)
        return matrix[selected].copy(), frame.iloc[selected].reset_index(drop=True)

    def fit_transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
    ) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float32)
        frame = _metadata_frame(studies, metadata, len(values))
        self.batch_key = _resolve_batch_key(frame, self.parameters)
        self.batch_levels = sorted(set(_batch_values(frame, self.batch_key)))
        if len(self.batch_levels) < 2:
            raise ValueError(f"{self.method_id} requires at least two batches")
        self.anchor_matrix, self.anchor_metadata = self._select_anchors(values, frame)
        return self._run(values, frame)

    def transform(
        self,
        matrix: np.ndarray,
        studies: Iterable[object],
        metadata: pd.DataFrame | None,
        *,
        allow_transductive: bool,
    ) -> np.ndarray:
        if not allow_transductive:
            raise ValueError(
                f"{self.method_id} requires explicit transductive preprocessing"
            )
        if self.anchor_matrix is None or self.anchor_metadata is None:
            raise RuntimeError(f"{self.method_id} is not fitted")
        values = np.asarray(matrix, dtype=np.float32)
        frame = _metadata_frame(studies, metadata, len(values))
        combined = np.concatenate([self.anchor_matrix, values])
        combined_frame = pd.concat(
            [self.anchor_metadata, frame], ignore_index=True, sort=False
        ).fillna("__missing__")
        adjusted = self._run(combined, combined_frame)
        return adjusted[len(self.anchor_matrix) :]

    def save(self, directory: Path) -> dict[str, object]:
        if self.anchor_matrix is None or self.anchor_metadata is None:
            raise RuntimeError(f"Cannot save an unfitted {self.method_id}")
        arrays_path = directory / f"{self.method_id}_harmonizer.npz"
        metadata_path = directory / f"{self.method_id}_anchor_metadata.tsv.gz"
        np.savez_compressed(arrays_path, anchor_matrix=self.anchor_matrix)
        self.anchor_metadata.to_csv(
            metadata_path, sep="\t", index=False, compression="gzip"
        )
        return {
            **self.audit(),
            "artifact": arrays_path.name,
            "anchor_metadata": metadata_path.name,
        }

    @classmethod
    def load(cls, directory: Path, payload: dict[str, object]) -> "MBatchHarmonizer":
        harmonizer = cls(
            method_id=str(payload["method"]),
            covariates=tuple(payload.get("covariates", ())),
            parameters=dict(payload.get("parameters", {})),
            device_spec="cpu",
            seed=int(payload.get("seed", 0)),
        )
        harmonizer.batch_key = str(payload["batch_key"])
        harmonizer.batch_levels = list(map(str, payload["batch_levels"]))
        harmonizer.r_runtime = dict(payload.get("r_runtime", {}))
        harmonizer.source_identity = dict(payload.get("source_identity", {}))
        harmonizer.nonfinite_audit = list(payload.get("nonfinite_audit", []))
        arrays = np.load(directory / str(payload["artifact"]))
        harmonizer.anchor_matrix = np.asarray(
            arrays["anchor_matrix"], dtype=np.float32
        )
        harmonizer.anchor_metadata = pd.read_csv(
            directory / str(payload["anchor_metadata"]), sep="\t"
        ).fillna("__missing__")
        return harmonizer

    def audit(self) -> dict[str, object]:
        return {
            "method": self.method_id,
            "fold_behavior": "transductive_with_training_anchors",
            "batch_key": self.batch_key,
            "batch_levels": self.batch_levels,
            "covariates": list(self.covariates),
            "outcome_informed": False,
            "parameters": self.parameters,
            "seed": self.seed,
            "anchor_samples": (
                0 if self.anchor_matrix is None else len(self.anchor_matrix)
            ),
            "r_runtime": self.r_runtime,
            "source_identity": self.source_identity,
            "nonfinite_audit": self.nonfinite_audit,
        }
