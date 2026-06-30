"""Classifier backend adapters for TabPFN3 OSDR runs."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


@dataclass
class BackendStatus:
    backend: str
    available: bool
    device: str
    model_version: str
    message: str


def detect_device(requested: str = "auto") -> str:
    torch = None
    try:
        import torch as torch_module

        torch = torch_module
    except Exception:
        torch = None
    if requested and requested != "auto":
        return requested
    if torch is not None and bool(torch.cuda.is_available()):
        return "cuda"
    return "cpu"


def make_classifier(
    *,
    backend: str,
    device: str = "auto",
    n_estimators: int = 8,
    random_state: int = 0,
):
    backend = backend.strip().lower()
    if backend == "tabpfn":
        from tabpfn import TabPFNClassifier
        from tabpfn.constants import ModelVersion

        return TabPFNClassifier.create_default_for_version(
            ModelVersion.V3,
            n_estimators=int(n_estimators),
            device=detect_device(device),
            random_state=random_state,
            show_progress_bar=False,
        )
    if backend == "sklearn_logreg":
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        return make_pipeline(
            StandardScaler(),
            LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                random_state=random_state,
                solver="liblinear",
            ),
        )
    raise ValueError(f"Unsupported backend: {backend}")


def backend_status(
    *,
    backend: str,
    device: str = "auto",
    n_estimators: int = 1,
    random_state: int = 0,
) -> BackendStatus:
    """Check whether a backend can fit a tiny problem.

    This intentionally exercises TabPFN model loading because local TabPFN3
    inference requires accepted Prior Labs terms and a `TABPFN_TOKEN` when model
    weights are not already cached.
    """

    np = __import__("numpy")
    selected_device = detect_device(device)
    try:
        model = make_classifier(
            backend=backend,
            device=selected_device,
            n_estimators=n_estimators,
            random_state=random_state,
        )
        x = np.random.default_rng(random_state).normal(size=(20, 8)).astype("float32")
        y = np.array([0, 1] * 10, dtype="int64")
        model.fit(x, y)
        model.predict_proba(x[:2])
    except Exception as exc:  # noqa: BLE001 - surfaced in run manifest.
        token_hint = ""
        if backend == "tabpfn" and not os.environ.get("TABPFN_TOKEN"):
            token_hint = " TABPFN_TOKEN is not set."
        return BackendStatus(
            backend=backend,
            available=False,
            device=selected_device,
            model_version="v3" if backend == "tabpfn" else "",
            message=f"{type(exc).__name__}: {exc}{token_hint}",
        )
    return BackendStatus(
        backend=backend,
        available=True,
        device=selected_device,
        model_version="v3" if backend == "tabpfn" else "",
        message="backend smoke fit succeeded",
    )


def local_model_cache_status() -> str:
    candidates = [
        Path.home() / ".cache" / "tabpfn",
        Path.home() / ".cache" / "huggingface",
    ]
    existing = [str(path) for path in candidates if path.exists()]
    if not existing:
        return "no local TabPFN/Hugging Face cache detected"
    return "; ".join(existing)
