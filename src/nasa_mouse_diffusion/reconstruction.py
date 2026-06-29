"""Landmark-to-target gene reconstruction models."""

from __future__ import annotations

from dataclasses import dataclass

from nasa_mouse_glare.io import require_import


@dataclass
class LinearReconstructor:
    model: object
    landmark_indices: object
    target_indices: object

    def reconstruct_full(self, landmark_matrix):
        np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
        full = np.zeros((landmark_matrix.shape[0], len(self.landmark_indices) + len(self.target_indices)), dtype="float32")
        full[:, self.landmark_indices] = landmark_matrix
        if len(self.target_indices):
            full[:, self.target_indices] = self.model.predict(landmark_matrix).astype("float32")
        return full


def train_linear_reconstructor(full_matrix, landmark_indices, target_indices, *, alpha: float = 1.0):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_linear = require_import("sklearn.linear_model", "pip install -r requirements-nasa-mouse-glare.txt")
    x = np.asarray(full_matrix[:, landmark_indices], dtype=np.float32)
    y = np.asarray(full_matrix[:, target_indices], dtype=np.float32)
    model = sklearn_linear.Ridge(alpha=float(alpha), fit_intercept=True)
    if y.shape[1] > 0:
        model.fit(x, y)
    return LinearReconstructor(model=model, landmark_indices=landmark_indices, target_indices=target_indices)


def reconstruction_metrics(reconstructor: LinearReconstructor, full_matrix):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    x = np.asarray(full_matrix[:, reconstructor.landmark_indices], dtype=np.float32)
    truth = np.asarray(full_matrix, dtype=np.float32)
    pred = reconstructor.reconstruct_full(x)
    err = pred[:, reconstructor.target_indices] - truth[:, reconstructor.target_indices]
    if err.size == 0:
        return {"target_mae": 0.0, "target_rmse": 0.0}
    return {
        "target_mae": float(np.mean(np.abs(err))),
        "target_rmse": float(np.sqrt(np.mean(err ** 2))),
    }
