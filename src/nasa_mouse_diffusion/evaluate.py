"""Generated-expression quality metrics for diffusion runs."""

from __future__ import annotations

from nasa_mouse_glare.io import require_import


def corrcoef_safe(x, y):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def generated_quality(real, fake, *, max_pr_samples: int = 2000):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    real = np.asarray(real, dtype=np.float32)
    fake = np.asarray(fake, dtype=np.float32)
    n = min(real.shape[0], fake.shape[0])
    real = real[:n]
    fake = fake[:n]
    real_mean = real.mean(axis=0)
    fake_mean = fake.mean(axis=0)
    real_std = real.std(axis=0)
    fake_std = fake.std(axis=0)
    metrics = {
        "n_real": int(real.shape[0]),
        "n_fake": int(fake.shape[0]),
        "genes": int(real.shape[1]),
        "gene_mean_correlation": corrcoef_safe(real_mean, fake_mean),
        "gene_std_correlation": corrcoef_safe(real_std, fake_std),
        "mean_rmse": float(np.sqrt(np.mean((real_mean - fake_mean) ** 2))),
        "std_rmse": float(np.sqrt(np.mean((real_std - fake_std) ** 2))),
        "real_global_mean": float(real.mean()),
        "fake_global_mean": float(fake.mean()),
        "real_global_std": float(real.std()),
        "fake_global_std": float(fake.std()),
    }
    metrics.update(manifold_metrics(real, fake, max_samples=max_pr_samples))
    return metrics


def manifold_metrics(real, fake, *, max_samples: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_decomposition = require_import("sklearn.decomposition", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_neighbors = require_import("sklearn.neighbors", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_linalg = require_import("scipy.linalg", "pip install -r requirements-nasa-mouse-glare.txt")
    n = min(real.shape[0], fake.shape[0], int(max_samples))
    if n < 5:
        return {"precision": float("nan"), "recall": float("nan"), "frechet_pca": float("nan"), "adversarial_accuracy": float("nan")}
    rng = np.random.default_rng(2026)
    real_idx = rng.choice(real.shape[0], n, replace=False)
    fake_idx = rng.choice(fake.shape[0], n, replace=False)
    real = real[real_idx]
    fake = fake[fake_idx]
    dim = min(50, real.shape[1], n - 1)
    pca = sklearn_decomposition.PCA(n_components=dim, random_state=0)
    real_z = pca.fit_transform(real)
    fake_z = pca.transform(fake)
    k = min(10, n - 1)
    real_nn = sklearn_neighbors.NearestNeighbors(n_neighbors=k + 1).fit(real_z)
    fake_nn = sklearn_neighbors.NearestNeighbors(n_neighbors=k + 1).fit(fake_z)
    real_radius = real_nn.kneighbors(real_z, return_distance=True)[0][:, -1]
    fake_radius = fake_nn.kneighbors(fake_z, return_distance=True)[0][:, -1]
    d_fake_real, i_fake_real = real_nn.kneighbors(fake_z, n_neighbors=1, return_distance=True)
    d_real_fake, i_real_fake = fake_nn.kneighbors(real_z, n_neighbors=1, return_distance=True)
    precision = float((d_fake_real[:, 0] <= real_radius[i_fake_real[:, 0]]).mean())
    recall = float((d_real_fake[:, 0] <= fake_radius[i_real_fake[:, 0]]).mean())
    mu1, mu2 = real_z.mean(axis=0), fake_z.mean(axis=0)
    cov1, cov2 = np.cov(real_z, rowvar=False), np.cov(fake_z, rowvar=False)
    covmean = scipy_linalg.sqrtm(cov1.dot(cov2))
    if np.iscomplexobj(covmean):
        covmean = covmean.real
    frechet = float((mu1 - mu2).dot(mu1 - mu2) + np.trace(cov1 + cov2 - 2 * covmean))
    real_self = real_nn.kneighbors(real_z, n_neighbors=2, return_distance=True)[0][:, 1]
    fake_self = fake_nn.kneighbors(fake_z, n_neighbors=2, return_distance=True)[0][:, 1]
    aa = 0.5 * float((d_real_fake[:, 0] > real_self).mean() + (d_fake_real[:, 0] > fake_self).mean())
    return {
        "precision": precision,
        "recall": recall,
        "f1": float(2 * precision * recall / (precision + recall + 1e-12)),
        "frechet_pca": frechet,
        "adversarial_accuracy": aa,
    }


def reverse_validation(real, fake, labels):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_model_selection = require_import("sklearn.model_selection", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_pipeline = require_import("sklearn.pipeline", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_preprocessing = require_import("sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_neural = require_import("sklearn.neural_network", "pip install -r requirements-nasa-mouse-glare.txt")
    real = np.asarray(real, dtype=np.float32)
    fake = np.asarray(fake, dtype=np.float32)
    labels = np.asarray(labels, dtype=str)
    if len(set(labels.tolist())) < 2 or real.shape[0] < 10:
        return {"real_train_real_test_accuracy": float("nan"), "synthetic_train_real_test_accuracy": float("nan")}
    train_idx, test_idx = sklearn_model_selection.train_test_split(
        np.arange(real.shape[0]), test_size=0.3, random_state=0, stratify=labels
    )
    def clf():
        return sklearn_pipeline.make_pipeline(
            sklearn_preprocessing.StandardScaler(),
            sklearn_neural.MLPClassifier(hidden_layer_sizes=(64,), max_iter=300, random_state=0),
        )
    real_model = clf().fit(real[train_idx], labels[train_idx])
    synth_model = clf().fit(fake[train_idx], labels[train_idx])
    return {
        "real_train_real_test_accuracy": float(real_model.score(real[test_idx], labels[test_idx])),
        "synthetic_train_real_test_accuracy": float(synth_model.score(real[test_idx], labels[test_idx])),
    }
