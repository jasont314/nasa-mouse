"""Cross-validation, metrics, and permutation importance for TabPFN3 runs."""

from __future__ import annotations

from dataclasses import dataclass
import math

from nasa_mouse_glare.io import require_import

from .features import select_features, univariate_feature_rank
from .models import make_classifier


@dataclass
class FoldSpec:
    cv_scheme: str
    fold_id: str
    train_index: object
    test_index: object
    heldout_accession: str = ""


def make_folds(y, groups, *, cv_schemes: tuple[str, ...], random_state: int = 0):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    skms = require_import("sklearn.model_selection", "pip install scikit-learn")

    y = np.asarray(y, dtype="int64")
    groups = np.asarray(groups, dtype="object")
    folds: list[FoldSpec] = []
    min_class = int(min((y == 0).sum(), (y == 1).sum()))
    unique_groups = np.unique(groups)

    if "random" in cv_schemes and min_class >= 2:
        n_splits = min(5, min_class)
        splitter = skms.StratifiedKFold(
            n_splits=n_splits, shuffle=True, random_state=random_state
        )
        for fold, (train, test) in enumerate(splitter.split(y, y), start=1):
            folds.append(FoldSpec("random", str(fold), train, test))

    if "grouped" in cv_schemes and min_class >= 2 and unique_groups.size >= 2:
        n_splits = min(5, int(unique_groups.size), min_class)
        if n_splits >= 2:
            splitter = skms.StratifiedGroupKFold(
                n_splits=n_splits, shuffle=True, random_state=random_state
            )
            for fold, (train, test) in enumerate(splitter.split(y, y, groups), start=1):
                if np.unique(y[train]).size < 2:
                    continue
                folds.append(FoldSpec("grouped", str(fold), train, test))

    if "loo_accession" in cv_schemes and unique_groups.size >= 2:
        all_index = np.arange(y.size)
        for accession in sorted(unique_groups.astype(str).tolist()):
            test = np.flatnonzero(groups.astype(str) == accession)
            train = np.setdiff1d(all_index, test, assume_unique=False)
            if np.unique(y[train]).size < 2:
                continue
            folds.append(
                FoldSpec(
                    "loo_accession",
                    str(accession),
                    train,
                    test,
                    heldout_accession=str(accession),
                )
            )
    return folds


def metric_dict(y_true, y_prob):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    metrics = require_import("sklearn.metrics", "pip install scikit-learn")

    y_true = np.asarray(y_true, dtype="int64")
    y_prob = np.asarray(y_prob, dtype="float64")
    y_pred = (y_prob >= 0.5).astype("int64")
    out = {
        "n_test": int(y_true.size),
        "n_flight": int((y_true == 1).sum()),
        "n_ground_control": int((y_true == 0).sum()),
        "accuracy": float(metrics.accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(metrics.balanced_accuracy_score(y_true, y_pred)),
        "f1": float(metrics.f1_score(y_true, y_pred, zero_division=0)),
        "brier": float(metrics.brier_score_loss(y_true, y_prob)),
        "tn": 0,
        "fp": 0,
        "fn": 0,
        "tp": 0,
    }
    labels = [0, 1]
    cm = metrics.confusion_matrix(y_true, y_pred, labels=labels)
    out["tn"], out["fp"], out["fn"], out["tp"] = [int(v) for v in cm.ravel()]
    if np.unique(y_true).size >= 2:
        out["auroc"] = float(metrics.roc_auc_score(y_true, y_prob))
        out["auprc"] = float(metrics.average_precision_score(y_true, y_prob))
    else:
        out["auroc"] = math.nan
        out["auprc"] = math.nan
    return out


def permutation_importance(
    model,
    x_test,
    y_test,
    *,
    candidate_indices,
    genes: list[str],
    repeats: int = 3,
    random_state: int = 0,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    x_test = np.asarray(x_test, dtype="float32")
    y_test = np.asarray(y_test, dtype="int64")
    rng = np.random.default_rng(random_state)
    baseline = metric_dict(y_test, model.predict_proba(x_test)[:, 1])[
        "balanced_accuracy"
    ]
    rows = []
    for local_idx in candidate_indices:
        drops = []
        for _ in range(int(repeats)):
            permuted = x_test.copy()
            permuted[:, int(local_idx)] = rng.permutation(permuted[:, int(local_idx)])
            score = metric_dict(y_test, model.predict_proba(permuted)[:, 1])[
                "balanced_accuracy"
            ]
            drops.append(float(baseline - score))
        rows.append(
            {
                "gene_id": genes[int(local_idx)],
                "local_feature_index": int(local_idx),
                "baseline_balanced_accuracy": float(baseline),
                "permutation_repeats": int(repeats),
                "mean_decrease_balanced_accuracy": float(np.mean(drops)),
                "std_decrease_balanced_accuracy": float(np.std(drops)),
            }
        )
    return rows


def _clean_covariate_value(value) -> str:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return "unknown"
    return text.replace("\t", " ").replace("\n", " ")


def encode_covariates(obs, train_index, test_index, columns: tuple[str, ...]):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    train_parts = []
    test_parts = []
    names: list[str] = []
    train_obs = obs.iloc[train_index]
    test_obs = obs.iloc[test_index]
    for column in columns:
        if column not in obs.columns:
            continue
        train_values = train_obs[column].map(_clean_covariate_value).astype(str)
        test_values = test_obs[column].map(_clean_covariate_value).astype(str)
        categories = sorted(pd.unique(train_values).tolist())
        if len(categories) <= 1:
            continue
        for category in categories:
            feature_name = f"covariate:{column}={category}"
            names.append(feature_name)
            train_parts.append((train_values.to_numpy() == category).astype("float32")[:, None])
            test_parts.append((test_values.to_numpy() == category).astype("float32")[:, None])

    if not train_parts:
        return (
            np.zeros((len(train_index), 0), dtype="float32"),
            np.zeros((len(test_index), 0), dtype="float32"),
            [],
        )
    return np.concatenate(train_parts, axis=1), np.concatenate(test_parts, axis=1), names


def run_cv(
    *,
    dataset,
    expression,
    feature_mode: str,
    cv_schemes: tuple[str, ...],
    backend: str,
    device: str,
    n_estimators: int,
    hvg_top_n: int,
    min_expr_fraction: float,
    max_features: int,
    importance_candidates: int,
    permutation_repeats: int,
    random_state: int,
    covariate_columns: tuple[str, ...] = (),
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    y = dataset.obs["tabpfn_label"].to_numpy(dtype="int64")
    groups = dataset.obs["id.accession"].astype(str).to_numpy()
    folds = make_folds(y, groups, cv_schemes=cv_schemes, random_state=random_state)
    if not folds:
        raise ValueError(f"{dataset.dataset_id}: no valid CV folds")

    metric_rows = []
    prediction_rows = []
    importance_rows = []
    for fold_number, fold in enumerate(folds, start=1):
        selection = select_features(
            expression[fold.train_index],
            y[fold.train_index],
            dataset.genes,
            feature_mode=feature_mode,
            min_expr_fraction=min_expr_fraction,
            hvg_top_n=hvg_top_n,
            max_features=max_features,
        )
        x_train = expression[fold.train_index][:, selection.indices]
        x_test = expression[fold.test_index][:, selection.indices]
        cov_train, cov_test, covariate_names = encode_covariates(
            dataset.obs, fold.train_index, fold.test_index, covariate_columns
        )
        if covariate_names:
            x_train = np.concatenate([x_train, cov_train], axis=1)
            x_test = np.concatenate([x_test, cov_test], axis=1)
        combined_features = [*selection.genes, *covariate_names]
        combined_rank = univariate_feature_rank(x_train, y[fold.train_index])
        model = make_classifier(
            backend=backend,
            device=device,
            n_estimators=n_estimators,
            random_state=random_state + fold_number,
        )
        model.fit(x_train, y[fold.train_index])
        prob = model.predict_proba(x_test)[:, 1]

        metrics = metric_dict(y[fold.test_index], prob)
        metric_rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "tissue": dataset.tissue,
                "split_group": dataset.split_group,
                "feature_mode": feature_mode,
                "cv_scheme": fold.cv_scheme,
                "fold_id": fold.fold_id,
                "heldout_accession": fold.heldout_accession,
                "backend": backend,
                "n_selected_features": int(len(combined_features)),
                "n_selected_genes": int(len(selection.genes)),
                "n_covariate_features": int(len(covariate_names)),
                **metrics,
            }
        )

        test_obs = dataset.obs.iloc[fold.test_index].copy()
        fold_predictions = pd.DataFrame(
            {
                "dataset_id": dataset.dataset_id,
                "feature_mode": feature_mode,
                "cv_scheme": fold.cv_scheme,
                "fold_id": fold.fold_id,
                "profile_id": test_obs["profile_id"].astype(str).to_numpy(),
                "id.accession": test_obs["id.accession"].astype(str).to_numpy(),
                "condition_inferred": test_obs["condition_inferred"].astype(str).to_numpy(),
                "y_true": y[fold.test_index],
                "p_flight": prob,
                "y_pred": (prob >= 0.5).astype("int64"),
            }
        )
        prediction_rows.append(fold_predictions)

        if importance_candidates > 0 and permutation_repeats > 0:
            candidate_count = min(int(importance_candidates), len(combined_features))
            candidates = combined_rank[:candidate_count]
            for row in permutation_importance(
                model,
                x_test,
                y[fold.test_index],
                candidate_indices=candidates,
                genes=combined_features,
                repeats=permutation_repeats,
                random_state=random_state + fold_number,
            ):
                row.update(
                    {
                        "dataset_id": dataset.dataset_id,
                        "tissue": dataset.tissue,
                        "split_group": dataset.split_group,
                        "feature_mode": feature_mode,
                        "cv_scheme": fold.cv_scheme,
                        "fold_id": fold.fold_id,
                        "heldout_accession": fold.heldout_accession,
                        "backend": backend,
                    }
                )
                importance_rows.append(row)

    predictions = pd.concat(prediction_rows, ignore_index=True) if prediction_rows else pd.DataFrame()
    metrics_frame = pd.DataFrame(metric_rows)
    importance = pd.DataFrame(importance_rows)
    return metrics_frame, predictions, importance

