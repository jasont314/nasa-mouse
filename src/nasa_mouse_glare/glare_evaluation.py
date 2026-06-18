"""GLARE-style evaluation and verification for the mouse workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import dense_matrix, load_matrix_bundle, require_import


DEFAULT_REPRESENTATION = "outputs/glare_hpt_tms_facs_osdr/FTSAE_representation.npy"
DEFAULT_TARGET_MANIFEST = "data/processed/tms_facs_osdr_aligned.target.manifest.json"
DEFAULT_POST_DIR = "outputs/glare_hpt_tms_facs_osdr/post_finetune"
DEFAULT_OUTPUT_DIR = "outputs/glare_hpt_tms_facs_osdr/post_finetune/evaluation"


def stable_sample_indices(n_rows: int, max_rows: int, seed: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    if max_rows <= 0 or n_rows <= max_rows:
        return np.arange(n_rows)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_rows, size=max_rows, replace=False))


def scaled_array(values, scale: bool = True):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    if not scale:
        return np.asarray(values, dtype="float32")
    StandardScaler = require_import(
        "sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt"
    ).StandardScaler
    return StandardScaler().fit_transform(values).astype("float32", copy=False)


def compute_pca(values, n_components: int, seed: int):
    PCA = require_import(
        "sklearn.decomposition", "pip install -r requirements-nasa-mouse-glare.txt"
    ).PCA
    n_components = min(n_components, values.shape[0], values.shape[1])
    model = PCA(n_components=n_components, svd_solver="randomized", random_state=seed)
    return model.fit_transform(values).astype("float32", copy=False), model


def kmeans_labels(values, n_clusters: int, seed: int):
    KMeans = require_import(
        "sklearn.cluster", "pip install -r requirements-nasa-mouse-glare.txt"
    ).KMeans
    n_clusters = min(n_clusters, values.shape[0])
    if n_clusters < 2:
        raise SystemExit("Need at least two rows for KMeans evaluation")
    model = KMeans(n_clusters=n_clusters, random_state=seed, n_init=20)
    return model.fit_predict(values), model


def silhouette(values, labels, sample_size: int, seed: int) -> float | None:
    silhouette_score = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    ).silhouette_score
    if len(set(labels.tolist())) < 2 or values.shape[0] <= len(set(labels.tolist())):
        return None
    kwargs = {"metric": "euclidean", "random_state": seed}
    if sample_size and values.shape[0] > sample_size:
        kwargs["sample_size"] = sample_size
    return float(silhouette_score(values, labels, **kwargs))


def knn_cluster_accuracy(values, labels, folds: int, neighbors: int, seed: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    KFold = require_import(
        "sklearn.model_selection", "pip install -r requirements-nasa-mouse-glare.txt"
    ).KFold
    KNeighborsClassifier = require_import(
        "sklearn.neighbors", "pip install -r requirements-nasa-mouse-glare.txt"
    ).KNeighborsClassifier
    accuracy_score = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    ).accuracy_score

    folds = min(folds, values.shape[0])
    if folds < 2:
        return None, None
    scores = []
    splitter = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for train_idx, test_idx in splitter.split(values):
        n_neighbors = min(neighbors, len(train_idx))
        if n_neighbors < 1:
            continue
        model = KNeighborsClassifier(
            n_neighbors=n_neighbors,
            metric="cosine",
            n_jobs=-1,
        )
        model.fit(values[train_idx], labels[train_idx])
        pred = model.predict(values[test_idx])
        scores.append(accuracy_score(labels[test_idx], pred))
    if not scores:
        return None, None
    return float(np.mean(scores)), float(np.std(scores))


def trustworthiness_value(source_values, representation, neighbors: int) -> float | None:
    trustworthiness = require_import(
        "sklearn.manifold", "pip install -r requirements-nasa-mouse-glare.txt"
    ).trustworthiness
    if source_values.shape[0] <= 2:
        return None
    n_neighbors = min(neighbors, max(1, (source_values.shape[0] - 1) // 2))
    return float(trustworthiness(source_values, representation, n_neighbors=n_neighbors))


def evaluate_representations(args, output_dir: Path) -> dict[str, str]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    bundle = load_matrix_bundle(args.target_manifest)
    raw = dense_matrix(bundle.matrix, max_dense_gb=args.max_dense_gb)
    latent = np.load(args.representation)
    if latent.shape[0] != raw.shape[0]:
        raise SystemExit(
            "Representation rows do not match target genes: "
            f"{latent.shape[0]} vs {raw.shape[0]}"
        )

    eval_idx = stable_sample_indices(raw.shape[0], args.max_eval_rows, args.seed)
    raw_eval = raw[eval_idx]
    latent_eval = latent[eval_idx]
    raw_eval_scaled = scaled_array(raw_eval, scale=not args.no_scale)
    latent_eval_scaled = scaled_array(latent_eval, scale=not args.no_scale)

    max_raw_pca = max(args.raw_pca_components)
    raw_pca, raw_pca_model = compute_pca(raw_eval_scaled, max_raw_pca, args.seed)

    max_latent_pca = max(args.latent_pca_components)
    latent_pca, latent_pca_model = compute_pca(latent_eval_scaled, max_latent_pca, args.seed)

    representations = {
        "raw_pca_2d": raw_pca[:, : min(2, raw_pca.shape[1])],
        "ftsae_pca_2d": latent_pca[:, : min(2, latent_pca.shape[1])],
        "ftsae_latent": latent_eval_scaled,
    }
    for n_components in args.raw_pca_components:
        if n_components != 2:
            representations[f"raw_pca_{n_components}d"] = raw_pca[:, :n_components]
    for n_components in args.latent_pca_components:
        if n_components != 2:
            representations[f"ftsae_pca_{n_components}d"] = latent_pca[:, :n_components]

    rows = []
    cluster_rows = [{"gene": bundle.genes[int(idx)]} for idx in eval_idx]
    for name, values in representations.items():
        labels, kmeans_model = kmeans_labels(values, args.n_clusters, args.seed)
        sil = silhouette(values, labels, args.silhouette_sample_size, args.seed)
        knn_mean, knn_std = knn_cluster_accuracy(
            values,
            labels,
            args.cv_folds,
            args.knn_neighbors,
            args.seed,
        )
        trust = trustworthiness_value(raw_eval_scaled, values, args.trust_neighbors)
        rows.append(
            {
                "representation": name,
                "n_rows_evaluated": int(values.shape[0]),
                "n_features": int(values.shape[1]),
                "n_clusters": int(len(set(labels.tolist()))),
                "kmeans_inertia": float(kmeans_model.inertia_),
                "silhouette": sil,
                "knn_cv_accuracy_mean": knn_mean,
                "knn_cv_accuracy_std": knn_std,
                "trustworthiness_vs_raw_expression": trust,
            }
        )
        for row, label in zip(cluster_rows, labels):
            row[f"{name}_kmeans_cluster"] = int(label)

    eval_df = pd.DataFrame(rows).sort_values("representation")
    eval_path = output_dir / "representation_evaluation.tsv"
    eval_df.to_csv(eval_path, sep="\t", index=False)

    cluster_df = pd.DataFrame(cluster_rows)
    cluster_path = output_dir / "representation_evaluation_clusters.tsv"
    cluster_df.to_csv(cluster_path, sep="\t", index=False)

    variance_rows = []
    for prefix, model in [("raw_pca", raw_pca_model), ("ftsae_pca", latent_pca_model)]:
        for idx, value in enumerate(model.explained_variance_ratio_, start=1):
            variance_rows.append(
                {
                    "representation": prefix,
                    "component": idx,
                    "explained_variance_ratio": float(value),
                }
            )
    variance_path = output_dir / "representation_pca_variance.tsv"
    pd.DataFrame(variance_rows).to_csv(variance_path, sep="\t", index=False)

    return {
        "representation_evaluation": str(eval_path),
        "representation_evaluation_clusters": str(cluster_path),
        "representation_pca_variance": str(variance_path),
    }


def make_verification_model(seed: int):
    make_pipeline = require_import(
        "sklearn.pipeline", "pip install -r requirements-nasa-mouse-glare.txt"
    ).make_pipeline
    StandardScaler = require_import(
        "sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt"
    ).StandardScaler
    SGDClassifier = require_import(
        "sklearn.linear_model", "pip install -r requirements-nasa-mouse-glare.txt"
    ).SGDClassifier
    return make_pipeline(
        StandardScaler(),
        SGDClassifier(
            loss="log_loss",
            penalty="elasticnet",
            alpha=1e-4,
            l1_ratio=0.15,
            max_iter=2000,
            tol=1e-3,
            early_stopping=True,
            class_weight="balanced",
            random_state=seed,
        ),
    )


def classifier_scores(model, values):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(values)[:, 1]
    if hasattr(model, "decision_function"):
        return model.decision_function(values)
    return model.predict(values)


def cv_metrics(values, labels, splitter, splitter_name: str, groups=None):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    metrics = require_import(
        "sklearn.metrics", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    rows = []
    split_iter = splitter.split(values, labels, groups) if groups is not None else splitter.split(values, labels)
    for fold, (train_idx, test_idx) in enumerate(split_iter, start=1):
        y_test = labels[test_idx]
        if len(set(y_test.tolist())) < 2:
            rows.append(
                {
                    "cv": splitter_name,
                    "fold": fold,
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    "accuracy": None,
                    "f1": None,
                    "roc_auc": None,
                    "status": "skipped_single_class_test_fold",
                }
            )
            continue
        model = make_verification_model(seed=1996 + fold)
        model.fit(values[train_idx], labels[train_idx])
        pred = model.predict(values[test_idx])
        scores = classifier_scores(model, values[test_idx])
        rows.append(
            {
                "cv": splitter_name,
                "fold": fold,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "accuracy": float(metrics.accuracy_score(y_test, pred)),
                "f1": float(metrics.f1_score(y_test, pred)),
                "roc_auc": float(metrics.roc_auc_score(y_test, scores)),
                "status": "ok",
            }
        )
    return rows


def write_verification(args, output_dir: Path) -> dict[str, str | None]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    model_selection = require_import(
        "sklearn.model_selection", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    metadata_path = Path(args.post_dir) / "profile_metadata.tsv"
    if not metadata_path.exists():
        raise SystemExit(f"Missing profile metadata: {metadata_path}")
    metadata = pd.read_csv(metadata_path, sep="\t", keep_default_na=False)
    if args.verification_label_col not in metadata:
        raise SystemExit(
            f"Missing verification label column '{args.verification_label_col}' "
            f"in {metadata_path}"
        )

    bundle = load_matrix_bundle(args.target_manifest)
    raw = dense_matrix(bundle.matrix, max_dense_gb=args.max_dense_gb)
    profile_values = metadata["profile"].astype(str).tolist()
    if profile_values != list(map(str, bundle.profiles)):
        raise SystemExit(
            "Profile metadata order does not match target manifest profiles. "
            "Rerun post_finetune.py before verification."
        )

    label_values = metadata[args.verification_label_col].astype(str)
    positive_mask = label_values.eq(args.positive_label)
    negative_mask = label_values.eq(args.negative_label)
    sample_mask = (positive_mask | negative_mask).to_numpy()
    if not sample_mask.any():
        raise SystemExit("No samples matched the requested verification labels")

    values = raw.T[sample_mask].astype("float32", copy=False)
    labels = positive_mask[sample_mask].astype(int).to_numpy()
    selected_metadata = metadata.loc[sample_mask].reset_index(drop=True)
    if min(np.bincount(labels)) < args.cv_folds:
        raise SystemExit(
            "Not enough samples per class for requested CV folds: "
            f"{np.bincount(labels).tolist()} samples"
        )

    rows = []
    random_cv = model_selection.StratifiedKFold(
        n_splits=args.cv_folds,
        shuffle=True,
        random_state=args.seed,
    )
    rows.extend(cv_metrics(values, labels, random_cv, "stratified_random_cv"))

    group_col = args.verification_group_col
    group_status = "not_run"
    if group_col and group_col in selected_metadata:
        groups = selected_metadata[group_col].replace("", "NA").astype(str).to_numpy()
        if len(set(groups.tolist())) >= args.cv_folds and hasattr(
            model_selection, "StratifiedGroupKFold"
        ):
            group_cv = model_selection.StratifiedGroupKFold(
                n_splits=args.cv_folds,
                shuffle=True,
                random_state=args.seed,
            )
            try:
                rows.extend(
                    cv_metrics(
                        values,
                        labels,
                        group_cv,
                        f"stratified_group_cv_by_{group_col}",
                        groups=groups,
                    )
                )
                group_status = "ok"
            except ValueError as exc:
                group_status = f"skipped: {exc}"
        else:
            group_status = (
                "skipped: not enough groups or StratifiedGroupKFold unavailable"
            )

    fold_df = pd.DataFrame(rows)
    fold_path = output_dir / "verification_folds.tsv"
    fold_df.to_csv(fold_path, sep="\t", index=False)

    summary_df = (
        fold_df[fold_df["status"].eq("ok")]
        .groupby("cv", as_index=False)
        .agg(
            n_folds=("fold", "count"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            f1_mean=("f1", "mean"),
            f1_std=("f1", "std"),
            roc_auc_mean=("roc_auc", "mean"),
            roc_auc_std=("roc_auc", "std"),
        )
    )
    summary_path = output_dir / "verification_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    final_model = make_verification_model(seed=args.seed)
    final_model.fit(values, labels)
    classifier = final_model.steps[-1][1]
    coef = getattr(classifier, "coef_", None)
    importance_path = None
    if coef is not None:
        coef = coef.ravel()
        order = np.argsort(np.abs(coef))[::-1][: args.top_features]
        importance = pd.DataFrame(
            {
                "gene": [bundle.genes[int(idx)] for idx in order],
                "coefficient": coef[order],
                "abs_coefficient": np.abs(coef[order]),
                "direction": [
                    args.positive_label if value > 0 else args.negative_label
                    for value in coef[order]
                ],
                "rank": np.arange(1, len(order) + 1),
            }
        )
        importance_path = output_dir / "verification_feature_importance.tsv"
        importance.to_csv(importance_path, sep="\t", index=False)

    shap_status = {
        "run": False,
        "reason": (
            "GLARE used XGBoost SHAP in post_pipeline.py. This environment does "
            "not have xgboost/shap installed, so this mouse evaluation writes "
            "linear classifier coefficients instead."
        ),
        "xgboost_available": bool(__import__("importlib").util.find_spec("xgboost")),
        "shap_available": bool(__import__("importlib").util.find_spec("shap")),
    }
    shap_path = output_dir / "shap_status.json"
    shap_path.write_text(json.dumps(shap_status, indent=2) + "\n", encoding="utf-8")

    label_counts = pd.DataFrame(
        {
            "label": [args.negative_label, args.positive_label],
            "encoded": [0, 1],
            "n_samples": [int((labels == 0).sum()), int((labels == 1).sum())],
        }
    )
    label_counts_path = output_dir / "verification_label_counts.tsv"
    label_counts.to_csv(label_counts_path, sep="\t", index=False)

    return {
        "verification_folds": str(fold_path),
        "verification_summary": str(summary_path),
        "verification_feature_importance": str(importance_path) if importance_path else None,
        "verification_label_counts": str(label_counts_path),
        "shap_status": str(shap_path),
        "group_cv_status": group_status,
    }


def run(args) -> Path:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}
    outputs.update(evaluate_representations(args, output_dir))
    if not args.skip_verification:
        outputs.update(write_verification(args, output_dir))

    summary = {
        "representation": str(args.representation),
        "target_manifest": str(args.target_manifest),
        "post_dir": str(args.post_dir),
        "output_dir": str(output_dir),
        "seed": args.seed,
        "n_clusters": args.n_clusters,
        "max_eval_rows": args.max_eval_rows,
        "knn_neighbors": args.knn_neighbors,
        "trust_neighbors": args.trust_neighbors,
        "verification": {
            "run": not args.skip_verification,
            "label_col": args.verification_label_col,
            "positive_label": args.positive_label,
            "negative_label": args.negative_label,
            "group_col": args.verification_group_col,
            "classifier": "sklearn SGDClassifier(log_loss), GLARE-compatible fallback",
        },
        "outputs": outputs,
    }
    summary_path = output_dir / "glare_evaluation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"representation_evaluation={outputs['representation_evaluation']}")
    if not args.skip_verification:
        print(f"verification_summary={outputs['verification_summary']}")
    print(f"summary={summary_path}")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run GLARE-style evaluation on mouse GLARE outputs."
    )
    parser.add_argument("--representation", default=DEFAULT_REPRESENTATION)
    parser.add_argument("--target-manifest", default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--post-dir", default=DEFAULT_POST_DIR)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-clusters", type=int, default=15)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--max-eval-rows", type=int, default=5000)
    parser.add_argument("--max-dense-gb", type=float, default=8.0)
    parser.add_argument("--silhouette-sample-size", type=int, default=5000)
    parser.add_argument("--knn-neighbors", type=int, default=500)
    parser.add_argument("--trust-neighbors", type=int, default=500)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--raw-pca-components", nargs="+", type=int, default=[2, 10])
    parser.add_argument("--latent-pca-components", nargs="+", type=int, default=[2, 10])
    parser.add_argument("--no-scale", action="store_true")
    parser.add_argument("--skip-verification", action="store_true")
    parser.add_argument("--verification-label-col", default="condition_inferred")
    parser.add_argument("--positive-label", default="flight")
    parser.add_argument("--negative-label", default="ground_control")
    parser.add_argument("--verification-group-col", default="id.accession")
    parser.add_argument("--top-features", type=int, default=200)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
