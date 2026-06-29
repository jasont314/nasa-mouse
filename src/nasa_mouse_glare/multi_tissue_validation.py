"""Paper-style validation stack for multi-tissue API GLARE outputs.

This module validates existing GLARE outputs without retraining. It mirrors
the GLARE manuscript validation layers where the mouse data support them:

- melted FLT-vs-GC XGBoost verification;
- representation QC against raw-expression PCA;
- consensus/base-clustering QC, including sampled average-linkage EAC;
- DEG-enrichment proportion by base method and consensus clusters;
- recurring intersection versus GLARE-only module-score validation;
- Metascape-ready gene-list export for validated modules.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Iterable

os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(os.cpu_count() or 1))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    f1_score,
    normalized_mutual_info_score,
    roc_auc_score,
    silhouette_score,
)
from sklearn.model_selection import GroupKFold, KFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

from .cluster_enrichment import bh_fdr, read_gmt
from .multi_tissue_reports import clean_reactome_term


DEFAULT_ROOT = "outputs/glare_multi_tissue_api"
DEFAULT_OUTPUT_DIR = "outputs/glare_multi_tissue_api/validation_stack"
DEFAULT_PANGLAO_GMT = "src/expiMap_reproducibility/metadata/PanglaoDB_markers_27_Mar_2020_mouseEID.gmt"

ARTIFACT_RE = re.compile(
    r"influenza|hiv|viral|vif|rhodopsin|"
    r"^signaling by gpcr$|^gpcr downstream signaling$|"
    r"^gpcr ligand binding$|class a1|peptide ligand binding receptors|"
    r"defensins|^immune system$",
    re.IGNORECASE,
)
OLFACTORY_RE = re.compile(r"olfactory|odorant", re.IGNORECASE)


def is_excluded_candidate(tissue: str, clean_term: str) -> bool:
    """Filter broad labels from candidate-module selection.

    Olfactory receptor biology can be relevant in liver, so liver olfactory
    terms are kept but marked as high-caution instead of being dropped.
    """
    if OLFACTORY_RE.search(clean_term):
        return tissue != "liver"
    return bool(ARTIFACT_RE.search(clean_term))


def interpretation_note(tissue: str, clean_term: str) -> str:
    if OLFACTORY_RE.search(clean_term):
        if tissue == "liver":
            return (
                "liver olfactory/chemosensory candidate; interpret with caution "
                "because large receptor gene families can dominate enrichment"
            )
        return "excluded olfactory/chemosensory label outside liver"
    if ARTIFACT_RE.search(clean_term):
        return "excluded broad or artifact-prone label"
    return ""


@dataclass(frozen=True)
class Scope:
    tissue: str
    scope: str
    run_dir: Path

    @property
    def label(self) -> str:
        return f"{self.tissue}/{self.scope}"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None, max_rows: int = 24) -> list[str]:
    if frame.empty:
        return ["No rows."]
    display = frame.copy()
    if columns is not None:
        display = display[[col for col in columns if col in display.columns]]
    display = display.head(max_rows).fillna("")
    headers = list(display.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in display.iterrows():
        values = []
        for value in row.tolist():
            if isinstance(value, float):
                values.append(f"{value:.4g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def safe_float(value) -> float:
    try:
        value = float(value)
    except Exception:
        return math.nan
    return value if math.isfinite(value) else math.nan


def discover_tissues(root: Path) -> list[str]:
    tissues = []
    for path in sorted(root.iterdir()):
        if not path.is_dir() or path.name in {"input_audit", "retina", "validation_stack"}:
            continue
        if (path / "aggregate" / "controlled_target.npz").exists():
            tissues.append(path.name)
    return tissues


def discover_scopes(root: Path, include_per_study: bool, include_mober: bool) -> list[Scope]:
    scopes: list[Scope] = []
    for tissue in discover_tissues(root):
        tissue_dir = root / tissue
        scopes.append(Scope(tissue, "aggregate", tissue_dir / "aggregate"))
        if include_mober and (tissue_dir / "aggregate_mober" / "controlled_target.npz").exists():
            scopes.append(Scope(tissue, "aggregate_mober", tissue_dir / "aggregate_mober"))
        if include_per_study:
            per_dir = tissue_dir / "per_study"
            for run_dir in sorted(per_dir.iterdir()) if per_dir.exists() else []:
                if run_dir.is_dir() and (run_dir / "controlled_target.npz").exists():
                    scopes.append(Scope(tissue, f"per_study/{run_dir.name}", run_dir))
    return scopes


def load_target(run_dir: Path) -> dict[str, object]:
    target = np.load(run_dir / "controlled_target.npz")
    retained_path = run_dir / "retained_profile_features.tsv"
    retained = read_tsv(retained_path) if retained_path.exists() else pd.DataFrame()
    return {
        "flt": target["flt"].astype(np.float32, copy=False),
        "gc": target["gc"].astype(np.float32, copy=False),
        "genes": target["genes"].astype(str),
        "flt_features": target["flt_features"].astype(str),
        "gc_features": target["gc_features"].astype(str),
        "retained": retained,
    }


def balanced_feature_indices(run: dict[str, object]) -> tuple[list[int], list[int], list[str], pd.DataFrame]:
    flt_features = [str(x) for x in run["flt_features"]]
    gc_features = [str(x) for x in run["gc_features"]]
    flt_lookup = {feature: idx for idx, feature in enumerate(flt_features)}
    gc_lookup = {feature: idx for idx, feature in enumerate(gc_features)}
    retained = run["retained"]

    rows = []
    flt_indices: list[int] = []
    gc_indices: list[int] = []
    balanced_names: list[str] = []

    if isinstance(retained, pd.DataFrame) and not retained.empty and "location" in retained:
        accession_col = "h5_accession" if "h5_accession" in retained else "id.accession"
        feature_col = "feature"
        retained = retained.copy()
        retained[feature_col] = retained[feature_col].astype(str)
        for accession, group in retained.groupby(accession_col, sort=True):
            flt_rows = group.loc[group["location"].eq("FLT")].copy()
            gc_rows = group.loc[group["location"].eq("GC")].copy()
            flt_rows = flt_rows[flt_rows[feature_col].isin(flt_lookup)].sort_values(feature_col)
            gc_rows = gc_rows[gc_rows[feature_col].isin(gc_lookup)].sort_values(feature_col)
            n_keep = min(len(flt_rows), len(gc_rows))
            for within, (_, flt_row), (_, gc_row) in zip(
                range(1, n_keep + 1),
                flt_rows.head(n_keep).iterrows(),
                gc_rows.head(n_keep).iterrows(),
            ):
                flt_feature = str(flt_row[feature_col])
                gc_feature = str(gc_row[feature_col])
                flt_indices.append(flt_lookup[flt_feature])
                gc_indices.append(gc_lookup[gc_feature])
                name = f"{accession}_balanced_{within:03d}"
                balanced_names.append(name)
                rows.append(
                    {
                        "accession": accession,
                        "balanced_feature": name,
                        "flt_feature": flt_feature,
                        "gc_feature": gc_feature,
                        "selection": "kept",
                    }
                )

    if not flt_indices:
        n_keep = min(len(flt_features), len(gc_features))
        flt_indices = list(range(n_keep))
        gc_indices = list(range(n_keep))
        balanced_names = [f"balanced_{idx + 1:03d}" for idx in range(n_keep)]
        rows = [
            {
                "accession": "",
                "balanced_feature": balanced_names[idx],
                "flt_feature": flt_features[idx],
                "gc_feature": gc_features[idx],
                "selection": "fallback_order_balanced",
            }
            for idx in range(n_keep)
        ]

    return flt_indices, gc_indices, balanced_names, pd.DataFrame(rows)


def xgb_classifier(seed: int, n_estimators: int, max_depth: int):
    import xgboost

    return xgboost.XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=0.08,
        subsample=0.9,
        colsample_bytree=0.9,
        n_jobs=4,
        random_state=seed,
        tree_method="hist",
        eval_metric="logloss",
    )


def verification_for_scope(
    scope: Scope,
    output_dir: Path,
    n_estimators: int,
    max_depth: int,
    folds: int,
    seed: int,
    run_shap: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    run = load_target(scope.run_dir)
    flt_idx, gc_idx, balanced_features, selection = balanced_feature_indices(run)
    scope_out = output_dir / "verification" / scope.tissue / scope.scope.replace("/", "__")
    scope_out.mkdir(parents=True, exist_ok=True)
    selection.to_csv(scope_out / "balanced_verification_features.tsv", sep="\t", index=False)

    flt = run["flt"][:, flt_idx]
    gc = run["gc"][:, gc_idx]
    genes = np.asarray(run["genes"], dtype=str)
    x = np.vstack([flt, gc]).astype(np.float32, copy=False)
    y = np.concatenate([np.ones(len(genes), dtype=np.int8), np.zeros(len(genes), dtype=np.int8)])
    groups = np.concatenate([np.arange(len(genes)), np.arange(len(genes))])

    rows = []
    n_splits = min(folds, len(genes))
    schemes = [
        ("random_kfold_glare", KFold(n_splits=n_splits, shuffle=True, random_state=seed).split(x)),
        ("gene_grouped_kfold_audit", GroupKFold(n_splits=n_splits).split(x, y, groups)),
    ]
    for scheme, splits in schemes:
        for fold, (train, test) in enumerate(splits, start=1):
            model = xgb_classifier(seed + fold, n_estimators, max_depth)
            model.fit(x[train], y[train])
            prob = model.predict_proba(x[test])[:, 1]
            pred = (prob >= 0.5).astype(np.int8)
            rows.append(
                {
                    "tissue": scope.tissue,
                    "scope": scope.scope,
                    "scheme": scheme,
                    "fold": fold,
                    "n_train": int(len(train)),
                    "n_test": int(len(test)),
                    "n_genes": int(len(genes)),
                    "balanced_features": int(len(balanced_features)),
                    "accuracy": float(accuracy_score(y[test], pred)),
                    "f1": float(f1_score(y[test], pred)),
                    "roc_auc": float(roc_auc_score(y[test], prob)),
                }
            )
    folds_df = pd.DataFrame(rows)
    folds_df.to_csv(scope_out / "xgboost_verification_folds.tsv", sep="\t", index=False)
    summary = (
        folds_df.groupby(["tissue", "scope", "scheme"], as_index=False)
        .agg(
            n_folds=("fold", "count"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            f1_mean=("f1", "mean"),
            f1_std=("f1", "std"),
            roc_auc_mean=("roc_auc", "mean"),
            roc_auc_std=("roc_auc", "std"),
            n_genes=("n_genes", "first"),
            balanced_features=("balanced_features", "first"),
        )
    )

    final = xgb_classifier(seed, n_estimators, max_depth)
    final.fit(x, y)
    feature_importance = pd.DataFrame(
        {
            "tissue": scope.tissue,
            "scope": scope.scope,
            "balanced_feature": balanced_features,
            "xgboost_importance": final.feature_importances_,
        }
    ).sort_values("xgboost_importance", ascending=False)
    feature_importance.to_csv(scope_out / "xgboost_feature_importance.tsv", sep="\t", index=False)

    if run_shap:
        try:
            import shap

            scaled = StandardScaler().fit_transform(x)
            explainer = shap.TreeExplainer(final)
            shap_values = explainer(scaled)
            values = np.asarray(shap_values.values, dtype=np.float32)
            np.savez_compressed(
                scope_out / "shap_values.npz",
                values=values,
                labels=y,
                genes=np.concatenate([genes, genes]),
                features=np.asarray(balanced_features, dtype=str),
            )
            pd.DataFrame(
                {
                    "tissue": scope.tissue,
                    "scope": scope.scope,
                    "balanced_feature": balanced_features,
                    "mean_absolute_shap": np.mean(np.abs(values), axis=0),
                    "mean_signed_shap": np.mean(values, axis=0),
                    "xgboost_importance": final.feature_importances_,
                }
            ).sort_values("mean_absolute_shap", ascending=False).to_csv(
                scope_out / "shap_feature_importance.tsv", sep="\t", index=False
            )
            summary["shap_status"] = "ok"
        except Exception as exc:
            summary["shap_status"] = f"skipped: {exc}"
    else:
        summary["shap_status"] = "not_requested"

    return folds_df, summary


def sample_indices(n: int, max_n: int, seed: int) -> np.ndarray:
    if max_n <= 0 or n <= max_n:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max_n, replace=False))


def safe_silhouette(x: np.ndarray, labels: np.ndarray, seed: int) -> float:
    labels = np.asarray(labels)
    if len(set(labels.tolist())) < 2 or len(labels) <= len(set(labels.tolist())):
        return math.nan
    try:
        return float(silhouette_score(x, labels, sample_size=min(5000, len(labels)), random_state=seed))
    except Exception:
        return math.nan


def knn_accuracy(x: np.ndarray, labels: np.ndarray, seed: int) -> float:
    labels = np.asarray(labels)
    counts = pd.Series(labels).value_counts()
    valid_classes = counts[counts >= 2].index
    mask = np.isin(labels, valid_classes)
    if mask.sum() < 20 or len(valid_classes) < 2:
        return math.nan
    x = x[mask]
    labels = labels[mask]
    n_splits = min(5, int(pd.Series(labels).value_counts().min()))
    if n_splits < 2:
        return math.nan
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = []
    for train, test in splitter.split(x):
        model = KNeighborsClassifier(n_neighbors=min(15, len(train)), metric="cosine")
        model.fit(x[train], labels[train])
        scores.append(model.score(x[test], labels[test]))
    return float(np.mean(scores)) if scores else math.nan


def representation_qc_for_scope(scope: Scope, max_genes: int, seed: int) -> pd.DataFrame:
    run = load_target(scope.run_dir)
    rows = []
    for location, matrix_key, latent_name in [
        ("FLT", "flt", "FLT_FTSAE_representation.npy"),
        ("GC", "gc", "GC_FTSAE_representation.npy"),
    ]:
        latent_path = scope.run_dir / latent_name
        clusters_path = scope.run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        if not latent_path.exists() or not clusters_path.exists():
            continue
        raw = np.asarray(run[matrix_key], dtype=np.float32)
        latent = np.load(latent_path).astype(np.float32, copy=False)
        clusters = pd.read_csv(clusters_path, sep="\t")
        labels = clusters["consensus"].to_numpy()
        idx = sample_indices(raw.shape[0], max_genes, seed)
        raw_scaled = StandardScaler().fit_transform(raw[idx])
        latent_scaled = StandardScaler().fit_transform(latent[idx])
        n_raw_pca = min(10, raw_scaled.shape[0], raw_scaled.shape[1])
        raw_pca = PCA(n_components=n_raw_pca, random_state=seed).fit_transform(raw_scaled)
        n_latent_pca = min(10, latent_scaled.shape[0], latent_scaled.shape[1])
        latent_pca = PCA(n_components=n_latent_pca, random_state=seed).fit_transform(latent_scaled)
        labels_eval = labels[idx]

        trustworthiness = math.nan
        try:
            from sklearn.manifold import trustworthiness as trustworthiness_fn

            neighbors = min(15, max(1, (len(idx) - 1) // 2))
            trustworthiness = float(trustworthiness_fn(raw_scaled, latent_scaled, n_neighbors=neighbors))
        except Exception:
            pass

        for name, values in [
            ("raw_pca_10d", raw_pca),
            ("ftsae_latent", latent_scaled),
            ("ftsae_pca_10d", latent_pca),
        ]:
            rows.append(
                {
                    "tissue": scope.tissue,
                    "scope": scope.scope,
                    "location": location,
                    "representation": name,
                    "n_genes_evaluated": int(len(idx)),
                    "n_features": int(values.shape[1]),
                    "consensus_silhouette": safe_silhouette(values, labels_eval, seed),
                    "knn_consensus_accuracy": knn_accuracy(values, labels_eval, seed),
                    "trustworthiness_vs_raw_expression": trustworthiness if name.startswith("ftsae") else math.nan,
                }
            )
    return pd.DataFrame(rows)


def sampled_eac(labels_frame: pd.DataFrame, n_clusters: int, max_genes: int, seed: int) -> tuple[float, float]:
    cols = [col for col in ["gmm", "hdbscan", "spectral"] if col in labels_frame]
    if not cols:
        return math.nan, math.nan
    idx = sample_indices(len(labels_frame), max_genes, seed)
    base = labels_frame.iloc[idx][cols].to_numpy(dtype=int)
    consensus = labels_frame.iloc[idx]["consensus"].to_numpy(dtype=int)
    n = len(idx)
    if n < 20 or n_clusters < 2:
        return math.nan, math.nan

    coassoc = np.zeros((n, n), dtype=np.float32)
    used = np.zeros((n, n), dtype=np.float32)
    for col in range(base.shape[1]):
        labels = base[:, col]
        valid = labels >= 0
        valid_idx = np.flatnonzero(valid)
        if len(valid_idx) < 2:
            continue
        used[np.ix_(valid_idx, valid_idx)] += 1
        for cluster in np.unique(labels[valid]):
            members = np.flatnonzero(labels == cluster)
            coassoc[np.ix_(members, members)] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        similarity = np.divide(coassoc, used, out=np.zeros_like(coassoc), where=used > 0)
    np.fill_diagonal(similarity, 1.0)
    distance = 1.0 - similarity
    condensed = squareform(distance, checks=False)
    z = linkage(condensed, method="average")
    eac_labels = fcluster(z, t=n_clusters, criterion="maxclust") - 1
    return (
        float(adjusted_rand_score(consensus, eac_labels)),
        float(normalized_mutual_info_score(consensus, eac_labels)),
    )


def clustering_qc_for_scope(scope: Scope, max_eac_genes: int, seed: int) -> pd.DataFrame:
    rows = []
    for location in ("FLT", "GC"):
        clusters_path = scope.run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        latent_path = scope.run_dir / f"{location}_FTSAE_representation.npy"
        if not clusters_path.exists() or not latent_path.exists():
            continue
        clusters = pd.read_csv(clusters_path, sep="\t")
        labels = clusters["consensus"].to_numpy(dtype=int)
        latent = StandardScaler().fit_transform(np.load(latent_path).astype(np.float32, copy=False))
        n_clusters = len(set(labels.tolist()))
        eac_ari, eac_nmi = sampled_eac(clusters, n_clusters, max_eac_genes, seed)
        base_metrics = {}
        for alg in ["gmm", "hdbscan", "spectral"]:
            if alg not in clusters:
                continue
            alg_labels = clusters[alg].to_numpy(dtype=int)
            valid = alg_labels >= 0
            base_metrics[f"{alg}_clusters"] = len(set(alg_labels[valid].tolist()))
            base_metrics[f"{alg}_noise_genes"] = int((~valid).sum())
            if valid.sum() > 0:
                base_metrics[f"{alg}_ari_vs_consensus"] = float(adjusted_rand_score(labels[valid], alg_labels[valid]))
                base_metrics[f"{alg}_nmi_vs_consensus"] = float(normalized_mutual_info_score(labels[valid], alg_labels[valid]))
            else:
                base_metrics[f"{alg}_ari_vs_consensus"] = math.nan
                base_metrics[f"{alg}_nmi_vs_consensus"] = math.nan
        rows.append(
            {
                "tissue": scope.tissue,
                "scope": scope.scope,
                "location": location,
                "n_genes": int(len(clusters)),
                "consensus_clusters": int(n_clusters),
                "consensus_silhouette_latent": safe_silhouette(latent, labels, seed),
                "sampled_eac_ari_vs_consensus": eac_ari,
                "sampled_eac_nmi_vs_consensus": eac_nmi,
                **base_metrics,
            }
        )
    return pd.DataFrame(rows)


def base_cluster_dgea_for_tissue(root: Path, tissue: str, alpha: float) -> pd.DataFrame:
    dgea_dir = root / tissue / "dgea_comparison"
    rows = []
    for gene_path in sorted(dgea_dir.glob("*_gene_level_glare_dgea.tsv")):
        accession = gene_path.name.replace("_gene_level_glare_dgea.tsv", "")
        run_dir = root / tissue / "per_study" / accession
        if not run_dir.exists():
            continue
        gene_table = pd.read_csv(gene_path, sep="\t")
        gene_table["significant_padj05"] = gene_table["significant_padj05"].fillna(False).astype(bool)
        gene_table["tested_dgea"] = gene_table["tested_dgea"].fillna(False).astype(bool)
        for location in ("FLT", "GC"):
            clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
            if not clusters_path.exists():
                continue
            clusters = pd.read_csv(clusters_path, sep="\t")
            merged = clusters.merge(gene_table[["gene_id", "tested_dgea", "significant_padj05"]], on="gene_id", how="left")
            tested = merged[merged["tested_dgea"].fillna(False)]
            total_sig = int(tested["significant_padj05"].sum())
            total_not = int(len(tested) - total_sig)
            for alg in [col for col in ["gmm", "hdbscan", "spectral", "consensus"] if col in clusters]:
                alg_tested = tested[tested[alg] >= 0]
                pvals = []
                alg_rows = []
                for cluster, group in alg_tested.groupby(alg):
                    sig = int(group["significant_padj05"].sum())
                    non = int(len(group) - sig)
                    if total_sig and len(group):
                        pval = float(stats.fisher_exact([[sig, non], [total_sig - sig, total_not - non]], alternative="greater").pvalue)
                    else:
                        pval = math.nan
                    pvals.append(pval)
                    alg_rows.append(
                        {
                            "tissue": tissue,
                            "accession": accession,
                            "location": location,
                            "algorithm": alg,
                            "cluster": int(cluster),
                            "tested_dgea_genes": int(len(group)),
                            "significant_padj05_genes": sig,
                            "significant_fraction": sig / len(group) if len(group) else math.nan,
                            "fisher_p_value": pval,
                        }
                    )
                fdrs = bh_fdr(np.asarray([p if math.isfinite(p) else 1.0 for p in pvals], dtype=float)) if pvals else []
                for row, fdr in zip(alg_rows, fdrs):
                    row["fisher_fdr_bh"] = float(fdr)
                    row["significant_cluster_fdr05"] = bool(fdr < alpha)
                    rows.append(row)
    return pd.DataFrame(rows)


def class_term_candidates(root: Path, tissue: str, terms_per_class: int) -> pd.DataFrame:
    dgea_dir = root / tissue / "dgea_comparison"
    rec_path = dgea_dir / "recurring_dgea_glare_pathway_overlap.tsv"
    sig_path = dgea_dir / "significant_glare_reactome_terms_by_study.tsv"
    if not sig_path.exists():
        return pd.DataFrame()

    supported_terms: set[str] = set()
    intersection_rows = []
    if rec_path.exists():
        rec = pd.read_csv(rec_path, sep="\t")
        for row in rec.itertuples(index=False):
            clean = str(row.clean_term)
            supported_terms.add(clean)
            if not is_excluded_candidate(tissue, clean):
                intersection_rows.append(
                    {
                        "tissue": tissue,
                        "module_class": "intersection",
                        "term": row.term,
                        "clean_term": clean,
                        "interpretation_note": interpretation_note(tissue, clean),
                        "study_count": int(row.study_count),
                        "best_dgea_fdr_bh": safe_float(row.best_dgea_fdr_bh),
                        "best_glare_fdr_bh": safe_float(row.best_glare_fdr_bh),
                        "mean_wald_stat_shift_mean": safe_float(row.mean_wald_stat_shift_mean),
                    }
                )
    intersection = pd.DataFrame(intersection_rows)
    if not intersection.empty:
        intersection = intersection.sort_values(
            ["study_count", "best_dgea_fdr_bh", "best_glare_fdr_bh"],
            ascending=[False, True, True],
        ).head(terms_per_class)

    sig = pd.read_csv(sig_path, sep="\t")
    sig["clean_term"] = sig["term"].map(clean_reactome_term)
    hidden_rows = []
    for clean, group in sig.groupby("clean_term"):
        if clean in supported_terms or is_excluded_candidate(tissue, clean):
            continue
        study_count = group["accession"].nunique()
        if study_count < min(3, max(2, sig["accession"].nunique())):
            continue
        hidden_rows.append(
            {
                "tissue": tissue,
                "module_class": "glare_only",
                "term": group.sort_values("glare_best_fdr_bh").iloc[0]["term"],
                "clean_term": clean,
                "interpretation_note": interpretation_note(tissue, clean),
                "study_count": int(study_count),
                "best_dgea_fdr_bh": math.nan,
                "best_glare_fdr_bh": float(group["glare_best_fdr_bh"].min()),
                "mean_wald_stat_shift_mean": math.nan,
            }
        )
    hidden = pd.DataFrame(hidden_rows)
    if not hidden.empty:
        hidden = hidden.sort_values(["study_count", "best_glare_fdr_bh"], ascending=[False, True]).head(terms_per_class)
    return pd.concat([intersection, hidden], ignore_index=True)


def term_module_genes(enrichment: pd.DataFrame, term: str) -> dict[str, set[str]]:
    subset = enrichment[enrichment["term"].astype(str).eq(term)]
    result: dict[str, set[str]] = {}
    for accession, group in subset.groupby("accession"):
        genes: set[str] = set()
        for values in group["overlap_genes"].fillna("").astype(str):
            genes.update(gene for gene in values.split(",") if gene)
        result[str(accession)] = genes
    return result


def module_score_validation(
    root: Path,
    candidates: pd.DataFrame,
    random_sets: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_rows = []
    random_rows = []
    rng = np.random.default_rng(seed)

    for row in candidates.itertuples(index=False):
        tissue = row.tissue
        term = row.term
        enrich_path = root / tissue / "per_study" / "glare_cluster_reactome_enrichment.tsv"
        if not enrich_path.exists():
            continue
        enrichment = pd.read_csv(enrich_path, sep="\t")
        genes_by_study = term_module_genes(enrichment, term)
        for accession, module_genes in genes_by_study.items():
            run_dir = root / tissue / "per_study" / accession
            target_path = run_dir / "controlled_target.npz"
            if not target_path.exists() or len(module_genes) < 5:
                continue
            run = load_target(run_dir)
            genes = np.asarray(run["genes"], dtype=str)
            gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
            idx = np.asarray([gene_to_idx[g] for g in module_genes if g in gene_to_idx], dtype=int)
            if len(idx) < 5:
                continue
            flt_scores = np.asarray(run["flt"][idx, :].mean(axis=0), dtype=float)
            gc_scores = np.asarray(run["gc"][idx, :].mean(axis=0), dtype=float)
            effect = float(flt_scores.mean() - gc_scores.mean())
            try:
                pvalue = float(stats.ttest_ind(flt_scores, gc_scores, equal_var=False).pvalue)
            except Exception:
                pvalue = math.nan
            try:
                mw_pvalue = float(stats.mannwhitneyu(flt_scores, gc_scores, alternative="two-sided").pvalue)
            except Exception:
                mw_pvalue = math.nan
            random_effects = []
            universe = np.arange(len(genes))
            if random_sets > 0 and len(idx) < len(universe):
                for _ in range(random_sets):
                    ridx = rng.choice(universe, size=len(idx), replace=False)
                    random_effects.append(float(run["flt"][ridx, :].mean() - run["gc"][ridx, :].mean()))
            if random_effects:
                random_effects_array = np.asarray(random_effects)
                empirical_abs_p = float((np.abs(random_effects_array) >= abs(effect)).mean())
                random_mean = float(random_effects_array.mean())
                random_sd = float(random_effects_array.std(ddof=1)) if len(random_effects_array) > 1 else math.nan
                random_z = (effect - random_mean) / random_sd if random_sd and math.isfinite(random_sd) and random_sd > 0 else math.nan
            else:
                empirical_abs_p = random_mean = random_sd = random_z = math.nan
            score_rows.append(
                {
                    "tissue": tissue,
                    "module_class": row.module_class,
                    "term": term,
                    "clean_term": row.clean_term,
                    "accession": accession,
                    "module_genes": int(len(idx)),
                    "n_flight": int(run["flt"].shape[1]),
                    "n_ground": int(run["gc"].shape[1]),
                    "flight_mean_score": float(flt_scores.mean()),
                    "ground_mean_score": float(gc_scores.mean()),
                    "flight_minus_ground": effect,
                    "welch_p_value": pvalue,
                    "mannwhitney_p_value": mw_pvalue,
                    "random_sets": int(random_sets),
                    "random_effect_mean": random_mean,
                    "random_effect_sd": random_sd,
                    "random_effect_z": random_z,
                    "empirical_abs_p": empirical_abs_p,
                }
            )
            random_rows.extend(
                {
                    "tissue": tissue,
                    "module_class": row.module_class,
                    "term": term,
                    "clean_term": row.clean_term,
                    "accession": accession,
                    "random_effect": value,
                }
                for value in random_effects
            )
    scores = pd.DataFrame(score_rows)
    if not scores.empty:
        scores["welch_fdr_bh"] = bh_fdr(scores["welch_p_value"].fillna(1.0).to_numpy(dtype=float))
        scores["mannwhitney_fdr_bh"] = bh_fdr(scores["mannwhitney_p_value"].fillna(1.0).to_numpy(dtype=float))
    return scores, pd.DataFrame(random_rows)


def meta_analyze_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    rows = []
    for keys, group in scores.groupby(["tissue", "module_class", "term", "clean_term"]):
        pvals = group["welch_p_value"].fillna(1.0).clip(lower=1e-300).to_numpy(dtype=float)
        stat, combined_p = stats.combine_pvalues(pvals, method="fisher")
        effects = group["flight_minus_ground"].to_numpy(dtype=float)
        rows.append(
            {
                "tissue": keys[0],
                "module_class": keys[1],
                "term": keys[2],
                "clean_term": keys[3],
                "studies_tested": int(group["accession"].nunique()),
                "combined_welch_p_fisher": float(combined_p),
                "mean_flight_minus_ground": float(np.nanmean(effects)),
                "median_flight_minus_ground": float(np.nanmedian(effects)),
                "direction_consistency": float(np.mean(np.sign(effects) == np.sign(np.nanmean(effects)))) if len(effects) else math.nan,
                "min_empirical_abs_p": float(group["empirical_abs_p"].min()),
                "median_empirical_abs_p": float(group["empirical_abs_p"].median()),
                "mean_random_effect_z": float(group["random_effect_z"].mean()),
                "total_module_genes_median": float(group["module_genes"].median()),
            }
        )
    result = pd.DataFrame(rows)
    result["combined_welch_fdr_bh"] = bh_fdr(result["combined_welch_p_fisher"].fillna(1.0).to_numpy(dtype=float))
    return result.sort_values(["module_class", "combined_welch_fdr_bh", "median_empirical_abs_p"])


def panglao_enrichment(candidates: pd.DataFrame, scores: pd.DataFrame, root: Path, gmt_path: Path, output_dir: Path) -> pd.DataFrame:
    if candidates.empty or not gmt_path.exists():
        return pd.DataFrame()
    gene_sets = read_gmt(gmt_path)
    rows = []
    for row in candidates.itertuples(index=False):
        tissue = row.tissue
        enrich_path = root / tissue / "per_study" / "glare_cluster_reactome_enrichment.tsv"
        if not enrich_path.exists():
            continue
        enrichment = pd.read_csv(enrich_path, sep="\t")
        genes_by_study = term_module_genes(enrichment, row.term)
        genes = set().union(*genes_by_study.values()) if genes_by_study else set()
        if len(genes) < 5:
            continue
        universe = set(np.load(root / tissue / "aggregate" / "controlled_target.npz")["genes"].astype(str))
        query = genes & universe
        for gene_set in gene_sets:
            term_genes = gene_set["genes"] & universe
            overlap = query & term_genes
            if len(overlap) < 3:
                continue
            pvalue = stats.hypergeom.sf(len(overlap) - 1, len(universe), len(term_genes), len(query))
            rows.append(
                {
                    "tissue": tissue,
                    "module_class": row.module_class,
                    "term": row.term,
                    "clean_term": row.clean_term,
                    "panglao_term": gene_set["term"],
                    "module_genes": len(query),
                    "overlap": len(overlap),
                    "panglao_genes": len(term_genes),
                    "p_value": float(pvalue),
                    "overlap_genes": ",".join(sorted(overlap)),
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_bh"] = result.groupby(["tissue", "module_class", "term"], group_keys=False)["p_value"].transform(bh_fdr)
    result = result.sort_values(["tissue", "module_class", "term", "fdr_bh", "p_value"])
    result.to_csv(output_dir / "candidate_module_panglao_enrichment.tsv", sep="\t", index=False)
    return result


def write_metascape_inputs(root: Path, candidates: pd.DataFrame, output_dir: Path) -> None:
    meta_dir = output_dir / "metascape_gene_lists"
    meta_dir.mkdir(parents=True, exist_ok=True)
    lists: dict[str, pd.Series] = {}
    background: set[str] = set()
    manifest_rows = []
    for row in candidates.itertuples(index=False):
        enrich_path = root / row.tissue / "per_study" / "glare_cluster_reactome_enrichment.tsv"
        if not enrich_path.exists():
            continue
        enrichment = pd.read_csv(enrich_path, sep="\t")
        genes_by_study = term_module_genes(enrichment, row.term)
        genes = sorted(set().union(*genes_by_study.values()) if genes_by_study else set())
        if len(genes) < 5:
            continue
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", f"{row.tissue}_{row.module_class}_{row.clean_term}").strip("_")[:90]
        lists[safe_name] = pd.Series(genes)
        background.update(np.load(root / row.tissue / "aggregate" / "controlled_target.npz")["genes"].astype(str).tolist())
        manifest_rows.append(
            {
                "list": safe_name,
                "tissue": row.tissue,
                "module_class": row.module_class,
                "term": row.term,
                "clean_term": row.clean_term,
                "genes": len(genes),
            }
        )
        (meta_dir / f"{safe_name}.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    if lists:
        pd.DataFrame(lists).to_csv(meta_dir / "candidate_module_gene_lists.csv", index=False)
    pd.DataFrame(manifest_rows).to_csv(meta_dir / "candidate_module_gene_list_manifest.tsv", sep="\t", index=False)
    (meta_dir / "background_all_candidate_tissue_genes.txt").write_text(
        "\n".join(sorted(background)) + "\n",
        encoding="utf-8",
    )


def summarize_validation(
    output_dir: Path,
    verification: pd.DataFrame,
    rep: pd.DataFrame,
    clustering: pd.DataFrame,
    base_summary: pd.DataFrame,
    candidates: pd.DataFrame,
    module_meta: pd.DataFrame,
) -> None:
    lines = [
        "# Multi-Tissue GLARE Validation Stack",
        "",
        "This report validates existing API-derived multi-tissue GLARE outputs using a paper-style stack.",
        "",
        "## What Was Run",
        "",
        "- XGBoost melted FLT-vs-GC verification with both GLARE-like random folds and gene-grouped audit folds.",
        "- Representation QC for raw PCA, FT-SAE latent, and FT-SAE PCA representations.",
        "- Consensus/base-clustering QC, including sampled average-linkage EAC agreement.",
        "- DEG-enrichment proportion by GMM, HDBSCAN, Spectral, and consensus labels.",
        "- Candidate module-score validation for DGEA-intersection and GLARE-only recurring modules.",
        "- Random-gene-set controls and Metascape-ready gene-list export.",
        "- PanglaoDB marker enrichment as a cell-type-marker proxy for the paper's cell-type follow-up.",
        "- TF/stress-network validation is documented as unavailable here because the repo does not contain a curated mouse spaceflight stress/TF network.",
        "",
        "## Key Counts",
        "",
        f"- Verification summaries: {len(verification)} scheme rows.",
        f"- Representation QC rows: {len(rep)}.",
        f"- Clustering QC rows: {len(clustering)}.",
        f"- Candidate modules tested: {len(candidates)}.",
        f"- Module-score meta rows: {len(module_meta)}.",
        "",
        "## Verification Snapshot",
        "",
    ]
    if verification.empty:
        lines.append("No verification rows.")
    else:
        snap = verification.sort_values(["scope", "tissue", "scheme"]).head(24)
        lines.extend(markdown_table(snap))
    lines.extend(["", "## Strongest Module-Score Meta Results", ""])
    if module_meta.empty:
        lines.append("No module-score meta rows.")
    else:
        cols = [
            "tissue",
            "module_class",
            "clean_term",
            "studies_tested",
            "combined_welch_fdr_bh",
            "mean_flight_minus_ground",
            "direction_consistency",
            "median_empirical_abs_p",
        ]
        top = module_meta.sort_values(["module_class", "combined_welch_fdr_bh"]).groupby("module_class").head(12)
        lines.extend(markdown_table(top, cols))
    lines.extend(["", "## Interpretation Rules", ""])
    lines.extend(
        [
            "- Intersection modules have the strongest support because they recur in both per-study DGEA and GLARE cluster enrichment.",
            "- GLARE-only modules are candidate hidden modules only when module-score tests are consistent across studies and stronger than random gene sets.",
            "- Generic GPCR, defensin, viral, and similarly broad labels are excluded from candidate hidden-module selection.",
            "- Liver olfactory/chemosensory labels are retained as high-caution candidates because liver expression can be biologically relevant, but large receptor gene families can also dominate enrichment.",
            "- Sampled EAC is a scalability audit of the GLARE paper's average-linkage EAC idea, not a full dense 20k-gene co-association matrix.",
        ]
    )
    (output_dir / "VALIDATION_STACK_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> Path:
    root = Path(args.root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scopes = discover_scopes(root, args.include_per_study, args.include_mober)
    if args.scope_filter:
        pattern = re.compile(args.scope_filter)
        scopes = [scope for scope in scopes if pattern.search(scope.label)]
    log(f"Validating {len(scopes)} scopes")

    verification_frames = []
    verification_summary_frames = []
    rep_frames = []
    clustering_frames = []
    for idx, scope in enumerate(scopes, start=1):
        log(f"[{idx}/{len(scopes)}] {scope.label}")
        run_shap = args.shap_aggregate and scope.scope == "aggregate"
        try:
            folds, summary = verification_for_scope(
                scope,
                output_dir,
                args.verification_estimators,
                args.verification_max_depth,
                args.verification_folds,
                args.seed,
                run_shap,
            )
            verification_frames.append(folds)
            verification_summary_frames.append(summary)
        except Exception as exc:
            verification_summary_frames.append(
                pd.DataFrame(
                    [
                        {
                            "tissue": scope.tissue,
                            "scope": scope.scope,
                            "scheme": "verification",
                            "status": f"failed: {exc}",
                        }
                    ]
                )
            )
        try:
            rep_frames.append(representation_qc_for_scope(scope, args.max_eval_genes, args.seed))
        except Exception as exc:
            rep_frames.append(pd.DataFrame([{"tissue": scope.tissue, "scope": scope.scope, "status": f"failed: {exc}"}]))
        try:
            clustering_frames.append(clustering_qc_for_scope(scope, args.max_eac_genes, args.seed))
        except Exception as exc:
            clustering_frames.append(pd.DataFrame([{"tissue": scope.tissue, "scope": scope.scope, "status": f"failed: {exc}"}]))

    verification = pd.concat(verification_frames, ignore_index=True) if verification_frames else pd.DataFrame()
    verification_summary = pd.concat(verification_summary_frames, ignore_index=True) if verification_summary_frames else pd.DataFrame()
    rep = pd.concat(rep_frames, ignore_index=True) if rep_frames else pd.DataFrame()
    clustering = pd.concat(clustering_frames, ignore_index=True) if clustering_frames else pd.DataFrame()

    verification.to_csv(output_dir / "xgboost_verification_folds.tsv", sep="\t", index=False)
    verification_summary.to_csv(output_dir / "xgboost_verification_summary.tsv", sep="\t", index=False)
    rep.to_csv(output_dir / "representation_qc.tsv", sep="\t", index=False)
    clustering.to_csv(output_dir / "clustering_qc.tsv", sep="\t", index=False)

    base_frames = []
    candidate_frames = []
    for tissue in discover_tissues(root):
        log(f"Summarizing DGEA and candidate modules for {tissue}")
        base = base_cluster_dgea_for_tissue(root, tissue, args.alpha)
        if not base.empty:
            base_frames.append(base)
        candidate = class_term_candidates(root, tissue, args.terms_per_class)
        if not candidate.empty:
            candidate_frames.append(candidate)
    base_clusters = pd.concat(base_frames, ignore_index=True) if base_frames else pd.DataFrame()
    base_clusters.to_csv(output_dir / "base_cluster_dgea_enrichment.tsv", sep="\t", index=False)
    if base_clusters.empty:
        base_summary = pd.DataFrame()
    else:
        base_summary = (
            base_clusters.groupby(["tissue", "algorithm"], as_index=False)
            .agg(
                clusters_tested=("cluster", "count"),
                significant_clusters_fdr05=("significant_cluster_fdr05", "sum"),
                mean_significant_fraction=("significant_fraction", "mean"),
            )
            .sort_values(["tissue", "significant_clusters_fdr05"], ascending=[True, False])
        )
    base_summary.to_csv(output_dir / "base_cluster_dgea_summary.tsv", sep="\t", index=False)

    candidates = pd.concat(candidate_frames, ignore_index=True) if candidate_frames else pd.DataFrame()
    candidates.to_csv(output_dir / "candidate_modules.tsv", sep="\t", index=False)
    write_metascape_inputs(root, candidates, output_dir)
    scores, random_scores = module_score_validation(root, candidates, args.random_sets, args.seed)
    scores.to_csv(output_dir / "candidate_module_score_validation.tsv", sep="\t", index=False)
    random_scores.to_csv(output_dir / "candidate_module_random_set_effects.tsv", sep="\t", index=False)
    module_meta = meta_analyze_scores(scores)
    module_meta.to_csv(output_dir / "candidate_module_score_meta.tsv", sep="\t", index=False)
    panglao_enrichment(candidates, scores, root, Path(args.panglao_gmt), output_dir)

    summarize_validation(output_dir, verification_summary, rep, clustering, base_summary, candidates, module_meta)
    manifest = {
        "root": str(root),
        "output_dir": str(output_dir),
        "scopes": [scope.label for scope in scopes],
        "args": vars(args),
    }
    (output_dir / "validation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    log(f"Saved validation stack to {output_dir}")
    return output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--include-per-study", action="store_true")
    parser.add_argument("--include-mober", action="store_true")
    parser.add_argument("--scope-filter", default="")
    parser.add_argument("--verification-estimators", type=int, default=80)
    parser.add_argument("--verification-max-depth", type=int, default=6)
    parser.add_argument("--verification-folds", type=int, default=5)
    parser.add_argument("--max-eval-genes", type=int, default=4000)
    parser.add_argument("--max-eac-genes", type=int, default=1200)
    parser.add_argument("--terms-per-class", type=int, default=6)
    parser.add_argument("--random-sets", type=int, default=100)
    parser.add_argument("--panglao-gmt", default=DEFAULT_PANGLAO_GMT)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--shap-aggregate", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
