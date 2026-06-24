"""Run GLARE-original-style post-analysis for aggregate liver runs.

This intentionally avoids the later cross-condition cluster split analysis
(for example FLT14_not_GC8). It follows the GLARE manuscript workflow more
closely: verify FLT-vs-GC separability with melted gene rows, run SHAP on the
verification classifier, then interpret direct FLT and GC consensus clusters.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
import shap
import xgboost
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_RUN_DIR = "outputs/glare_tms_liver_mober_ribo6_osdr_12_muscle_outliers"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.fillna(False).map(
        lambda value: str(value).strip().lower() in {"true", "1", "yes"}
    )


def clean_label(row: pd.Series) -> str:
    symbol = row.get("gene_symbol")
    if isinstance(symbol, str) and symbol and symbol != "NA":
        return symbol
    return str(row["gene_id"])


def safe_mean(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.mean())


def load_target(run_dir: Path) -> dict[str, Any]:
    target = np.load(run_dir / "controlled_target.npz")
    retained = pd.read_csv(run_dir / "retained_profile_features.tsv", sep="\t")
    return {
        "flt": target["flt"].astype(np.float32),
        "gc": target["gc"].astype(np.float32),
        "genes": target["genes"].astype(str),
        "flt_features": target["flt_features"].astype(str),
        "gc_features": target["gc_features"].astype(str),
        "retained": retained,
    }


def select_balanced_features(run: dict[str, Any], output_dir: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    retained = run["retained"].copy()
    flt_feature_to_index = {
        feature: index for index, feature in enumerate(run["flt_features"].tolist())
    }
    gc_feature_to_index = {
        feature: index for index, feature in enumerate(run["gc_features"].tolist())
    }
    rows = []
    flt_indices: list[int] = []
    gc_indices: list[int] = []
    balanced_feature_names: list[str] = []

    for accession, accession_frame in retained.groupby("h5_accession", sort=True):
        flt_rows = accession_frame.loc[accession_frame["location"].eq("FLT")].copy()
        gc_rows = accession_frame.loc[accession_frame["location"].eq("GC")].copy()
        flt_rows = flt_rows.sort_values(["accession_order", "feature"])
        gc_rows = gc_rows.sort_values(["accession_order", "feature"])
        n_keep = min(len(flt_rows), len(gc_rows))
        for within, (_, flt_row), (_, gc_row) in zip(
            range(1, n_keep + 1),
            flt_rows.head(n_keep).iterrows(),
            gc_rows.head(n_keep).iterrows(),
        ):
            flt_feature = str(flt_row["feature"])
            gc_feature = str(gc_row["feature"])
            flt_indices.append(flt_feature_to_index[flt_feature])
            gc_indices.append(gc_feature_to_index[gc_feature])
            balanced_name = f"{accession}_balanced_{within:02d}"
            balanced_feature_names.append(balanced_name)
            rows.append(
                {
                    "h5_accession": accession,
                    "balanced_feature": balanced_name,
                    "flt_feature": flt_feature,
                    "gc_feature": gc_feature,
                    "flt_profile": flt_row.get("profile", ""),
                    "gc_profile": gc_row.get("profile", ""),
                    "selection": "kept",
                }
            )
        for _, row in pd.concat([flt_rows.iloc[n_keep:], gc_rows.iloc[n_keep:]]).iterrows():
            rows.append(
                {
                    "h5_accession": accession,
                    "balanced_feature": "",
                    "flt_feature": str(row["feature"]) if row["location"] == "FLT" else "",
                    "gc_feature": str(row["feature"]) if row["location"] == "GC" else "",
                    "flt_profile": row.get("profile", "") if row["location"] == "FLT" else "",
                    "gc_profile": row.get("profile", "") if row["location"] == "GC" else "",
                    "selection": "dropped_to_balance_within_accession",
                }
            )

    selection = pd.DataFrame(rows)
    selection.to_csv(output_dir / "balanced_verification_features.tsv", sep="\t", index=False)
    flt_balanced = run["flt"][:, flt_indices]
    gc_balanced = run["gc"][:, gc_indices]
    return flt_balanced, gc_balanced, balanced_feature_names


def new_classifier(seed: int, n_estimators: int) -> xgboost.XGBClassifier:
    return xgboost.XGBClassifier(
        n_estimators=n_estimators,
        max_depth=10,
        n_jobs=-1,
        random_state=seed,
        tree_method="hist",
        eval_metric="logloss",
    )


def run_cv(expression: np.ndarray, labels: np.ndarray, groups: np.ndarray, seed: int, n_estimators: int) -> pd.DataFrame:
    records = []
    schemes = [
        ("random_kfold_glare", KFold(n_splits=5, shuffle=True, random_state=2023).split(expression)),
        ("gene_grouped_kfold_audit", GroupKFold(n_splits=5).split(expression, labels, groups)),
    ]
    for scheme, splits in schemes:
        for fold, (train, test) in enumerate(splits, start=1):
            model = new_classifier(seed, n_estimators)
            model.fit(expression[train], labels[train])
            probability = model.predict_proba(expression[test])[:, 1]
            prediction = (probability >= 0.5).astype(np.int8)
            row = {
                "scheme": scheme,
                "fold": fold,
                "train_rows": int(len(train)),
                "test_rows": int(len(test)),
                "accuracy": float(accuracy_score(labels[test], prediction)),
                "f1": float(f1_score(labels[test], prediction)),
                "roc_auc": float(roc_auc_score(labels[test], probability)),
            }
            records.append(row)
            log(f"{scheme} fold {fold}: AUC={row['roc_auc']:.4f} F1={row['f1']:.4f}")
    return pd.DataFrame(records)


def run_verification(run: dict[str, Any], output_dir: Path, seed: int, n_estimators: int) -> dict[str, Any]:
    verification_dir = output_dir / "verification"
    verification_dir.mkdir(parents=True, exist_ok=True)
    flt, gc, features = select_balanced_features(run, verification_dir)
    genes = run["genes"].astype(str).tolist()
    expression = np.vstack([flt, gc])
    labels = np.concatenate(
        [np.ones(len(genes), dtype=np.int8), np.zeros(len(genes), dtype=np.int8)]
    )
    gene_ids = np.asarray(genes + genes, dtype=str)
    groups = np.concatenate([np.arange(len(genes)), np.arange(len(genes))])

    log("Running GLARE-style XGBoost verification")
    folds = run_cv(expression, labels, groups, seed, n_estimators)
    folds.to_csv(verification_dir / "xgboost_verification_folds.tsv", sep="\t", index=False)
    summary = (
        folds.groupby("scheme")[["accuracy", "f1", "roc_auc"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary.columns = [
        "scheme" if col[0] == "scheme" else f"{col[0]}_{col[1]}"
        for col in summary.columns
    ]
    summary.to_csv(verification_dir / "xgboost_verification_summary.tsv", sep="\t", index=False)

    log("Running SHAP on final standardized XGBoost model")
    scaler = StandardScaler()
    scaled = scaler.fit_transform(expression)
    model = new_classifier(seed, n_estimators)
    model.fit(scaled, labels)
    probability = model.predict_proba(scaled)[:, 1]
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(scaled)
    values = np.asarray(shap_values.values, dtype=np.float32)
    np.savez_compressed(
        verification_dir / "shap_values.npz",
        values=values,
        base_values=np.asarray(shap_values.base_values, dtype=np.float32),
        data=np.asarray(shap_values.data, dtype=np.float32),
        gene_ids=gene_ids,
        labels=labels,
        features=np.asarray(features, dtype=str),
    )
    pd.DataFrame(
        {
            "feature": features,
            "mean_absolute_shap": np.mean(np.abs(values), axis=0),
            "mean_signed_shap": np.mean(values, axis=0),
            "xgboost_gain_importance": model.feature_importances_,
        }
    ).sort_values("mean_absolute_shap", ascending=False).to_csv(
        verification_dir / "shap_feature_importance.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "gene_id": gene_ids,
            "location": np.where(labels == 1, "FLT", "GC"),
            "flight_probability": probability,
            "mean_absolute_shap": np.mean(np.abs(values), axis=1),
            "sum_signed_shap": np.sum(values, axis=1),
        }
    ).sort_values("mean_absolute_shap", ascending=False).to_csv(
        verification_dir / "shap_gene_condition.tsv", sep="\t", index=False
    )
    shap.summary_plot(values, scaled, feature_names=features, max_display=20, show=False)
    plt.tight_layout()
    plt.savefig(verification_dir / "shap_beeswarm.png", dpi=180, bbox_inches="tight")
    plt.close()

    result = {
        "method": "GLARE-style melted-data XGBoost verification with within-accession balanced FLT/GC feature columns",
        "genes": int(len(genes)),
        "balanced_features": int(len(features)),
        "rows": int(len(expression)),
        "classifier": {
            "n_estimators": int(n_estimators),
            "max_depth": 10,
            "random_state": int(seed),
        },
        "metrics": summary.to_dict(orient="records"),
    }
    (verification_dir / "verification_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def load_gene_table(run_dir: Path) -> pd.DataFrame:
    candidates = [
        run_dir / "post_analysis" / "raw_count_deseq2_glare_overlap" / "deseq2_genes_with_glare_clusters.tsv",
        run_dir / "post_analysis" / "paired_cluster_report" / "paired_gene_level_table.tsv",
        run_dir / "post_analysis" / "gene_cluster_comparison.tsv",
    ]
    for path in candidates:
        if path.exists():
            table = pd.read_csv(path, sep="\t")
            break
    else:
        raise FileNotFoundError("Could not find a gene-level cluster table")
    for column in ("eligible_meta", "significant_fdr05_abs_log2fc1"):
        if column in table:
            table[column] = to_bool(table[column])
        else:
            table[column] = False
    if "gene_symbol" not in table:
        table["gene_symbol"] = ""
    table["gene_label"] = table.apply(clean_label, axis=1)
    return table


def direct_cluster_summary(location: str, clusters: pd.DataFrame, gene_table: pd.DataFrame) -> pd.DataFrame:
    cluster_col = f"{location.lower()}_cluster"
    if cluster_col not in gene_table.columns:
        cluster_col = "consensus"
    rows = []
    for cluster, cluster_genes in clusters.groupby("consensus", sort=True):
        merged = cluster_genes[["gene_id", "consensus"]].merge(
            gene_table.drop(columns=["consensus"], errors="ignore"),
            on="gene_id",
            how="left",
        )
        sig = merged.loc[merged["significant_fdr05_abs_log2fc1"]]
        top_sig = (
            sig.assign(abs_lfc=sig.get("meta_log2_fold_change", pd.Series(dtype=float)).abs())
            .sort_values("abs_lfc", ascending=False)
            if "meta_log2_fold_change" in sig
            else sig
        )
        rows.append(
            {
                "location": location,
                "cluster": int(cluster),
                "gene_count": int(len(merged)),
                "eligible_dgea_genes": int(merged["eligible_meta"].sum()),
                "significant_dgea_genes": int(len(sig)),
                "significant_up_genes": int((sig.get("meta_log2_fold_change", pd.Series(dtype=float)) > 0).sum()),
                "significant_down_genes": int((sig.get("meta_log2_fold_change", pd.Series(dtype=float)) < 0).sum()),
                "mean_meta_log2fc_sig": safe_mean(sig.get("meta_log2_fold_change", pd.Series(dtype=float))),
                "included_for_metascape": bool(10 <= len(merged) <= 3000),
                "top_significant_genes": ",".join(top_sig["gene_label"].head(20).astype(str)),
            }
        )
    return pd.DataFrame(rows)


def write_direct_cluster_analysis(run_dir: Path, output_dir: Path) -> dict[str, Any]:
    cluster_dir = output_dir / "direct_clusters"
    metascape_dir = cluster_dir / "metascape_gene_lists"
    metascape_dir.mkdir(parents=True, exist_ok=True)
    gene_table = load_gene_table(run_dir)
    summaries = []
    lists: dict[str, pd.Series] = {}

    for location in ("FLT", "GC"):
        clusters = pd.read_csv(run_dir / "clustering" / f"{location}_gene_clusters.tsv", sep="\t")
        annotated = clusters.merge(
            gene_table[
                [
                    "gene_id",
                    "gene_symbol",
                    "gene_label",
                    "eligible_meta",
                    "significant_fdr05_abs_log2fc1",
                    "meta_log2_fold_change",
                    "consistent_direction",
                ]
            ],
            on="gene_id",
            how="left",
        )
        annotated.to_csv(cluster_dir / f"{location}_gene_clusters_annotated.tsv", sep="\t", index=False, na_rep="NA")
        summaries.append(direct_cluster_summary(location, clusters, gene_table))
        for cluster, group in annotated.groupby("consensus", sort=True):
            if 10 <= len(group) <= 3000:
                name = f"{location}_cluster_{int(cluster):02d}"
                values = group["gene_label"].dropna().astype(str)
                lists[name] = values.reset_index(drop=True)
                (metascape_dir / f"{name}.txt").write_text(
                    "\n".join(values.tolist()) + "\n", encoding="utf-8"
                )

    summary = pd.concat(summaries, ignore_index=True).sort_values(
        ["significant_dgea_genes", "gene_count"], ascending=[False, False]
    )
    summary.to_csv(cluster_dir / "direct_cluster_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(lists).to_csv(metascape_dir / "metascape_direct_cluster_gene_lists.csv", index=False)
    pd.DataFrame(
        [{"list": name, "genes": int(values.nunique())} for name, values in lists.items()]
    ).to_csv(metascape_dir / "metascape_direct_cluster_manifest.tsv", sep="\t", index=False)
    background = gene_table["gene_label"].dropna().astype(str).drop_duplicates()
    (metascape_dir / "metascape_background_all_glare_genes.txt").write_text(
        "\n".join(background.tolist()) + "\n", encoding="utf-8"
    )
    return {
        "direct_cluster_count": int(len(summary)),
        "metascape_list_count": int(len(lists)),
        "top_direct_clusters": summary.head(10).to_dict(orient="records"),
    }


def write_report(output_dir: Path, verification: dict[str, Any], cluster_result: dict[str, Any]) -> None:
    summary_path = output_dir / "verification" / "xgboost_verification_summary.tsv"
    direct_path = output_dir / "direct_clusters" / "direct_cluster_summary.tsv"
    summary = pd.read_csv(summary_path, sep="\t") if summary_path.exists() else pd.DataFrame()
    direct = pd.read_csv(direct_path, sep="\t")
    metric_lines = (
        summary.to_csv(sep="\t", index=False).strip().splitlines()
        if not summary.empty
        else ["verification_not_run"]
    )
    top_lines = direct.head(12).to_csv(sep="\t", index=False).strip().splitlines()
    lines = [
        "# GLARE Original-Style Analysis",
        "",
        "This run intentionally follows the original GLARE post-analysis style:",
        "direct FLT and GC representations are analyzed separately, direct",
        "consensus clusters are interpreted, and XGBoost/SHAP is used as a",
        "verification study. It does not use the later paired-cluster split",
        "analysis such as `FLT14_not_GC8`.",
        "",
        "Because the 12-filter aggregate liver target has unequal profile counts",
        "after filtering (`73` FLT and `71` GC), the verification feature matrix",
        "was balanced within each OSD accession before running the GLARE-style",
        "melted-data classifier. This keeps the method usable without mixing",
        "studies as feature positions.",
        "",
        "## Verification",
        "",
        f"- Genes: {verification['genes']:,}",
        f"- Balanced feature columns: {verification['balanced_features']:,}",
        f"- Rows in melted classifier table: {verification['rows']:,}",
        "",
        "```tsv",
        *metric_lines,
        "```",
        "",
        "## Direct Cluster Interpretation",
        "",
        f"- Direct FLT/GC clusters summarized: {cluster_result['direct_cluster_count']}",
        f"- Metascape-eligible direct cluster lists: {cluster_result['metascape_list_count']}",
        "",
        "```tsv",
        *top_lines,
        "```",
        "",
        "## Outputs",
        "",
        "- `verification/xgboost_verification_summary.tsv`",
        "- `verification/shap_feature_importance.tsv`",
        "- `verification/shap_gene_condition.tsv`",
        "- `verification/shap_beeswarm.png`",
        "- `direct_clusters/direct_cluster_summary.tsv`",
        "- `direct_clusters/metascape_gene_lists/metascape_direct_cluster_gene_lists.csv`",
    ]
    (output_dir / "GLARE_ORIGINAL_STYLE_ANALYSIS.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GLARE-original-style verification and direct cluster analysis."
    )
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--output-dir",
        help="Defaults to <run-dir>/post_analysis/glare_original_style.",
    )
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument(
        "--skip-verification",
        action="store_true",
        help="Only export direct FLT/GC cluster interpretation files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else run_dir / "post_analysis" / "glare_original_style"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    run = load_target(run_dir)

    if args.skip_verification:
        verification = {
            "genes": int(len(run["genes"])),
            "balanced_features": 0,
            "rows": 0,
            "metrics": [],
        }
    else:
        verification = run_verification(run, output_dir, args.seed, args.n_estimators)
    cluster_result = write_direct_cluster_analysis(run_dir, output_dir)
    write_report(output_dir, verification, cluster_result)
    log(f"Saved GLARE original-style analysis to {output_dir}")


if __name__ == "__main__":
    main()
