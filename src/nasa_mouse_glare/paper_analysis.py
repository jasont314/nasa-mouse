"""GLARE verification, SHAP, DEG, and cluster interpretation for OSD-379."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import shap
import xgboost
from scipy.stats import ttest_ind
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    f1_score,
    normalized_mutual_info_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def load_target(run_dir: Path) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    target = np.load(run_dir / "controlled_target.npz")
    return (
        target["flt"].astype(np.float32),
        target["gc"].astype(np.float32),
        target["genes"].astype(str).tolist(),
        target["features"].astype(str).tolist(),
    )


def melted_data(
    flt: np.ndarray, gc: np.ndarray, genes: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    expression = np.vstack([flt, gc])
    labels = np.concatenate(
        [np.ones(len(genes), dtype=np.int8), np.zeros(len(genes), dtype=np.int8)]
    )
    gene_ids = np.asarray(genes + genes, dtype=str)
    groups = np.concatenate([np.arange(len(genes)), np.arange(len(genes))])
    return expression, labels, gene_ids, groups


def new_classifier(seed: int, n_estimators: int) -> xgboost.XGBClassifier:
    return xgboost.XGBClassifier(
        n_estimators=n_estimators,
        max_depth=10,
        n_jobs=-1,
        random_state=seed,
        tree_method="hist",
        eval_metric="logloss",
    )


def evaluate_splits(
    expression: np.ndarray,
    labels: np.ndarray,
    splits,
    scheme: str,
    seed: int,
    n_estimators: int,
) -> list[dict]:
    records = []
    for fold, (train, test) in enumerate(splits, start=1):
        model = new_classifier(seed, n_estimators)
        model.fit(expression[train], labels[train])
        probability = model.predict_proba(expression[test])[:, 1]
        prediction = (probability >= 0.5).astype(np.int8)
        records.append(
            {
                "scheme": scheme,
                "fold": fold,
                "train_rows": int(len(train)),
                "test_rows": int(len(test)),
                "accuracy": float(accuracy_score(labels[test], prediction)),
                "f1": float(f1_score(labels[test], prediction)),
                "roc_auc": float(roc_auc_score(labels[test], probability)),
            }
        )
        log(
            f"verification {scheme} fold {fold}: "
            f"AUC={records[-1]['roc_auc']:.4f} F1={records[-1]['f1']:.4f}"
        )
    return records


def run_verification(run_dir: Path, seed: int, n_estimators: int) -> dict:
    output_dir = run_dir / "verification"
    output_dir.mkdir(parents=True, exist_ok=True)
    flt, gc, genes, features = load_target(run_dir)
    expression, labels, gene_ids, groups = melted_data(flt, gc, genes)

    random_cv = KFold(n_splits=5, shuffle=True, random_state=2023)
    grouped_cv = GroupKFold(n_splits=5)
    records = evaluate_splits(
        expression,
        labels,
        random_cv.split(expression),
        "random_kfold_glare",
        seed,
        n_estimators,
    )
    records.extend(
        evaluate_splits(
            expression,
            labels,
            grouped_cv.split(expression, labels, groups),
            "gene_grouped_kfold",
            seed,
            n_estimators,
        )
    )
    folds = pd.DataFrame(records)
    folds.to_csv(output_dir / "xgboost_verification_folds.tsv", sep="\t", index=False)
    summary_table = (
        folds.groupby("scheme")[["accuracy", "f1", "roc_auc"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary_table.columns = [
        "scheme" if column[0] == "scheme" else f"{column[0]}_{column[1]}"
        for column in summary_table.columns
    ]
    summary_table.to_csv(
        output_dir / "xgboost_verification_summary.tsv", sep="\t", index=False
    )

    log("Fitting final standardized XGBoost model for SHAP")
    scaler = StandardScaler()
    scaled = scaler.fit_transform(expression)
    model = new_classifier(seed, n_estimators)
    model.fit(scaled, labels)
    probabilities = model.predict_proba(scaled)[:, 1]
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(scaled)
    values = np.asarray(shap_values.values, dtype=np.float32)
    np.savez_compressed(
        output_dir / "shap_values.npz",
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
        output_dir / "shap_feature_importance.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "gene_id": gene_ids,
            "location": np.where(labels == 1, "FLT", "GC"),
            "flight_probability": probabilities,
            "mean_absolute_shap": np.mean(np.abs(values), axis=1),
            "sum_signed_shap": np.sum(values, axis=1),
        }
    ).sort_values("mean_absolute_shap", ascending=False).to_csv(
        output_dir / "shap_gene_condition.tsv", sep="\t", index=False
    )
    shap.summary_plot(
        values,
        scaled,
        feature_names=features,
        max_display=20,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(output_dir / "shap_beeswarm.png", dpi=180, bbox_inches="tight")
    plt.close()

    summary = {
        "method": "GLARE melted-data XGBoost verification and SHAP",
        "rows": int(len(expression)),
        "genes": int(len(genes)),
        "features": int(len(features)),
        "classifier": {
            "n_estimators": n_estimators,
            "max_depth": 10,
            "random_state": seed,
        },
        "random_cv_note": (
            "Matches GLARE's shuffled KFold but allows paired FLT/GC rows from "
            "one gene to appear across train and test folds"
        ),
        "grouped_cv_note": "Keeps both rows for every gene in the same fold",
        "metrics": summary_table.to_dict(orient="records"),
        "outputs": {
            "folds": str(output_dir / "xgboost_verification_folds.tsv"),
            "summary": str(output_dir / "xgboost_verification_summary.tsv"),
            "shap_values": str(output_dir / "shap_values.npz"),
            "shap_features": str(output_dir / "shap_feature_importance.tsv"),
            "shap_genes": str(output_dir / "shap_gene_condition.tsv"),
            "shap_plot": str(output_dir / "shap_beeswarm.png"),
        },
    }
    (output_dir / "verification_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def pooled_welch_differential_expression(
    flt: np.ndarray, gc: np.ndarray, genes: list[str]
) -> pd.DataFrame:
    log_flt = np.log2(flt + 1.0)
    log_gc = np.log2(gc + 1.0)
    statistic, p_value = ttest_ind(
        log_flt, log_gc, axis=1, equal_var=False, nan_policy="omit"
    )
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)
    fdr = multipletests(p_value, method="fdr_bh")[1]
    log2_fc = np.mean(log_flt, axis=1) - np.mean(log_gc, axis=1)
    significant = (fdr < 0.05) & (np.abs(log2_fc) >= 1.0)
    return pd.DataFrame(
        {
            "gene_id": genes,
            "mean_log2_flt": np.mean(log_flt, axis=1),
            "mean_log2_gc": np.mean(log_gc, axis=1),
            "log2_fold_change": log2_fc,
            "welch_t": statistic,
            "p_value": p_value,
            "fdr_bh": fdr,
            "significant_fdr05_abs_log2fc1": significant,
            "direction": np.where(significant, np.where(log2_fc > 0, "up", "down"), "not_deg"),
        }
    )


def official_differential_expression(
    path: Path, genes: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load the four age/euthanasia-matched NASA flight-vs-ground contrasts."""
    columns = pd.read_csv(path, nrows=0).columns.tolist()
    log2_columns = []
    for column in columns:
        if not column.startswith("Log2fc_(Space Flight") or "v(Ground Control" not in column:
            continue
        same_age = (
            column.count("10 to 12 week") == 2 or column.count("32 week") == 2
        )
        same_collection = (
            column.count("Carcass") == 2
            or column.count("Upon euthanasia") == 2
        )
        if same_age and same_collection:
            log2_columns.append(column)
    if len(log2_columns) != 4:
        raise ValueError(
            f"Expected four matched flight-vs-ground contrasts, found {len(log2_columns)}"
        )
    contrast_names = [column.removeprefix("Log2fc_") for column in log2_columns]
    p_columns = [f"P.value_{name}" for name in contrast_names]
    adjusted_columns = [f"Adj.p.value_{name}" for name in contrast_names]
    usecols = ["ENSEMBL", "SYMBOL"] + log2_columns + p_columns + adjusted_columns
    source = pd.read_csv(path, usecols=usecols)
    source["ENSEMBL"] = source["ENSEMBL"].astype(str)
    source = source.drop_duplicates("ENSEMBL").set_index("ENSEMBL")
    missing = [gene for gene in genes if gene not in source.index]
    if missing:
        raise ValueError(
            f"Official differential-expression file is missing {len(missing)} genes"
        )
    source = source.loc[genes]

    long_frames = []
    for name, log2_column, p_column, adjusted_column in zip(
        contrast_names, log2_columns, p_columns, adjusted_columns
    ):
        frame = pd.DataFrame(
            {
                "gene_id": genes,
                "gene_symbol": source["SYMBOL"].fillna("").astype(str).to_numpy(),
                "contrast": name,
                "log2_fold_change": pd.to_numeric(
                    source[log2_column], errors="coerce"
                ).to_numpy(),
                "p_value": pd.to_numeric(source[p_column], errors="coerce").to_numpy(),
                "fdr_bh": pd.to_numeric(
                    source[adjusted_column], errors="coerce"
                ).to_numpy(),
            }
        )
        frame["significant_fdr05_abs_log2fc1"] = (
            frame["fdr_bh"].lt(0.05) & frame["log2_fold_change"].abs().ge(1.0)
        )
        long_frames.append(frame)
    long = pd.concat(long_frames, ignore_index=True)

    aggregate_rows = []
    for gene_id, group in long.groupby("gene_id", sort=False):
        significant = group.loc[group["significant_fdr05_abs_log2fc1"]]
        strongest = group.loc[group["log2_fold_change"].abs().idxmax()]
        signs = set(np.sign(significant["log2_fold_change"]).astype(int))
        if not signs:
            direction = "not_deg"
        elif signs == {1}:
            direction = "up"
        elif signs == {-1}:
            direction = "down"
        else:
            direction = "mixed"
        aggregate_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": group["gene_symbol"].iloc[0],
                "log2_fold_change": float(strongest["log2_fold_change"]),
                "p_value": float(group["p_value"].min()),
                "fdr_bh": float(group["fdr_bh"].min()),
                "significant_fdr05_abs_log2fc1": bool(len(significant)),
                "significant_contrast_count": int(len(significant)),
                "direction": direction,
            }
        )
    aggregate = pd.DataFrame(aggregate_rows)
    return aggregate, long, contrast_names


def export_cluster_analysis(
    run_dir: Path, official_de_path: Path | None = None
) -> dict:
    output_dir = run_dir / "biological_analysis"
    metascape_dir = output_dir / "metascape_gene_lists"
    output_dir.mkdir(parents=True, exist_ok=True)
    metascape_dir.mkdir(parents=True, exist_ok=True)
    flt, gc, genes, _ = load_target(run_dir)
    if official_de_path:
        deg, official_long, contrast_names = official_differential_expression(
            official_de_path, genes
        )
        official_long.to_csv(
            output_dir / "official_matched_contrast_deg_long.tsv",
            sep="\t",
            index=False,
        )
        deg_method = (
            "NASA OSD-379 GLbulkRNAseq differential-expression results; union of "
            "four age- and collection-matched Space Flight vs Ground Control "
            "contrasts at adjusted p < 0.05 and absolute log2 fold change >= 1"
        )
    else:
        deg = pooled_welch_differential_expression(flt, gc, genes)
        deg.insert(1, "gene_symbol", "")
        contrast_names = []
        deg_method = (
            "Fallback Welch t-test on log2(expression + 1), Benjamini-Hochberg "
            "FDR; threshold FDR < 0.05 and absolute log2 fold change >= 1"
        )
    deg.to_csv(output_dir / "osd379_flt_vs_gc_deg.tsv", sep="\t", index=False)

    cluster_tables = {}
    proportion_rows = []
    metascape_rows = []
    metascape_lists: dict[str, list[str]] = {}
    for location in ("FLT", "GC"):
        clusters = pd.read_csv(
            run_dir / "clustering" / f"{location}_gene_clusters.tsv", sep="\t"
        )
        clusters = clusters.merge(
            deg[["gene_id", "gene_symbol"]], on="gene_id", validate="one_to_one"
        )
        clusters.to_csv(
            output_dir / f"{location}_gene_clusters_annotated.tsv",
            sep="\t",
            index=False,
        )
        cluster_tables[location] = clusters
        merged = clusters[["gene_id", "consensus"]].merge(
            deg[[
                "gene_id",
                "gene_symbol",
                "log2_fold_change",
                "fdr_bh",
                "significant_fdr05_abs_log2fc1",
                "direction",
            ]],
            on="gene_id",
            validate="one_to_one",
        )
        for cluster, group in merged.groupby("consensus"):
            significant = group["significant_fdr05_abs_log2fc1"]
            proportion_rows.append(
                {
                    "location": location,
                    "cluster": int(cluster),
                    "gene_count": int(len(group)),
                    "deg_count": int(significant.sum()),
                    "deg_proportion": float(significant.mean()),
                    "up_deg_count": int((group["direction"] == "up").sum()),
                    "down_deg_count": int((group["direction"] == "down").sum()),
                    "median_log2_fold_change": float(group["log2_fold_change"].median()),
                }
            )
            include = 10 <= len(group) <= 3000
            list_path = metascape_dir / f"{location}_cluster_{int(cluster):02d}.txt"
            if include:
                identifiers = group["gene_symbol"].where(
                    group["gene_symbol"].astype(str).str.len().gt(0),
                    group["gene_id"],
                )
                list_name = f"{location}_cluster_{int(cluster):02d}"
                metascape_lists[list_name] = identifiers.astype(str).tolist()
                list_path.write_text(
                    "\n".join(identifiers.astype(str)) + "\n", encoding="utf-8"
                )
            metascape_rows.append(
                {
                    "location": location,
                    "cluster": int(cluster),
                    "gene_count": int(len(group)),
                    "included": include,
                    "exclusion_reason": (
                        "" if include else "fewer_than_10" if len(group) < 10 else "more_than_3000"
                    ),
                    "gene_list": str(list_path) if include else "",
                }
            )

    proportions = pd.DataFrame(proportion_rows).sort_values(
        ["location", "deg_proportion"], ascending=[True, False]
    )
    proportions.to_csv(output_dir / "cluster_deg_proportions.tsv", sep="\t", index=False)
    metascape = pd.DataFrame(metascape_rows)
    metascape.to_csv(metascape_dir / "manifest.tsv", sep="\t", index=False)
    background = deg["gene_symbol"].where(
        deg["gene_symbol"].astype(str).str.len().gt(0), deg["gene_id"]
    ).astype(str).tolist()
    pd.DataFrame(
        {name: pd.Series(values) for name, values in metascape_lists.items()}
    ).to_csv(metascape_dir / "metascape_multiple_gene_lists.csv", index=False)
    (metascape_dir / "metascape_background.txt").write_text(
        "\n".join(background) + "\n", encoding="utf-8"
    )
    for location in ("FLT", "GC"):
        location_lists = {
            name: values
            for name, values in metascape_lists.items()
            if name.startswith(f"{location}_")
        }
        pd.DataFrame(
            {name: pd.Series(values) for name, values in location_lists.items()}
        ).to_csv(
            metascape_dir / f"metascape_{location.lower()}_gene_lists.csv",
            index=False,
        )
    (metascape_dir / "README.txt").write_text(
        "Upload metascape_multiple_gene_lists.csv (or one location-specific CSV) "
        "to Metascape after selecting Multiple Gene Lists. Confirm that the first "
        "row is detected as the column header, then choose Custom Analysis. Select "
        "Mus musculus as both input and analysis species. In the enrichment settings, "
        "paste metascape_background.txt into the custom background dialog. Each CSV "
        "column is one foreground list. Lists over 3,000 genes and lists under 10 "
        "genes are excluded to match GLARE and Metascape limits. Metascape has no "
        "public API, so web submission and result download are manual.\n",
        encoding="utf-8",
    )

    comparison = cluster_tables["FLT"][["gene_id", "consensus"]].merge(
        cluster_tables["GC"][["gene_id", "consensus"]],
        on="gene_id",
        suffixes=("_flt", "_gc"),
        validate="one_to_one",
    )
    contingency = pd.crosstab(
        comparison["consensus_flt"], comparison["consensus_gc"]
    )
    contingency.to_csv(output_dir / "flt_gc_cluster_contingency.tsv", sep="\t")
    significant_count = int(deg["significant_fdr05_abs_log2fc1"].sum())
    summary = {
        "deg_method": deg_method,
        "official_de_source": str(official_de_path) if official_de_path else "",
        "official_matched_contrasts": contrast_names,
        "genes_tested": int(len(deg)),
        "significant_degs": significant_count,
        "significant_up": int((deg["direction"] == "up").sum()),
        "significant_down": int((deg["direction"] == "down").sum()),
        "significant_mixed_direction": int((deg["direction"] == "mixed").sum()),
        "significant_by_official_contrast": (
            official_long.groupby("contrast")["significant_fdr05_abs_log2fc1"]
            .sum()
            .astype(int)
            .to_dict()
            if official_de_path
            else {}
        ),
        "flt_gc_cluster_adjusted_rand": float(
            adjusted_rand_score(comparison["consensus_flt"], comparison["consensus_gc"])
        ),
        "flt_gc_cluster_normalized_mutual_information": float(
            normalized_mutual_info_score(
                comparison["consensus_flt"], comparison["consensus_gc"]
            )
        ),
        "metascape": (
            "Gene lists and manifest exported; Metascape is an external web "
            "submission and is not executed locally"
        ),
        "outputs": {
            "deg": str(output_dir / "osd379_flt_vs_gc_deg.tsv"),
            "deg_proportions": str(output_dir / "cluster_deg_proportions.tsv"),
            "cluster_contingency": str(output_dir / "flt_gc_cluster_contingency.tsv"),
            "metascape_manifest": str(metascape_dir / "manifest.tsv"),
            "metascape_multiple_lists": str(
                metascape_dir / "metascape_multiple_gene_lists.csv"
            ),
            "metascape_background": str(
                metascape_dir / "metascape_background.txt"
            ),
        },
    }
    (output_dir / "biological_analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GLARE paper-style OSD-379 analyses.")
    parser.add_argument(
        "stage", choices=["verification", "post", "all"], help="Analysis stage to run"
    )
    parser.add_argument(
        "--run-dir", default="outputs/glare_paper_tms_liver_osd379"
    )
    parser.add_argument("--seed", type=int, default=1996)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument(
        "--official-de",
        default="assets/osdr/GLDS-379_rna_seq_differential_expression_GLbulkRNAseq.csv",
        help="Official NASA differential-expression table; pass an empty value for Welch fallback.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    if args.stage in {"verification", "all"}:
        run_verification(run_dir, args.seed, args.n_estimators)
    if args.stage in {"post", "all"}:
        official_de = Path(args.official_de) if args.official_de else None
        export_cluster_analysis(run_dir, official_de)


if __name__ == "__main__":
    main()
