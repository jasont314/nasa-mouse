"""Cluster expiMap pathway scores separately by condition and compare clusters."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import entropy
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


PATHWAY_PREFIX = "R-MMU"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def clean_term(term: str) -> str:
    return term.removeprefix("R-MMU-").replace("_", " ")


def canonical_condition(value: object) -> str:
    text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if text in {"flight", "flt", "space_flight", "spaceflight"}:
        return "flight"
    if text in {"ground_control", "gc", "ground"}:
        return "ground_control"
    return text


def bh_fdr(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(out, index=p_values.index)
    pv = values[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    valid_indices = np.flatnonzero(valid)
    out[valid_indices[order]] = adjusted
    return pd.Series(out, index=p_values.index)


def choose_kmeans(
    x: np.ndarray,
    *,
    min_k: int,
    max_k: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    n = x.shape[0]
    upper = min(max_k, n - 1)
    lower = min(min_k, upper)
    if upper < 2:
        return np.zeros(n, dtype=int), {
            "selected_k": 1,
            "silhouette": math.nan,
            "scores": [],
            "reason": "too_few_samples",
        }

    scores: list[dict[str, Any]] = []
    best_labels: np.ndarray | None = None
    best_score = -np.inf
    best_k = lower
    for k in range(lower, upper + 1):
        labels = KMeans(n_clusters=k, random_state=seed, n_init=50).fit_predict(x)
        if np.unique(labels).size < 2:
            score = math.nan
        else:
            score = float(silhouette_score(x, labels))
        scores.append({"k": k, "silhouette": score})
        if np.isfinite(score) and score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    if best_labels is None:
        best_labels = KMeans(n_clusters=lower, random_state=seed, n_init=50).fit_predict(x)
        best_score = math.nan
        best_k = lower
    return best_labels.astype(int), {
        "selected_k": int(best_k),
        "silhouette": None if not np.isfinite(best_score) else float(best_score),
        "scores": scores,
        "reason": "max_silhouette",
    }


def top_counts(series: pd.Series, n: int = 8) -> str:
    counts = series.fillna("NA").astype(str).value_counts()
    return ",".join(f"{idx}:{int(value)}" for idx, value in counts.head(n).items())


def cluster_summaries(
    assignments: pd.DataFrame,
    *,
    condition: str,
    accession_col: str,
) -> pd.DataFrame:
    rows = []
    subset = assignments.loc[assignments["condition"] == condition]
    for cluster, group in subset.groupby("cluster", sort=True):
        accession_counts = group[accession_col].fillna("NA").astype(str).value_counts()
        probs = accession_counts / accession_counts.sum()
        rows.append(
            {
                "condition": condition,
                "cluster": int(cluster),
                "n_samples": int(len(group)),
                "n_accessions": int(accession_counts.size),
                "accession_entropy": float(entropy(probs, base=2)),
                "top_accessions": top_counts(group[accession_col]),
            }
        )
    return pd.DataFrame(rows)


def summarize_pathway_shifts(
    centroid_delta: np.ndarray,
    pathway_cols: list[str],
    *,
    top_n: int,
) -> pd.DataFrame:
    order = np.argsort(np.abs(centroid_delta))[::-1][:top_n]
    return pd.DataFrame(
        {
            "pathway": [pathway_cols[i] for i in order],
            "pathway_label": [clean_term(pathway_cols[i]) for i in order],
            "delta_flight_minus_ground_z": centroid_delta[order],
        }
    )


def compare_clusters(
    z_scores: np.ndarray,
    assignments: pd.DataFrame,
    pathway_cols: list[str],
    *,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    centroids: dict[tuple[str, int], np.ndarray] = {}
    for (condition, cluster), index in assignments.groupby(
        ["condition", "cluster"], sort=True
    ).groups.items():
        centroids[(condition, int(cluster))] = z_scores[list(index)].mean(axis=0)

    gc_keys = sorted(key for key in centroids if key[0] == "ground_control")
    flt_keys = sorted(key for key in centroids if key[0] == "flight")
    gc_matrix = np.vstack([centroids[key] for key in gc_keys])
    flt_matrix = np.vstack([centroids[key] for key in flt_keys])
    cosine_distance = cdist(flt_matrix, gc_matrix, metric="cosine")
    euclidean_distance = cdist(flt_matrix, gc_matrix, metric="euclidean")

    pair_rows = []
    shift_frames = []
    for flt_i, flt_key in enumerate(flt_keys):
        best_j = int(np.argmin(cosine_distance[flt_i]))
        gc_key = gc_keys[best_j]
        delta = centroids[flt_key] - centroids[gc_key]
        pair_rows.append(
            {
                "flight_cluster": flt_key[1],
                "matched_ground_cluster": gc_key[1],
                "cosine_similarity": float(1.0 - cosine_distance[flt_i, best_j]),
                "euclidean_distance": float(euclidean_distance[flt_i, best_j]),
                "mean_abs_pathway_delta_z": float(np.mean(np.abs(delta))),
                "max_abs_pathway_delta_z": float(np.max(np.abs(delta))),
            }
        )
        shifts = summarize_pathway_shifts(delta, pathway_cols, top_n=top_n)
        shifts.insert(0, "matched_ground_cluster", gc_key[1])
        shifts.insert(0, "flight_cluster", flt_key[1])
        shift_frames.append(shifts)

    return pd.DataFrame(pair_rows), pd.concat(shift_frames, ignore_index=True)


def pathway_cluster_tests(
    z_df: pd.DataFrame,
    assignments: pd.DataFrame,
    pathway_cols: list[str],
) -> pd.DataFrame:
    """One-vs-rest pathway contrasts within each condition."""
    rows = []
    for condition in ["ground_control", "flight"]:
        condition_index = assignments.index[assignments["condition"] == condition]
        condition_scores = z_df.loc[condition_index, pathway_cols]
        condition_assignments = assignments.loc[condition_index]
        for cluster in sorted(condition_assignments["cluster"].unique()):
            in_cluster = condition_assignments["cluster"].eq(cluster)
            outside = ~in_cluster
            if in_cluster.sum() < 2 or outside.sum() < 2:
                continue
            cluster_mean = condition_scores.loc[in_cluster].mean(axis=0)
            outside_mean = condition_scores.loc[outside].mean(axis=0)
            cluster_var = condition_scores.loc[in_cluster].var(axis=0, ddof=1)
            outside_var = condition_scores.loc[outside].var(axis=0, ddof=1)
            se = np.sqrt(
                cluster_var / int(in_cluster.sum())
                + outside_var / int(outside.sum())
            )
            delta = cluster_mean - outside_mean
            z_stat = delta / se.replace(0, np.nan)
            # Normal approximation is enough here; these are descriptive
            # cluster signatures, not final biological inference.
            p_values = 2.0 * (
                1.0
                - pd.Series(np.abs(z_stat), index=pathway_cols).map(
                    lambda value: 0.5 * (1.0 + math.erf(float(value) / math.sqrt(2.0)))
                    if np.isfinite(value)
                    else math.nan
                )
            )
            frame = pd.DataFrame(
                {
                    "condition": condition,
                    "cluster": int(cluster),
                    "pathway": pathway_cols,
                    "pathway_label": [clean_term(col) for col in pathway_cols],
                    "cluster_mean_z": cluster_mean.to_numpy(),
                    "rest_mean_z": outside_mean.to_numpy(),
                    "delta_cluster_minus_rest_z": delta.to_numpy(),
                    "z_stat": z_stat.to_numpy(),
                    "p_value": p_values.to_numpy(),
                }
            )
            frame["fdr_bh"] = bh_fdr(frame["p_value"])
            rows.append(frame)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    return out.sort_values(
        ["condition", "cluster", "fdr_bh", "delta_cluster_minus_rest_z"],
        ascending=[True, True, True, False],
    )


def write_readme(
    output_dir: Path,
    *,
    score_tsv: Path,
    summary: dict[str, Any],
    pairs: pd.DataFrame,
    cluster_summary: pd.DataFrame,
) -> None:
    lines = [
        "# expiMap Condition-Specific Clustering",
        "",
        f"Input scores: `{score_tsv}`",
        "",
        "This analysis clusters posterior mean expiMap pathway scores separately",
        "within ground control and flight samples, then compares flight cluster",
        "centroids to their nearest ground-control centroids by pathway signature.",
        "The cluster labels are condition-specific; use matched centroids and",
        "pathway shifts rather than assuming equal numeric cluster IDs correspond.",
        "",
        "## Selected cluster counts",
        "",
        "| condition | selected_k | silhouette |",
        "|---|---:|---:|",
    ]
    for condition, info in summary["clustering"].items():
        silhouette = info["silhouette"]
        silhouette_text = "" if silhouette is None else f"{silhouette:.4f}"
        lines.append(
            f"| {condition} | {info['selected_k']} | {silhouette_text} |"
        )
    lines.extend(["", "## Cluster sizes", "", "| condition | cluster | samples | top accessions |", "|---|---:|---:|---|"])
    for row in cluster_summary.sort_values(["condition", "cluster"]).itertuples(index=False):
        lines.append(
            f"| {row.condition} | {int(row.cluster)} | {int(row.n_samples)} | `{row.top_accessions}` |"
        )
    lines.extend(["", "## Nearest centroid matches", "", "| flight cluster | matched GC cluster | cosine similarity | mean abs delta z |", "|---:|---:|---:|---:|"])
    for row in pairs.sort_values("flight_cluster").itertuples(index=False):
        lines.append(
            f"| {int(row.flight_cluster)} | {int(row.matched_ground_cluster)} | "
            f"{float(row.cosine_similarity):.4f} | {float(row.mean_abs_pathway_delta_z):.4f} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `condition_cluster_assignments.tsv`: sample-level condition-specific labels.",
            "- `condition_cluster_summary.tsv`: cluster sizes and accession composition.",
            "- `flight_to_ground_cluster_matches.tsv`: nearest GC centroid for each FLT cluster.",
            "- `flight_vs_matched_ground_pathway_shifts.tsv`: strongest pathway deltas for matched pairs.",
            "- `cluster_pathway_signatures.tsv`: one-vs-rest pathway signatures within each condition.",
            "- `analysis_summary.json`: run parameters and selected-k diagnostics.",
            "",
            "Interpretation note: this is exploratory structure discovery. Any pathway",
            "called from these clusters still needs accession-aware validation, because",
            "OSDR accession/batch structure is strong enough to create apparent clusters.",
        ]
    )
    output_dir.joinpath("README.md").write_text("\n".join(lines) + "\n")


def run(args: argparse.Namespace) -> None:
    score_tsv = Path(args.score_tsv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"Loading {score_tsv}")
    df = pd.read_csv(score_tsv, sep="\t")
    pathway_cols = [col for col in df.columns if col.startswith(PATHWAY_PREFIX)]
    if not pathway_cols:
        raise ValueError(f"No {PATHWAY_PREFIX} pathway columns found in {score_tsv}")
    if args.condition_col not in df.columns:
        raise ValueError(f"Missing condition column: {args.condition_col}")
    if args.accession_col not in df.columns:
        raise ValueError(f"Missing accession column: {args.accession_col}")
    if args.sample_id_col not in df.columns:
        raise ValueError(f"Missing sample id column: {args.sample_id_col}")

    df = df.reset_index(drop=True)
    conditions = df[args.condition_col].map(canonical_condition)
    keep = conditions.isin(["flight", "ground_control"])
    if not keep.all():
        log(f"Dropping {int((~keep).sum())} samples without FLT/GC condition")
    df = df.loc[keep].reset_index(drop=True)
    conditions = conditions.loc[keep].reset_index(drop=True)

    scores = df[pathway_cols].apply(pd.to_numeric, errors="coerce")
    scores = scores.fillna(scores.median(axis=0))
    scaler = StandardScaler()
    z_scores = scaler.fit_transform(scores.to_numpy(dtype=float))
    z_df = pd.DataFrame(z_scores, columns=pathway_cols)

    n_components = min(args.n_pcs, z_scores.shape[0] - 1, z_scores.shape[1])
    pca = PCA(n_components=n_components, random_state=args.seed)
    pcs = pca.fit_transform(z_scores)

    assignments = pd.DataFrame(
        {
            "sample_id": df[args.sample_id_col].astype(str),
            args.accession_col: df[args.accession_col].astype(str),
            "condition": conditions,
        }
    )
    clustering_summary: dict[str, Any] = {}
    assignments["cluster"] = -1
    for condition in ["ground_control", "flight"]:
        idx = assignments.index[assignments["condition"] == condition].to_numpy()
        labels, info = choose_kmeans(
            pcs[idx],
            min_k=args.min_k,
            max_k=args.max_k,
            seed=args.seed,
        )
        assignments.loc[idx, "cluster"] = labels
        clustering_summary[condition] = info
        log(
            f"{condition}: selected k={info['selected_k']} "
            f"silhouette={info['silhouette']}"
        )

    for component in range(min(2, pcs.shape[1])):
        assignments[f"pc{component + 1}"] = pcs[:, component]

    cluster_summary = cluster_summaries(
        assignments,
        condition="ground_control",
        accession_col=args.accession_col,
    )
    cluster_summary = pd.concat(
        [
            cluster_summary,
            cluster_summaries(
                assignments,
                condition="flight",
                accession_col=args.accession_col,
            ),
        ],
        ignore_index=True,
    )

    pairs, matched_shifts = compare_clusters(
        z_scores,
        assignments,
        pathway_cols,
        top_n=args.top_pathways,
    )
    signatures = pathway_cluster_tests(z_df, assignments, pathway_cols)

    assignments.to_csv(output_dir / "condition_cluster_assignments.tsv", sep="\t", index=False)
    cluster_summary.to_csv(output_dir / "condition_cluster_summary.tsv", sep="\t", index=False)
    pairs.to_csv(output_dir / "flight_to_ground_cluster_matches.tsv", sep="\t", index=False)
    matched_shifts.to_csv(
        output_dir / "flight_vs_matched_ground_pathway_shifts.tsv",
        sep="\t",
        index=False,
    )
    signatures.to_csv(output_dir / "cluster_pathway_signatures.tsv", sep="\t", index=False)

    summary = {
        "score_tsv": str(score_tsv),
        "n_samples": int(len(df)),
        "n_pathways": int(len(pathway_cols)),
        "condition_counts": conditions.value_counts().to_dict(),
        "n_pcs": int(n_components),
        "pca_explained_variance_ratio_sum": float(pca.explained_variance_ratio_.sum()),
        "min_k": int(args.min_k),
        "max_k": int(args.max_k),
        "seed": int(args.seed),
        "clustering": clustering_summary,
        "outputs": {
            "assignments": str(output_dir / "condition_cluster_assignments.tsv"),
            "cluster_summary": str(output_dir / "condition_cluster_summary.tsv"),
            "matches": str(output_dir / "flight_to_ground_cluster_matches.tsv"),
            "matched_pathway_shifts": str(
                output_dir / "flight_vs_matched_ground_pathway_shifts.tsv"
            ),
            "cluster_pathway_signatures": str(output_dir / "cluster_pathway_signatures.tsv"),
        },
    }
    (output_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_readme(
        output_dir,
        score_tsv=score_tsv,
        summary=summary,
        pairs=pairs,
        cluster_summary=cluster_summary,
    )
    log(f"Wrote {output_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster GC and FLT expiMap pathway scores separately and compare "
            "condition-specific pathway centroids."
        )
    )
    parser.add_argument("--score-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--condition-col", default="condition_inferred")
    parser.add_argument("--accession-col", default="id.accession")
    parser.add_argument("--sample-id-col", default="obs_name")
    parser.add_argument("--min-k", type=int, default=2)
    parser.add_argument("--max-k", type=int, default=6)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--top-pathways", type=int, default=25)
    parser.add_argument("--seed", type=int, default=2024)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
