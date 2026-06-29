"""Analyze WGAN critic feature scores for FLT vs GC shifts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.cluster_enrichment import bh_fdr
from nasa_mouse_glare.io import require_import
from nasa_mouse_glare.validate_expimap_accession_effects import random_effects


def feature_columns(frame) -> list[str]:
    return [
        column
        for column in frame.columns
        if str(column).startswith("WGAN_FEATURE_") or column == "CRITIC_SCORE"
    ]


def compare_features(scores, features):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")
    flight = scores.loc[scores["condition_inferred"].eq("flight")]
    ground = scores.loc[scores["condition_inferred"].eq("ground_control")]
    rows = []
    for feature in features:
        x = flight[feature].astype(float).to_numpy()
        y = ground[feature].astype(float).to_numpy()
        if len(x) < 2 or len(y) < 2:
            welch_p = 1.0
            mann_p = 1.0
        else:
            welch_p = float(scipy_stats.ttest_ind(x, y, equal_var=False, nan_policy="omit").pvalue)
            mann_p = float(scipy_stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
        rows.append(
            {
                "feature": feature,
                "n_flight": int(len(x)),
                "n_ground_control": int(len(y)),
                "flight_mean": float(np.nanmean(x)) if len(x) else float("nan"),
                "ground_control_mean": float(np.nanmean(y)) if len(y) else float("nan"),
                "flight_minus_ground": float(np.nanmean(x) - np.nanmean(y))
                if len(x) and len(y)
                else float("nan"),
                "welch_p": welch_p if np.isfinite(welch_p) else 1.0,
                "mannwhitney_p": mann_p if np.isfinite(mann_p) else 1.0,
            }
        )
    result = pd.DataFrame(rows)
    result["welch_fdr"] = bh_fdr(result["welch_p"].fillna(1.0).to_numpy())
    result["mannwhitney_fdr"] = bh_fdr(result["mannwhitney_p"].fillna(1.0).to_numpy())
    return result.sort_values(["welch_fdr", "welch_p"], kind="stable")


def accession_effects(scores, features):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    for accession, frame in scores.groupby("id.accession", dropna=False):
        flight = frame.loc[frame["condition_inferred"].eq("flight")]
        ground = frame.loc[frame["condition_inferred"].eq("ground_control")]
        if flight.empty or ground.empty:
            continue
        for feature in features:
            x = flight[feature].astype(float).to_numpy()
            y = ground[feature].astype(float).to_numpy()
            variance = float(np.var(x, ddof=1) / len(x) + np.var(y, ddof=1) / len(y))
            rows.append(
                {
                    "id.accession": accession,
                    "feature": feature,
                    "n_flight": int(len(x)),
                    "n_ground_control": int(len(y)),
                    "flight_minus_ground": float(np.mean(x) - np.mean(y)),
                    "effect_variance": variance,
                }
            )
    return pd.DataFrame(rows)


def meta_analysis(effects):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    for feature, frame in effects.groupby("feature", sort=False):
        result = random_effects(
            frame["flight_minus_ground"].to_numpy(),
            frame["effect_variance"].to_numpy(),
            scipy_stats,
        )
        rows.append({"feature": feature, **result})
    result = pd.DataFrame(rows)
    direction_rows = []
    for feature, frame in effects.groupby("feature", sort=False):
        direction = float(result.loc[result["feature"].eq(feature), "meta_effect"].iloc[0])
        individual = np.sign(frame["flight_minus_ground"].astype(float).to_numpy())
        direction_rows.append(
            {
                "feature": feature,
                "n_accession_same_direction": int((individual == np.sign(direction)).sum()),
                "n_accession_opposite_direction": int((individual == -np.sign(direction)).sum()),
            }
        )
    result = result.merge(pd.DataFrame(direction_rows), on="feature", how="left")
    result["meta_fdr"] = bh_fdr(result["meta_p"].fillna(1.0).to_numpy())
    return result.sort_values(["meta_fdr", "meta_p"], kind="stable")


def leave_one_accession_out(effects):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")
    rows = []
    accessions = sorted(effects["id.accession"].astype(str).unique())
    for held_out in accessions:
        retained = effects.loc[effects["id.accession"].astype(str).ne(held_out)]
        for feature, frame in retained.groupby("feature", sort=False):
            result = random_effects(
                frame["flight_minus_ground"].to_numpy(),
                frame["effect_variance"].to_numpy(),
                scipy_stats,
            )
            rows.append({"held_out_accession": held_out, "feature": feature, **result})
    result = pd.DataFrame(rows)
    if result.empty:
        return result, pd.DataFrame()
    result["meta_fdr"] = result.groupby("held_out_accession", sort=False)["meta_p"].transform(
        lambda values: bh_fdr(values.fillna(1.0).to_numpy())
    )
    primary = meta_analysis(effects).set_index("feature")["meta_effect"]
    summary_rows = []
    for feature, frame in result.groupby("feature", sort=False):
        direction = float(np.sign(primary.get(feature, 0.0)))
        loo_direction = np.sign(frame["meta_effect"].astype(float).to_numpy())
        summary_rows.append(
            {
                "feature": feature,
                "n_leave_one_out": int(len(frame)),
                "n_same_direction": int((loo_direction == direction).sum()),
                "minimum_leave_one_out_fdr": float(frame["meta_fdr"].min()),
                "maximum_leave_one_out_fdr": float(frame["meta_fdr"].max()),
            }
        )
    return result, pd.DataFrame(summary_rows)


def plot_pca(scores, features, output_dir: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_decomposition = require_import(
        "sklearn.decomposition", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    x = scores[features].astype(float).to_numpy()
    x = np.nan_to_num(x)
    n_components = min(2, x.shape[0], x.shape[1])
    if n_components < 2:
        return {}
    pca = sklearn_decomposition.PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(x)
    frame = pd.DataFrame(
        {
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
            "condition_inferred": scores["condition_inferred"].astype(str).to_numpy(),
            "id.accession": scores["id.accession"].astype(str).to_numpy(),
        }
    )
    coord_path = output_dir / "wgan_feature_pca.tsv"
    frame.to_csv(coord_path, sep="\t", index=False)
    paths = {"pca_coordinates": str(coord_path)}
    for color_col, name in [
        ("condition_inferred", "wgan_feature_pca.png"),
        ("id.accession", "wgan_feature_pca_by_accession.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 5))
        for label, group in frame.groupby(color_col, sort=True):
            ax.scatter(group["PC1"], group["PC2"], s=25, label=label, alpha=0.85)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("WGAN critic features")
        ax.legend(loc="best", fontsize=7)
        fig.tight_layout()
        path = output_dir / name
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths[name.removesuffix(".png")] = str(path)
    return paths


def plot_umap(scores, features, output_dir: Path):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")
    try:
        umap_module = require_import("umap", "pip install umap-learn")
    except Exception:
        return {}
    x = np.nan_to_num(scores[features].astype(float).to_numpy())
    if x.shape[0] < 5 or x.shape[1] < 2:
        return {}
    reducer = umap_module.UMAP(n_components=2, random_state=0, n_neighbors=min(15, x.shape[0] - 1))
    coords = reducer.fit_transform(x)
    frame = pd.DataFrame(
        {
            "UMAP1": coords[:, 0],
            "UMAP2": coords[:, 1],
            "condition_inferred": scores["condition_inferred"].astype(str).to_numpy(),
            "id.accession": scores["id.accession"].astype(str).to_numpy(),
        }
    )
    coord_path = output_dir / "wgan_feature_umap.tsv"
    frame.to_csv(coord_path, sep="\t", index=False)
    paths = {"umap_coordinates": str(coord_path)}
    for color_col, name in [
        ("condition_inferred", "wgan_feature_umap.png"),
        ("id.accession", "wgan_feature_umap_by_accession.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 5))
        for label, group in frame.groupby(color_col, sort=True):
            ax.scatter(group["UMAP1"], group["UMAP2"], s=25, label=label, alpha=0.85)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("WGAN critic features")
        ax.legend(loc="best", fontsize=7)
        fig.tight_layout()
        path = output_dir / name
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths[name.removesuffix(".png")] = str(path)
    return paths


def plot_heatmap(scores, comparison, features, output_dir: Path, *, top_n: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")
    selected = comparison.sort_values(["welch_fdr", "welch_p"], kind="stable")["feature"].head(top_n)
    selected = [feature for feature in selected if feature in features]
    if not selected:
        return {}
    order = scores.sort_values(["condition_inferred", "id.accession"]).index
    matrix = scores.loc[order, selected].astype(float).to_numpy()
    matrix = np.nan_to_num(matrix)
    fig, ax = plt.subplots(figsize=(max(6, len(selected) * 0.35), 6))
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="coolwarm")
    ax.set_xticks(range(len(selected)))
    ax.set_xticklabels(selected, rotation=90, fontsize=7)
    ax.set_yticks([])
    ax.set_title("Top WGAN FLT vs GC feature shifts")
    fig.colorbar(im, ax=ax, shrink=0.75)
    fig.tight_layout()
    path = output_dir / "top_wgan_feature_shift_heatmap.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    terms_path = output_dir / "top_wgan_feature_shift_heatmap_features.tsv"
    comparison.loc[comparison["feature"].isin(selected)].to_csv(terms_path, sep="\t", index=False)
    return {"top_feature_heatmap": str(path), "top_feature_heatmap_features": str(terms_path)}


def write_readme(output_dir: Path, summary: dict, meta, loo_summary) -> None:
    merged = meta.merge(loo_summary, on="feature", how="left") if not loo_summary.empty else meta
    lines = [
        "# WGAN Feature Analysis",
        "",
        f"- Scores: `{summary['scores']}`",
        f"- Samples: {summary['n_samples']}",
        f"- Features: {summary['n_features']}",
        f"- Ordinary Welch FDR < 0.05: {summary['ordinary_welch_fdr_lt_005']}",
        f"- Random-effects FDR < 0.05: {summary['random_effects_fdr_lt_005']}",
        f"- Strict LOO-stable FDR < 0.05: {summary['loo_stable_fdr_lt_005']}",
        "",
        "| feature | meta effect | meta FDR | same-direction accessions | max LOO FDR |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for _, row in merged.head(10).iterrows():
        max_loo = row.get("maximum_leave_one_out_fdr", float("nan"))
        same = row.get("n_accession_same_direction", 0)
        n_acc = row.get("n_accessions", 0)
        lines.append(
            f"| {row.feature} | {row.meta_effect:.4g} | {row.meta_fdr:.4g} | "
            f"{int(same)}/{int(n_acc)} | {max_loo:.4g} |"
        )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.scores, sep="\t")
    required = {"id.accession", "condition_inferred"}
    missing = required.difference(scores.columns)
    if missing:
        raise SystemExit(f"Scores missing required columns: {sorted(missing)}")
    features = feature_columns(scores)
    if not features:
        raise SystemExit("No WGAN feature columns found.")
    comparison = compare_features(scores, features)
    effects = accession_effects(scores, features)
    meta = meta_analysis(effects) if not effects.empty else pd.DataFrame()
    loo, loo_summary = leave_one_accession_out(effects) if not effects.empty else (pd.DataFrame(), pd.DataFrame())

    paths = {
        "comparison": output_dir / "flt_vs_gc_wgan_feature_comparison.tsv",
        "per_accession_effects": output_dir / "per_accession_wgan_feature_effects.tsv",
        "random_effects_meta_analysis": output_dir / "random_effects_meta_analysis.tsv",
        "leave_one_accession_out": output_dir / "leave_one_accession_out.tsv",
        "leave_one_out_summary": output_dir / "leave_one_out_summary.tsv",
    }
    comparison.to_csv(paths["comparison"], sep="\t", index=False)
    effects.to_csv(paths["per_accession_effects"], sep="\t", index=False)
    meta.to_csv(paths["random_effects_meta_analysis"], sep="\t", index=False)
    loo.to_csv(paths["leave_one_accession_out"], sep="\t", index=False)
    loo_summary.to_csv(paths["leave_one_out_summary"], sep="\t", index=False)
    plot_paths = {}
    plot_paths.update(plot_pca(scores, features, output_dir))
    plot_paths.update(plot_umap(scores, features, output_dir))
    plot_paths.update(plot_heatmap(scores, comparison, features, output_dir, top_n=args.top_features))

    if not meta.empty and not loo_summary.empty:
        stable = meta.merge(loo_summary, on="feature", how="inner")
        stable = stable.loc[
            (stable["meta_fdr"] < 0.05)
            & (stable["maximum_leave_one_out_fdr"] < 0.05)
            & (stable["n_same_direction"] == stable["n_leave_one_out"])
        ]
    else:
        stable = pd.DataFrame()

    summary = {
        "scores": str(args.scores),
        "n_samples": int(len(scores)),
        "n_features": int(len(features)),
        "ordinary_welch_fdr_lt_005": int((comparison["welch_fdr"] < 0.05).sum()),
        "ordinary_mannwhitney_fdr_lt_005": int((comparison["mannwhitney_fdr"] < 0.05).sum()),
        "random_effects_fdr_lt_005": int((meta["meta_fdr"] < 0.05).sum()) if not meta.empty else 0,
        "loo_stable_fdr_lt_005": int(len(stable)),
        "outputs": {key: str(path) for key, path in paths.items()} | plot_paths,
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_readme(output_dir, summary, meta, loo_summary)
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-features", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
