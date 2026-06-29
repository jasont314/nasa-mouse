"""Analyze expiMap pathway scores for FLT-vs-GC shifts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cluster_enrichment import bh_fdr
from .io import require_import


def pathway_columns(frame, include_de_novo: bool = False) -> list[str]:
    metadata_prefixes = {
        "obs_name",
        "profile_id",
        "profile",
        "id.accession",
        "id.assay name",
        "id.sample name",
        "condition_inferred",
        "tissue_final",
    }
    annotated_prefixes = ("R-MMU-", "MUSCLE_")
    prefixes = (
        (*annotated_prefixes, "unconstrained_")
        if include_de_novo
        else annotated_prefixes
    )
    return [
        column
        for column in frame.columns
        if str(column).startswith(prefixes) and column not in metadata_prefixes
    ]


def compare_pathways(scores, terms: list[str]):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for term in terms:
        flight = scores.loc[scores["condition_inferred"].eq("flight"), term].astype(float)
        ground = scores.loc[
            scores["condition_inferred"].eq("ground_control"),
            term,
        ].astype(float)
        if len(flight) == 0 or len(ground) == 0:
            continue
        t_p = scipy_stats.ttest_ind(
            flight,
            ground,
            equal_var=False,
            nan_policy="omit",
        ).pvalue
        try:
            mw_p = scipy_stats.mannwhitneyu(
                flight,
                ground,
                alternative="two-sided",
            ).pvalue
        except ValueError:
            mw_p = float("nan")
        rows.append(
            {
                "term": term,
                "n_flight": int(len(flight)),
                "n_ground_control": int(len(ground)),
                "mean_flight": float(np.nanmean(flight)),
                "mean_ground_control": float(np.nanmean(ground)),
                "flight_minus_ground": float(np.nanmean(flight) - np.nanmean(ground)),
                "welch_p": float(t_p),
                "mannwhitney_p": float(mw_p),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["welch_fdr"] = bh_fdr(result["welch_p"].fillna(1.0).to_numpy())
        result["mannwhitney_fdr"] = bh_fdr(
            result["mannwhitney_p"].fillna(1.0).to_numpy()
        )
        result = result.sort_values(
            ["welch_fdr", "mannwhitney_fdr", "flight_minus_ground"],
            ascending=[True, True, False],
        )
    return result


def accession_effects(scores, terms: list[str]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_stats = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    for accession, group in scores.groupby("id.accession", dropna=False):
        if not {
            "flight",
            "ground_control",
        }.issubset(set(group["condition_inferred"].astype(str))):
            continue
        for term in terms:
            flight = group.loc[group["condition_inferred"].eq("flight"), term].astype(float)
            ground = group.loc[
                group["condition_inferred"].eq("ground_control"),
                term,
            ].astype(float)
            rows.append(
                {
                    "id.accession": accession,
                    "term": term,
                    "n_flight": int(len(flight)),
                    "n_ground_control": int(len(ground)),
                    "mean_flight": float(flight.mean()),
                    "mean_ground_control": float(ground.mean()),
                    "flight_minus_ground": float(flight.mean() - ground.mean()),
                }
            )
    effects = pd.DataFrame(rows)
    tests = []
    if not effects.empty:
        for term, group in effects.groupby("term"):
            values = group["flight_minus_ground"].astype(float)
            if len(values) < 2:
                p_value = float("nan")
            else:
                try:
                    p_value = scipy_stats.wilcoxon(values).pvalue
                except ValueError:
                    p_value = 1.0
            tests.append(
                {
                    "term": term,
                    "n_accessions": int(len(values)),
                    "mean_accession_effect": float(values.mean()),
                    "median_accession_effect": float(values.median()),
                    "wilcoxon_p": float(p_value),
                }
            )
    test_frame = pd.DataFrame(tests)
    if not test_frame.empty:
        test_frame["wilcoxon_fdr"] = bh_fdr(
            test_frame["wilcoxon_p"].fillna(1.0).to_numpy()
        )
        test_frame = test_frame.sort_values(
            ["wilcoxon_fdr", "mean_accession_effect"],
            ascending=[True, False],
        )
    return effects, test_frame


def metadata_frame(scores):
    return scores[
        [
            column
            for column in [
                "obs_name",
                "profile_id",
                "id.accession",
                "condition_inferred",
                "tissue_final",
                "study.characteristics.sex",
                "study.characteristics.strain",
                "study.characteristics.genotype",
            ]
            if column in scores
        ]
    ].copy()


def write_pca(scores, terms: list[str], output_dir: Path) -> dict[str, str]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn_decomp = require_import("sklearn.decomposition", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    x = scores[terms].astype(float).to_numpy()
    x = np.nan_to_num(x, copy=False)
    n_components = min(2, x.shape[0], x.shape[1])
    pca = sklearn_decomp.PCA(n_components=n_components)
    coords = pca.fit_transform(x)
    frame = metadata_frame(scores)
    frame["PC1"] = coords[:, 0]
    frame["PC2"] = coords[:, 1] if n_components > 1 else 0.0
    path = output_dir / "pathway_score_pca.tsv"
    frame.to_csv(path, sep="\t", index=False)

    colors = scores["condition_inferred"].map(
        {"flight": "#c43c39", "ground_control": "#2878b5"}
    ).fillna("#777777")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(frame["PC1"], frame["PC2"], c=colors, s=18, alpha=0.8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
        if n_components > 1
        else "PC2"
    )
    ax.set_title("expiMap pathway scores")
    fig.tight_layout()
    png = output_dir / "pathway_score_pca.png"
    fig.savefig(png, dpi=180)
    plt.close(fig)

    accession_png = output_dir / "pathway_score_pca_by_accession.png"
    if "id.accession" in scores:
        codes, labels = pd.factorize(scores["id.accession"].astype(str))
        fig, ax = plt.subplots(figsize=(7, 5))
        scatter = ax.scatter(
            frame["PC1"],
            frame["PC2"],
            c=codes,
            cmap="tab20",
            s=18,
            alpha=0.85,
        )
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
        ax.set_ylabel(
            f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
            if n_components > 1
            else "PC2"
        )
        ax.set_title("expiMap pathway scores by OSD accession")
        if len(labels) <= 15:
            handles = [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    color=scatter.cmap(scatter.norm(code)),
                    label=label,
                )
                for code, label in enumerate(labels)
            ]
            ax.legend(handles=handles, fontsize=7, frameon=False, loc="best")
        fig.tight_layout()
        fig.savefig(accession_png, dpi=180)
        plt.close(fig)
    return {
        "pca_coordinates": str(path),
        "pca_by_condition": str(png),
        "pca_by_accession": str(accession_png),
    }


def write_umap(scores, terms: list[str], output_dir: Path) -> dict[str, str]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    umap_module = require_import("umap", "pip install umap-learn")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    x = scores[terms].astype(float).to_numpy()
    x = np.nan_to_num(x, copy=False)
    n_neighbors = max(2, min(15, x.shape[0] - 1))
    reducer = umap_module.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.3,
        metric="euclidean",
        random_state=2020,
    )
    coords = reducer.fit_transform(x)
    frame = metadata_frame(scores)
    frame["UMAP1"] = coords[:, 0]
    frame["UMAP2"] = coords[:, 1]
    path = output_dir / "pathway_score_umap.tsv"
    frame.to_csv(path, sep="\t", index=False)

    colors = scores["condition_inferred"].map(
        {"flight": "#c43c39", "ground_control": "#2878b5"}
    ).fillna("#777777")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(frame["UMAP1"], frame["UMAP2"], c=colors, s=18, alpha=0.8)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("expiMap pathway scores")
    fig.tight_layout()
    condition_png = output_dir / "pathway_score_umap.png"
    fig.savefig(condition_png, dpi=180)
    plt.close(fig)

    accession_png = output_dir / "pathway_score_umap_by_accession.png"
    if "id.accession" in scores:
        codes, labels = pd.factorize(scores["id.accession"].astype(str))
        fig, ax = plt.subplots(figsize=(7, 5))
        scatter = ax.scatter(
            frame["UMAP1"],
            frame["UMAP2"],
            c=codes,
            cmap="tab20",
            s=18,
            alpha=0.85,
        )
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("expiMap pathway scores by OSD accession")
        if len(labels) <= 15:
            handles = [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    color=scatter.cmap(scatter.norm(code)),
                    label=label,
                )
                for code, label in enumerate(labels)
            ]
            ax.legend(handles=handles, fontsize=7, frameon=False, loc="best")
        fig.tight_layout()
        fig.savefig(accession_png, dpi=180)
        plt.close(fig)

    return {
        "umap_coordinates": str(path),
        "umap_by_condition": str(condition_png),
        "umap_by_accession": str(accession_png),
    }


def write_top_heatmap(scores, comparison, terms: list[str], output_dir: Path, top_n: int = 30) -> dict[str, str]:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    plt = require_import("matplotlib.pyplot", "pip install -r requirements-nasa-mouse-glare.txt")

    if comparison.empty:
        return {}
    top_terms = (
        comparison.assign(abs_effect=comparison["flight_minus_ground"].abs())
        .sort_values(["abs_effect", "welch_fdr"], ascending=[False, True])
        .head(top_n)["term"]
        .tolist()
    )
    top_terms = [term for term in top_terms if term in terms]
    if not top_terms:
        return {}
    ordered = scores.copy()
    sort_cols = [
        column
        for column in ["condition_inferred", "id.accession", "obs_name"]
        if column in ordered
    ]
    ordered = ordered.sort_values(sort_cols, kind="stable")
    x = ordered[top_terms].astype(float).to_numpy()
    mean = np.nanmean(x, axis=0, keepdims=True)
    std = np.nanstd(x, axis=0, keepdims=True)
    std[std == 0] = 1.0
    z = (x - mean) / std

    fig_height = max(5, min(14, 0.22 * len(top_terms)))
    fig, ax = plt.subplots(figsize=(9, fig_height))
    im = ax.imshow(z.T, aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels(top_terms, fontsize=6)
    ax.set_xticks([])
    ax.set_xlabel("samples sorted by condition/accession")
    ax.set_title("Top expiMap FLT-vs-GC pathway shifts")
    fig.colorbar(im, ax=ax, label="pathway score z-score", fraction=0.025, pad=0.02)
    fig.tight_layout()
    heatmap_path = output_dir / "top_pathway_shift_heatmap.png"
    fig.savefig(heatmap_path, dpi=180)
    plt.close(fig)

    table_path = output_dir / "top_pathway_shift_heatmap_terms.tsv"
    comparison.loc[comparison["term"].isin(top_terms)].to_csv(
        table_path,
        sep="\t",
        index=False,
    )
    return {
        "top_pathway_heatmap": str(heatmap_path),
        "top_pathway_heatmap_terms": str(table_path),
    }


def write_summary(output_dir: Path, comparison, study_tests, scores_path: str) -> Path:
    top = comparison.head(10) if not comparison.empty else comparison
    lines = [
        "# expiMap Pathway Analysis",
        "",
        f"- Pathway scores: `{scores_path}`",
        f"- Tested pathways: {len(comparison)}",
        f"- Study-aware tested pathways: {len(study_tests)}",
        "",
        "## Top Aggregate FLT-vs-GC Pathways",
        "",
    ]
    if top.empty:
        lines.append("No pathway comparisons were available.")
    else:
        columns = [
            "term",
            "flight_minus_ground",
            "welch_p",
            "welch_fdr",
            "mannwhitney_fdr",
        ]
        lines.append("| " + " | ".join(columns) + " |")
        lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
        for _, row in top[columns].iterrows():
            values = []
            for column in columns:
                value = row[column]
                if isinstance(value, float):
                    values.append(f"{value:.6g}")
                else:
                    values.append(str(value))
            lines.append("| " + " | ".join(values) + " |")
    path = output_dir / "README.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(args.scores, sep="\t", keep_default_na=False)
    terms = pathway_columns(scores, include_de_novo=args.include_de_novo)
    if not terms:
        raise SystemExit("No requested expiMap score columns found.")

    comparison = compare_pathways(scores, terms)
    comparison_path = output_dir / "flt_vs_gc_pathway_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)

    effects, study_tests = accession_effects(scores, terms)
    effects_path = output_dir / "flight_ground_effects_by_accession.tsv"
    tests_path = output_dir / "flight_ground_study_aware_tests.tsv"
    effects.to_csv(effects_path, sep="\t", index=False)
    study_tests.to_csv(tests_path, sep="\t", index=False)
    pca_outputs = write_pca(scores, terms, output_dir)
    umap_outputs = write_umap(scores, terms, output_dir)
    heatmap_outputs = write_top_heatmap(scores, comparison, terms, output_dir)
    readme_path = write_summary(output_dir, comparison, study_tests, args.scores)

    summary = {
        "scores": str(args.scores),
        "n_samples": int(len(scores)),
        "n_terms": int(len(terms)),
        "includes_de_novo_programs": bool(args.include_de_novo),
        "outputs": {
            "comparison": str(comparison_path),
            "accession_effects": str(effects_path),
            "study_aware_tests": str(tests_path),
            **pca_outputs,
            **umap_outputs,
            **heatmap_outputs,
            "readme": str(readme_path),
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze expiMap pathway scores by FLT/GC condition."
    )
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--include-de-novo",
        action="store_true",
        help="Include unconstrained_ expiMap programs alongside Reactome pathways.",
    )
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
