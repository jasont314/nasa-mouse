"""Triangulate the corrected kidney and spleen expiMap reassessment models."""

from __future__ import annotations

import json
from pathlib import Path
import re

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

from .build_asgsr_paper import ROOT
from .reviewer_robustness_analysis import (
    TMS_PATH,
    TissueData,
    accession_and_project_effects,
    atlas_marker_signatures,
    composition_adjusted_effects,
    composition_marker_scores,
    project_heldout_folds,
    run_preranked_gsea,
    run_ssgsea,
)
from .run_kidney_spleen_seed_sensitivity import (
    CONFIGS,
    OUTPUT_DIR,
    consensus_table,
    latent_directions,
)


MINIMUM_HELDOUT_CONCORDANCE = 2 / 3
MINIMUM_COMPOSITION_EFFECT_RATIO = 0.25
TOP_EFFECT_PERCENTILE = 0.9
COLORS = {"kidney": "#287b8e", "spleen": "#9b4f65"}
CONFOUNDED_ACCESSIONS = {"kidney": (), "spleen": ("OSD-288",)}


def _dense(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _log2cpm(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    return np.log2(counts / totals * 1_000_000.0 + 1.0)


def display_term(term: str) -> str:
    label = re.sub(r"^R-MMU-\d+_", "", str(term)).replace("_", " ").strip()
    return label.title()


def load_data(config) -> TissueData:
    query = ad.read_h5ad(config.query_input)
    counts = query.layers["counts"] if "counts" in query.layers else query.X
    counts = _dense(counts)
    genes = query.var_names.astype(str).to_numpy()
    terms = list(map(str, query.uns["terms"]))
    mask = np.asarray(query.varm["I"])
    gene_sets = {
        term: genes[np.flatnonzero(mask[:, index] != 0)].tolist()
        for index, term in enumerate(terms)
    }

    query_dir = config.query_dir(2020)
    score_frame = pd.read_csv(query_dir / "query_pathway_scores.tsv", sep="\t")
    direction = latent_directions(query_dir)
    active = [term for term in terms if term in score_frame and direction.get(term, 0) != 0]
    scores = score_frame[active].astype(float).mul(direction.loc[active], axis=1)
    scores.index = score_frame["obs_name"].astype(str)
    obs = query.obs.copy()
    obs.index = obs.index.astype(str)
    scores = scores.reindex(obs.index)
    excluded = CONFOUNDED_ACCESSIONS[config.tissue]
    retained = ~obs["id.accession"].astype(str).isin(excluded)
    obs = obs.loc[retained].copy()
    scores = scores.loc[retained].copy()
    counts = counts[retained.to_numpy()]
    if scores.isna().any().any():
        raise RuntimeError(f"Could not align expiMap scores for {config.tissue}")
    return TissueData(
        tissue=config.tissue,
        obs=obs,
        genes=genes,
        log2cpm=_log2cpm(counts),
        gene_sets=gene_sets,
        expimap_scores=scores,
    )


def method_benchmark(data: TissueData) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expi_accession, expi_project, expi_summary = accession_and_project_effects(
        data.expimap_scores, data.obs, data.tissue, "expimap"
    )
    print(f"ssGSEA {data.tissue}", flush=True)
    ss_scores = run_ssgsea(data)
    ss_accession, _, ss_summary = accession_and_project_effects(
        ss_scores, data.obs, data.tissue, "ssgsea"
    )
    print(f"preranked GSEA {data.tissue}", flush=True)
    gsea = run_preranked_gsea(data)

    expi = expi_summary.rename(
        columns={
            column: f"expimap_{column}"
            for column in expi_summary.columns
            if column not in {"tissue", "method", "term"}
        }
    ).drop(columns="method")
    ss = ss_summary.rename(
        columns={
            column: f"ssgsea_{column}"
            for column in ss_summary.columns
            if column not in {"tissue", "method", "term"}
        }
    ).drop(columns="method")
    benchmark = expi.merge(ss, on=["tissue", "term"], how="left").merge(
        gsea, on=["tissue", "term"], how="left"
    )
    benchmark["expimap_ssgsea_direction_match"] = (
        np.sign(benchmark["expimap_accession_balanced_effect"])
        == np.sign(benchmark["ssgsea_accession_balanced_effect"])
    )
    benchmark["expimap_gsea_direction_match"] = (
        np.sign(benchmark["expimap_accession_balanced_effect"])
        == np.sign(benchmark["gsea_nes"])
    )

    empty = pd.DataFrame(columns=["tissue", "term"])
    heldout = project_heldout_folds(expi_project, empty, empty)
    heldout_summary = (
        heldout.groupby(["tissue", "term"])
        .agg(
            heldout_project_direction_concordance=("direction_concordant", "mean"),
            heldout_projects=("heldout_project", "size"),
            heldout_top_decile_selection_fraction=(
                "selected_in_training_top_decile",
                "mean",
            ),
        )
        .reset_index()
    )
    selected = heldout.loc[heldout["selected_in_training_top_decile"]]
    selected_summary = (
        selected.groupby(["tissue", "term"])["direction_concordant"]
        .agg(
            heldout_selected_direction_concordance="mean",
            heldout_selected_folds="size",
        )
        .reset_index()
    )
    heldout_summary = heldout_summary.merge(
        selected_summary, on=["tissue", "term"], how="left"
    )
    accession = pd.concat([expi_accession, ss_accession], ignore_index=True)
    return benchmark, heldout, heldout_summary


def composition_sensitivity(
    data: TissueData, atlas_obs: pd.DataFrame, atlas_genes: pd.Index
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print(f"composition proxies {data.tissue}", flush=True)
    _, markers = atlas_marker_signatures(data, atlas_obs, atlas_genes)
    marker_scores = composition_marker_scores(data, markers)
    effects, sample_scores = composition_adjusted_effects(data, marker_scores)
    return effects, markers, sample_scores


def old_model_comparison(config, new_effects: pd.DataFrame) -> dict[str, object]:
    old_dir = (
        config.base_dir
        / "tutorial_hvg_2000"
        / "query_nb_250epoch_seed2020"
    )
    if not (old_dir / "query_pathway_scores.tsv").exists():
        return {"tissue": config.tissue, "old_model_available": False}
    scores = pd.read_csv(old_dir / "query_pathway_scores.tsv", sep="\t")
    directions = latent_directions(old_dir)
    terms = [term for term in directions.index if term in scores and directions[term] != 0]
    oriented = scores[terms].astype(float).mul(directions.loc[terms], axis=1)
    condition = scores["condition_inferred"].astype(str)
    excluded = set(CONFOUNDED_ACCESSIONS[config.tissue])
    rows = []
    for accession, indexes in scores.groupby("id.accession", observed=True).indices.items():
        if str(accession) in excluded:
            continue
        indexes = np.asarray(indexes)
        local = condition.iloc[indexes]
        flight = local.eq("flight").to_numpy()
        ground = local.eq("ground_control").to_numpy()
        if flight.any() and ground.any():
            effect = (
                oriented.iloc[indexes].iloc[flight].mean(axis=0)
                - oriented.iloc[indexes].iloc[ground].mean(axis=0)
            )
            rows.append(effect.rename(str(accession)))
    old = pd.DataFrame(rows).mean(axis=0).rename("old_effect")
    new = new_effects.set_index("term")["expimap_accession_balanced_effect"].rename(
        "new_effect"
    )
    common = pd.concat([old, new], axis=1, join="inner").dropna()
    old_top = set(common.index[common["old_effect"].abs().rank(pct=True).ge(0.9)])
    new_top = set(common.index[common["new_effect"].abs().rank(pct=True).ge(0.9)])
    return {
        "tissue": config.tissue,
        "old_model_available": True,
        "common_active_pathways": int(len(common)),
        "old_new_spearman_rho": float(
            stats.spearmanr(common["old_effect"], common["new_effect"]).statistic
        ),
        "old_new_direction_agreement": float(
            (np.sign(common["old_effect"]) == np.sign(common["new_effect"])).mean()
        ),
        "old_new_top_decile_jaccard": float(
            len(old_top & new_top) / len(old_top | new_top)
        ),
    }


def evidence_matrix(
    benchmark: pd.DataFrame,
    heldout: pd.DataFrame,
    composition: pd.DataFrame,
    seeds: pd.DataFrame,
) -> pd.DataFrame:
    seed_columns = [
        "tissue",
        "term",
        "effect_seed2020",
        "effect_seed2021",
        "effect_seed2022",
        "seed_effect_median",
        "seed_effect_minimum",
        "seed_effect_maximum",
        "primary_absolute_percentile",
        "all_three_seeds_available",
        "all_three_direction_concordant",
    ]
    composition_columns = [
        "tissue",
        "term",
        "within_accession_unadjusted_effect",
        "composition_proxy_adjusted_effect",
        "adjusted_direction_matches_unadjusted",
        "absolute_effect_ratio_adjusted_to_unadjusted",
    ]
    matrix = (
        benchmark.merge(heldout, on=["tissue", "term"], how="left")
        .merge(seeds[seed_columns], on=["tissue", "term"], how="left")
        .merge(composition[composition_columns], on=["tissue", "term"], how="left")
    )
    matrix["display_term"] = matrix["term"].map(display_term)
    matrix["primary_top_decile"] = matrix["primary_absolute_percentile"].ge(
        TOP_EFFECT_PERCENTILE
    )
    matrix["ssgsea_direction_support"] = matrix[
        "expimap_ssgsea_direction_match"
    ].fillna(False)
    matrix["preranked_gsea_direction_support"] = matrix[
        "expimap_gsea_direction_match"
    ].fillna(False)
    matrix["heldout_direction_support"] = matrix[
        "heldout_project_direction_concordance"
    ].ge(MINIMUM_HELDOUT_CONCORDANCE)
    matrix["seed_direction_support"] = (
        matrix["all_three_seeds_available"].fillna(False)
        & matrix["all_three_direction_concordant"].fillna(False)
    )
    matrix["composition_proxy_support"] = (
        matrix["adjusted_direction_matches_unadjusted"].fillna(False)
        & matrix["absolute_effect_ratio_adjusted_to_unadjusted"].ge(
            MINIMUM_COMPOSITION_EFFECT_RATIO
        )
    )
    support_columns = [
        "ssgsea_direction_support",
        "preranked_gsea_direction_support",
        "heldout_direction_support",
        "seed_direction_support",
        "composition_proxy_support",
    ]
    matrix["robustness_support_count"] = matrix[support_columns].sum(axis=1)
    conventional = matrix[
        ["ssgsea_direction_support", "preranked_gsea_direction_support"]
    ].sum(axis=1)
    internal = (
        matrix["heldout_direction_support"]
        & matrix["seed_direction_support"]
        & matrix["composition_proxy_support"]
    )
    matrix["robustness_status"] = "sensitivity-dependent"
    matrix.loc[
        internal & conventional.lt(2), "robustness_status"
    ] = "internally robust, incomplete conventional support"
    matrix.loc[internal & conventional.eq(2), "robustness_status"] = "triangulated"
    matrix.loc[
        conventional.eq(2) & matrix["heldout_direction_support"] & ~internal,
        "robustness_status",
    ] = "method-supported, model-sensitive"
    return matrix.sort_values(
        ["tissue", "primary_top_decile", "robustness_support_count", "seed_effect_median"],
        ascending=[True, False, False, False],
    )


def plot_top_pathways(matrix: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(17, 10), constrained_layout=True)
    for ax, config in zip(axes, CONFIGS):
        frame = matrix.loc[matrix["tissue"].eq(config.tissue)].copy()
        frame = frame.assign(abs_effect=frame["seed_effect_median"].abs()).nlargest(
            20, "abs_effect"
        )
        frame = frame.sort_values("seed_effect_median")
        y = np.arange(len(frame))
        colors = [
            "#2f7d4a" if status == "triangulated" else
            "#5685a8" if status == "internally robust, incomplete conventional support" else
            "#d18b32" if status == "method-supported, model-sensitive" else
            "#a5aaad"
            for status in frame["robustness_status"]
        ]
        ax.hlines(
            y,
            frame["seed_effect_minimum"],
            frame["seed_effect_maximum"],
            color=colors,
            lw=2,
            alpha=0.8,
        )
        ax.scatter(frame["seed_effect_median"], y, c=colors, s=48, zorder=3)
        ax.axvline(0, color="#575d60", lw=0.8)
        labels = [
            label if len(label) <= 54 else label[:51] + "..."
            for label in frame["display_term"]
        ]
        ax.set_yticks(y, labels, fontsize=8)
        ax.set_xlabel("Study-balanced flight minus ground expiMap shift")
        ax.set_title(config.tissue.title(), fontweight="bold")
        ax.grid(axis="x", color="#e1e4e5", linewidth=0.7)
    fig.suptitle(
        "Kidney and spleen: largest corrected-HVG pathway shifts across three seeds",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(OUTPUT_DIR / "top_pathway_seed_stability.png", dpi=300)
    fig.savefig(OUTPUT_DIR / "top_pathway_seed_stability.pdf")
    plt.close(fig)


def write_report(
    matrix: pd.DataFrame,
    seed_summary: pd.DataFrame,
    model_comparison: pd.DataFrame,
    manifests: dict[str, dict],
) -> None:
    lines = [
        "# Kidney and spleen expiMap reassessment",
        "",
        "This reassessment corrects the two main comparability issues in the historical screen: kidney now uses all eligible ARCHS4 reference samples, and spleen uses batch-aware HVG selection after excluding singleton ARCHS4 series from HVG ranking only. The primary spleen contrast also excludes OSD-288 because recorded strain is disjoint between flight and ground groups; the full query remains available as a sensitivity. It applies three full training seeds, ssGSEA, preranked GSEA, held-out-project direction checks, and broad atlas-derived composition-proxy adjustment.",
        "",
        "The pathway ranking below uses relative effect magnitude rather than a hard FDR gate. A top-decile label means the pathway is among the largest absolute seed-2020 expiMap shifts within that tissue.",
        "",
    ]
    for config in CONFIGS:
        tissue = config.tissue
        manifest = manifests[tissue]
        seed = seed_summary.loc[
            seed_summary["tissue"].eq(tissue)
            & seed_summary["pathway_set"].eq("primary_top_decile")
        ].iloc[0]
        comparison = model_comparison.loc[model_comparison["tissue"].eq(tissue)].iloc[0]
        top = matrix.loc[
            matrix["tissue"].eq(tissue) & matrix["primary_top_decile"]
        ].sort_values(
            ["robustness_support_count", "seed_effect_median"],
            ascending=[False, False],
        )
        lines.extend(
            [
                f"## {tissue.title()}",
                "",
                f"The corrected model uses {manifest['n_reference_samples']:,} ARCHS4 samples, {manifest['n_terms_after_hvg_filter']} Reactome programs, and {manifest['n_query_samples']} mapped OSDR samples. The primary effect summary uses {5 if tissue == 'spleen' else 6} unconfounded projects. Across the primary top-decile pathways, {seed.all_three_seed_direction_agreement:.0%} retained one direction across all three complete training runs.",
                f"Compared with the historical HVG model, the corrected model has Spearman r={comparison.old_new_spearman_rho:.2f}, {comparison.old_new_direction_agreement:.0%} directional agreement, and top-decile Jaccard overlap {comparison.old_new_top_decile_jaccard:.2f}.",
                "",
                "### Largest pathways after triangulation",
                "",
            ]
        )
        for row in top.head(20).itertuples(index=False):
            lines.append(
                f"- **{row.display_term}:** median shift {row.seed_effect_median:+.3f}; "
                f"{row.robustness_support_count}/5 directional checks; "
                f"{row.robustness_status}; held-out projects "
                f"{row.heldout_project_direction_concordance:.0%} concordant."
            )
        lines.append("")
    lines.extend(
        [
            "## Interpretation rule",
            "",
            "A pathway is suitable for biological follow-up when it combines a large relative effect with repeatable direction across seeds and projects and support from at least one conventional enrichment method. Off-tissue labels, tiny effects, and pathways that reverse across complete training runs are not promoted even if an isolated statistic is favorable.",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    seed_accession = pd.read_csv(OUTPUT_DIR / "seed_accession_effects.tsv.gz", sep="\t")
    retained_seed_rows = []
    for tissue, frame in seed_accession.groupby("tissue", observed=True):
        excluded = CONFOUNDED_ACCESSIONS[str(tissue)]
        retained_seed_rows.append(
            frame.loc[~frame["accession"].astype(str).isin(excluded)]
        )
    primary_seed_accession = pd.concat(retained_seed_rows, ignore_index=True)
    primary_seed_effects = (
        primary_seed_accession.groupby(["tissue", "seed", "term"])["effect"]
        .agg(
            accession_balanced_effect="mean",
            n_accessions="size",
            accessions_positive=lambda values: int((values > 0).sum()),
            accessions_negative=lambda values: int((values < 0).sum()),
        )
        .reset_index()
    )
    seeds, seed_summary = consensus_table(primary_seed_effects)
    primary_seed_effects.to_csv(
        OUTPUT_DIR / "primary_seed_pathway_effects.tsv", sep="\t", index=False
    )
    seeds.to_csv(OUTPUT_DIR / "primary_seed_consensus.tsv", sep="\t", index=False)
    seed_summary.to_csv(
        OUTPUT_DIR / "primary_seed_summary.tsv", sep="\t", index=False
    )

    atlas = ad.read_h5ad(TMS_PATH, backed="r")
    atlas_obs = atlas.obs.copy()
    atlas_genes = atlas.var_names.astype(str)
    try:
        atlas.file.close()
    except AttributeError:
        pass

    benchmarks = []
    heldout_folds = []
    heldout_summaries = []
    compositions = []
    markers = []
    marker_scores = []
    comparisons = []
    manifests = {}
    for config in CONFIGS:
        print(f"loading {config.tissue}", flush=True)
        data = load_data(config)
        benchmark, folds, heldout = method_benchmark(data)
        composition, tissue_markers, tissue_marker_scores = composition_sensitivity(
            data, atlas_obs, atlas_genes
        )
        benchmarks.append(benchmark)
        heldout_folds.append(folds)
        heldout_summaries.append(heldout)
        compositions.append(composition)
        markers.append(tissue_markers)
        marker_scores.append(tissue_marker_scores)
        comparisons.append(old_model_comparison(config, benchmark))
        manifests[config.tissue] = json.loads(
            (
                config.base_dir
                / "reassessment_hvg_2000/input/tutorial_hvg_input_manifest.json"
            ).read_text(encoding="utf-8")
        )

    benchmark = pd.concat(benchmarks, ignore_index=True)
    heldout_fold = pd.concat(heldout_folds, ignore_index=True)
    heldout = pd.concat(heldout_summaries, ignore_index=True)
    composition = pd.concat(compositions, ignore_index=True)
    marker = pd.concat(markers, ignore_index=True)
    sample_markers = pd.concat(marker_scores, ignore_index=True)
    comparison = pd.DataFrame(comparisons)
    matrix = evidence_matrix(benchmark, heldout, composition, seeds)

    benchmark.to_csv(OUTPUT_DIR / "conventional_method_benchmark.tsv", sep="\t", index=False)
    heldout_fold.to_csv(
        OUTPUT_DIR / "heldout_project_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    heldout.to_csv(OUTPUT_DIR / "heldout_project_summary.tsv", sep="\t", index=False)
    composition.to_csv(
        OUTPUT_DIR / "composition_proxy_adjusted_effects.tsv", sep="\t", index=False
    )
    marker.to_csv(OUTPUT_DIR / "atlas_compartment_markers.tsv", sep="\t", index=False)
    sample_markers.to_csv(
        OUTPUT_DIR / "sample_composition_proxy_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    comparison.to_csv(OUTPUT_DIR / "historical_model_comparison.tsv", sep="\t", index=False)
    matrix.to_csv(OUTPUT_DIR / "pathway_evidence_matrix.tsv", sep="\t", index=False)
    plot_top_pathways(matrix)
    write_report(matrix, seed_summary, comparison, manifests)
    print(
        matrix.loc[matrix["primary_top_decile"]]
        .groupby(["tissue", "robustness_status"], observed=True)
        .size()
        .rename("pathways")
        .to_string(),
        flush=True,
    )


if __name__ == "__main__":
    main()
