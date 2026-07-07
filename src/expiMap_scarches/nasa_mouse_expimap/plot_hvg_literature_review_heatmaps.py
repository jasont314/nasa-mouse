"""Plot HVG expiMap heatmaps using pathway-by-pathway literature review labels.

Unlike ``plot_hvg_interpretation_heatmaps.py``, this script does not infer
label colors from keyword patterns. It reads curated review TSVs and only keeps
gray as the quantitative low-effect category from the original label tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class ModelConfig:
    name: str
    run_dir: Path


MODELS = (
    ModelConfig(
        name="liver_hvg",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_liver/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    ),
    ModelConfig(
        name="skin_hvg",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_skin/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    ),
    ModelConfig(
        name="thymus_hvg",
        run_dir=ROOT
        / "outputs/expimap_archs4_reference_osdr_query_thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020",
    ),
    ModelConfig(
        name="soleus_hvg",
        run_dir=ROOT
        / "outputs/expimap_muscle_targeted_combined_min8/tutorial_hvg_soleus_2000/query_nb_250epoch_seed2020",
    ),
)


REVIEW_DIR = ROOT / "presentation/expimap/literature_review_manual"
PRESENTATION_DIR = ROOT / "presentation/expimap/literature_reviewed_hvg"
BASE_LABEL_DIR = ROOT / "presentation/expimap/annotated_hvg"


CATEGORY_COLORS = {
    "green_prior_aligned": "#16803c",
    "blue_plausible_complementary": "#1f63b5",
    "orange_uncertain_no_direct_support": "#a6611a",
    "red_conflict_or_artifact": "#c51b29",
    "gray_low_effect": "#8f8f8f",
}

CATEGORY_LABELS = {
    "green_prior_aligned": "prior literature aligned",
    "blue_plausible_complementary": "plausible/complementary",
    "orange_uncertain_no_direct_support": "uncertain/no direct support",
    "red_conflict_or_artifact": "conflict/artifact",
    "gray_low_effect": "low/no relative effect",
}


def read_matrix(config: ModelConfig) -> pd.DataFrame:
    matrix_path = config.run_dir / "analysis/all_program_accession_flt_minus_gc_matrix_signed_order.tsv"
    return pd.read_csv(matrix_path, sep="\t").set_index("term")


def read_reviewed_labels(config: ModelConfig, matrix: pd.DataFrame) -> pd.DataFrame:
    base_path = BASE_LABEL_DIR / f"{config.name}_pathway_interpretation_labels.tsv"
    review_path = REVIEW_DIR / f"{config.name}_manual_literature_review.tsv"
    if not base_path.exists():
        raise FileNotFoundError(base_path)
    if not review_path.exists():
        raise FileNotFoundError(review_path)

    base = pd.read_csv(base_path, sep="\t")
    review = pd.read_csv(review_path, sep="\t")
    required = {
        "review_row",
        "term",
        "observed_direction",
        "reviewed_category",
        "literature_alignment",
        "direction_assessment",
        "confidence",
        "rationale",
        "citations",
    }
    missing = required.difference(review.columns)
    if missing:
        raise ValueError(f"{review_path} is missing columns: {sorted(missing)}")

    valid = set(CATEGORY_COLORS).difference({"gray_low_effect"})
    invalid_categories = sorted(set(review["reviewed_category"].dropna()) - valid)
    if invalid_categories:
        raise ValueError(f"{review_path} has invalid categories: {invalid_categories}")

    base = base.set_index("term", drop=False)
    review = review.set_index("term", drop=False)

    non_low_terms = set(base.loc[~base["category"].eq("low_or_no_effect"), "term"])
    reviewed_terms = set(review.index)
    missing_reviews = sorted(non_low_terms - reviewed_terms)
    extra_reviews = sorted(reviewed_terms - set(base.index))
    if missing_reviews:
        raise ValueError(
            f"{review_path} does not review {len(missing_reviews)} non-low terms; "
            f"first missing: {missing_reviews[:5]}"
        )
    if extra_reviews:
        raise ValueError(
            f"{review_path} has {len(extra_reviews)} terms absent from base labels; "
            f"first extra: {extra_reviews[:5]}"
        )

    records = []
    for term in matrix.index.astype(str):
        if term in review.index:
            row = review.loc[term]
            category = row["reviewed_category"]
            records.append(
                {
                    "term": term,
                    "reviewed_category": category,
                    "reviewed_category_label": CATEGORY_LABELS[category],
                    "observed_direction": row["observed_direction"],
                    "literature_alignment": row["literature_alignment"],
                    "direction_assessment": row["direction_assessment"],
                    "confidence": row["confidence"],
                    "review_rationale": row["rationale"],
                    "citations": row["citations"],
                }
            )
        else:
            base_row = base.loc[term]
            if not base_row["category"] == "low_or_no_effect":
                raise ValueError(f"Unexpected unreviewed non-low term: {term}")
            records.append(
                {
                    "term": term,
                    "reviewed_category": "gray_low_effect",
                    "reviewed_category_label": CATEGORY_LABELS["gray_low_effect"],
                    "observed_direction": "",
                    "literature_alignment": "not reviewed; low relative effect",
                    "direction_assessment": "not assessed",
                    "confidence": "",
                    "review_rationale": "Kept gray because original quantitative effect threshold marked this row low/no relative effect.",
                    "citations": "",
                }
            )

    reviewed = pd.DataFrame(records)
    stats_cols = [
        "term",
        "flight_minus_ground",
        "welch_fdr",
        "abs_effect",
        "mean_abs_accession_effect",
        "abs_mean_accession_effect",
        "low_effect_threshold",
        "category",
        "category_label",
    ]
    # base keeps term both as index and column for lookup convenience;
    # reset the index before merging to avoid pandas treating it ambiguously.
    return reviewed.merge(base[stats_cols].reset_index(drop=True), on="term", how="left")


def plot_one(
    config: ModelConfig,
    matrix: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    mode: str,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels_by_term = labels.set_index("term")
    if mode == "no_gray_red":
        keep = labels_by_term.loc[
            matrix.index,
            "reviewed_category",
        ].ne("gray_low_effect") & labels_by_term.loc[
            matrix.index,
            "reviewed_category",
        ].ne(
            "red_conflict_or_artifact"
        )
        plot_matrix = matrix.loc[keep.to_numpy()]
        suffix = "literature_review_no_gray_red"
        title = "expiMap FLT-GC pathway shifts, literature-reviewed labels; gray/red removed"
        footnote = (
            "Rows ordered by signed study-mean FLT-GC effect; literature-reviewed gray low/no-effect "
            "and red conflict/artifact rows removed"
        )
        row_fontsize = 5.0
        row_scale = 0.11
        max_height = 90.0
    elif mode == "blue_only":
        keep = labels_by_term.loc[matrix.index, "reviewed_category"].eq(
            "blue_plausible_complementary"
        )
        plot_matrix = matrix.loc[keep.to_numpy()]
        suffix = "literature_review_blue_only"
        title = "expiMap FLT-GC pathway shifts, blue plausible/complementary labels only"
        footnote = (
            "Rows ordered by signed study-mean FLT-GC effect; only literature-reviewed "
            "blue plausible/complementary rows shown"
        )
        row_fontsize = 5.5
        row_scale = 0.13
        max_height = 90.0
    elif mode == "all_labels":
        plot_matrix = matrix
        suffix = "literature_review_all_labels"
        title = "All expiMap FLT-GC pathway shifts, literature-reviewed labels"
        footnote = "Rows ordered by signed study-mean FLT-GC effect (FLT-up to FLT-down)"
        row_fontsize = 4.0
        row_scale = 0.08
        max_height = 140.0
    else:
        raise ValueError(f"Unknown plot mode: {mode}")

    ordered_terms = plot_matrix.index.astype(str).tolist()
    categories = labels_by_term.loc[ordered_terms, "reviewed_category"].astype(str).tolist()
    row_colors = [CATEGORY_COLORS[category] for category in categories]
    data = plot_matrix.to_numpy(dtype=float)

    n_rows, n_cols = plot_matrix.shape
    fig_width = max(18.0, min(32.0, 15.0 + 1.0 * n_cols))
    fig_height = min(max(8.0, n_rows * row_scale + 2.0), max_height)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#f2f2f2")
    image = ax.imshow(
        np.ma.masked_invalid(data),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-2.5,
        vmax=2.5,
    )

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(ordered_terms, fontsize=row_fontsize)
    for tick, color in zip(ax.get_yticklabels(), row_colors):
        tick.set_color(color)
    ax.tick_params(axis="y", length=0, pad=2)

    accessions = plot_matrix.columns.astype(str).tolist()
    ax.set_xticks(np.arange(len(accessions)))
    ax.set_xticklabels(accessions, rotation=90, ha="center", va="top", fontsize=7.0)
    ax.set_xlabel("OSD accession / study", fontsize=8)
    ax.set_title(title, fontsize=10, pad=10)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("mean FLT - mean GC pathway score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.text(0.63, 0.012, footnote, ha="center", va="bottom", fontsize=6)
    fig.subplots_adjust(left=0.56, right=0.89, top=0.985, bottom=0.055)

    output = PRESENTATION_DIR / f"{config.name}_{suffix}.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def main() -> None:
    PRESENTATION_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    outputs = []
    for config in MODELS:
        matrix = read_matrix(config)
        labels = read_reviewed_labels(config, matrix)
        labels_path = PRESENTATION_DIR / f"{config.name}_literature_review_labels.tsv"
        labels.to_csv(labels_path, sep="\t", index=False)
        blue_labels = labels[
            labels["reviewed_category"].eq("blue_plausible_complementary")
        ].copy()
        blue_labels.to_csv(
            PRESENTATION_DIR / f"{config.name}_literature_review_blue_only_labels.tsv",
            sep="\t",
            index=False,
        )

        source_path = REVIEW_DIR / f"{config.name}_sources.md"
        if source_path.exists():
            shutil.copy2(source_path, PRESENTATION_DIR / source_path.name)

        outputs.append(plot_one(config, matrix, labels, mode="all_labels"))
        outputs.append(plot_one(config, matrix, labels, mode="no_gray_red"))
        outputs.append(plot_one(config, matrix, labels, mode="blue_only"))

        counts = labels["reviewed_category"].value_counts().to_dict()
        rows.append(
            {
                "model": config.name,
                "total_rows": len(labels),
                "green_prior_aligned": counts.get("green_prior_aligned", 0),
                "blue_plausible_complementary": counts.get("blue_plausible_complementary", 0),
                "orange_uncertain_no_direct_support": counts.get(
                    "orange_uncertain_no_direct_support", 0
                ),
                "red_conflict_or_artifact": counts.get("red_conflict_or_artifact", 0),
                "gray_low_effect": counts.get("gray_low_effect", 0),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(PRESENTATION_DIR / "literature_review_color_summary.tsv", sep="\t", index=False)
    for output in outputs:
        print(output.relative_to(ROOT))
    print((PRESENTATION_DIR / "literature_review_color_summary.tsv").relative_to(ROOT))


if __name__ == "__main__":
    main()
