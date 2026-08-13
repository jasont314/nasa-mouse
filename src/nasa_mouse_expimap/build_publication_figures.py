"""Build final-size, publication-oriented figures for the ASGSR expiMap paper."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import textwrap

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from .build_asgsr_paper import ROOT
from .integrate_reassessed_tissues_paper import (
    DISPLAY_LABELS,
    MAIN_TISSUES,
    RETAINED_TERMS,
    load_retained_evidence,
    project_effects,
)


PAPER_DIR = ROOT / "paper/asgsr_expimap_hvg"
SOURCE_DIR = PAPER_DIR / "source_data"
FIGURE_DIR = PAPER_DIR / "figures"
PROCESS_SUMMARY_DIR = ROOT / "outputs/expimap/analyses/publication_process_summary"
REASSESSMENT_DIR = ROOT / "outputs/expimap/analyses/kidney_spleen_reassessment"

FIGURE_WIDTH = 7.2

GROUND_COLOR = "#0072B2"
FLIGHT_COLOR = "#D55E00"
REFERENCE_COLOR = "#C3C8CB"
TEXT_COLOR = "#202629"
GRID_COLOR = "#DCE1E3"
ROLE_COLORS = {
    "aligned": "#009E73",
    "complementary": "#0072B2",
    "context_sensitive": "#D55E00",
}
ROLE_MARKERS = {
    "aligned": "o",
    "complementary": "s",
    "context_sensitive": "^",
}
TISSUE_COLORS = {
    "thymus": "#6B4C9A",
    "skin": "#C05A2B",
    "liver": "#16837A",
    "spleen": "#9B4F65",
    "kidney": "#287B8E",
    "soleus": "#5A6F8A",
}
CHECK_COLUMNS = (
    ("ssgsea_direction_support", "ssGSEA"),
    ("preranked_gsea_direction_support", "GSEA"),
    ("heldout_direction_support", "Held-out"),
    ("seed_direction_support", "3 seeds"),
    ("composition_proxy_support", "Composition"),
)


@dataclass(frozen=True)
class LatentConfig:
    tissue: str
    reference_path: Path
    query_path: Path
    representative_term: str


LATENT_CONFIGS = (
    LatentConfig(
        "thymus",
        ROOT
        / "outputs/expimap/runs/reference_query/thymus/tutorial_hvg_2000/reference_nb_400epoch_seed2020/trained_input_with_scores.h5ad",
        ROOT
        / "outputs/expimap/runs/reference_query/thymus/tutorial_hvg_2000/query_nb_250epoch_seed2020/mapped_query_with_scores.h5ad",
        "R-MMU-73894_DNA_REPAIR",
    ),
    LatentConfig(
        "skin",
        ROOT
        / "outputs/expimap/runs/reference_query/skin/tutorial_hvg_2000/reference_nb_400epoch_seed2020/trained_input_with_scores.h5ad",
        ROOT
        / "outputs/expimap/runs/reference_query/skin/tutorial_hvg_2000/query_nb_250epoch_seed2020/mapped_query_with_scores.h5ad",
        "R-MMU-421270_CELL_CELL_JUNCTION_ORGANIZATION",
    ),
    LatentConfig(
        "liver",
        ROOT
        / "outputs/expimap/runs/reference_query/liver/tutorial_hvg_2000/reference_nb_400epoch_seed2020/trained_input_with_scores.h5ad",
        ROOT
        / "outputs/expimap/runs/reference_query/liver/tutorial_hvg_2000/query_nb_250epoch_seed2020_primary_deduplicated/mapped_query_with_scores.h5ad",
        "R-MMU-2132295_MHC_CLASS_II_ANTIGEN_PRESENTATION",
    ),
    LatentConfig(
        "spleen",
        ROOT
        / "outputs/expimap/runs/reference_query/spleen/reassessment_hvg_2000/reference_nb_400epoch_seed2020/pathway_scores.tsv",
        ROOT
        / "outputs/expimap/runs/reference_query/spleen/reassessment_hvg_2000/query_nb_250epoch_seed2020/query_pathway_scores.tsv",
        "R-MMU-202403_TCR_SIGNALING",
    ),
)


OBSOLETE_FIGURE_STEMS = (
    "figure_1_workflow",
    "figure_2_tissue_pathway_shifts",
    "figure_3_evidence_map",
    "figure_4_kidney_spleen_reassessment",
    "figure_4_primary_analysis_sensitivity",
    "figure_5_complementary_process_model",
    "figure_6_generated_biological_processes",
    "figure_6_skin_protocol_context",
    "figure_7_skin_protocol_context",
    "figure_s8_generated_biological_processes",
    "figure_s8_process_summary",
    "figure_s9_original_tissue_sensitivity",
    "figure_s10_program_score_distributions",
)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.0,
            "axes.titlesize": 9.0,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.edgecolor": "#5B6468",
            "axes.labelcolor": TEXT_COLOR,
            "text.color": TEXT_COLOR,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, name: str, output_dir: Path = FIGURE_DIR) -> None:
    """Preserve the authored final dimensions in both raster and vector output."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{name}.png", dpi=300, facecolor="white")
    fig.savefig(output_dir / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=9.0,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def clean_axis(ax: plt.Axes, *, keep_left: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if not keep_left:
        ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID_COLOR, linewidth=0.6)
    ax.set_axisbelow(True)


def format_fdr(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def plot_workflow(scope: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(FIGURE_WIDTH, 4.65), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(0.88, 1.12), width_ratios=(1, 1))
    workflow_ax = fig.add_subplot(grid[0, :])
    architecture_ax = fig.add_subplot(grid[1, 0])
    cohort_ax = fig.add_subplot(grid[1, 1])

    workflow_ax.set_xlim(0, 1)
    workflow_ax.set_ylim(0, 1)
    workflow_ax.axis("off")
    workflow_boxes = (
        (0.02, 0.55, 0.22, 0.33, "ARCHS4", "Tissue-matched non-spaceflight RNA-seq", "#E8F1F4"),
        (0.02, 0.12, 0.22, 0.33, "Reactome", "Mouse gene-program mask", "#E8F2EA"),
        (0.39, 0.24, 0.25, 0.52, "expiMap", "Reference training\n2,000 selected HVGs", "#ECE9F4"),
        (0.76, 0.24, 0.22, 0.52, "NASA OSDR", "Flight and ground\nquery mapping", "#F6ECE5"),
    )
    for x, y, width, height, title, body, color in workflow_boxes:
        workflow_ax.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.014",
                facecolor=color,
                edgecolor="#586267",
                linewidth=0.8,
            )
        )
        workflow_ax.text(
            x + width / 2,
            y + height * 0.68,
            title,
            ha="center",
            va="center",
            fontsize=8.3,
            fontweight="bold",
        )
        workflow_ax.text(
            x + width / 2,
            y + height * 0.32,
            textwrap.fill(body, width=26) if "\n" not in body else body,
            ha="center",
            va="center",
            fontsize=7.1,
            linespacing=1.25,
        )
    for source_y in (0.715, 0.285):
        workflow_ax.add_patch(
            FancyArrowPatch(
                (0.252, source_y),
                (0.378, 0.50),
                arrowstyle="-|>",
                mutation_scale=10,
                color="#505A5F",
                linewidth=0.9,
                connectionstyle="arc3,rad=0.08" if source_y > 0.5 else "arc3,rad=-0.08",
            )
        )
    workflow_ax.add_patch(
        FancyArrowPatch(
            (0.652, 0.50),
            (0.748, 0.50),
            arrowstyle="-|>",
            mutation_scale=10,
            color="#505A5F",
            linewidth=0.9,
        )
    )
    panel_label(workflow_ax, "a", x=-0.015, y=0.96)

    architecture_ax.set_xlim(-0.8, 7.2)
    architecture_ax.set_ylim(-0.8, 6.5)
    architecture_ax.axis("off")
    mask = np.array(
        [
            [1, 0, 0, 0, 1],
            [1, 1, 0, 0, 0],
            [0, 1, 0, 1, 0],
            [0, 1, 1, 0, 0],
            [0, 0, 1, 0, 1],
            [0, 0, 1, 1, 0],
            [0, 0, 0, 1, 1],
        ],
        dtype=int,
    )
    architecture_ax.imshow(
        mask,
        cmap=LinearSegmentedColormap.from_list("mask", ["#F1F3F4", "#286F9E"]),
        vmin=0,
        vmax=1,
        extent=(0, 4.0, 0, 5.4),
        origin="lower",
        aspect="auto",
    )
    architecture_ax.set_xticks(np.linspace(0.4, 3.6, 5))
    architecture_ax.set_xticklabels(["GP1", "GP2", "GP3", "GP4", "GP5"], fontsize=6.8)
    architecture_ax.set_yticks(np.linspace(0.39, 5.01, 7))
    architecture_ax.set_yticklabels([f"gene {index}" for index in range(1, 8)], fontsize=6.8)
    architecture_ax.tick_params(length=0)
    architecture_ax.text(2.0, 5.95, "Sparse Reactome mask", ha="center", fontsize=8.1, fontweight="bold")
    architecture_ax.add_patch(
        FancyArrowPatch((4.2, 2.7), (5.05, 2.7), arrowstyle="-|>", mutation_scale=10, color="#505A5F")
    )
    architecture_ax.text(6.1, 3.2, "Interpretable\nlatent programs", ha="center", va="center", fontsize=7.2)
    architecture_ax.text(6.1, 1.8, "Accession-conditioned\nquery mapping", ha="center", va="center", fontsize=7.2)
    panel_label(architecture_ax, "b", x=-0.08, y=1.02)

    plot_scope = scope.loc[scope["tissue"].isin((*MAIN_TISSUES, "kidney"))].copy()
    order = ["thymus", "skin", "liver", "spleen", "kidney"]
    plot_scope["order"] = plot_scope["tissue"].map({name: i for i, name in enumerate(order)})
    plot_scope = plot_scope.sort_values("order")
    y = np.arange(len(plot_scope))[::-1]
    cohort_ax.scatter(
        plot_scope["reference_samples"],
        y + 0.11,
        marker="s",
        s=38,
        color="#4C78A8",
        label="ARCHS4 reference",
        zorder=3,
    )
    cohort_ax.scatter(
        plot_scope["primary_effect_samples"],
        y - 0.11,
        marker="o",
        s=38,
        color="#E17C3A",
        label="OSDR query",
        zorder=3,
    )
    for position, row in zip(y, plot_scope.itertuples(index=False)):
        cohort_ax.text(float(row.reference_samples) * 1.12, position + 0.11, f"{int(row.reference_samples):,}", va="center", fontsize=6.8)
        cohort_ax.text(float(row.primary_effect_samples) * 1.12, position - 0.11, f"{int(row.primary_effect_samples):,}", va="center", fontsize=6.8)
    cohort_ax.set_xscale("log")
    cohort_ax.set_xlim(60, 12000)
    cohort_ax.set_yticks(y)
    cohort_ax.set_yticklabels([name.title() for name in plot_scope["tissue"]], fontsize=7.2)
    cohort_ax.set_xlabel("Samples, log scale")
    cohort_ax.set_title(
        "Analysis scope",
        loc="left",
        y=1.15,
        fontsize=8.1,
        fontweight="bold",
    )
    cohort_ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.06),
        ncol=2,
        fontsize=6.8,
        borderaxespad=0,
    )
    clean_axis(cohort_ax)
    panel_label(cohort_ax, "c", x=-0.18, y=1.15)
    save_figure(fig, "figure_1_workflow_architecture")


def _load_latent_scores(
    path: Path,
    *,
    query: bool,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, np.ndarray]:
    if path.suffix == ".h5ad":
        adata = ad.read_h5ad(path, backed="r")
        terms = np.asarray(adata.uns["terms"], dtype=str)
        score_key = "X_expimap_query" if query else "X_expimap"
        scores = np.asarray(adata.obsm[score_key], dtype=float)
        obs = adata.obs.copy()
        sample_ids = adata.obs_names.astype(str).to_numpy()
        adata.file.close()
        return terms, scores, obs, sample_ids

    frame = pd.read_csv(path, sep="\t")
    terms = np.asarray(
        [column for column in frame.columns if column.startswith("R-MMU-")],
        dtype=str,
    )
    if not len(terms):
        raise RuntimeError(f"No Reactome score columns found in {path}")
    sample_column = "obs_name" if "obs_name" in frame else "profile_id"
    return (
        terms,
        frame.loc[:, terms].to_numpy(dtype=float),
        frame,
        frame[sample_column].astype(str).to_numpy(),
    )


def build_latent_mapping_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    coordinate_frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    for config in LATENT_CONFIGS:
        reference_terms, reference_scores, reference_obs, reference_ids = (
            _load_latent_scores(config.reference_path, query=False)
        )
        query_terms, query_scores, query_obs, query_ids = _load_latent_scores(
            config.query_path,
            query=True,
        )
        if not np.array_equal(reference_terms, query_terms):
            raise RuntimeError(f"Reference-query term mismatch for {config.tissue}")

        scaler = StandardScaler().fit(reference_scores)
        scaled_reference = scaler.transform(reference_scores)
        scaled_query = scaler.transform(query_scores)
        n_components = min(20, scaled_reference.shape[0] - 1, scaled_reference.shape[1])
        pca = PCA(n_components=n_components, random_state=2026).fit(scaled_reference)
        reference_pcs = pca.transform(scaled_reference)
        query_pcs = pca.transform(scaled_query)

        neighbors = NearestNeighbors(n_neighbors=2).fit(reference_pcs)
        reference_nn = neighbors.kneighbors(reference_pcs, return_distance=True)[0][:, 1]
        query_nn = neighbors.kneighbors(query_pcs, n_neighbors=1, return_distance=True)[0][:, 0]
        threshold = float(np.quantile(reference_nn, 0.95))
        query_within = query_nn <= threshold

        reference_frame = pd.DataFrame(
            {
                "tissue": config.tissue,
                "source": "ARCHS4_reference",
                "sample_id": reference_ids,
                "project": reference_obs["series_id"].astype(str).to_numpy(),
                "condition": "reference",
                "PC1": reference_pcs[:, 0],
                "PC2": reference_pcs[:, 1],
                "nearest_reference_distance_20pc": reference_nn,
                "within_reference_95pct_nn_distance": True,
            }
        )
        query_frame = pd.DataFrame(
            {
                "tissue": config.tissue,
                "source": "OSDR_query",
                "sample_id": query_ids,
                "project": query_obs["id.accession"].astype(str).to_numpy(),
                "condition": query_obs["condition_inferred"].astype(str).to_numpy(),
                "PC1": query_pcs[:, 0],
                "PC2": query_pcs[:, 1],
                "nearest_reference_distance_20pc": query_nn,
                "within_reference_95pct_nn_distance": query_within,
            }
        )
        coordinate_frames.extend((reference_frame, query_frame))
        qc_rows.append(
            {
                "tissue": config.tissue,
                "reference_samples": len(reference_frame),
                "query_samples": len(query_frame),
                "programs": len(reference_terms),
                "pc1_variance_fraction": pca.explained_variance_ratio_[0],
                "pc2_variance_fraction": pca.explained_variance_ratio_[1],
                "reference_nearest_neighbor_95pct_threshold_20pc": threshold,
                "query_within_reference_95pct_nn_fraction": float(query_within.mean()),
                "query_median_to_reference_nn_distance_20pc": float(np.median(query_nn)),
                "reference_median_nn_distance_20pc": float(np.median(reference_nn)),
            }
        )

    coordinates = pd.concat(coordinate_frames, ignore_index=True)
    qc = pd.DataFrame(qc_rows)
    coordinates.to_csv(
        SOURCE_DIR / "table_s31_latent_mapping_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    qc.to_csv(SOURCE_DIR / "table_s32_latent_mapping_qc.tsv", sep="\t", index=False)
    return coordinates, qc


def plot_latent_mapping(coordinates: pd.DataFrame, qc: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH, 5.85), layout="constrained")
    rng = np.random.default_rng(2026)
    for label, ax, config in zip("abcd", axes.flat, LATENT_CONFIGS):
        frame = coordinates.loc[coordinates["tissue"].eq(config.tissue)]
        reference = frame.loc[frame["source"].eq("ARCHS4_reference")]
        query = frame.loc[frame["source"].eq("OSDR_query")]
        if len(reference) > 2200:
            reference = reference.iloc[rng.choice(len(reference), 2200, replace=False)]
        ax.scatter(
            reference["PC1"],
            reference["PC2"],
            s=5,
            color=REFERENCE_COLOR,
            alpha=0.38,
            linewidth=0,
            rasterized=True,
            zorder=1,
        )
        for _, project in query.groupby("project", observed=True):
            means = project.groupby("condition")[["PC1", "PC2"]].mean()
            if {"flight", "ground_control"}.issubset(means.index):
                ax.plot(
                    [means.loc["ground_control", "PC1"], means.loc["flight", "PC1"]],
                    [means.loc["ground_control", "PC2"], means.loc["flight", "PC2"]],
                    color="#7B858A",
                    linewidth=0.7,
                    alpha=0.7,
                    zorder=2,
                )
        condition_specs = (
            ("ground_control", GROUND_COLOR, "o"),
            ("flight", FLIGHT_COLOR, "^"),
        )
        for condition, color, marker in condition_specs:
            subset = query.loc[query["condition"].eq(condition)]
            ax.scatter(
                subset["PC1"],
                subset["PC2"],
                s=18,
                color=color,
                marker=marker,
                edgecolor="white",
                linewidth=0.35,
                alpha=0.82,
                zorder=3,
            )
        row = qc.loc[qc["tissue"].eq(config.tissue)].iloc[0]
        ax.set_title(
            f"{config.tissue.title()}  |  {int(row.query_samples)} OSDR, {int(row.reference_samples):,} ARCHS4",
            loc="left",
            fontweight="bold",
        )
        ax.text(
            0.02,
            0.03,
            f"20-PC NN coverage: {row.query_within_reference_95pct_nn_fraction:.0%}",
            transform=ax.transAxes,
            fontsize=6.9,
            ha="left",
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
        )
        ax.set_xlabel(f"Reference PC1 ({row.pc1_variance_fraction:.1%})")
        ax.set_ylabel(f"Reference PC2 ({row.pc2_variance_fraction:.1%})")
        clean_axis(ax, keep_left=True)
        panel_label(ax, label)

    handles = (
        Line2D([0], [0], marker=".", linestyle="none", color=REFERENCE_COLOR, markersize=8, label="ARCHS4 reference"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=GROUND_COLOR, markeredgecolor="white", markersize=6, label="Ground control"),
        Line2D([0], [0], marker="^", linestyle="none", markerfacecolor=FLIGHT_COLOR, markeredgecolor="white", markersize=6, label="Flight"),
        Line2D([0], [0], color="#7B858A", linewidth=1, label="Project centroid shift"),
    )
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    save_figure(fig, "figure_2_latent_mapping")


def plot_pathway_shifts(evidence: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH, 6.25), layout="constrained")
    for label, ax, tissue in zip("abcd", axes.flat, MAIN_TISSUES):
        terms = RETAINED_TERMS[tissue]
        subset = evidence.loc[evidence["tissue"].eq(tissue)].set_index("term").loc[
            list(terms)
        ].reset_index()
        points = project_effects(tissue, terms)
        positions = np.arange(len(subset))[::-1]
        for position, row in zip(positions, subset.itertuples(index=False)):
            local = points.loc[points["term"].eq(row.term), "project_effect"].to_numpy(dtype=float)
            ax.scatter(
                local,
                np.full(len(local), position),
                s=22,
                facecolor="white",
                edgecolor="#707A7F",
                linewidth=0.7,
                alpha=0.95,
                zorder=2,
            )
            role = str(row.evidence_role)
            ax.hlines(
                position,
                float(row.seed_effect_minimum),
                float(row.seed_effect_maximum),
                color=ROLE_COLORS[role],
                linewidth=1.4,
                zorder=3,
            )
            ax.scatter(
                float(row.seed_effect_median),
                position,
                marker=ROLE_MARKERS[role],
                s=48,
                color=ROLE_COLORS[role],
                edgecolor="white",
                linewidth=0.6,
                zorder=4,
            )
        ax.axvline(0, color="#3F494E", linewidth=0.8)
        ax.set_yticks(positions)
        ax.set_yticklabels(subset["display_label"], fontsize=7.2, color=TEXT_COLOR)
        project_count = int(subset["expimap_n_projects"].max())
        ax.set_title(
            f"{tissue.title()} ({project_count} projects)",
            loc="left",
            fontweight="bold",
            fontsize=8.5,
        )
        if label in "cd":
            ax.set_xlabel("Flight - ground pathway shift")
        clean_axis(ax)
        panel_label(ax, label)

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#707A7F",
            markersize=5,
            label="OSDR project",
        ),
        *[
            Line2D(
                [0],
                [0],
                marker=ROLE_MARKERS[role],
                color=ROLE_COLORS[role],
                linewidth=1.4,
                markersize=5,
                label=label_text,
            )
            for role, label_text in (
                ("aligned", "Literature aligned"),
                ("complementary", "Complementary"),
                ("context_sensitive", "Context sensitive"),
            )
        ],
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False, fontsize=6.7)
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.08, wspace=0.08, hspace=0.08)
    save_figure(fig, "figure_3_tissue_pathway_shifts")


def _derive_latent_orientation(
    scores: np.ndarray,
    obs: pd.DataFrame,
    target_effect: float,
) -> float:
    frame = pd.DataFrame(
        {
            "score": scores,
            "project": obs["id.accession"].astype(str).to_numpy(),
            "condition": obs["condition_inferred"].astype(str).to_numpy(),
        }
    )
    effects = []
    for _, project in frame.groupby("project", observed=True):
        means = project.groupby("condition")["score"].mean()
        if {"flight", "ground_control"}.issubset(means.index):
            effects.append(float(means["flight"] - means["ground_control"]))
    raw_effect = float(np.mean(effects))
    if raw_effect == 0 or target_effect == 0:
        return 1.0
    return 1.0 if np.sign(raw_effect) == np.sign(target_effect) else -1.0


def build_program_score_table(evidence: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for config in LATENT_CONFIGS:
        terms, query_scores, query_obs, query_ids = _load_latent_scores(
            config.query_path,
            query=True,
        )
        matches = np.flatnonzero(terms == config.representative_term)
        if len(matches) != 1:
            raise RuntimeError(f"Representative term missing for {config.tissue}")
        scores = query_scores[:, matches[0]]
        target = evidence.loc[
            evidence["tissue"].eq(config.tissue)
            & evidence["term"].eq(config.representative_term),
            "effect_seed2020",
        ]
        if len(target) != 1:
            raise RuntimeError(f"Primary pathway effect missing for {config.tissue}")
        orientation = _derive_latent_orientation(scores, query_obs, float(target.iloc[0]))
        oriented = scores * orientation
        frame = pd.DataFrame(
            {
                "tissue": config.tissue,
                "term": config.representative_term,
                "display_label": DISPLAY_LABELS[config.representative_term],
                "sample_id": query_ids,
                "project": query_obs["id.accession"].astype(str).to_numpy(),
                "condition": query_obs["condition_inferred"].astype(str).to_numpy(),
                "oriented_pathway_score": oriented,
                "latent_orientation": orientation,
            }
        )
        if config.tissue == "spleen":
            frame = frame.loc[~frame["project"].eq("OSD-288")].copy()
        frame["project_centered_pathway_score"] = frame["oriented_pathway_score"] - frame.groupby(
            "project", observed=True
        )["oriented_pathway_score"].transform("mean")
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    result.to_csv(
        SOURCE_DIR / "table_s33_representative_program_sample_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return result


def plot_program_score_distributions(scores: pd.DataFrame, evidence: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH, 5.25), layout="constrained")
    rng = np.random.default_rng(2026)
    for label, ax, config in zip("abcd", axes.flat, LATENT_CONFIGS):
        subset = scores.loc[scores["tissue"].eq(config.tissue)].copy()
        x_lookup = {"ground_control": 0.0, "flight": 1.0}
        for _, project in subset.groupby("project", observed=True):
            means = project.groupby("condition")["project_centered_pathway_score"].mean()
            if {"flight", "ground_control"}.issubset(means.index):
                ax.plot(
                    [0, 1],
                    [means["ground_control"], means["flight"]],
                    color="#9AA2A6",
                    linewidth=0.7,
                    alpha=0.65,
                    zorder=1,
                )
        for condition, color, marker in (
            ("ground_control", GROUND_COLOR, "o"),
            ("flight", FLIGHT_COLOR, "^"),
        ):
            local = subset.loc[subset["condition"].eq(condition)]
            x = x_lookup[condition] + rng.uniform(-0.085, 0.085, len(local))
            ax.scatter(
                x,
                local["project_centered_pathway_score"],
                s=15,
                marker=marker,
                color=color,
                edgecolor="white",
                linewidth=0.25,
                alpha=0.52,
                zorder=2,
            )
            values = local["project_centered_pathway_score"].to_numpy(dtype=float)
            quartiles = np.quantile(values, [0.25, 0.5, 0.75])
            ax.vlines(x_lookup[condition], quartiles[0], quartiles[2], color=TEXT_COLOR, linewidth=2.0, zorder=3)
            ax.scatter(x_lookup[condition], quartiles[1], marker="D", s=26, color=TEXT_COLOR, edgecolor="white", linewidth=0.4, zorder=4)
        row = evidence.loc[
            evidence["tissue"].eq(config.tissue)
            & evidence["term"].eq(config.representative_term)
        ].iloc[0]
        ax.set_title(
            f"{config.tissue.title()}: {DISPLAY_LABELS[config.representative_term]}",
            loc="left",
            fontweight="bold",
        )
        ax.text(
            0.98,
            0.96,
            f"project-balanced shift {float(row.effect_seed2020):+.2f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=6.8,
        )
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Ground", "Flight"])
        ax.set_ylabel("Project-centered pathway score")
        clean_axis(ax, keep_left=True)
        panel_label(ax, label)
    save_figure(fig, "figure_s9_program_score_distributions")


def build_member_gene_tables(evidence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_results = pd.read_csv(SOURCE_DIR / "table_s4_gene_level_results.tsv.gz", sep="\t")
    detail_rows: list[dict[str, object]] = []
    config_lookup = {config.tissue: config for config in LATENT_CONFIGS}

    for tissue in ("thymus", "skin", "liver"):
        config = config_lookup[tissue]
        query = ad.read_h5ad(config.query_path, backed="r")
        terms = np.asarray(query.uns["terms"], dtype=str)
        membership = np.asarray(query.varm["I"], dtype=float)
        gene_ids = query.var_names.astype(str)
        local_gene_results = gene_results.loc[gene_results["tissue"].eq(tissue)].set_index("gene_id")
        for term in RETAINED_TERMS[tissue]:
            term_index = np.flatnonzero(terms == term)
            if len(term_index) != 1:
                raise RuntimeError(f"Pathway membership missing for {tissue}: {term}")
            member_ids = gene_ids[membership[:, term_index[0]] > 0]
            pathway_effect = float(
                evidence.loc[
                    evidence["tissue"].eq(tissue) & evidence["term"].eq(term),
                    "effect_seed2020",
                ].iloc[0]
            )
            for gene_id in member_ids:
                if gene_id not in local_gene_results.index:
                    continue
                row = local_gene_results.loc[gene_id]
                gene_effect = float(row["study_balanced_log2cpm_flight_minus_ground"])
                detail_rows.append(
                    {
                        "tissue": tissue,
                        "term": term,
                        "display_label": DISPLAY_LABELS[term],
                        "gene_id": gene_id,
                        "gene_symbol": row["gene_symbol"],
                        "pathway_effect_seed2020": pathway_effect,
                        "project_balanced_gene_log2cpm_effect": gene_effect,
                        "gene_fdr": row["study_t_fdr"],
                        "same_direction_as_pathway_score": bool(
                            gene_effect != 0 and np.sign(gene_effect) == np.sign(pathway_effect)
                        ),
                        "oriented_decoder_weight": np.nan,
                        "matches_decoder_predicted_direction": pd.NA,
                        "evidence_source": "OSDR project-balanced gene expression",
                    }
                )
        query.file.close()

    reassessed = pd.read_csv(
        REASSESSMENT_DIR / "main_candidate_member_gene_effects.tsv", sep="\t"
    )
    reassessed = reassessed.loc[reassessed["seed"].eq(2020)].copy()
    retained_pairs = set(zip(evidence["tissue"], evidence["term"]))
    reassessed = reassessed.loc[
        [(tissue, term) in retained_pairs for tissue, term in zip(reassessed["tissue"], reassessed["term"])]
    ]
    evidence_effect = evidence.set_index(["tissue", "term"])["effect_seed2020"]
    for row in reassessed.itertuples(index=False):
        detail_rows.append(
            {
                "tissue": row.tissue,
                "term": row.term,
                "display_label": DISPLAY_LABELS[row.term],
                "gene_id": row.gene_id,
                "gene_symbol": row.gene_symbol,
                "pathway_effect_seed2020": float(evidence_effect.loc[(row.tissue, row.term)]),
                "project_balanced_gene_log2cpm_effect": float(row.project_balanced_log2cpm_effect),
                "gene_fdr": np.nan,
                "same_direction_as_pathway_score": bool(row.same_direction_as_pathway_score),
                "oriented_decoder_weight": float(row.oriented_decoder_weight),
                "matches_decoder_predicted_direction": bool(row.matches_decoder_predicted_direction),
                "evidence_source": "OSDR project-balanced gene expression and seed-2020 decoder",
            }
        )

    detail = pd.DataFrame(detail_rows)
    detail["absolute_gene_effect"] = detail["project_balanced_gene_log2cpm_effect"].abs()
    summary_rows = []
    decoder_summary = pd.read_csv(
        SOURCE_DIR / "table_s29_kidney_spleen_member_gene_support.tsv", sep="\t"
    ).set_index(["tissue", "term"])
    for (tissue, term), group in detail.groupby(["tissue", "term"], observed=True):
        concordant = group.loc[group["same_direction_as_pathway_score"]].sort_values(
            "absolute_gene_effect", ascending=False
        )
        top_genes = "; ".join(
            f"{row.gene_symbol} ({row.project_balanced_gene_log2cpm_effect:+.2f})"
            for row in concordant.head(5).itertuples(index=False)
        )
        decoder_fraction = np.nan
        decoder_minimum = np.nan
        if (tissue, term) in decoder_summary.index:
            decoder_row = decoder_summary.loc[(tissue, term)]
            decoder_fraction = float(
                decoder_row["decoder_abs_weight_direction_match_fraction_median"]
            )
            decoder_minimum = float(
                decoder_row["decoder_abs_weight_direction_match_fraction_minimum"]
            )
        finite_fdr = pd.to_numeric(group["gene_fdr"], errors="coerce")
        summary_rows.append(
            {
                "tissue": tissue,
                "term": term,
                "display_label": DISPLAY_LABELS[term],
                "pathway_effect_seed2020": float(group["pathway_effect_seed2020"].iloc[0]),
                "member_genes_in_hvg_model": int(len(group)),
                "member_gene_same_direction_fraction": float(
                    group["same_direction_as_pathway_score"].mean()
                ),
                "median_member_gene_log2cpm_effect": float(
                    group["project_balanced_gene_log2cpm_effect"].median()
                ),
                "member_gene_fdr_lt_005_fraction": (
                    float((finite_fdr < 0.05).mean()) if finite_fdr.notna().any() else np.nan
                ),
                "decoder_abs_weight_direction_match_fraction_median": decoder_fraction,
                "decoder_abs_weight_direction_match_fraction_minimum": decoder_minimum,
                "top_concordant_member_genes": top_genes,
            }
        )
    summary = pd.DataFrame(summary_rows)
    order = {tissue: index for index, tissue in enumerate((*MAIN_TISSUES, "kidney"))}
    term_order = {
        (tissue, term): index
        for tissue, terms in RETAINED_TERMS.items()
        for index, term in enumerate(terms)
    }
    summary["tissue_order"] = summary["tissue"].map(order)
    summary["term_order"] = [
        term_order[(tissue, term)] for tissue, term in zip(summary["tissue"], summary["term"])
    ]
    summary = summary.sort_values(["tissue_order", "term_order"]).drop(
        columns=["tissue_order", "term_order"]
    )
    detail = detail.merge(
        summary[["tissue", "term"]].assign(retained=True),
        on=["tissue", "term"],
        how="inner",
    ).drop(columns="retained")
    summary.to_csv(
        SOURCE_DIR / "table_s34_retained_pathway_member_gene_support.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    detail.to_csv(
        SOURCE_DIR / "table_s35_retained_pathway_member_gene_effects.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        na_rep="NA",
    )
    return summary, detail


def plot_evidence_and_gene_support(evidence: pd.DataFrame, gene_summary: pd.DataFrame) -> None:
    order = {tissue: index for index, tissue in enumerate((*MAIN_TISSUES, "kidney"))}
    term_order = {
        (tissue, term): index
        for tissue, terms in RETAINED_TERMS.items()
        for index, term in enumerate(terms)
    }
    frame = evidence.copy()
    frame["tissue_order"] = frame["tissue"].map(order)
    frame["term_order"] = [
        term_order[(tissue, term)] for tissue, term in zip(frame["tissue"], frame["term"])
    ]
    frame = frame.sort_values(["tissue_order", "term_order"]).drop(
        columns=["tissue_order", "term_order"]
    )
    frame = frame.merge(
        gene_summary[
            [
                "tissue",
                "term",
                "member_genes_in_hvg_model",
                "member_gene_same_direction_fraction",
            ]
        ],
        on=["tissue", "term"],
        how="left",
        validate="one_to_one",
    )
    matrix = frame[[column for column, _ in CHECK_COLUMNS]].fillna(False).astype(bool).to_numpy()
    labels = [
        f"{row.tissue.title()}: {row.display_label}"
        + (" [secondary]" if row.tissue == "kidney" else "")
        for row in frame.itertuples(index=False)
    ]
    y = np.arange(len(frame))

    fig = plt.figure(figsize=(FIGURE_WIDTH, 6.75), layout="constrained")
    grid = fig.add_gridspec(1, 2, width_ratios=(2.35, 1.0), wspace=0.05)
    matrix_ax = fig.add_subplot(grid[0, 0])
    gene_ax = fig.add_subplot(grid[0, 1], sharey=matrix_ax)
    matrix_ax.imshow(
        matrix,
        aspect="auto",
        cmap=LinearSegmentedColormap.from_list("support", ["#EFF1F2", "#2675A6"]),
        vmin=0,
        vmax=1,
    )
    matrix_ax.set_xticks(np.arange(len(CHECK_COLUMNS)))
    matrix_ax.set_xticklabels([label for _, label in CHECK_COLUMNS], rotation=35, ha="left")
    matrix_ax.xaxis.tick_top()
    matrix_ax.set_yticks(y)
    matrix_ax.set_yticklabels(labels, fontsize=6.9, color=TEXT_COLOR)
    matrix_ax.set_xticks(np.arange(-0.5, len(CHECK_COLUMNS), 1), minor=True)
    matrix_ax.set_yticks(np.arange(-0.5, len(frame), 1), minor=True)
    matrix_ax.grid(which="minor", color="white", linewidth=1.0)
    matrix_ax.tick_params(which="minor", bottom=False, left=False)
    matrix_ax.tick_params(axis="both", length=0)
    for row_index in range(len(frame)):
        for column_index in range(len(CHECK_COLUMNS)):
            supported = bool(matrix[row_index, column_index])
            matrix_ax.text(
                column_index,
                row_index,
                "+" if supported else "-",
                ha="center",
                va="center",
                fontsize=7.0,
                color="white" if supported else "#687176",
                fontweight="bold" if supported else "normal",
            )
    fdr_x = len(CHECK_COLUMNS) + 0.18
    matrix_ax.set_xlim(-0.5, fdr_x + 1.05)
    matrix_ax.text(
        fdr_x,
        -0.95,
        "GSEA FDR",
        ha="left",
        va="bottom",
        fontsize=7.0,
        fontweight="bold",
    )
    for index, row in enumerate(frame.itertuples(index=False)):
        matrix_ax.text(
            fdr_x,
            index,
            format_fdr(float(row.gsea_fdr)),
            ha="left",
            va="center",
            fontsize=6.8,
        )
    matrix_ax.set_title("Directional support", loc="left", fontsize=9.0, fontweight="bold", pad=34)
    panel_label(matrix_ax, "a", x=-0.45, y=1.05)

    for index, row in enumerate(frame.itertuples(index=False)):
        role = str(row.evidence_role)
        gene_ax.scatter(
            float(row.member_gene_same_direction_fraction),
            index,
            s=18 + 0.20 * np.sqrt(float(row.member_genes_in_hvg_model)) * 12,
            marker=ROLE_MARKERS[role],
            color=ROLE_COLORS[role],
            edgecolor="white",
            linewidth=0.55,
            zorder=3,
        )
    gene_ax.axvline(0.5, color="#6C757A", linestyle="--", linewidth=0.8)
    gene_ax.set_xlim(0, 1.02)
    gene_ax.set_xlabel("Member genes moving\nwith pathway direction")
    gene_ax.set_title("Gene-level support", loc="left", fontsize=9.0, fontweight="bold", pad=34)
    gene_ax.tick_params(axis="y", left=False, labelleft=False)
    clean_axis(gene_ax, keep_left=True)
    panel_label(gene_ax, "b", x=-0.08, y=1.05)

    for index in range(1, len(frame)):
        if frame.iloc[index]["tissue"] != frame.iloc[index - 1]["tissue"]:
            for ax in (matrix_ax, gene_ax):
                ax.axhline(index - 0.5, color="#586267", linewidth=1.0)
    handles = [
        Line2D(
            [0],
            [0],
            marker=ROLE_MARKERS[role],
            linestyle="none",
            markerfacecolor=ROLE_COLORS[role],
            markeredgecolor="white",
            markersize=5,
            label=label,
        )
        for role, label in (
            ("aligned", "Literature aligned"),
            ("complementary", "Complementary"),
            ("context_sensitive", "Context sensitive"),
        )
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    save_figure(fig, "figure_4_evidence_gene_support")


def plot_skin_protocol_context() -> None:
    context = pd.read_csv(
        SOURCE_DIR / "table_s8_skin_protocol_context_effects.tsv", sep="\t"
    )
    groups = (
        ("MHU-2 gravity and site", range(0, 4)),
        ("RR-5 recovery and RR-6 endpoint", range(4, 8)),
        ("RR-7 strain and duration", range(8, 12)),
    )
    short_labels = {
        0: "Dorsal\nmicrogravity",
        1: "Dorsal\nartificial 1 g",
        2: "Femoral\nmicrogravity",
        3: "Femoral\nartificial 1 g",
        4: "RR-5 dorsal\n30-day recovery",
        5: "RR-5 femoral\n30-day recovery",
        6: "RR-6 live return\nabout 30 days",
        7: "RR-6 terminal\nabout 60 days",
        8: "C3H/HeJ\n25 days",
        9: "C3H/HeJ\n75 days",
        10: "C57BL/6J\n25 days",
        11: "C57BL/6J\n75 days",
    }
    cmap = LinearSegmentedColormap.from_list(
        "condition_shift", [GROUND_COLOR, "#FAFAFA", FLIGHT_COLOR]
    )
    vmax = float(np.ceil(context["flight_minus_ground"].abs().max() * 2) / 2)
    fig, axes = plt.subplots(3, 1, figsize=(FIGURE_WIDTH, 6.35), layout="constrained")
    image = None
    for label, ax, (title, contrast_orders) in zip("abc", axes, groups):
        subset = context.loc[context["contrast_order"].isin(contrast_orders)].copy()
        matrix = subset.pivot(
            index="pathway_order", columns="contrast_order", values="flight_minus_ground"
        ).sort_index(axis=0).sort_index(axis=1)
        pathways = (
            subset[["pathway_order", "short_label"]]
            .drop_duplicates()
            .sort_values("pathway_order")["short_label"]
            .tolist()
        )
        metadata = (
            subset[["contrast_order", "n_flight", "n_ground_control"]]
            .drop_duplicates()
            .set_index("contrast_order")
        )
        image = ax.imshow(matrix.to_numpy(), cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        xlabels = []
        for contrast_order in matrix.columns:
            row = metadata.loc[contrast_order]
            xlabels.append(
                f"{short_labels[int(contrast_order)]}\nn={int(row.n_flight)} F/{int(row.n_ground_control)} G"
            )
        ax.set_xticks(np.arange(len(xlabels)))
        ax.set_xticklabels(xlabels, fontsize=6.7)
        ax.set_yticks(np.arange(len(pathways)))
        ax.set_yticklabels(pathways, fontsize=6.8)
        ax.set_title(title, loc="left", fontsize=8.3, fontweight="bold")
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = float(matrix.iloc[row_index, column_index])
                ax.text(
                    column_index,
                    row_index,
                    f"{value:+.2f}",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="white" if abs(value) > 0.58 * vmax else TEXT_COLOR,
                )
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel_label(ax, label, x=-0.23, y=1.04)
    if image is None:
        raise RuntimeError("Skin protocol context table is empty")
    colorbar = fig.colorbar(image, ax=axes, orientation="horizontal", fraction=0.035, pad=0.03, aspect=45)
    colorbar.set_label("Flight minus matched-ground pathway shift", fontsize=7.3)
    colorbar.ax.tick_params(labelsize=6.8)
    save_figure(fig, "figure_5_skin_protocol_context")


def plot_broad_pathway_screen() -> None:
    screen = pd.read_csv(SOURCE_DIR / "table_s9_systematic_pathway_screen.tsv", sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(FIGURE_WIDTH, 8.8), layout="constrained")
    for label, ax, tissue in zip("abcd", axes.flat, ("thymus", "skin", "liver", "soleus")):
        subset = screen.loc[
            screen["tissue"].eq(tissue) & screen["top_20_absolute_effect_active"].astype(bool)
        ].copy()
        subset = subset.sort_values("mean_accession_effect")
        y = np.arange(len(subset))
        positive = subset["mean_accession_effect"].ge(0).to_numpy()
        for position, row, is_positive in zip(y, subset.itertuples(index=False), positive):
            curated = bool(row.curated_for_main_figures)
            ax.scatter(
                float(row.mean_accession_effect),
                position,
                marker="^" if is_positive else "o",
                s=28 if curated else 18,
                color=FLIGHT_COLOR if is_positive else GROUND_COLOR,
                edgecolor=TEXT_COLOR if curated else "white",
                linewidth=0.8 if curated else 0.35,
                zorder=3,
            )
        ax.axvline(0, color="#3F494E", linewidth=0.75)
        labels = [
            textwrap.shorten(
                f"{row.screen_display_label}{' *' if bool(row.curated_for_main_figures) else ''}",
                width=43,
                placeholder="...",
            )
            for row in subset.itertuples(index=False)
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.3, color=TEXT_COLOR)
        active = int(screen.loc[screen["tissue"].eq(tissue), "active_latent_program"].sum())
        ax.set_title(f"{tissue.title()}: top 20 of {active}", loc="left", fontsize=8.2, fontweight="bold")
        if label in "cd":
            ax.set_xlabel("Flight - ground score")
        clean_axis(ax)
        panel_label(ax, label)
    handles = (
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=GROUND_COLOR, markeredgecolor="white", markersize=5, label="Lower in flight"),
        Line2D([0], [0], marker="^", linestyle="none", markerfacecolor=FLIGHT_COLOR, markeredgecolor="white", markersize=5, label="Higher in flight"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor=TEXT_COLOR, markersize=5, label="Individually reviewed (*)"),
    )
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=6.7)
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.08, wspace=0.10, hspace=0.08)
    save_figure(fig, "figure_s1_broad_pathway_screen")


def plot_original_robustness_matrix() -> None:
    frame = pd.read_csv(
        SOURCE_DIR / "table_s24_pathway_robustness_evidence.tsv", sep="\t"
    ).copy()
    order = {tissue: index for index, tissue in enumerate(("thymus", "skin", "liver", "soleus"))}
    frame["tissue_order"] = frame["tissue"].map(order)
    frame = frame.sort_values(["tissue_order", "short_label"]).drop(columns="tissue_order")
    matrix = frame[[column for column, _ in CHECK_COLUMNS]].fillna(False).astype(bool).to_numpy()
    labels = [f"{row.tissue.title()}: {row.short_label}" for row in frame.itertuples(index=False)]
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, 7.45), layout="constrained")
    ax.imshow(
        matrix,
        aspect="auto",
        cmap=LinearSegmentedColormap.from_list("support", ["#EFF1F2", "#2675A6"]),
        vmin=0,
        vmax=1,
    )
    ax.set_xticks(np.arange(len(CHECK_COLUMNS)))
    ax.set_xticklabels([label for _, label in CHECK_COLUMNS], rotation=30, ha="left")
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(frame)))
    ax.set_yticklabels(labels, fontsize=6.9, color=TEXT_COLOR)
    ax.set_xticks(np.arange(-0.5, len(CHECK_COLUMNS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(frame), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="both", length=0)
    for row_index in range(len(frame)):
        for column_index in range(len(CHECK_COLUMNS)):
            supported = bool(matrix[row_index, column_index])
            ax.text(
                column_index,
                row_index,
                "+" if supported else "-",
                ha="center",
                va="center",
                fontsize=7.0,
                color="white" if supported else "#687176",
                fontweight="bold" if supported else "normal",
            )
    for index in range(1, len(frame)):
        if frame.iloc[index]["tissue"] != frame.iloc[index - 1]["tissue"]:
            ax.axhline(index - 0.5, color="#586267", linewidth=1.0)
    ax.set_title("Original four-model pathway robustness checks", loc="left", fontsize=9.2, fontweight="bold", pad=30)
    save_figure(fig, "figure_s7_pathway_robustness_matrix")


def plot_process_summary() -> None:
    rows = (
        (
            "Thymus",
            "DNA repair, cytoskeletal control,\nand lymphoid-stromal\ninteraction lower",
            "Adaptive immune decline includes\nweaker cell coordination and\nniche-associated programs.",
            "Bulk tissue cannot separate\nthymocytes from stromal-cell\nabundance or state.",
        ),
        (
            "Skin",
            "Chromatin, DNA repair, Hedgehog,\nsphingolipid, and junction\nprograms lower",
            "Barrier stress is accompanied by\na broad maintenance and\ncommunication decrease.",
            "Gravity, site, recovery, strain,\nand duration modify the\npooled response.",
        ),
        (
            "Liver",
            "MHC class II presentation and\nT-cell receptor signaling\nlower",
            "Metabolic disruption is\naccompanied by lower adaptive\nimmune communication.",
            "Immune-cell abundance and\nsignaling state are not separable\nin bulk RNA-seq.",
        ),
        (
            "Spleen",
            "T-cell receptor, neutrophil\ndegranulation, and C-type\nlectin programs lower",
            "Adaptive suppression extends to\nlower innate effector and\npathogen-sensing programs.",
            "Transcriptomic programs are not\ndirect functional assays of\nimmune activity.",
        ),
    )
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, 4.65), layout="constrained")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    columns = ((0.27, "Observed program direction"), (0.57, "Complementary perspective"), (0.86, "Interpretation boundary"))
    for x, heading in columns:
        ax.text(x, 0.95, heading, ha="center", va="top", fontsize=7.5, fontweight="bold")
    for index, (tissue, observed, interpretation, boundary) in enumerate(rows):
        y_top = 0.84 - index * 0.205
        y_mid = y_top - 0.07
        ax.add_patch(
            FancyBboxPatch(
                (0.01, y_top - 0.13),
                0.11,
                0.15,
                boxstyle="round,pad=0.008,rounding_size=0.01",
                facecolor=TISSUE_COLORS[tissue.lower()],
                edgecolor="none",
            )
        )
        ax.text(0.065, y_top - 0.055, tissue, ha="center", va="center", fontsize=6.7, color="white", fontweight="bold")
        ax.text(0.27, y_mid, observed, ha="center", va="center", fontsize=6.8, linespacing=1.25)
        ax.text(0.57, y_mid, interpretation, ha="center", va="center", fontsize=6.8, linespacing=1.25)
        ax.text(0.86, y_mid, boundary, ha="center", va="center", fontsize=6.8, linespacing=1.25, color="#4E575B")
        if index < len(rows) - 1:
            ax.plot([0.01, 0.99], [y_top - 0.16, y_top - 0.16], color="#E0E3E4", linewidth=0.8)
    ax.plot([0.42, 0.42], [0.08, 0.91], color="#E0E3E4", linewidth=0.8)
    ax.plot([0.71, 0.71], [0.08, 0.91], color="#E0E3E4", linewidth=0.8)
    ax.text(
        0.5,
        0.02,
        "Vector synthesis of retained pathway directions; columns organize interpretation and do not assert causality.",
        ha="center",
        va="bottom",
        fontsize=6.8,
        color="#596267",
    )
    save_figure(fig, "asgsr_process_summary", output_dir=PROCESS_SUMMARY_DIR)


def plot_tissue_state_hypotheses() -> None:
    """Separate prior evidence, observed bulk-tissue scores, and new hypotheses."""
    rows = (
        (
            "A",
            "Thymus",
            "Involution with reduced\nthymocyte proliferation and\nadaptive activity [6-9]",
            (
                "DNA repair",
                "RHOA cytoskeletal cycle",
                "Lymphoid-stromal interactions",
            ),
            "Known involution may also involve lower repair, thymocyte motility, and stromal-niche coordination.",
        ),
        (
            "B",
            "Skin",
            "Barrier injury, inflammation,\nand impaired regeneration\n[12-16]",
            (
                "Chromatin-modifying enzymes",
                "DNA repair",
                "Hedgehog signaling",
                "Sphingolipid metabolism",
                "Cell-cell junction organization",
            ),
            "Barrier injury may include a broader maintenance deficit spanning regulation, repair, barrier lipids, and cell coordination.",
        ),
        (
            "C",
            "Liver",
            "Metabolic and xenobiotic\ndysregulation [17-23]",
            (
                "MHC class II antigen presentation",
                "T-cell receptor signaling",
            ),
            "Metabolic heterogeneity may coexist with lower adaptive immune communication.",
        ),
        (
            "D",
            "Spleen",
            "Lower T-cell abundance,\nactivation, and responsiveness\n[35-37]",
            (
                "T-cell receptor signaling",
                "Neutrophil degranulation",
                "C-type lectin receptor signaling",
            ),
            "Reduced adaptive activity may extend to lower innate sensing and effector-related transcription.",
        ),
    )

    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, 6.35))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5,
        0.985,
        "Prior literature, observed pathway shifts, and hypotheses to test",
        ha="center",
        va="top",
        fontsize=10.0,
        fontweight="bold",
    )
    ax.text(
        0.17,
        0.925,
        "Established phenotype\nfrom prior literature",
        ha="center",
        va="center",
        fontsize=6.9,
        fontweight="bold",
        linespacing=1.15,
    )
    ax.text(
        0.49,
        0.925,
        "Observed lower expiMap\nscores in flight",
        ha="center",
        va="center",
        fontsize=6.9,
        fontweight="bold",
        linespacing=1.15,
    )
    ax.text(
        0.855,
        0.925,
        "Complementary tissue-state\nhypothesis to test",
        ha="center",
        va="center",
        fontsize=6.9,
        fontweight="bold",
        linespacing=1.15,
    )

    def draw_thymus(y: float) -> None:
        color = TISSUE_COLORS["thymus"]
        for x, angle in ((0.085, 25), (0.125, -25)):
            ax.add_patch(
                Ellipse(
                    (x, y),
                    0.065,
                    0.095,
                    angle=angle,
                    facecolor="#E7E0F2",
                    edgecolor=color,
                    linewidth=1.0,
                )
            )
        for x, offset in ((0.055, 0.008), (0.09, -0.018), (0.13, 0.018), (0.155, -0.012)):
            ax.add_patch(
                Circle(
                    (x, y + offset),
                    0.012,
                    facecolor="#9CC8E8",
                    edgecolor=GROUND_COLOR,
                    linewidth=0.7,
                )
            )
        ax.plot([0.045, 0.16], [y - 0.052, y - 0.052], color=color, linewidth=1.2)

    def draw_skin(y: float) -> None:
        colors = ("#F5C7B4", "#EFA88F", "#D9826B")
        for index, color in enumerate(colors):
            bottom = y - 0.055 + index * 0.037
            ax.add_patch(
                Rectangle(
                    (0.045, bottom),
                    0.12,
                    0.034,
                    facecolor=color,
                    edgecolor="#A85B49",
                    linewidth=0.6,
                )
            )
            for x in np.linspace(0.06, 0.15, 4):
                ax.add_patch(Circle((x, bottom + 0.017), 0.006, facecolor="#8F6681", edgecolor="none"))
        for x in (0.075, 0.105, 0.135):
            ax.plot([x, x], [y - 0.047, y + 0.05], color=GROUND_COLOR, linewidth=0.8)

    def draw_liver(y: float) -> None:
        color = TISSUE_COLORS["liver"]
        ax.add_patch(
            Polygon(
                (
                    (0.045, y - 0.035),
                    (0.065, y - 0.06),
                    (0.13, y - 0.055),
                    (0.16, y - 0.01),
                    (0.145, y + 0.05),
                    (0.075, y + 0.06),
                ),
                closed=True,
                facecolor="#F2C69E",
                edgecolor="#B86F3B",
                linewidth=0.9,
            )
        )
        ax.add_patch(Circle((0.095, y), 0.022, facecolor="#B97867", edgecolor="#8A5147", linewidth=0.6))
        ax.add_patch(Circle((0.145, y + 0.035), 0.014, facecolor="#DCD2EA", edgecolor=color, linewidth=0.8))
        ax.add_patch(Circle((0.165, y + 0.012), 0.012, facecolor="#9CC8E8", edgecolor=GROUND_COLOR, linewidth=0.7))
        ax.plot([0.145, 0.16], [y + 0.03, y + 0.016], color="#596267", linewidth=0.8)

    def draw_spleen(y: float) -> None:
        color = TISSUE_COLORS["spleen"]
        ax.add_patch(Ellipse((0.105, y), 0.12, 0.115, angle=-12, facecolor="#E9D4DB", edgecolor=color, linewidth=1.0))
        ax.add_patch(Circle((0.075, y + 0.02), 0.014, facecolor="#9CC8E8", edgecolor=GROUND_COLOR, linewidth=0.7))
        ax.add_patch(Circle((0.11, y), 0.019, facecolor="#C7A4D8", edgecolor="#6B4C9A", linewidth=0.7))
        for angle in np.linspace(0, 2 * np.pi, 5, endpoint=False):
            ax.add_patch(
                Circle(
                    (0.11 + 0.009 * np.cos(angle), y + 0.009 * np.sin(angle)),
                    0.006,
                    facecolor="#7A4B8C",
                    edgecolor="none",
                )
            )
        ax.add_patch(Circle((0.14, y - 0.027), 0.014, facecolor="#A9CFAE", edgecolor="#477A55", linewidth=0.7))

    icon_drawers = {
        "Thymus": draw_thymus,
        "Skin": draw_skin,
        "Liver": draw_liver,
        "Spleen": draw_spleen,
    }
    row_centers = (0.805, 0.61, 0.415, 0.22)
    for index, ((letter, tissue, prior_literature, programs, hypothesis), y) in enumerate(
        zip(rows, row_centers, strict=True)
    ):
        color = TISSUE_COLORS[tissue.lower()]
        ax.text(0.02, y + 0.07, letter, ha="left", va="top", fontsize=8.5, fontweight="bold")
        ax.text(0.045, y + 0.07, tissue.upper(), ha="left", va="top", fontsize=8.2, fontweight="bold", color=color)
        icon_drawers[tissue](y - 0.018)
        ax.text(
            0.18,
            y - 0.005,
            prior_literature,
            ha="left",
            va="center",
            fontsize=5.55,
            linespacing=1.2,
            color="#4E575B",
        )

        spread = min(0.125, 0.034 * (len(programs) - 1))
        program_ys = np.linspace(y + spread / 2, y - spread / 2, len(programs))
        for program, program_y in zip(programs, program_ys, strict=True):
            ax.add_patch(
                FancyArrowPatch(
                    (0.355, program_y + 0.012),
                    (0.355, program_y - 0.012),
                    arrowstyle="-|>",
                    mutation_scale=7.0,
                    linewidth=1.0,
                    color=GROUND_COLOR,
                )
            )
            ax.text(0.372, program_y, program, ha="left", va="center", fontsize=5.85)

        ax.add_patch(
            FancyArrowPatch(
                (0.64, y),
                (0.715, y),
                arrowstyle="-|>",
                mutation_scale=9.0,
                linewidth=1.0,
                linestyle=(0, (2, 2)),
                color="#6B7478",
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (0.735, y - 0.062),
                0.245,
                0.124,
                boxstyle="round,pad=0.008,rounding_size=0.01",
                facecolor="#F7F8F8",
                edgecolor=color,
                linewidth=1.0,
            )
        )
        ax.text(
            0.857,
            y,
            textwrap.fill(hypothesis, width=36),
            ha="center",
            va="center",
            fontsize=5.85,
            linespacing=1.2,
        )
        if index < len(rows) - 1:
            ax.plot([0.02, 0.98], [y - 0.098, y - 0.098], color="#E0E3E4", linewidth=0.8)

    ax.text(
        0.5,
        0.035,
        "Prior literature, this study's observations, and new hypotheses are shown separately; dotted arrows denote "
        "inference, not causality.\nLower scores do not prove pathway inhibition, and cell composition may contribute.",
        ha="center",
        va="center",
        fontsize=5.9,
        color="#596267",
    )
    save_figure(fig, "figure_6_tissue_state_hypotheses")


def clean_obsolete_assets() -> None:
    for stem in OBSOLETE_FIGURE_STEMS:
        for suffix in ("png", "pdf"):
            path = FIGURE_DIR / f"{stem}.{suffix}"
            if path.exists():
                path.unlink()
    for path in (
        FIGURE_DIR / "source/figure_6_gpt_image_2_process_art.png",
        FIGURE_DIR / "source/figure_6_gpt_image_2_prompt.md",
    ):
        if path.exists():
            path.unlink()
    source_directory = FIGURE_DIR / "source"
    if source_directory.exists() and not any(source_directory.iterdir()):
        source_directory.rmdir()


def write_figure_manifest() -> pd.DataFrame:
    rows = []
    for png_path in sorted(FIGURE_DIR.glob("*.png")):
        pdf_path = png_path.with_suffix(".pdf")
        with Image.open(png_path) as image:
            width, height = image.size
            rgb = np.asarray(image.convert("RGB"))
        border = np.concatenate(
            (
                rgb[:2].reshape(-1, 3),
                rgb[-2:].reshape(-1, 3),
                rgb[:, :2].reshape(-1, 3),
                rgb[:, -2:].reshape(-1, 3),
            ),
            axis=0,
        )
        rows.append(
            {
                "figure": png_path.stem,
                "png_width_pixels": width,
                "png_height_pixels": height,
                "authored_width_inches_at_300dpi": width / 300,
                "authored_height_inches_at_300dpi": height / 300,
                "png_bytes": png_path.stat().st_size,
                "pdf_exists": pdf_path.exists(),
                "pdf_bytes": pdf_path.stat().st_size if pdf_path.exists() else 0,
                "nonwhite_border_pixel_fraction": float((border < 248).any(axis=1).mean()),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(SOURCE_DIR / "figure_build_manifest.tsv", sep="\t", index=False)
    return manifest


def validate_new_figures(manifest: pd.DataFrame) -> None:
    expected = {
        "figure_1_workflow_architecture",
        "figure_2_latent_mapping",
        "figure_3_tissue_pathway_shifts",
        "figure_4_evidence_gene_support",
        "figure_5_skin_protocol_context",
        "figure_6_tissue_state_hypotheses",
        "figure_s1_broad_pathway_screen",
        "figure_s7_pathway_robustness_matrix",
        "figure_s9_program_score_distributions",
    }
    available = set(manifest["figure"])
    missing = expected - available
    if missing:
        raise RuntimeError(f"Publication figures were not generated: {sorted(missing)}")
    generated = manifest.loc[manifest["figure"].isin(expected)]
    if not generated["pdf_exists"].all():
        raise RuntimeError("One or more publication figures lack a vector PDF copy")
    if not np.allclose(generated["authored_width_inches_at_300dpi"], FIGURE_WIDTH):
        widths = generated.set_index("figure")["authored_width_inches_at_300dpi"].to_dict()
        raise RuntimeError(f"Unexpected authored figure widths: {widths}")
    if (generated["nonwhite_border_pixel_fraction"] > 0.001).any():
        failures = generated.loc[
            generated["nonwhite_border_pixel_fraction"] > 0.001,
            ["figure", "nonwhite_border_pixel_fraction"],
        ]
        raise RuntimeError(f"Figure content approaches the image boundary:\n{failures}")


def run() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_style()
    evidence = load_retained_evidence()
    scope = pd.read_csv(SOURCE_DIR / "table_s25_revised_model_scope.tsv", sep="\t")

    plot_workflow(scope)
    coordinates, qc = build_latent_mapping_tables()
    plot_latent_mapping(coordinates, qc)
    plot_pathway_shifts(evidence)
    scores = build_program_score_table(evidence)
    plot_program_score_distributions(scores, evidence)
    gene_summary, _ = build_member_gene_tables(evidence)
    plot_evidence_and_gene_support(evidence, gene_summary)
    plot_skin_protocol_context()
    plot_tissue_state_hypotheses()
    plot_broad_pathway_screen()
    plot_original_robustness_matrix()
    plot_process_summary()
    clean_obsolete_assets()
    manifest = write_figure_manifest()
    validate_new_figures(manifest)
    print(
        manifest.loc[
            manifest["figure"].str.match(r"figure_(?:[1-6]|s1|s7|s9)_"),
            ["figure", "png_width_pixels", "png_height_pixels", "pdf_exists"],
        ].to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    run()
