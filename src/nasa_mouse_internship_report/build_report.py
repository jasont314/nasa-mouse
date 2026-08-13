"""Build figures, source manifests, HTML, and PDF for the internship report."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[2]
PAPER_DIR = ROOT / "paper" / "slstp_internship_report"
FIGURE_DIR = PAPER_DIR / "figures"
SOURCE_DIR = PAPER_DIR / "source_data"

COPIED_FIGURES = {
    ROOT / "paper/asgsr_expimap_hvg/figures/figure_3_tissue_pathway_shifts.png": (
        FIGURE_DIR / "figure_6_expimap_pathway_shifts.png"
    ),
    ROOT / "paper/synthetic_guided_spaceflight/figures/figure_2b_locked_real_vs_synthetic_pca.png": (
        FIGURE_DIR / "figure_8_real_synthetic_pca.png"
    ),
    ROOT / "paper/synthetic_guided_spaceflight/figures/figure_3_thymus_biology.png": (
        FIGURE_DIR / "figure_10_synthetic_thymus_biology.png"
    ),
    ROOT / "paper/asgsr_expimap_hvg/figures/figure_4_evidence_gene_support.png": (
        FIGURE_DIR / "figure_s1_expimap_robustness.png"
    ),
    ROOT / "paper/synthetic_guided_spaceflight/figures/figure_2a_archs4_denoising_trajectory.png": (
        FIGURE_DIR / "figure_s2_ddim_trajectory.png"
    ),
}


def _required(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix in ("png", "pdf"):
        metadata = None
        if suffix == "pdf":
            metadata = {
                "Creator": "nasa-mouse-spaceflight",
                "CreationDate": None,
                "ModDate": None,
            }
        fig.savefig(
            FIGURE_DIR / f"{stem}.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
            metadata=metadata,
        )
    plt.close(fig)


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def figure_2_glare_batch_effects() -> None:
    source = _required(
        ROOT
        / "paper/slstp_internship_report/source_data/glare/"
        / "skeletal_muscle_aggregate_vs_mober_umap_by_accession.png"
    )
    image = plt.imread(source).copy()
    image[:176] = 1
    cropped = image[125:965]

    fig = plt.figure(figsize=(12.2, 4.35))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.80])
    ax.imshow(cropped)
    ax.axis("off")
    fig.text(0.19, 0.86, "A  Before correction", ha="center", fontsize=10.5, weight="bold")
    fig.text(0.49, 0.86, "B  After MOBER", ha="center", fontsize=10.5, weight="bold")
    fig.suptitle(
        "MOBER reduced study separation, but FLT and GC remained mixed",
        x=0.02,
        ha="left",
        fontsize=15,
        weight="bold",
        color="#213f52",
    )
    _save_figure(fig, "figure_2_glare_batch_effects")


def _diagram_arrow(ax, start, end, *, color="#68767d", width=1.4, style="-") -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=width,
            linestyle=style,
            color=color,
            shrinkA=2,
            shrinkB=2,
        )
    )


def figure_3_expimap_architecture() -> None:
    fig, ax = plt.subplots(figsize=(12.4, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    input_values = np.array([0.42, 0.76, 0.31, 0.63, 0.88, 0.54, 0.24, 0.69])
    recon_values = np.array([0.45, 0.72, 0.35, 0.60, 0.82, 0.57, 0.29, 0.65])
    bar_x = np.linspace(0.025, 0.108, len(input_values))
    for x, value in zip(bar_x, input_values):
        ax.add_patch(Rectangle((x, 0.23), 0.008, value * 0.52, facecolor="#3579a6", edgecolor="none"))
    ax.plot([0.018, 0.122], [0.23, 0.23], color="#60727b", linewidth=0.8)
    ax.text(0.07, 0.12, "Gene-expression\nprofile", ha="center", va="center", fontsize=8.5, weight="bold")
    ax.text(0.07, 0.06, "about 2,000 HVGs", ha="center", va="center", fontsize=7.2, color="#52616a")

    encoder = Polygon(
        [[0.16, 0.17], [0.16, 0.88], [0.31, 0.69], [0.31, 0.36]],
        closed=True,
        facecolor="#dce9ee",
        edgecolor="#50616a",
        linewidth=1.2,
    )
    ax.add_patch(encoder)
    ax.text(0.225, 0.57, "Encoder", ha="center", va="center", fontsize=10, weight="bold")
    ax.text(0.225, 0.46, "3 x 300", ha="center", va="center", fontsize=7.5, color="#52616a")
    _diagram_arrow(ax, (0.122, 0.53), (0.155, 0.53))

    latent_y = np.linspace(0.25, 0.78, 5)
    latent_labels = ["DNA repair", "TCR signaling", "Hedgehog", "RHOA cycle", "Other programs"]
    latent_colors = ["#5d7eb5", "#238985", "#cc7a3b", "#7c67a5", "#8f9da4"]
    for y, label, color in zip(latent_y[::-1], latent_labels, latent_colors):
        ax.add_patch(Circle((0.395, y), 0.019, facecolor=color, edgecolor="white", linewidth=0.8, zorder=4))
        ax.text(0.425, y, label, va="center", fontsize=7.6, color="#2c3b42")
    ax.text(0.395, 0.91, "Reactome latent space", ha="center", va="center", fontsize=10, weight="bold")
    ax.text(0.395, 0.855, "319 to 387 named programs per tissue", ha="center", va="center", fontsize=7.5, color="#52616a")
    _diagram_arrow(ax, (0.31, 0.53), (0.368, 0.53))

    gene_y = np.linspace(0.23, 0.81, 8)
    connection_map = {
        0: [0, 1, 4],
        1: [1, 3, 5],
        2: [2, 3, 6],
        3: [0, 5, 7],
        4: [2, 4, 6, 7],
    }
    for latent_index, targets in connection_map.items():
        source_y = latent_y[::-1][latent_index]
        for target in targets:
            ax.plot([0.535, 0.675], [source_y, gene_y[target]], color=latent_colors[latent_index], alpha=0.45, linewidth=0.85)
    for y in gene_y:
        ax.add_patch(Circle((0.68, y), 0.012, facecolor="#edf1f2", edgecolor="#64747c", linewidth=0.7, zorder=4))
    ax.text(0.61, 0.91, "Masked decoder", ha="center", va="center", fontsize=10, weight="bold")
    ax.text(0.61, 0.855, "Reactome defines preferred program-gene links", ha="center", va="center", fontsize=7.5, color="#52616a")

    decoder = Polygon(
        [[0.72, 0.36], [0.72, 0.69], [0.82, 0.88], [0.82, 0.17]],
        closed=True,
        facecolor="#eeeaf5",
        edgecolor="#625a75",
        linewidth=1.2,
    )
    ax.add_patch(decoder)
    ax.text(0.77, 0.57, "Decoder", ha="center", va="center", fontsize=9.5, weight="bold")
    _diagram_arrow(ax, (0.692, 0.53), (0.715, 0.53))

    output_x = np.linspace(0.87, 0.953, len(recon_values))
    for x, value in zip(output_x, recon_values):
        ax.add_patch(Rectangle((x, 0.23), 0.008, value * 0.52, facecolor="#8a6ea8", edgecolor="none"))
    ax.plot([0.863, 0.967], [0.23, 0.23], color="#60727b", linewidth=0.8)
    ax.text(0.915, 0.12, "Reconstructed\ngene profile", ha="center", va="center", fontsize=8.5, weight="bold")
    _diagram_arrow(ax, (0.825, 0.53), (0.858, 0.53))

    score_box = FancyBboxPatch(
        (0.33, 0.015),
        0.18,
        0.11,
        boxstyle="round,pad=0.009,rounding_size=0.014",
        facecolor="#e6f1ef",
        edgecolor="#287d7a",
        linewidth=1.0,
    )
    ax.add_patch(score_box)
    ax.text(0.42, 0.069, "Named program scores\nused for FLT vs GC", ha="center", va="center", fontsize=7.5, weight="bold")
    _diagram_arrow(ax, (0.395, 0.22), (0.402, 0.13), color="#287d7a", width=1.1)

    fig.suptitle("expiMap maps each sample into a constrained pathway space", x=0.02, ha="left", fontsize=15, weight="bold", color="#213f52", y=0.98)
    fig.subplots_adjust(top=0.88, bottom=0.06, left=0.03, right=0.98)
    _save_figure(fig, "figure_3_expimap_architecture")


def _expimap_heatmap_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    retained = pd.read_csv(
        _required(ROOT / "paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv"),
        sep="\t",
    )
    retained = retained.loc[
        retained["tissue"].isin(["thymus", "skin", "spleen", "kidney"]),
        ["tissue", "term", "display_label"],
    ].copy()

    primary = pd.read_csv(
        _required(ROOT / "paper/asgsr_expimap_hvg/source_data/table_s5_accession_pathway_effects.tsv.gz"),
        sep="\t",
    ).rename(columns={"id.accession": "accession", "flight_minus_ground": "effect"})
    primary = primary.loc[
        primary["tissue"].isin(["thymus", "skin"]),
        ["tissue", "accession", "term", "effect"],
    ]

    reassessed = pd.read_csv(
        _required(SOURCE_DIR / "expimap_kidney_spleen_seed_accession_effects.tsv.gz"),
        sep="\t",
    )
    reassessed = reassessed.loc[
        (reassessed["seed"] == 2020)
        & reassessed["tissue"].isin(["spleen", "kidney"])
        & ~((reassessed["tissue"] == "spleen") & (reassessed["accession"] == "OSD-288")),
        ["tissue", "accession", "term", "effect"],
    ]
    effects = pd.concat([primary, reassessed], ignore_index=True)
    effects = effects.merge(retained, on=["tissue", "term"], how="inner", validate="many_to_one")
    return retained, effects


def figure_4_expimap_pathway_heatmap() -> None:
    retained, effects = _expimap_heatmap_data()
    tissues = ["thymus", "skin", "spleen", "kidney"]
    panel_letters = ["A", "B", "C", "D"]
    short_labels = {
        "Chromatin-modifying enzymes": "Chromatin modifiers",
        "Cell-cell junction organization": "Cell-cell junctions",
        "Neutrophil degranulation program": "Neutrophil degranulation",
        "C-type lectin receptor signaling": "C-type lectin signaling",
    }
    fig = plt.figure(figsize=(12.4, 7.6))
    outer = fig.add_gridspec(
        2,
        2,
        left=0.025,
        right=0.985,
        top=0.84,
        bottom=0.08,
        hspace=0.52,
        wspace=0.16,
    )

    for index, (tissue, letter) in enumerate(zip(tissues, panel_letters)):
        row, column = divmod(index, 2)
        panel = outer[row, column].subgridspec(
            1,
            3,
            width_ratios=[1.65, 3.55, 0.13],
            wspace=0.06,
        )
        label_ax = fig.add_subplot(panel[0, 0])
        ax = fig.add_subplot(panel[0, 1])
        colorbar_ax = fig.add_subplot(panel[0, 2])
        term_order = retained.loc[retained["tissue"] == tissue, "display_label"].tolist()
        tissue_effects = effects.loc[effects["tissue"] == tissue].copy()
        accession_order = sorted(tissue_effects["accession"].unique())
        matrix = tissue_effects.pivot_table(
            index="display_label",
            columns="accession",
            values="effect",
            aggfunc="mean",
        ).reindex(index=term_order, columns=accession_order)
        values = matrix.to_numpy(dtype=float)
        limit = float(np.nanmax(np.abs(values)))
        image_handle = ax.imshow(values, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")

        label_ax.set_xlim(0, 1)
        label_ax.set_ylim(len(term_order) - 0.5, -0.5)
        label_ax.axis("off")
        label_ax.text(
            0.0,
            1.12,
            f"{letter}  {tissue.replace('_', ' ').title()}",
            transform=label_ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=10.5,
            weight="bold",
            color="#213f52",
        )
        for row_index, label in enumerate(term_order):
            display = short_labels.get(label, label)
            label_ax.text(
                0.98,
                row_index,
                textwrap.fill(display, width=25),
                ha="right",
                va="center",
                fontsize=7.6,
                color="#24333a",
                linespacing=1.05,
            )

        for row in range(values.shape[0]):
            for column in range(values.shape[1]):
                value = values[row, column]
                if np.isnan(value):
                    continue
                color = "white" if abs(value) > limit * 0.58 else "#24333a"
                ax.text(column, row, f"{value:+.2f}", ha="center", va="center", fontsize=6.8, color=color)
        ax.set_xticks(np.arange(len(accession_order)), accession_order, rotation=35, ha="right", fontsize=7.4)
        ax.set_yticks([])
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        colorbar = fig.colorbar(image_handle, cax=colorbar_ax)
        colorbar.ax.tick_params(labelsize=6.6, length=2)
        colorbar.set_label("FLT - GC", fontsize=7.0)

    fig.suptitle("Pathway shifts across OSDR studies", x=0.02, ha="left", fontsize=15, weight="bold", color="#213f52")
    fig.text(0.02, 0.945, "Blue cells are lower in flight; red cells are higher. Each tissue has its own symmetric color scale.", fontsize=9, color="#52616a")
    _save_figure(fig, "figure_4_expimap_pathway_heatmap")


def figure_5_expimap_annotation_workflow() -> None:
    fig = plt.figure(figsize=(12.2, 5.4))
    grid = fig.add_gridspec(2, 5, height_ratios=[1.0, 0.38], hspace=0.22, wspace=0.42)
    axes = [fig.add_subplot(grid[0, index]) for index in range(5)]

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    ax = axes[0]
    score_matrix = np.array(
        [
            [-0.6, -0.4, 0.2, 0.0, -0.3, -0.7],
            [0.3, 0.1, -0.2, -0.5, -0.4, -0.1],
            [-0.1, -0.4, -0.7, -0.5, 0.0, 0.2],
            [0.5, 0.2, 0.1, -0.2, -0.4, -0.5],
        ]
    )
    ax.imshow(score_matrix, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_title("1  Program scores", fontsize=9.2, weight="bold", color="#213f52")
    ax.text(0.5, -0.14, "samples", transform=ax.transAxes, ha="center", fontsize=7.5, color="#52616a")
    ax.text(-0.10, 0.5, "Reactome\nprograms", transform=ax.transAxes, ha="right", va="center", fontsize=7.5, color="#52616a")

    ax = axes[1]
    y = np.arange(5)
    values = np.array([-0.74, -0.42, -0.20, -0.61, -0.31])
    ax.axvline(0, color="#68767d", linewidth=0.9)
    ax.scatter(values, y, color="#2e78a4", s=35, zorder=3)
    ax.hlines(y, 0, values, color="#aebcc3", linewidth=1)
    ax.set_xlim(-0.9, 0.35)
    ax.set_ylim(-0.7, 4.7)
    ax.set_title("2  Study effects", fontsize=9.2, weight="bold", color="#213f52")
    ax.text(0.5, -0.14, "FLT - GC within each study", transform=ax.transAxes, ha="center", fontsize=7.5, color="#52616a")

    ax = axes[2]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    check_labels = [
        "ssGSEA direction",
        "GSEA direction",
        "held-out projects",
        "repeat training",
        "cell-mixture check",
    ]
    check_colors = ["#2e78a4", "#2e78a4", "#238985", "#8065a5", "#cc7a3b"]
    for y, label, color in zip(np.linspace(0.82, 0.18, 5), check_labels, check_colors):
        ax.add_patch(
            FancyBboxPatch(
                (0.07, y - 0.055),
                0.86,
                0.11,
                boxstyle="round,pad=0.008,rounding_size=0.02",
                facecolor="#f4f6f7",
                edgecolor="#a6b2b8",
                linewidth=0.8,
            )
        )
        ax.add_patch(Circle((0.14, y), 0.018, facecolor=color, edgecolor="none"))
        ax.text(0.20, y, label, ha="left", va="center", fontsize=6.7, color="#34444c")
    ax.set_title("3  Agreement checks", fontsize=9.2, weight="bold", color="#213f52")

    ax = axes[3]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    gene_steps = [
        (0.73, "FLT - GC for each gene", "#e7eff3", "#2e78a4"),
        (0.48, "Compare directions", "#e7f1eb", "#238985"),
        (0.23, "Review supporting genes", "#f4e9df", "#cc7a3b"),
    ]
    for y, label, fill, edge in gene_steps:
        ax.add_patch(
            FancyBboxPatch(
                (0.07, y - 0.08),
                0.86,
                0.16,
                boxstyle="round,pad=0.008,rounding_size=0.02",
                facecolor=fill,
                edgecolor=edge,
                linewidth=0.9,
            )
        )
        ax.text(0.50, y, label, ha="center", va="center", fontsize=6.8, color="#34444c", weight="bold")
    ax.set_title("4  Member-gene review", fontsize=9.2, weight="bold", color="#213f52")
    ax.text(0.5, -0.14, "support or contradict the pathway", transform=ax.transAxes, ha="center", fontsize=7.2, color="#52616a")

    ax = axes[4]
    for index, (x, y, color) in enumerate(
        [(0.10, 0.18, "#e7eff3"), (0.23, 0.30, "#e7f1eb"), (0.36, 0.42, "#f4e9df")]
    ):
        page = Rectangle((x, y), 0.48, 0.42, transform=ax.transAxes, facecolor=color, edgecolor="#68767d", linewidth=0.9)
        ax.add_patch(page)
        for line in range(3):
            ax.plot([x + 0.07, x + 0.40], [y + 0.30 - line * 0.085] * 2, transform=ax.transAxes, color="#7d8a90", linewidth=0.8)
    ax.set_title("5  Literature review", fontsize=9.2, weight="bold", color="#213f52")
    ax.text(0.5, -0.14, "compare with prior findings", transform=ax.transAxes, ha="center", fontsize=7.5, color="#52616a")

    for left_ax, right_ax in zip(axes[:-1], axes[1:]):
        start = left_ax.get_position().bounds
        end = right_ax.get_position().bounds
        arrow = FancyArrowPatch(
            (start[0] + start[2] + 0.005, start[1] + start[3] * 0.52),
            (end[0] - 0.005, end[1] + end[3] * 0.52),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=9,
            color="#849198",
            linewidth=1.1,
        )
        fig.add_artist(arrow)

    ax = fig.add_subplot(grid[1, :])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.20), 0.44, 0.58, boxstyle="round,pad=0.01,rounding_size=0.02", facecolor="#f0f3f4", edgecolor="#8d9aa0", linewidth=1.0))
    ax.text(0.045, 0.59, "Finding status", fontsize=8.8, weight="bold", color="#213f52")
    ax.text(0.045, 0.36, "supported in several analyses   |   mixed support   |   exploratory", fontsize=8.0, color="#3c4a51")
    ax.add_patch(FancyBboxPatch((0.50, 0.20), 0.48, 0.58, boxstyle="round,pad=0.01,rounding_size=0.02", facecolor="#f7f8f8", edgecolor="#8d9aa0", linewidth=1.0))
    ax.text(0.525, 0.59, "Literature relationship", fontsize=8.8, weight="bold", color="#213f52")
    literature_labels = [("aligned", "#238985"), ("complementary", "#3579a6"), ("context-sensitive", "#cc7a3b")]
    x = 0.525
    for label, color in literature_labels:
        ax.add_patch(Circle((x, 0.36), 0.012, facecolor=color, edgecolor="none"))
        ax.text(x + 0.019, 0.36, label, va="center", fontsize=8.0, color="#3c4a51")
        x += 0.145 if label != "complementary" else 0.18

    fig.suptitle("From pathway score to biological interpretation", x=0.02, ha="left", fontsize=15, weight="bold", color="#213f52", y=0.995)
    fig.text(0.02, 0.935, "Workflow schematic; measured results are shown in Figures 3 and 5.", fontsize=9, color="#52616a")
    fig.subplots_adjust(top=0.82, bottom=0.06, left=0.055, right=0.98)
    _save_figure(fig, "figure_5_expimap_annotation_workflow")


def figure_7_generator_validation() -> None:
    source_root = ROOT / "paper/synthetic_guided_spaceflight/source_data"
    arch = pd.read_csv(
        _required(source_root / "table_s1_archs4_ddim_metrics.tsv"),
        sep="\t",
    ).set_index("metric")["value"]
    locked = pd.read_csv(
        _required(source_root / "table_s2_locked_ddim_repeats.tsv"),
        sep="\t",
    )
    model_screen = pd.read_csv(
        _required(source_root / "table_4_generator_model_selection.tsv"),
        sep="\t",
    ).set_index("model")

    teal = "#238985"
    coral = "#d76552"
    gray = "#8d969b"
    gold = "#d99b28"
    light = "#edf1f2"
    fig = plt.figure(figsize=(12.2, 7.3))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.05, 0.82],
        width_ratios=[0.86, 1.14],
        hspace=0.54,
        wspace=0.34,
    )

    ax = fig.add_subplot(grid[0, 0])
    values = [
        arch["Real train to real test tissue BA"],
        arch["Synthetic train to real test tissue BA"],
    ]
    bars = ax.bar(["Real", "Synthetic"], values, color=[gray, teal], width=0.62)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Balanced accuracy")
    ax.set_title("A  ARCHS4 tissue recovery", loc="left")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.3f}", ha="center", fontsize=8.5)
    ax.grid(axis="y", color="#dce3e6", linewidth=0.7)

    ax = fig.add_subplot(grid[0, 1])
    metric_columns = [
        ("Corr", "correlation"),
        ("Precision", "precision"),
        ("Recall", "recall"),
        ("F1", "f1"),
        ("AA", "adversarial_accuracy"),
        ("FD ratio", "frechet_ratio"),
    ]
    labels = [item[0] for item in metric_columns]
    x = np.arange(len(labels))
    width = 0.36
    wgan_values = [
        model_screen.loc["Study-conditioned WGAN-GP", column]
        for _, column in metric_columns
    ]
    ddim_values = [
        model_screen.loc["Factorized DDIM", column]
        for _, column in metric_columns
    ]
    ax.add_patch(
        Rectangle(
            (3.55, 0.40),
            0.90,
            0.20,
            facecolor=light,
            edgecolor="none",
            zorder=0,
        )
    )
    ax.bar(x - width / 2, wgan_values, width, label="WGAN", color=coral, zorder=2)
    ax.bar(x + width / 2, ddim_values, width, label="DDIM", color=teal, zorder=2)
    ax.set_xticks(x, labels, rotation=25, ha="right", fontsize=7.5)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Metric value")
    ax.set_title("B  Generator comparison", loc="left", pad=24)
    ax.legend(
        frameon=False,
        fontsize=7.5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        borderaxespad=0,
        columnspacing=1.4,
        handlelength=1.5,
    )
    ax.grid(axis="y", color="#dce3e6", linewidth=0.7)

    ax = fig.add_subplot(grid[1, :])
    gate_labels = ["Corr", "Precision", "Recall", "F1", "AA", "FD", "FLT/GC"]
    gate_passes = [
        (locked["correlation"] >= locked["correlation_minimum"]).mean(),
        (locked["precision"] >= 0.95).mean(),
        (locked["recall"] >= 0.85).mean(),
        (locked["f1"] >= 0.90).mean(),
        locked["adversarial_accuracy"].between(0.40, 0.60).mean(),
        (locked["frechet_ratio"] <= 1.0).mean(),
        locked["condition_effect_pass"].astype(bool).mean(),
    ]
    y = np.arange(len(gate_labels))
    bars = ax.barh(y, gate_passes, color=[teal] * 6 + [gold])
    ax.set_yticks(y, gate_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Fraction of four runs meeting target")
    ax.set_title("C  DDIM consistency across runs", loc="left")
    for bar, value in zip(bars, gate_passes):
        ax.text(value + 0.025, bar.get_y() + bar.get_height() / 2, f"{int(value * 4)}/4", va="center", fontsize=7.5)
    ax.grid(axis="x", color="#dce3e6", linewidth=0.7)

    fig.suptitle(
        "DDIM best matched the real expression distribution",
        x=0.02,
        ha="left",
        fontsize=15,
        weight="bold",
        color="#213f52",
    )
    fig.subplots_adjust(top=0.81, bottom=0.08, left=0.07, right=0.98)
    _save_figure(fig, "figure_7_generator_validation")


def figure_9_synthetic_tissue_evidence() -> None:
    source_dir = ROOT / "paper" / "synthetic_guided_spaceflight" / "source_data"
    utility = pd.read_csv(
        _required(source_dir / "table_s18_matched_all_gene_utility.tsv"),
        sep="\t",
    )
    candidates = pd.read_csv(
        _required(source_dir / "table_s19_matched_all_gene_candidates.tsv"),
        sep="\t",
    )
    comparison = pd.read_csv(
        _required(source_dir / "table_s21_matched_consensus_comparison.tsv"),
        sep="\t",
    ).iloc[0]

    dark = "#213f52"
    teal = "#238985"
    coral = "#d86652"
    blue = "#3579a6"
    gray = "#8f9da4"

    fig = plt.figure(figsize=(12.4, 8.2))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.9, 0.95],
        width_ratios=[1.15, 0.85],
        hspace=0.48,
        wspace=0.28,
    )

    ax = fig.add_subplot(grid[0, :])
    performance = utility.loc[
        utility["arm"].eq("real_plus_generated")
    ].sort_values("mean_delta_balanced_accuracy", ascending=False, kind="stable")
    y = np.arange(len(performance))
    real_scores = performance["real_mean_balanced_accuracy"].to_numpy(float)
    synthetic_scores = performance["arm_mean_balanced_accuracy"].to_numpy(float)
    improved = synthetic_scores >= real_scores

    for yi, real_score, synthetic_score, is_improved in zip(
        y,
        real_scores,
        synthetic_scores,
        improved,
    ):
        ax.plot(
            [real_score, synthetic_score],
            [yi, yi],
            color=teal if is_improved else coral,
            alpha=0.65,
            linewidth=1.4,
            zorder=1,
        )
    ax.scatter(
        real_scores,
        y,
        s=28,
        facecolor="white",
        edgecolor=dark,
        linewidth=0.9,
        label="Real only",
        zorder=3,
    )
    ax.scatter(
        synthetic_scores[improved],
        y[improved],
        s=30,
        color=teal,
        label="Real + synthetic, higher or tied",
        zorder=4,
    )
    ax.scatter(
        synthetic_scores[~improved],
        y[~improved],
        s=30,
        color=coral,
        label="Real + synthetic, lower",
        zorder=4,
    )

    aliases = {
        "edl": "EDL",
        "skeletal_muscle": "Skeletal muscle",
        "tibialis_anterior": "Tibialis anterior",
        "white_adipose_tissue": "White adipose tissue",
        "brown_adipose_tissue": "Brown adipose tissue",
        "adrenal_gland": "Adrenal gland",
        "mammary_gland": "Mammary gland",
        "optic_nerve": "Optic nerve",
    }
    labels = []
    for row in performance.itertuples(index=False):
        label = aliases.get(row.tissue, str(row.tissue).replace("_", " ").title())
        if row.scope == "muscle_group":
            label = f"{label} (muscle group)"
        labels.append(label)
    ax.set_yticks(y, labels, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlim(0.42, 1.02)
    ax.set_xlabel("FLT/GC balanced accuracy on real test samples")
    ax.set_title("A  Classifier performance across tissues and muscle groups", loc="left", pad=12)
    ax.grid(axis="x", color="#dce3e6", linewidth=0.7)
    ax.legend(
        frameon=False,
        fontsize=7.2,
        ncol=3,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.005),
        handletextpad=0.4,
        columnspacing=1.2,
    )

    ax = fig.add_subplot(grid[1, 0])
    candidate_status = (
        candidates.groupby(["analysis_scope", "tissue", "gene"], observed=True)[
            "matched_status"
        ]
        .agg(
            lambda values: (
                "promoted" if "promoted" in set(values) else "shared_importance"
            )
        )
        .rename("matched_status")
        .reset_index()
    )
    gene_order = ["thymus", "liver", "skin", "spleen"]
    gene_counts = (
        candidate_status.groupby(["tissue", "matched_status"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(gene_order, fill_value=0)
    )
    promoted = gene_counts.get("promoted", pd.Series(0, index=gene_counts.index))
    reinforced = gene_counts.get(
        "shared_importance",
        pd.Series(0, index=gene_counts.index),
    )
    y_genes = np.arange(len(gene_counts))
    ax.barh(y_genes, reinforced, color=teal, label="Reinforced")
    ax.barh(y_genes, promoted, left=reinforced, color=coral, label="Promoted")
    ax.set_yticks(y_genes, [name.title() for name in gene_counts.index])
    ax.invert_yaxis()
    ax.set_xlabel("Genes identified by the analysis")
    ax.set_title("B  Synthetic-supported genes", loc="left")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.grid(axis="x", color="#dce3e6", linewidth=0.7)

    ax = fig.add_subplot(grid[1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("C  Overlap between feature analyses", loc="left")
    ax.add_patch(Circle((0.42, 0.55), 0.29, color=blue, alpha=0.24))
    ax.add_patch(Circle((0.61, 0.55), 0.29, color=coral, alpha=0.24))
    ax.text(0.24, 0.84, "Consensus panels", color=blue, fontsize=8.2, weight="bold", ha="center")
    ax.text(0.79, 0.84, "Matched all-gene", color=coral, fontsize=8.2, weight="bold", ha="center")
    for x, value, color in [
        (0.28, comparison.consensus_only, blue),
        (0.515, comparison.supported_by_both, dark),
        (0.76, comparison.matched_only, coral),
    ]:
        ax.text(
            x,
            0.54,
            f"{int(value)}",
            fontsize=18,
            color=color,
            weight="bold",
            ha="center",
            va="center",
        )
    ax.text(0.28, 0.22, "consensus only", fontsize=7.2, color=gray, ha="center")
    ax.text(0.515, 0.22, "both", fontsize=7.2, color=gray, ha="center")
    ax.text(0.76, 0.22, "matched only", fontsize=7.2, color=gray, ha="center")

    fig.suptitle(
        "Synthetic training changed classification and feature selection",
        x=0.02,
        ha="left",
        fontsize=15,
        weight="bold",
        color=dark,
    )
    fig.subplots_adjust(top=0.90, bottom=0.07, left=0.15, right=0.98)
    _save_figure(fig, "figure_9_synthetic_tissue_evidence")


def figure_s3_synthetic_soleus() -> None:
    genes = pd.read_csv(
        _required(
            ROOT
            / "paper/synthetic_guided_spaceflight/source_data/table_s7_soleus_genes.tsv"
        ),
        sep="\t",
    )
    genes = genes.set_index("symbol").loc[
        ["Bdh1", "Bnip3", "Ech1", "Decr1", "Tpm1"]
    ].reset_index()

    lower = "#d96653"
    higher = "#218b87"
    dark = "#213f52"
    gray = "#52616a"
    fig = plt.figure(figsize=(11.6, 4.5))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[0.92, 1.35],
        left=0.07,
        right=0.98,
        top=0.82,
        bottom=0.16,
        wspace=0.26,
    )

    ax = fig.add_subplot(grid[0, 0])
    effects = genes["real_meta_effect"].to_numpy(dtype=float)
    positions = np.arange(len(genes))
    colors = [lower if value < 0 else higher for value in effects]
    ax.barh(positions, effects, color=colors, height=0.66)
    ax.axvline(0, color="#75838a", linewidth=0.9)
    ax.set_yticks(positions, genes["symbol"], fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Flight - ground expression")
    ax.set_title("A  Soleus gene effects", loc="left", color=dark)
    limit = max(abs(effects.min()), abs(effects.max())) * 1.35
    ax.set_xlim(-limit, limit)
    for y, value in zip(positions, effects):
        offset = limit * 0.035
        ax.text(
            value + (offset if value >= 0 else -offset),
            y,
            f"{value:+.3f}",
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=8,
            color=dark,
        )
    ax.text(
        0.01,
        -0.19,
        "Coral: lower in flight    Teal: higher in flight",
        transform=ax.transAxes,
        fontsize=7.7,
        color=gray,
    )

    interpretation = fig.add_subplot(grid[0, 1])
    interpretation.set_xlim(0, 1)
    interpretation.set_ylim(0, 1)
    interpretation.axis("off")
    interpretation.set_title("B  Biological interpretation", loc="left", color=dark)
    rows = [
        (0.70, "Oxidative substrate use and lipid oxidation", "Bdh1, Ech1, Decr1", "#f7e4df", lower),
        (0.43, "Mitochondrial quality control", "Bnip3", "#eee8f5", "#8065a5"),
        (0.16, "Contractile remodeling", "Tpm1", "#def0ee", higher),
    ]
    for y, label, symbols, fill, edge in rows:
        box = FancyBboxPatch(
            (0.05, y),
            0.90,
            0.19,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            facecolor=fill,
            edgecolor=edge,
            linewidth=1.2,
        )
        interpretation.add_patch(box)
        interpretation.text(0.09, y + 0.125, label, fontsize=9.2, weight="bold", color=dark, va="center")
        interpretation.text(0.09, y + 0.058, symbols, fontsize=8.5, color=gray, va="center", fontstyle="italic")

    fig.suptitle(
        "A soleus consensus panel links metabolism with contractile remodeling",
        x=0.02,
        ha="left",
        fontsize=15,
        weight="bold",
        color=dark,
    )
    _save_figure(fig, "figure_s3_synthetic_soleus")


def copy_figures() -> None:
    for source, destination in COPIED_FIGURES.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(_required(source), destination)


def write_manifest() -> None:
    rows = []
    tracked_sources = [
        ROOT / "outputs/glare/study_effects/aggregate_vs_mober_study_effect_summary.tsv",
        ROOT
        / "paper/slstp_internship_report/source_data/glare/"
        / "skeletal_muscle_aggregate_vs_mober_umap_by_accession.png",
        ROOT / "paper/asgsr_expimap_hvg/manuscript.md",
        ROOT / "paper/asgsr_expimap_hvg/source_data/table_2_retained_pathway_evidence.tsv",
        ROOT / "paper/asgsr_expimap_hvg/source_data/table_s5_accession_pathway_effects.tsv.gz",
        SOURCE_DIR / "expimap_kidney_spleen_seed_accession_effects.tsv.gz",
        ROOT / "paper/synthetic_guided_spaceflight/manuscript.md",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_4_generator_model_selection.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s1_archs4_ddim_metrics.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s2_locked_ddim_repeats.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s18_matched_all_gene_utility.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s19_matched_all_gene_candidates.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s21_matched_consensus_comparison.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s22_matched_gene_literature_annotations.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s23_grouped_pathway_literature_annotations.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s24_importance_literature_sources.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s16_promoted_gene_literature_annotations.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s17_promoted_gene_literature_sources.tsv",
        ROOT / "paper/synthetic_guided_spaceflight/source_data/table_s7_soleus_genes.tsv",
        ROOT / "docs/distributed_response_hypotheses.md",
        *COPIED_FIGURES.keys(),
    ]
    for path in tracked_sources:
        path = _required(path)
        rows.append((str(path.relative_to(ROOT)), path.stat().st_size, _sha256(path)))
    with (SOURCE_DIR / "source_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("path\tbytes\tsha256\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def render_document(markdown_path: Path, title: str) -> tuple[Path, Path]:
    try:
        import markdown
        from weasyprint import HTML
    except ImportError as error:
        raise RuntimeError("Rendering requires markdown and weasyprint") from error

    body = markdown.markdown(
        markdown_path.read_text(encoding="utf-8"),
        extensions=["tables", "fenced_code", "sane_lists", "md_in_html"],
    )
    css = """
    @page {
      size: A4;
      margin: 14mm 15mm 15mm 15mm;
      @bottom-center { content: counter(page); color: #69767d; font-size: 7.4pt; }
    }
    @page:first { @bottom-center { content: none; } }
    html { font-family: "DejaVu Serif", Georgia, serif; color: #1f2c33; }
    body { font-size: 9pt; line-height: 1.31; }
    h1 { color: #203f52; font-family: "DejaVu Sans", Arial, sans-serif; font-size: 20pt; line-height: 1.12; margin: 0 0 3.5mm 0; letter-spacing: 0; }
    h2 { color: #203f52; font-family: "DejaVu Sans", Arial, sans-serif; font-size: 11.8pt; margin: 4.5mm 0 1.7mm 0; border-bottom: 0.45pt solid #aebcc3; padding-bottom: 0.8mm; letter-spacing: 0; }
    h3 { color: #176f70; font-family: "DejaVu Sans", Arial, sans-serif; font-size: 9.4pt; margin: 3.1mm 0 1mm 0; letter-spacing: 0; }
    p { margin: 0 0 1.8mm 0; text-align: justify; hyphens: auto; widows: 2; orphans: 2; }
    strong { color: #16252d; }
    a { color: #1f6685; text-decoration: none; }
    ul, ol { margin: 1mm 0 2mm 4.2mm; padding-left: 3.7mm; }
    li { margin-bottom: 0.6mm; }
    table { width: 100%; border-collapse: collapse; margin: 2.1mm 0 3.2mm 0; font-size: 6.9pt; line-height: 1.18; break-inside: avoid; }
    .table-block { break-inside: avoid; }
    th { background: #e8eef0; color: #203f52; font-family: "DejaVu Sans", Arial, sans-serif; text-align: left; padding: 1.05mm; border-bottom: 0.8pt solid #68767d; }
    td { padding: 0.95mm 1.05mm; border-bottom: 0.35pt solid #ced6da; vertical-align: top; }
    img { max-width: 100%; max-height: 170mm; display: block; margin: 1.5mm auto 1mm auto; }
    .figure { break-inside: avoid; margin: 2.2mm 0 3mm 0; }
    .figure img { width: 100%; max-height: 170mm; object-fit: contain; }
    .figure.compact img { width: 92%; max-height: 135mm; }
    .figure.narrow img { width: 82%; max-height: 130mm; }
    .caption { font-size: 7.15pt; line-height: 1.24; color: #3a484f; text-align: left; margin: 0 0 3.4mm 0; }
    .report-meta { font-family: "DejaVu Sans", Arial, sans-serif; color: #4b5d67; font-size: 8.5pt; text-align: left; margin-bottom: 3.5mm; }
    .abstract { background: #f3f6f7; border-left: 2.4pt solid #287d7a; padding: 2.6mm 3.2mm 1.8mm 3.2mm; margin: 2.5mm 0 3.5mm 0; }
    .abstract h2 { margin-top: 0; border: 0; padding: 0; }
    .keywords { font-size: 7.7pt; color: #4b5d67; text-align: left; }
    .callout { background: #f5f2e9; border-left: 2.4pt solid #c58c2b; padding: 2mm 3mm 1mm 3mm; margin: 2.2mm 0; break-inside: avoid; }
    .page-break { page-break-before: always; }
    .references { font-size: 7.1pt; line-height: 1.18; }
    .references p { margin-bottom: 0.72mm; text-align: left; }
    .appendix { font-size: 7.8pt; line-height: 1.22; }
    code { font-family: "DejaVu Sans Mono", monospace; font-size: 7pt; background: #f0f3f4; padding: 0.15mm 0.4mm; }
    """
    html_text = f"""<!doctype html>
<html><head><meta charset=\"utf-8\"><title>{title}</title><style>{css}</style></head>
<body>{body}</body></html>"""
    html_path = markdown_path.with_suffix(".html")
    pdf_path = markdown_path.with_suffix(".pdf")
    html_path.write_text(html_text, encoding="utf-8")
    HTML(string=html_text, base_url=str(PAPER_DIR)).write_pdf(str(pdf_path))
    return html_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-render", action="store_true")
    args = parser.parse_args()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    _style()
    for stale_name in ("figure_2_glare_validation.png", "figure_2_glare_validation.pdf"):
        (FIGURE_DIR / stale_name).unlink(missing_ok=True)
    for stale_name in ("glare_validation_summary.tsv", "glare_direction_summary.tsv"):
        (SOURCE_DIR / stale_name).unlink(missing_ok=True)
    figure_2_glare_batch_effects()
    figure_3_expimap_architecture()
    figure_4_expimap_pathway_heatmap()
    figure_5_expimap_annotation_workflow()
    figure_7_generator_validation()
    figure_9_synthetic_tissue_evidence()
    figure_s3_synthetic_soleus()
    copy_figures()
    write_manifest()

    if not args.skip_render:
        render_document(
            _required(PAPER_DIR / "manuscript.md"),
            "Interpretable and generative modeling of mouse spaceflight transcriptomes",
        )
    print(f"Built internship report: {PAPER_DIR}")


if __name__ == "__main__":
    main()
