"""Create UMAP copies with prior-work and hidden-module legend highlights."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

from .multi_tissue_reports import clean_reactome_term, cluster_colors


DEFAULT_ROOT = Path("outputs/glare_multi_tissue_api")
DEFAULT_VALIDATION_DIR = DEFAULT_ROOT / "validation_stack_terms15"

PRIOR_COLOR = "#2563eb"
HIDDEN_COLOR = "#f97316"
OTHER_COLOR = "#6b7280"
UNCLEAR_COLOR = "#9ca3af"

HIDDEN_PRIORITY_TERMS = {
    "skeletal_muscle": {
        "Bmal1 Clock Npas2 Activates Circadian Expression",
        "Cyclin E Associated Events During G1 S Transition",
    },
    "skeletal_muscle_soleus": {
        "Signalling By Ngf",
    },
    "kidney": {
        "Membrane Trafficking",
        "Signaling By Insulin Receptor",
    },
    "thymus": {
        "Response To Elevated Platelet Cytosolic Ca2",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy aggregate GLARE UMAPs with legend entries highlighted as "
            "prior-work aligned or hidden/novel."
        )
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument(
        "--scopes",
        nargs="+",
        default=["aggregate", "aggregate_mober"],
        help="Run scopes to scan under each tissue directory.",
    )
    return parser.parse_args()


def load_prior_terms(validation_dir: Path) -> dict[str, set[str]]:
    prior_terms: dict[str, set[str]] = {}
    candidate_modules = validation_dir / "candidate_modules.tsv"
    if candidate_modules.exists():
        table = pd.read_csv(candidate_modules, sep="\t")
        table = table.loc[table["module_class"].astype(str).eq("intersection")]
        for tissue, group in table.groupby("tissue"):
            prior_terms.setdefault(str(tissue), set()).update(group["clean_term"].dropna())

    candidate_summary = validation_dir / "novel_candidate_validation" / "candidate_validation_summary.tsv"
    if candidate_summary.exists():
        table = pd.read_csv(candidate_summary, sep="\t")
        known = table.loc[
            table["candidate_group"].astype(str).str.contains("validation_anchor", na=False)
            | table["theme"].astype(str).str.contains("known", case=False, na=False)
        ]
        for tissue, group in known.groupby("tissue"):
            prior_terms.setdefault(str(tissue), set()).update(group["clean_term"].dropna())
    return prior_terms


def tissue_from_run(run_dir: Path) -> str:
    if run_dir.name in {"aggregate", "aggregate_mober"}:
        return run_dir.parent.name
    if run_dir.parent.name == "per_study":
        return run_dir.parent.parent.name
    return run_dir.parent.name


def load_enrichment_terms(run_dir: Path, location: str, alpha: float) -> pd.DataFrame:
    path = run_dir / "clustering" / "reactome_enrichment.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table = pd.read_csv(path, sep="\t")
    if table.empty:
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table = table.loc[table["location"].astype(str).eq(location)].copy()
    table = table.loc[pd.to_numeric(table["fdr_bh"], errors="coerce").le(alpha)].copy()
    if table.empty:
        return pd.DataFrame(columns=["cluster", "clean_term", "fdr_bh"])
    table["cluster"] = table["cluster"].astype(int)
    table["clean_term"] = table["term"].map(clean_reactome_term)
    return table[["cluster", "clean_term", "fdr_bh"]]


def classify_clusters(
    annotations: pd.DataFrame,
    enrichment: pd.DataFrame,
    tissue: str,
    prior_terms: dict[str, set[str]],
) -> pd.DataFrame:
    hidden = HIDDEN_PRIORITY_TERMS.get(tissue, set())
    prior = prior_terms.get(tissue, set())
    rows = []
    for row in annotations.itertuples(index=False):
        cluster = int(row.cluster)
        terms = enrichment.loc[enrichment["cluster"].eq(cluster), "clean_term"].dropna().tolist()
        hidden_matches = sorted({term for term in terms if term in hidden})
        prior_matches = sorted({term for term in terms if term in prior})
        if hidden_matches:
            category = "hidden_novel"
            category_label = "Hidden/novel"
            category_color = HIDDEN_COLOR
            matched_terms = hidden_matches
        elif prior_matches:
            category = "prior_work"
            category_label = "Prior-work aligned"
            category_color = PRIOR_COLOR
            matched_terms = prior_matches
        elif str(getattr(row, "annotation_status", "")) == "ambiguous":
            category = "unclear"
            category_label = "Unclear/ambiguous"
            category_color = UNCLEAR_COLOR
            matched_terms = []
        else:
            category = "other"
            category_label = "Other"
            category_color = OTHER_COLOR
            matched_terms = []
        data = row._asdict()
        data.update(
            {
                "highlight_category": category,
                "highlight_label": category_label,
                "highlight_color": category_color,
                "matched_hidden_terms": "; ".join(hidden_matches),
                "matched_prior_terms": "; ".join(prior_matches),
                "matched_highlight_terms": "; ".join(matched_terms[:4]),
            }
        )
        rows.append(data)
    return pd.DataFrame(rows)


def plot_highlighted_umap(
    coords: pd.DataFrame,
    labels: np.ndarray,
    annotations: pd.DataFrame,
    title: str,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clusters = annotations["cluster"].astype(int).tolist()
    colors = cluster_colors(clusters)
    grouped = [
        ("hidden_novel", "Hidden/novel modules"),
        ("prior_work", "Prior-work aligned modules"),
        ("other", "Other modules"),
        ("unclear", "Unclear/ambiguous modules"),
    ]

    wrapped_groups = []
    total_sections = 0
    total_label_lines = 0
    for category, legend_title in grouped:
        subset = annotations.loc[annotations["highlight_category"].eq(category)].copy()
        if subset.empty:
            continue
        entries = []
        total_sections += 1
        for row in subset.itertuples(index=False):
            cluster = int(row.cluster)
            matched = str(row.matched_highlight_terms)
            has_matched = bool(matched and matched != "nan")
            primary_description = matched if has_matched else str(row.cluster_description)
            source_description = str(row.cluster_description)
            suffix = ""
            if has_matched and primary_description != source_description:
                suffix = f" [top cluster label: {source_description}]"
            label = (
                f"C{cluster} ({int(row.gene_count)} genes): "
                f"{primary_description}{suffix}"
            )
            wrapped = textwrap.wrap(label, width=62)
            entries.append((row, wrapped))
            total_label_lines += max(1, len(wrapped))
        wrapped_groups.append((category, legend_title, entries))

    legend_height = (
        1.0
        + 0.32 * total_sections
        + 0.14 * total_label_lines
        + 0.07 * sum(len(entries) for _, _, entries in wrapped_groups)
    )
    figure_height = max(8.6, legend_height)

    figure = plt.figure(figsize=(15.5, figure_height), dpi=170)
    grid = figure.add_gridspec(1, 2, width_ratios=[2.0, 1.35], wspace=0.04)
    axis = figure.add_subplot(grid[0, 0])
    legend_axis = figure.add_subplot(grid[0, 1])
    legend_axis.set_xlim(0, 1)
    legend_axis.set_ylim(0, 1)
    legend_axis.axis("off")

    for cluster in clusters:
        mask = labels == cluster
        axis.scatter(
            coords.loc[mask, "x"],
            coords.loc[mask, "y"],
            s=3,
            alpha=0.58,
            linewidths=0,
            color=colors[cluster],
        )

    axis.set_title(title, fontsize=11, pad=14)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")

    legend_axis.text(
        0.0,
        0.995,
        "Cluster Legend",
        fontsize=11,
        fontweight="bold",
        va="top",
        transform=legend_axis.transAxes,
    )
    legend_axis.text(
        0.0,
        0.958,
        "Orange text = hidden/novel\nBlue text = prior-work aligned\nDots keep cluster colors",
        fontsize=7.5,
        color="#374151",
        va="top",
        transform=legend_axis.transAxes,
    )

    line_step = 0.14 / figure_height
    entry_gap = 0.07 / figure_height
    section_step = 0.28 / figure_height
    section_gap = 0.12 / figure_height
    y = 0.89
    for category, legend_title, entries in wrapped_groups:
        title_color = {
            "hidden_novel": HIDDEN_COLOR,
            "prior_work": PRIOR_COLOR,
            "other": "#374151",
            "unclear": UNCLEAR_COLOR,
        }[category]
        legend_axis.text(
            0.0,
            max(y, 0.035),
            legend_title,
            fontsize=8.2,
            fontweight="bold",
            color=title_color,
            va="top",
            transform=legend_axis.transAxes,
        )
        y -= section_step
        for row, wrapped in entries:
            cluster = int(row.cluster)
            text_color = (
                HIDDEN_COLOR
                if category == "hidden_novel"
                else PRIOR_COLOR
                if category == "prior_work"
                else UNCLEAR_COLOR
                if category == "unclear"
                else "#374151"
            )
            legend_axis.scatter(
                [0.018],
                [y - 0.004],
                s=32,
                color=[colors[cluster]],
                edgecolors=[colors[cluster]],
                linewidths=0.7,
                transform=legend_axis.transAxes,
                clip_on=True,
            )
            legend_axis.text(
                0.052,
                y,
                "\n".join(wrapped),
                fontsize=6.7,
                color=text_color,
                va="top",
                transform=legend_axis.transAxes,
                clip_on=True,
            )
            y -= line_step * max(1, len(wrapped)) + entry_gap
        y -= section_gap

    figure.savefig(output_path)
    plt.close(figure)


def process_run(run_dir: Path, prior_terms: dict[str, set[str]], alpha: float) -> list[dict]:
    tissue = tissue_from_run(run_dir)
    output_rows: list[dict] = []
    for location in ("FLT", "GC"):
        coords_path = run_dir / "plots" / f"{location}_umap_coordinates.tsv"
        clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        annotations_path = run_dir / "plots" / f"{location}_cluster_annotations.tsv"
        if not coords_path.exists() or not clusters_path.exists() or not annotations_path.exists():
            continue
        coords = pd.read_csv(coords_path, sep="\t")
        clusters = pd.read_csv(clusters_path, sep="\t")
        labels = clusters["consensus"].astype(int).to_numpy()
        annotations = pd.read_csv(annotations_path, sep="\t")
        enrichment = load_enrichment_terms(run_dir, location, alpha)
        highlighted = classify_clusters(annotations, enrichment, tissue, prior_terms)

        output_annotations = run_dir / "plots" / f"{location}_cluster_highlight_annotations.tsv"
        highlighted.to_csv(output_annotations, sep="\t", index=False)

        output_plot = run_dir / "plots" / f"{location}_umap_prior_hidden_legend.png"
        plot_highlighted_umap(
            coords,
            labels,
            highlighted,
            f"{tissue} {run_dir.name} {location} GLARE clusters",
            output_plot,
        )

        counts = highlighted["highlight_category"].value_counts().to_dict()
        output_rows.append(
            {
                "tissue": tissue,
                "scope": run_dir.name,
                "location": location,
                "plot": str(output_plot),
                "annotations": str(output_annotations),
                "hidden_novel_clusters": int(counts.get("hidden_novel", 0)),
                "prior_work_clusters": int(counts.get("prior_work", 0)),
                "other_clusters": int(counts.get("other", 0)),
                "unclear_clusters": int(counts.get("unclear", 0)),
            }
        )
    return output_rows


def main() -> None:
    args = parse_args()
    prior_terms = load_prior_terms(args.validation_dir)
    manifest_rows: list[dict] = []
    for tissue_dir in sorted(args.root.iterdir()):
        if not tissue_dir.is_dir():
            continue
        for scope in args.scopes:
            run_dir = tissue_dir / scope
            if run_dir.exists():
                manifest_rows.extend(process_run(run_dir, prior_terms, args.alpha))

    manifest = pd.DataFrame(manifest_rows)
    output_path = args.root / "priority_umap_highlight_manifest.tsv"
    manifest.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {len(manifest)} highlighted UMAP rows to {output_path}")


if __name__ == "__main__":
    main()
