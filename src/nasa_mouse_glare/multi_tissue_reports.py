"""Reporting utilities for API-derived multi-tissue GLARE runs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .io import require_import
from .per_study_glare import DEFAULT_REACTOME_GMT, cluster_enrichment_for_run


DEFAULT_ROOT = "outputs/glare_multi_tissue_api"
DEFAULT_ALPHA = 0.05


def clean_reactome_term(term: str) -> str:
    text = str(term).replace("REACTOME_", "")
    text = text.replace("_", " ").strip()
    text = re.sub(r"\s+", " ", text)
    return text.title()


def short_text(text: str, max_chars: int = 70) -> str:
    text = re.sub(r"\s+", " ", str(text)).strip()
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 1].rstrip() + "..."


def find_tissue_dir(run_dir: Path) -> Path:
    run_dir = run_dir.resolve()
    if run_dir.name in {"aggregate", "aggregate_mober"}:
        return run_dir.parent
    if run_dir.parent.name == "per_study":
        return run_dir.parent.parent
    return run_dir.parent


def run_kind(run_dir: Path) -> str:
    if run_dir.name == "aggregate":
        return "aggregate"
    if run_dir.name == "aggregate_mober":
        return "aggregate_mober"
    if run_dir.parent.name == "per_study":
        return f"per_study/{run_dir.name}"
    return run_dir.name


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_input_summary(run_dir: Path) -> dict:
    summary = read_json(run_dir / "input_summary.json")
    if summary:
        return summary
    source_scope = read_json(run_dir / "input_summary.json").get("source_scope_dir")
    if source_scope:
        return read_json(Path(source_scope) / "input_summary.json")
    return {}


def location_profile_count(input_summary: dict, location: str) -> int | None:
    key = "flt_profiles" if location == "FLT" else "gc_profiles"
    if key in input_summary:
        return int(input_summary[key])
    counts = pd.DataFrame(input_summary.get("condition_counts", []))
    if location in counts:
        return int(counts[location].sum())
    return None


def load_or_compute_reactome(run_dir: Path, gmt_path: Path, min_overlap: int) -> pd.DataFrame:
    tissue_dir = find_tissue_dir(run_dir)
    accession = run_dir.name if run_dir.parent.name == "per_study" else None
    per_study_path = tissue_dir / "per_study" / "glare_cluster_reactome_enrichment.tsv"
    if accession and per_study_path.exists():
        table = pd.read_csv(per_study_path, sep="\t")
        return table.loc[table["accession"].astype(str).eq(accession)].copy()

    output_path = run_dir / "clustering" / "reactome_enrichment.tsv"
    if output_path.exists():
        return pd.read_csv(output_path, sep="\t")
    table = cluster_enrichment_for_run(run_dir, gmt_path, min_overlap)
    if not table.empty:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(output_path, sep="\t", index=False)
    return table


def load_dgea_cluster_flags(run_dir: Path) -> pd.DataFrame:
    if run_dir.parent.name != "per_study":
        return pd.DataFrame()
    tissue_dir = find_tissue_dir(run_dir)
    path = tissue_dir / "dgea_comparison" / f"{run_dir.name}_cluster_dgea_enrichment.tsv"
    if not path.exists():
        return pd.DataFrame()
    table = pd.read_csv(path, sep="\t")
    if "fisher_fdr_bh" not in table:
        return pd.DataFrame()
    return table


def build_cluster_annotations(
    run_dir: Path,
    location: str,
    alpha: float,
    gmt_path: Path,
    min_overlap: int,
) -> pd.DataFrame:
    clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
    clusters = pd.read_csv(clusters_path, sep="\t")
    counts = (
        clusters.groupby("consensus", as_index=False)
        .size()
        .rename(columns={"consensus": "cluster", "size": "gene_count"})
    )

    reactome = load_or_compute_reactome(run_dir, gmt_path, min_overlap)
    if reactome.empty:
        best_reactome = pd.DataFrame(
            columns=[
                "cluster",
                "best_reactome_term",
                "best_reactome_fdr_bh",
                "best_reactome_overlap",
                "term_genes_in_universe",
            ]
        )
    else:
        reactome = reactome.loc[reactome["location"].astype(str).eq(location)].copy()
        reactome["clean_term"] = reactome["term"].map(clean_reactome_term)
        best_reactome = (
            reactome.sort_values(["cluster", "fdr_bh", "p_value"])
            .groupby("cluster", as_index=False)
            .first()[
                ["cluster", "clean_term", "fdr_bh", "overlap", "term_genes_in_universe"]
            ]
            .rename(
                columns={
                    "clean_term": "best_reactome_term",
                    "fdr_bh": "best_reactome_fdr_bh",
                    "overlap": "best_reactome_overlap",
                }
            )
        )

    dgea = load_dgea_cluster_flags(run_dir)
    if dgea.empty:
        best_dgea = pd.DataFrame(
            columns=[
                "cluster",
                "dgea_fisher_fdr_bh",
                "significant_padj05_genes",
                "significant_padj05_abs_lfc_genes",
                "top_significant_genes",
            ]
        )
    else:
        dgea = dgea.loc[dgea["location"].astype(str).eq(location)].copy()
        best_dgea = (
            dgea.sort_values(["cluster", "fisher_fdr_bh"])
            .groupby("cluster", as_index=False)
            .first()[
                [
                    "cluster",
                    "fisher_fdr_bh",
                    "significant_padj05_genes",
                    "significant_padj05_abs_lfc_genes",
                    "top_significant_genes",
                ]
            ]
            .rename(columns={"fisher_fdr_bh": "dgea_fisher_fdr_bh"})
        )

    annotations = counts.merge(best_reactome, on="cluster", how="left").merge(
        best_dgea, on="cluster", how="left"
    )
    annotations["reactome_significant"] = (
        annotations["best_reactome_fdr_bh"].fillna(1.0).astype(float) < alpha
    )
    annotations["dgea_enriched"] = (
        annotations["dgea_fisher_fdr_bh"].fillna(1.0).astype(float) < alpha
    )
    descriptions = []
    statuses = []
    for row in annotations.itertuples(index=False):
        term = getattr(row, "best_reactome_term", "")
        term = "" if pd.isna(term) else str(term)
        if bool(row.dgea_enriched):
            status = "DGEA-enriched"
            detail = f"; {term}" if term else "; no significant Reactome label"
        elif bool(row.reactome_significant):
            status = "Reactome-enriched"
            detail = f"; {term}"
        else:
            status = "ambiguous"
            detail = "; no significant DGEA overlap or Reactome label"
        statuses.append(status)
        descriptions.append(f"{status}{detail}")
    annotations["annotation_status"] = statuses
    annotations["cluster_description"] = descriptions
    return annotations.sort_values("cluster")


def compute_coordinates(
    representation: np.ndarray,
    run_dir: Path,
    location: str,
    skip_umap: bool,
) -> dict[str, Path]:
    output_dir = run_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    scaled = StandardScaler().fit_transform(representation)
    pca = PCA(n_components=2, random_state=2024).fit_transform(scaled)
    pca_path = output_dir / f"{location}_pca_coordinates.tsv"
    pd.DataFrame({"x": pca[:, 0], "y": pca[:, 1]}).to_csv(
        pca_path, sep="\t", index=False
    )

    if skip_umap:
        return {"pca": pca_path}

    umap_path = output_dir / f"{location}_umap_coordinates.tsv"
    if not umap_path.exists():
        umap_module = require_import("umap", "pip install umap-learn")
        n_neighbors = min(15, max(2, representation.shape[0] - 1))
        reducer = umap_module.UMAP(
            n_components=2,
            n_neighbors=n_neighbors,
            min_dist=0.25,
            n_epochs=100,
            metric="euclidean",
            init="random",
            low_memory=True,
            n_jobs=-1,
        )
        umap = reducer.fit_transform(scaled)
        pd.DataFrame({"x": umap[:, 0], "y": umap[:, 1]}).to_csv(
            umap_path, sep="\t", index=False
        )
    return {"pca": pca_path, "umap": umap_path}


def cluster_colors(clusters: list[int]):
    plt = require_import("matplotlib.pyplot", "pip install matplotlib")
    cmap = plt.get_cmap("tab20")
    return {cluster: cmap(index % 20) for index, cluster in enumerate(clusters)}


def plot_coordinates(
    coords: pd.DataFrame,
    labels: np.ndarray,
    annotations: pd.DataFrame,
    title: str,
    subtitle: str,
    output_path: Path,
) -> None:
    matplotlib = require_import("matplotlib", "pip install matplotlib")
    matplotlib.use("Agg")
    plt = require_import("matplotlib.pyplot", "pip install matplotlib")
    from matplotlib.lines import Line2D

    clusters = annotations["cluster"].astype(int).tolist()
    colors = cluster_colors(clusters)
    figure, axis = plt.subplots(figsize=(13.5, 8), dpi=160)
    for cluster in clusters:
        mask = labels == cluster
        axis.scatter(
            coords.loc[mask, "x"],
            coords.loc[mask, "y"],
            s=3,
            alpha=0.65,
            linewidths=0,
            color=colors[cluster],
        )
    axis.set_title(f"{title}\n{subtitle}", fontsize=11, pad=14)
    axis.set_xlabel("Dimension 1")
    axis.set_ylabel("Dimension 2")

    interpretable_handles = []
    ambiguous_handles = []
    for row in annotations.itertuples(index=False):
        cluster = int(row.cluster)
        label = f"C{cluster} ({int(row.gene_count)} genes): {short_text(row.cluster_description)}"
        handle = Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=5,
            markerfacecolor=colors[cluster],
            markeredgecolor=colors[cluster],
            label=label,
        )
        if row.annotation_status == "ambiguous":
            ambiguous_handles.append(handle)
        else:
            interpretable_handles.append(handle)

    if interpretable_handles:
        first = axis.legend(
            handles=interpretable_handles,
            title="Interpretable modules",
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0,
            fontsize=7,
            title_fontsize=8,
            frameon=False,
        )
        axis.add_artist(first)
    if ambiguous_handles:
        y_anchor = 0.48 if interpretable_handles else 1.0
        axis.legend(
            handles=ambiguous_handles,
            title="Ambiguous/unclear modules",
            loc="upper left",
            bbox_to_anchor=(1.01, y_anchor),
            borderaxespad=0,
            fontsize=7,
            title_fontsize=8,
            frameon=False,
        )
    figure.tight_layout(rect=[0, 0, 0.72, 1])
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)


def plot_run(run_dir: Path, args: argparse.Namespace) -> Path:
    run_dir = Path(run_dir)
    input_summary = read_input_summary(run_dir)
    tissue_label = input_summary.get("tissue_label") or find_tissue_dir(run_dir).name
    output_dir = run_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_rows = []
    for location in ("FLT", "GC"):
        representation_path = run_dir / f"{location}_FTSAE_representation.npy"
        clusters_path = run_dir / "clustering" / f"{location}_gene_clusters.tsv"
        if not representation_path.exists() or not clusters_path.exists():
            continue
        representation = np.load(representation_path)
        clusters = pd.read_csv(clusters_path, sep="\t")
        labels = clusters["consensus"].astype(int).to_numpy()
        annotations = build_cluster_annotations(
            run_dir, location, args.alpha, Path(args.reactome_gmt), args.min_overlap
        )
        annotations.to_csv(
            output_dir / f"{location}_cluster_annotations.tsv", sep="\t", index=False
        )
        coord_paths = compute_coordinates(
            representation, run_dir, location, args.skip_umap
        )
        profile_count = location_profile_count(input_summary, location)
        profile_text = "profiles=NA" if profile_count is None else f"profiles={profile_count}"
        subtitle = (
            f"points=genes; genes={representation.shape[0]}; {profile_text}; "
            f"run={run_kind(run_dir)}"
        )
        for coord_name, coord_path in coord_paths.items():
            coords = pd.read_csv(coord_path, sep="\t")
            output_path = output_dir / f"{location}_{coord_name}_by_cluster_legend.png"
            plot_coordinates(
                coords,
                labels,
                annotations,
                f"{tissue_label} {location} GLARE clusters ({coord_name.upper()})",
                subtitle,
                output_path,
            )
            plot_rows.append(
                {
                    "location": location,
                    "projection": coord_name,
                    "coordinates": str(coord_path),
                    "plot": str(output_path),
                    "genes": int(representation.shape[0]),
                    "profiles": profile_count,
                    "clusters": int(annotations["cluster"].nunique()),
                    "interpretable_clusters": int(
                        annotations["annotation_status"].ne("ambiguous").sum()
                    ),
                    "ambiguous_clusters": int(
                        annotations["annotation_status"].eq("ambiguous").sum()
                    ),
                }
            )
    manifest = pd.DataFrame(plot_rows)
    manifest.to_csv(output_dir / "plot_manifest.tsv", sep="\t", index=False)
    return output_dir


def iter_run_dirs(tissue_dir: Path, include_per_study: bool, include_mober: bool) -> list[Path]:
    run_dirs = []
    aggregate = tissue_dir / "aggregate"
    if (aggregate / "clustering").exists():
        run_dirs.append(aggregate)
    mober = tissue_dir / "aggregate_mober"
    if include_mober and (mober / "clustering").exists():
        run_dirs.append(mober)
    if include_per_study:
        per_study = tissue_dir / "per_study"
        if per_study.exists():
            run_dirs.extend(
                sorted(
                    path
                    for path in per_study.iterdir()
                    if (path / "clustering").exists()
                )
            )
    return run_dirs


def plot_tissue(tissue_dir: Path, args: argparse.Namespace) -> Path:
    for run_dir in iter_run_dirs(tissue_dir, args.include_per_study, args.include_mober):
        print(f"Plotting {run_dir}", flush=True)
        plot_run(run_dir, args)
    return tissue_dir


def format_float(value) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.3g}"
    except (TypeError, ValueError):
        return str(value)


def tissue_summary_row(tissue_dir: Path) -> dict:
    prep = read_json(tissue_dir / "PREP_SUMMARY.json")
    aggregate_input = read_json(tissue_dir / "aggregate" / "input_summary.json")
    clustering = read_json(tissue_dir / "aggregate" / "clustering" / "clustering_summary.json")
    locations = {item["location"]: item for item in clustering.get("locations", [])}
    condition_counts = pd.DataFrame(aggregate_input.get("condition_counts", []))
    accessions = sorted(condition_counts["id.accession"].astype(str).tolist()) if not condition_counts.empty and "id.accession" in condition_counts else []
    per_study = tissue_dir / "per_study"
    per_study_total = len(list(per_study.glob("*/input_summary.json"))) if per_study.exists() else 0
    per_study_done = len(list(per_study.glob("*/finetune_summary.json"))) if per_study.exists() else 0
    dgea_summary = tissue_dir / "dgea_comparison" / "per_study_glare_dgea_summary.tsv"
    dgea_done = dgea_summary.exists()
    recurring = tissue_dir / "dgea_comparison" / "recurring_dgea_glare_pathway_overlap.tsv"
    top_terms = ""
    if recurring.exists():
        table = pd.read_csv(recurring, sep="\t")
        if not table.empty and "clean_term" in table:
            top_terms = "; ".join(
                table.sort_values(["study_count", "best_dgea_fdr_bh"], ascending=[False, True])
                .head(5)["clean_term"]
                .astype(str)
            )
    return {
        "tissue": tissue_dir.name,
        "tissue_label": aggregate_input.get("tissue_label", prep.get("tissue_label", tissue_dir.name)),
        "accessions": ",".join(accessions),
        "studies": len(accessions),
        "flt_profiles": aggregate_input.get("flt_profiles", ""),
        "gc_profiles": aggregate_input.get("gc_profiles", ""),
        "aggregate_flt_silhouette": locations.get("FLT", {}).get("consensus_silhouette", ""),
        "aggregate_gc_silhouette": locations.get("GC", {}).get("consensus_silhouette", ""),
        "per_study_done": per_study_done,
        "per_study_total": per_study_total,
        "dgea_done": dgea_done,
        "mober_done": (tissue_dir / "aggregate_mober" / "clustering" / "clustering_summary.json").exists(),
        "top_recurring_dgea_glare_terms": top_terms,
    }


def write_root_summary(root: Path) -> Path:
    rows = []
    for tissue_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        if (tissue_dir / "aggregate" / "input_summary.json").exists():
            rows.append(tissue_summary_row(tissue_dir))
    summary = pd.DataFrame(rows)
    summary_path = root / "MULTI_TISSUE_SUMMARY.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    lines = [
        "# Multi-Tissue API GLARE Summary",
        "",
        "All scopes use NASA OSDR API-derived expression/count inputs. GLARE uses log2(CPM+1) expression aligned to the matching Tabula Muris Senis FACS tissue where available; DESeq2 uses the matched raw-count inputs.",
        "",
        "## Run Status",
        "",
        "| tissue | FLT/GC profiles | studies | per-study GLARE | DGEA | MOBER | aggregate silhouette FLT/GC | top recurring DGEA/GLARE pathways |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        profiles = f"{row['flt_profiles']}/{row['gc_profiles']}"
        per_study = f"{row['per_study_done']}/{row['per_study_total']}"
        silhouette = (
            f"{format_float(row['aggregate_flt_silhouette'])}/"
            f"{format_float(row['aggregate_gc_silhouette'])}"
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["tissue_label"]),
                    profiles,
                    str(row["studies"]),
                    per_study,
                    "yes" if row["dgea_done"] else "no",
                    "yes" if row["mober_done"] else "no",
                    silhouette,
                    short_text(row["top_recurring_dgea_glare_terms"], 120),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- The PCA/UMAP points are genes, not samples. The plot subtitles list the number of source FLT or GC profiles used to learn that gene representation.",
            "- Cluster descriptions are evidence labels: DGEA-enriched clusters have significant same-study DEG over-representation; Reactome-enriched clusters have significant pathway over-representation; ambiguous clusters lack either label at FDR < 0.05.",
            "- Per-study DESeq2 recurrence is the primary FLT-vs-GC evidence. GLARE is used as module discovery that supports the DGEA signal when the same pathways recur in both layers.",
            "- Retina is skipped unless a matching TMS FACS pretraining tissue becomes available.",
        ]
    )
    md_path = root / "MULTI_TISSUE_SUMMARY.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return md_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=DEFAULT_ROOT)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--min-overlap", type=int, default=3)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("plot-run")
    run.add_argument("--run-dir", required=True)
    run.add_argument("--skip-umap", action="store_true")

    tissue = subparsers.add_parser("plot-tissue")
    tissue.add_argument("--tissue-dir", required=True)
    tissue.add_argument("--include-per-study", action="store_true")
    tissue.add_argument("--include-mober", action="store_true")
    tissue.add_argument("--skip-umap", action="store_true")

    root = subparsers.add_parser("plot-root")
    root.add_argument("--include-per-study", action="store_true")
    root.add_argument("--include-mober", action="store_true")
    root.add_argument("--skip-umap", action="store_true")

    subparsers.add_parser("summarize-root")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "plot-run":
        output = plot_run(Path(args.run_dir), args)
    elif args.command == "plot-tissue":
        output = plot_tissue(Path(args.tissue_dir), args)
    elif args.command == "plot-root":
        root = Path(args.root)
        for tissue_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            if (tissue_dir / "aggregate" / "input_summary.json").exists():
                plot_tissue(tissue_dir, args)
        output = root
    elif args.command == "summarize-root":
        output = write_root_summary(Path(args.root))
    else:
        raise ValueError(args.command)
    print(output)


if __name__ == "__main__":
    main()
