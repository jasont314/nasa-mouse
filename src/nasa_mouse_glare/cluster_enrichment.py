"""Gene-cluster enrichment and driver summaries for post-finetune GLARE output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import load_matrix_bundle, require_import


DEFAULT_POST_DIR = "outputs/glare_hpt_tms_facs_osdr/post_finetune"
DEFAULT_TARGET_MANIFEST = "data/processed/tms_facs_osdr_aligned.target.manifest.json"
DEFAULT_CLUSTERS = [13, 10, 8, 0, 6]
DEFAULT_REACTOME_GMT = (
    "src/expiMap_reproducibility/metadata/c2.cp.reactome.v4.0_mouseEID.gmt"
)
DEFAULT_PANGLAO_GMT = (
    "src/expiMap_reproducibility/metadata/"
    "PanglaoDB_markers_27_Mar_2020_mouseEID.gmt"
)


def read_gmt(path: str | Path) -> list[dict[str, object]]:
    gene_sets = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, description, *genes = parts
            gene_sets.append(
                {
                    "term": name,
                    "description": description,
                    "genes": {gene for gene in genes if gene},
                }
            )
    return gene_sets


def bh_fdr(pvalues):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    values = np.asarray(pvalues, dtype=float)
    if len(values) == 0:
        return values
    order = np.argsort(values)
    ordered = values[order]
    adjusted_ordered = np.empty_like(ordered)
    running_min = 1.0
    n_tests = len(values)
    for idx in range(n_tests - 1, -1, -1):
        rank = idx + 1
        running_min = min(running_min, ordered[idx] * n_tests / rank)
        adjusted_ordered[idx] = running_min
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(adjusted_ordered, 1.0)
    return adjusted


def load_focus_tables(post_dir: Path, clusters: list[int]):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    cluster_path = post_dir / "gene_clusters.tsv"
    shift_path = post_dir / "gene_cluster_flight_ground_summary.tsv"
    condition_path = post_dir / "gene_cluster_expression_by_condition_inferred.tsv"
    study_path = post_dir / "gene_cluster_expression_by_id.accession.tsv"
    profile_path = post_dir / "profile_metadata.tsv"

    for path in [cluster_path, shift_path, condition_path, study_path, profile_path]:
        if not path.exists():
            raise SystemExit(f"Missing required post-finetune output: {path}")

    gene_clusters = pd.read_csv(cluster_path, sep="\t")
    shift_summary = pd.read_csv(shift_path, sep="\t")
    condition_summary = pd.read_csv(condition_path, sep="\t")
    study_summary = pd.read_csv(study_path, sep="\t")
    profile_metadata = pd.read_csv(profile_path, sep="\t", keep_default_na=False)

    for frame in [gene_clusters, shift_summary, condition_summary, study_summary]:
        frame["gene_cluster"] = frame["gene_cluster"].astype(int)

    missing = sorted(set(clusters) - set(gene_clusters["gene_cluster"].unique()))
    if missing:
        raise SystemExit(f"Requested clusters not found in gene_clusters.tsv: {missing}")

    return gene_clusters, shift_summary, condition_summary, study_summary, profile_metadata


def write_gene_lists(gene_clusters, clusters: list[int], output_dir: Path):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    paths = {}
    for cluster in clusters:
        genes = (
            gene_clusters.loc[gene_clusters["gene_cluster"] == cluster, "gene"]
            .astype(str)
            .sort_values()
            .tolist()
        )
        path = output_dir / f"cluster_{cluster}_genes.txt"
        path.write_text("\n".join(genes) + "\n", encoding="utf-8")
        paths[str(cluster)] = str(path)
        rows.extend({"gene_cluster": cluster, "gene": gene} for gene in genes)

    combined = pd.DataFrame(rows)
    combined_path = output_dir / "focus_cluster_genes.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    return paths | {"combined": str(combined_path)}


def run_enrichment(
    gene_clusters,
    clusters: list[int],
    library_name: str,
    gmt_path: str | Path,
    min_overlap: int,
):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    hypergeom = require_import(
        "scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt"
    ).hypergeom

    library = read_gmt(gmt_path)
    all_genes = set(gene_clusters["gene"].astype(str))
    annotated_genes = set()
    for gene_set in library:
        annotated_genes.update(gene_set["genes"])
    universe = all_genes & annotated_genes
    universe_size = len(universe)
    if universe_size == 0:
        raise SystemExit(f"No overlap between clustered genes and {gmt_path}")

    rows = []
    for cluster in clusters:
        query_genes = set(
            gene_clusters.loc[gene_clusters["gene_cluster"] == cluster, "gene"].astype(str)
        )
        query = query_genes & universe
        query_size = len(query)
        if query_size == 0:
            continue

        for gene_set in library:
            term_genes = gene_set["genes"] & universe
            term_size = len(term_genes)
            if term_size == 0:
                continue
            overlap = query & term_genes
            overlap_size = len(overlap)
            if overlap_size < min_overlap:
                continue

            pvalue = float(
                hypergeom.sf(overlap_size - 1, universe_size, term_size, query_size)
            )
            enrichment_ratio = (overlap_size / query_size) / (term_size / universe_size)
            rows.append(
                {
                    "library": library_name,
                    "gene_cluster": cluster,
                    "term": gene_set["term"],
                    "description": gene_set["description"],
                    "overlap": overlap_size,
                    "cluster_genes_in_universe": query_size,
                    "term_genes_in_universe": term_size,
                    "universe_genes": universe_size,
                    "p_value": pvalue,
                    "enrichment_ratio": enrichment_ratio,
                    "overlap_genes": ",".join(sorted(overlap)),
                }
            )

    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(
            columns=[
                "library",
                "gene_cluster",
                "term",
                "description",
                "overlap",
                "cluster_genes_in_universe",
                "term_genes_in_universe",
                "universe_genes",
                "p_value",
                "fdr_bh",
                "enrichment_ratio",
                "overlap_genes",
            ]
        )

    result["fdr_bh"] = result.groupby(["library", "gene_cluster"], group_keys=False)[
        "p_value"
    ].transform(bh_fdr)
    return result.sort_values(
        ["gene_cluster", "fdr_bh", "p_value", "overlap"],
        ascending=[True, True, True, False],
    )


def write_cluster_shift_summary(shift_summary, clusters: list[int], output_dir: Path) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    order = pd.Series(range(len(clusters)), index=clusters)
    focus = shift_summary[shift_summary["gene_cluster"].isin(clusters)].copy()
    focus["requested_order"] = focus["gene_cluster"].map(order)
    focus = focus.sort_values("requested_order").drop(columns=["requested_order"])
    path = output_dir / "cluster_shift_summary.tsv"
    focus.to_csv(path, sep="\t", index=False)
    return str(path)


def write_condition_driver_summary(condition_summary, clusters: list[int], output_dir: Path) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    focus = condition_summary[condition_summary["gene_cluster"].isin(clusters)].copy()
    value_cols = [col for col in focus.columns if col not in {"gene_cluster", "n_genes"}]
    rows = []
    for _, row in focus.iterrows():
        values = row[value_cols].astype(float)
        for condition, value in values.sort_values(ascending=False).items():
            rows.append(
                {
                    "gene_cluster": int(row["gene_cluster"]),
                    "n_genes": int(row["n_genes"]),
                    "condition": condition,
                    "mean_expression": float(value),
                    "condition_rank_desc": int(values.rank(ascending=False)[condition]),
                    "delta_vs_cluster_condition_mean": float(value - values.mean()),
                }
            )
    result = pd.DataFrame(rows)
    path = output_dir / "condition_driver_summary.tsv"
    result.to_csv(path, sep="\t", index=False)
    return str(path)


def write_study_driver_summary(study_summary, clusters: list[int], top_n: int, output_dir: Path) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    study_cols = [col for col in study_summary.columns if str(col).startswith("OSD-")]
    rows = []
    for _, row in study_summary[study_summary["gene_cluster"].isin(clusters)].iterrows():
        values = row[study_cols].astype(float)
        std = float(values.std(ddof=0))
        zscores = (values - values.mean()) / std if std else values * 0.0
        high = zscores.sort_values(ascending=False).head(top_n)
        low = zscores.sort_values(ascending=True).head(top_n)
        for direction, selected in [("high", high), ("low", low)]:
            for rank, (accession, zscore) in enumerate(selected.items(), start=1):
                rows.append(
                    {
                        "gene_cluster": int(row["gene_cluster"]),
                        "n_genes": int(row["n_genes"]),
                        "direction": direction,
                        "rank": rank,
                        "id.accession": accession,
                        "mean_expression": float(values[accession]),
                        "z_score_across_studies": float(zscore),
                    }
                )
    result = pd.DataFrame(rows)
    path = output_dir / "study_driver_summary.tsv"
    result.to_csv(path, sep="\t", index=False)
    return str(path)


def write_source_driver_summary(
    target_manifest: str | Path,
    profile_metadata,
    gene_clusters,
    clusters: list[int],
    top_n: int,
    output_dir: Path,
) -> str | None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    if "study.source name" not in profile_metadata:
        return None
    target_manifest = Path(target_manifest)
    if not target_manifest.exists():
        return None

    bundle = load_matrix_bundle(target_manifest)
    if list(map(str, bundle.genes)) != gene_clusters["gene"].astype(str).tolist():
        raise SystemExit(
            "Target manifest genes do not match gene_clusters.tsv order. "
            "Rerun post_finetune.py before source-level driver analysis."
        )
    if list(map(str, bundle.profiles)) != profile_metadata["profile"].astype(str).tolist():
        raise SystemExit(
            "Target manifest profiles do not match profile_metadata.tsv order. "
            "Rerun post_finetune.py before source-level driver analysis."
        )

    source = (
        profile_metadata["study.source name"]
        .fillna("NA")
        .astype(str)
        .replace("", "NA")
    )
    source_values = sorted(source.unique().tolist())
    rows = []
    for cluster in clusters:
        gene_mask = gene_clusters["gene_cluster"].to_numpy() == cluster
        cluster_matrix = bundle.matrix[gene_mask, :]
        sample_means = np.asarray(cluster_matrix.mean(axis=0)).ravel()
        if scipy_sparse.issparse(cluster_matrix):
            sample_means = np.asarray(cluster_matrix.mean(axis=0)).ravel()

        grouped = pd.Series(sample_means).groupby(source, sort=True).mean()
        grouped = grouped.reindex(source_values).dropna()
        std = float(grouped.std(ddof=0))
        zscores = (grouped - grouped.mean()) / std if std else grouped * 0.0
        counts = source.value_counts()

        high = zscores.sort_values(ascending=False).head(top_n)
        low = zscores.sort_values(ascending=True).head(top_n)
        for direction, selected in [("high", high), ("low", low)]:
            for rank, (source_name, zscore) in enumerate(selected.items(), start=1):
                rows.append(
                    {
                        "gene_cluster": cluster,
                        "n_genes": int(gene_mask.sum()),
                        "direction": direction,
                        "rank": rank,
                        "study.source name": source_name,
                        "n_profiles": int(counts[source_name]),
                        "mean_expression": float(grouped[source_name]),
                        "z_score_across_sources": float(zscore),
                    }
                )

    result = pd.DataFrame(rows)
    path = output_dir / "source_driver_summary.tsv"
    result.to_csv(path, sep="\t", index=False)
    return str(path)


def write_top_enrichment_summary(enrichment_tables, top_n: int, output_dir: Path) -> str:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    combined = pd.concat(enrichment_tables, ignore_index=True)
    if combined.empty:
        top = combined
    else:
        top = (
            combined.sort_values(
                ["library", "gene_cluster", "fdr_bh", "p_value", "overlap"],
                ascending=[True, True, True, True, False],
            )
            .groupby(["library", "gene_cluster"], group_keys=False)
            .head(top_n)
        )
    path = output_dir / "top_enrichment_summary.tsv"
    top.to_csv(path, sep="\t", index=False)
    return str(path)


def run(args) -> Path:
    output_dir = Path(args.output_dir) if args.output_dir else Path(args.post_dir) / "enrichment"
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters = [int(cluster) for cluster in args.clusters]
    (
        gene_clusters,
        shift_summary,
        condition_summary,
        study_summary,
        profile_metadata,
    ) = load_focus_tables(Path(args.post_dir), clusters)

    gene_list_paths = write_gene_lists(gene_clusters, clusters, output_dir)
    cluster_shift_path = write_cluster_shift_summary(shift_summary, clusters, output_dir)
    condition_driver_path = write_condition_driver_summary(
        condition_summary, clusters, output_dir
    )
    study_driver_path = write_study_driver_summary(
        study_summary, clusters, args.top_n, output_dir
    )
    source_driver_path = write_source_driver_summary(
        args.target_manifest,
        profile_metadata,
        gene_clusters,
        clusters,
        args.top_n,
        output_dir,
    )

    reactome = run_enrichment(
        gene_clusters,
        clusters,
        "reactome",
        args.reactome_gmt,
        args.min_overlap,
    )
    reactome_path = output_dir / "reactome_enrichment.tsv"
    reactome.to_csv(reactome_path, sep="\t", index=False)

    panglao = run_enrichment(
        gene_clusters,
        clusters,
        "panglao",
        args.panglao_gmt,
        args.min_overlap,
    )
    panglao_path = output_dir / "panglao_enrichment.tsv"
    panglao.to_csv(panglao_path, sep="\t", index=False)

    top_enrichment_path = write_top_enrichment_summary(
        [reactome, panglao], args.top_n, output_dir
    )

    summary = {
        "post_dir": str(args.post_dir),
        "target_manifest": str(args.target_manifest),
        "clusters": clusters,
        "min_overlap": args.min_overlap,
        "top_n": args.top_n,
        "gmt_files": {
            "reactome": str(args.reactome_gmt),
            "panglao": str(args.panglao_gmt),
        },
        "outputs": {
            "gene_lists": gene_list_paths,
            "cluster_shift_summary": cluster_shift_path,
            "condition_driver_summary": condition_driver_path,
            "study_driver_summary": study_driver_path,
            "source_driver_summary": source_driver_path,
            "reactome_enrichment": str(reactome_path),
            "panglao_enrichment": str(panglao_path),
            "top_enrichment_summary": top_enrichment_path,
        },
    }
    summary_path = output_dir / "cluster_enrichment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"clusters={','.join(map(str, clusters))}")
    print(f"reactome_terms={len(reactome)}")
    print(f"panglao_terms={len(panglao)}")
    print(f"summary={summary_path}")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run enrichment and driver summaries for GLARE gene clusters."
    )
    parser.add_argument("--post-dir", default=DEFAULT_POST_DIR)
    parser.add_argument("--output-dir", default="")
    parser.add_argument(
        "--target-manifest",
        default=DEFAULT_TARGET_MANIFEST,
        help="Aligned target manifest; used only for study.source name driver summary.",
    )
    parser.add_argument("--clusters", nargs="+", type=int, default=DEFAULT_CLUSTERS)
    parser.add_argument("--reactome-gmt", default=DEFAULT_REACTOME_GMT)
    parser.add_argument("--panglao-gmt", default=DEFAULT_PANGLAO_GMT)
    parser.add_argument("--min-overlap", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=20)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
