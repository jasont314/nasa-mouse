"""Summarize query-specific expiMap programs and their relation to Reactome terms."""

from __future__ import annotations

import argparse
from pathlib import Path

from .cluster_enrichment import bh_fdr
from .inspect_archs4_mouse import DEFAULT_ARCHS4
from .io import require_import


def decode(values) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def load_gene_symbols(path: str | Path) -> dict[str, str]:
    h5py = require_import("h5py", "pip install -r requirements-nasa-mouse-glare.txt")
    with h5py.File(path, "r") as handle:
        ensembl = decode(handle["/meta/genes/ensembl_gene"][:])
        symbols = decode(handle["/meta/genes/symbol"][:])
    return {
        gene: symbol
        for gene, symbol in zip(ensembl, symbols)
        if gene and symbol and symbol != "nan"
    }


def program_enrichment(adata, program_genes: list[str], top_n: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    hypergeom = require_import("scipy.stats", "pip install -r requirements-nasa-mouse-glare.txt").hypergeom

    mask = np.asarray(adata.varm["I"], dtype=bool)
    term_labels = [str(term) for term in adata.uns["terms"][: mask.shape[1]]]
    if mask.shape[1] != len(term_labels):
        raise SystemExit("Reactome mask columns and term labels are not aligned.")
    reactome_indices = [
        index for index, term in enumerate(term_labels) if term.startswith("R-MMU-")
    ]
    reactome_terms = [term_labels[index] for index in reactome_indices]
    mask = mask[:, reactome_indices]
    gene_to_index = {gene: index for index, gene in enumerate(adata.var_names.astype(str))}
    indices = [gene_to_index[gene] for gene in program_genes[:top_n] if gene in gene_to_index]
    if not indices:
        return pd.DataFrame()

    n_genes = int(mask.shape[0])
    sample_size = int(len(indices))
    overlaps = mask[indices].sum(axis=0)
    pathway_sizes = mask.sum(axis=0)
    p_values = hypergeom.sf(overlaps - 1, n_genes, pathway_sizes, sample_size)
    result = pd.DataFrame(
        {
            "reactome_term": reactome_terms,
            "overlap_genes": overlaps.astype(int),
            "pathway_genes": pathway_sizes.astype(int),
            "enrichment_p": p_values,
        }
    )
    result["enrichment_fdr"] = bh_fdr(result["enrichment_p"].to_numpy())
    if not int(result["overlap_genes"].max()):
        return pd.DataFrame(columns=result.columns)
    return result.sort_values(
        ["enrichment_fdr", "overlap_genes", "reactome_term"],
        ascending=[True, False, True],
        kind="stable",
    )


def run(args) -> Path:
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.mapped_h5ad)
    scores = pd.read_csv(args.scores, sep="\t")
    comparison = pd.read_csv(args.comparison, sep="\t")
    study_tests = pd.read_csv(args.study_tests, sep="\t")
    programs = pd.read_csv(args.programs, sep="\t")
    loadings = pd.read_csv(args.gene_loadings, sep="\t")
    symbols = load_gene_symbols(args.archs4_h5)

    de_novo_terms = [column for column in scores if str(column).startswith("unconstrained_")]
    reactome_terms = [column for column in scores if str(column).startswith("R-MMU-")]
    loading_rows = []
    summary_rows = []
    enrichment_rows = []
    comparison_by_term = comparison.set_index("term")
    study_by_term = study_tests.set_index("term")
    program_sizes = programs.set_index("term")["n_nonzero_genes"].to_dict()

    for term in de_novo_terms:
        genes = loadings.loc[loadings["term"].eq(term)].copy()
        genes["gene_symbol"] = genes["gene"].map(symbols).fillna("")
        genes.to_csv(
            output_dir / f"{term}_gene_loadings.tsv",
            sep="\t",
            index=False,
        )
        loading_rows.append(genes)
        enrichment = program_enrichment(
            adata,
            genes["gene"].astype(str).tolist(),
            args.top_genes,
        )
        if not enrichment.empty:
            enrichment.insert(0, "term", term)
            enrichment_rows.append(enrichment)
            best_enrichment = enrichment.iloc[0]
        else:
            best_enrichment = {}

        correlations = scores[reactome_terms].corrwith(scores[term].astype(float))
        correlations = correlations.dropna()
        if correlations.empty:
            strongest_term, strongest_corr = "", float("nan")
        else:
            strongest_term = str(correlations.abs().idxmax())
            strongest_corr = float(correlations.loc[strongest_term])
        test = comparison_by_term.loc[term] if term in comparison_by_term.index else {}
        study = study_by_term.loc[term] if term in study_by_term.index else {}
        top_gene_labels = [
            f"{row.gene_symbol or row.gene} ({row.decoder_weight:.3g})"
            for row in genes.head(args.top_genes).itertuples(index=False)
        ]
        summary_rows.append(
            {
                "term": term,
                "score_std": float(np.nanstd(scores[term].astype(float))),
                "n_nonzero_gene_weights": int(program_sizes.get(term, len(genes))),
                "flight_minus_ground": test.get("flight_minus_ground", float("nan")),
                "welch_p": test.get("welch_p", float("nan")),
                "welch_fdr": test.get("welch_fdr", float("nan")),
                "study_mean_accession_effect": study.get(
                    "mean_accession_effect", float("nan")
                ),
                "study_wilcoxon_fdr": study.get("wilcoxon_fdr", float("nan")),
                "top_reactome_gene_set": best_enrichment.get("reactome_term", ""),
                "top_reactome_overlap_genes": best_enrichment.get("overlap_genes", 0),
                "top_reactome_enrichment_fdr": best_enrichment.get(
                    "enrichment_fdr", float("nan")
                ),
                "strongest_score_correlated_reactome_term": strongest_term,
                "strongest_score_correlation": strongest_corr,
                "top_weighted_genes": "; ".join(top_gene_labels),
            }
        )

    annotated_loadings = pd.concat(loading_rows, ignore_index=True) if loading_rows else pd.DataFrame()
    enrichment_table = pd.concat(enrichment_rows, ignore_index=True) if enrichment_rows else pd.DataFrame()
    summary = pd.DataFrame(summary_rows).sort_values("welch_fdr", kind="stable")
    summary_path = output_dir / "de_novo_program_summary.tsv"
    loadings_path = output_dir / "de_novo_program_gene_loadings_annotated.tsv"
    enrichment_path = output_dir / "de_novo_reactome_gene_set_enrichment.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    annotated_loadings.to_csv(loadings_path, sep="\t", index=False)
    enrichment_table.to_csv(enrichment_path, sep="\t", index=False)

    lines = [
        "# De Novo expiMap Programs",
        "",
        "The program labels are arbitrary latent-dimension identifiers. Gene weights are decoder loadings, and the Reactome association is a post-hoc enrichment of the highest-magnitude weighted genes.",
        "",
        f"- Programs tested: {len(summary)}",
        f"- Programs with aggregate Welch FDR < 0.05: {int((summary['welch_fdr'].astype(float) < 0.05).sum()) if not summary.empty else 0}",
        f"- Programs with study-aware FDR < 0.05: {int((summary['study_wilcoxon_fdr'].astype(float) < 0.05).sum()) if not summary.empty else 0}",
        "",
        "| program | score SD | FLT-GC | Welch FDR | study FDR | strongest Reactome score correlation | top weighted genes |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "| {term} | {score_std:.3g} | {effect:.3g} | {welch_fdr:.3g} | {study_fdr:.3g} | {corr_term} ({corr:.3g}) | {genes} |".format(
                term=row["term"],
                score_std=float(row["score_std"]),
                effect=float(row["flight_minus_ground"]),
                welch_fdr=float(row["welch_fdr"]),
                study_fdr=float(row["study_wilcoxon_fdr"]),
                corr_term=row["strongest_score_correlated_reactome_term"],
                corr=float(row["strongest_score_correlation"]),
                genes=row["top_weighted_genes"],
            )
        )
    readme_path = output_dir / "README.md"
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_csv(sep="\t", index=False))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize de novo expiMap query programs and Reactome associations."
    )
    parser.add_argument("--mapped-h5ad", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--study-tests", required=True)
    parser.add_argument("--programs", required=True)
    parser.add_argument("--gene-loadings", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--archs4-h5", default=DEFAULT_ARCHS4)
    parser.add_argument("--top-genes", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
