"""Map OSDR tissue data into a trained ARCHS4 expiMap reference model."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from .io import require_import


def patch_stable_expimap_hsic() -> None:
    """Patch scArches expiMap HSIC for high-dimensional latent spaces.

    scArches 0.6.1 computes 2 * gamma((d + 1) / 2) / gamma(d / 2) directly.
    That overflows for the 300+ dimensional latent used by the expiMap
    tutorial-style model here. The log-gamma form is mathematically equivalent
    and stays finite.
    """
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    losses = require_import(
        "scarches.models.expimap.losses",
        "pip install -r requirements-nasa-mouse-glare.txt",
    )
    expimap_module = require_import(
        "scarches.models.expimap.expimap",
        "pip install -r requirements-nasa-mouse-glare.txt",
    )

    def stable_bandwidth(d):
        gz = math.exp(math.log(2.0) + math.lgamma(0.5 * (d + 1)) - math.lgamma(0.5 * d))
        return 1.0 / (2.0 * gz**2)

    def stable_hsic(z, s):
        zz = losses.kernel_matrix(z, stable_bandwidth(z.shape[1]))
        ss = losses.kernel_matrix(s, stable_bandwidth(s.shape[1]))
        h = (zz * ss).mean() + zz.mean() * ss.mean() - 2 * torch.mean(
            zz.mean(1) * ss.mean(1)
        )
        return torch.sqrt(torch.clamp(h, min=0.0))

    losses.bandwidth = stable_bandwidth
    losses.hsic = stable_hsic
    expimap_module.hsic = stable_hsic


def update_terms_for_extensions(model, adata) -> tuple[list[str], list[str]]:
    """Append stable labels and descriptions for unconstrained query programs."""
    model.update_terms("terms")
    terms = list(map(str, adata.uns["terms"]))
    descriptions = list(map(str, adata.uns.get("term_descriptions", terms)))
    if len(descriptions) > len(terms):
        descriptions = descriptions[: len(terms)]
    if len(descriptions) < len(terms):
        descriptions.extend(
            f"Unconstrained de novo expiMap program {term.rsplit('_', 1)[-1]}"
            if term.startswith("unconstrained_")
            else term
            for term in terms[len(descriptions) :]
        )
    adata.uns["term_descriptions"] = descriptions
    return terms, descriptions


def write_de_novo_gene_loadings(model, terms: list[str], output_dir: Path, top_n: int):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    rows = []
    summary_rows = []
    de_novo_terms = [term for term in terms if term.startswith("unconstrained_")]
    for term in de_novo_terms:
        genes = model.term_genes(term, terms=terms).copy()
        genes["absolute_weight"] = genes["weights"].abs()
        genes = genes.sort_values("absolute_weight", ascending=False, kind="stable")
        summary_rows.append(
            {
                "term": term,
                "n_nonzero_genes": int(len(genes)),
                "top_genes": ";".join(genes["genes"].head(top_n).astype(str)),
            }
        )
        for rank, (_, row) in enumerate(genes.head(top_n).iterrows(), start=1):
            rows.append(
                {
                    "term": term,
                    "rank": rank,
                    "gene": str(row["genes"]),
                    "decoder_weight": float(row["weights"]),
                    "absolute_weight": float(row["absolute_weight"]),
                }
            )

    loadings_path = output_dir / "de_novo_program_gene_loadings.tsv"
    summary_path = output_dir / "de_novo_programs.tsv"
    pd.DataFrame(rows).to_csv(loadings_path, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    return de_novo_terms, {
        "de_novo_programs": str(summary_path),
        "de_novo_gene_loadings": str(loadings_path),
    }


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    sca = require_import("scarches", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")

    if args.use_hsic and not args.no_stable_hsic_patch:
        patch_stable_expimap_hsic()

    query = ad.read_h5ad(args.query_h5ad)
    if args.query_condition_source not in query.obs:
        raise SystemExit(
            f"query_condition_source {args.query_condition_source!r} not in query obs."
        )
    query.obs[args.reference_condition_key] = query.obs[
        args.query_condition_source
    ].astype(str)
    if args.recon_loss == "nb":
        if "counts" not in query.layers:
            raise SystemExit("NB query mapping requires query.layers['counts'].")
        query.X = query.layers["counts"].copy()

    new_n_ext = args.n_de_novo_programs or None
    query_load_kwargs = {}
    if args.use_hsic:
        query_load_kwargs["hsic_one_vs_all"] = args.hsic_one_vs_all
    model = sca.models.EXPIMAP.load_query_data(
        query,
        reference_model=args.reference_model,
        freeze=not args.no_freeze,
        freeze_expression=not args.no_freeze_expression,
        remove_dropout=True,
        new_n_ext=new_n_ext,
        use_hsic=args.use_hsic,
        **query_load_kwargs,
    )
    model.train(
        n_epochs=args.epochs,
        lr=args.learning_rate,
        alpha=None if args.no_alpha else args.alpha,
        alpha_kl=args.alpha_kl,
        alpha_epoch_anneal=args.alpha_epoch_anneal,
        weight_decay=args.weight_decay,
        gamma_ext=args.gamma_ext if new_n_ext else None,
        gamma_epoch_anneal=args.gamma_epoch_anneal if new_n_ext else None,
        beta=args.beta if args.use_hsic and new_n_ext else None,
        seed=args.seed,
        use_early_stopping=args.early_stopping,
        monitor=not args.no_monitor,
    )

    latent = model.get_latent(mean=args.mean_latent)
    if isinstance(latent, tuple):
        latent = latent[0]
    query_mapped = model.adata
    terms, descriptions = update_terms_for_extensions(model, query_mapped)
    if len(terms) != latent.shape[1]:
        raise SystemExit(
            f"Latent dimension ({latent.shape[1]}) does not match term labels ({len(terms)})."
        )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    query_mapped.obsm["X_expimap_query"] = latent
    mapped_h5ad = output_dir / "mapped_query_with_scores.h5ad"
    query_mapped.write_h5ad(mapped_h5ad)

    scores = pd.concat(
        [
            query_mapped.obs.reset_index(names="obs_name"),
            pd.DataFrame(latent, columns=terms),
        ],
        axis=1,
    )
    scores_path = output_dir / "query_pathway_scores.tsv"
    scores.to_csv(scores_path, sep="\t", index=False)
    de_novo_terms, de_novo_outputs = write_de_novo_gene_loadings(
        model,
        terms,
        output_dir,
        args.de_novo_top_genes,
    )

    terms_path = output_dir / "terms.tsv"
    pd.DataFrame({"term": terms, "description": descriptions}).to_csv(
        terms_path,
        sep="\t",
        index=False,
    )

    model_dir = output_dir / "query_model"
    model.save(str(model_dir), overwrite=True, save_anndata=False)
    summary = {
        "reference_model": str(args.reference_model),
        "query_h5ad": str(args.query_h5ad),
        "reference_condition_key": args.reference_condition_key,
        "query_condition_source": args.query_condition_source,
        "epochs": args.epochs,
        "n_query_samples": int(query_mapped.n_obs),
        "n_query_genes": int(query_mapped.n_vars),
        "n_score_terms": int(latent.shape[1]),
        "posterior_mean_latent": bool(args.mean_latent),
        "torch": {
            "version": str(torch.__version__),
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_count": int(torch.cuda.device_count()),
            "cuda_device_name": (
                torch.cuda.get_device_name(0) if torch.cuda.is_available() else ""
            ),
            "trainer_device": str(getattr(model.trainer, "device", "")),
            "model_parameter_device": str(next(model.model.parameters()).device),
        },
        "alpha": None if args.no_alpha else args.alpha,
        "n_de_novo_programs": int(len(de_novo_terms)),
        "de_novo_programs": de_novo_terms,
        "use_hsic": bool(args.use_hsic and new_n_ext),
        "hsic_one_vs_all": bool(args.hsic_one_vs_all if args.use_hsic and new_n_ext else False),
        "stable_hsic_patch": bool(args.use_hsic and new_n_ext and not args.no_stable_hsic_patch),
        "gamma_ext": args.gamma_ext if new_n_ext else None,
        "outputs": {
            "query_model": str(model_dir),
            "query_scores": str(scores_path),
            "mapped_h5ad": str(mapped_h5ad),
            "terms": str(terms_path),
            **de_novo_outputs,
        },
    }
    summary_path = output_dir / "query_mapping_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map an OSDR query AnnData into a trained ARCHS4 expiMap reference."
    )
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--reference-condition-key", default="archs4_condition")
    parser.add_argument("--query-condition-source", default="id.accession")
    parser.add_argument("--recon-loss", choices=["nb"], default="nb")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--alpha", type=float, default=0.7)
    parser.add_argument(
        "--no-alpha",
        action="store_true",
        help="Do not pass group-lasso alpha during query training.",
    )
    parser.add_argument("--alpha-kl", type=float, default=0.35)
    parser.add_argument("--alpha-epoch-anneal", type=int, default=10)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument(
        "--n-de-novo-programs",
        type=int,
        default=0,
        help="Add this many unconstrained, query-specific latent programs.",
    )
    parser.add_argument(
        "--gamma-ext",
        type=float,
        default=0.7,
        help="L1 sparsity coefficient for unconstrained program gene weights.",
    )
    parser.add_argument(
        "--gamma-epoch-anneal",
        type=int,
        default=50,
        help="Epochs over which unconstrained-program sparsity is annealed.",
    )
    parser.add_argument(
        "--use-hsic",
        action="store_true",
        help="Encourage de novo programs to be independent of annotated programs.",
    )
    parser.add_argument(
        "--hsic-one-vs-all",
        action="store_true",
        help="Use the tutorial-style one-vs-all HSIC penalty for extension nodes.",
    )
    parser.add_argument(
        "--no-stable-hsic-patch",
        action="store_true",
        help="Disable the repo runtime patch for scArches' high-dimensional HSIC overflow.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=3.0,
        help="HSIC regularization coefficient when --use-hsic is enabled.",
    )
    parser.add_argument(
        "--de-novo-top-genes",
        type=int,
        default=30,
        help="Number of largest-magnitude gene weights to export per de novo program.",
    )
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--no-monitor", action="store_true")
    parser.set_defaults(mean_latent=True)
    parser.add_argument(
        "--mean-latent",
        dest="mean_latent",
        action="store_true",
        help="Write posterior-mean latent scores (the default).",
    )
    parser.add_argument(
        "--sample-latent",
        dest="mean_latent",
        action="store_false",
        help="Write one stochastic latent sample; use only for diagnostics.",
    )
    parser.add_argument("--no-freeze", action="store_true")
    parser.add_argument("--no-freeze-expression", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
