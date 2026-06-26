"""Map OSDR tissue data into a trained ARCHS4 expiMap reference model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    sca = require_import("scarches", "pip install -r requirements-nasa-mouse-glare.txt")

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

    model = sca.models.EXPIMAP.load_query_data(
        query,
        reference_model=args.reference_model,
        freeze=not args.no_freeze,
        freeze_expression=not args.no_freeze_expression,
        remove_dropout=True,
    )
    model.train(
        n_epochs=args.epochs,
        lr=args.learning_rate,
        alpha=args.alpha,
        alpha_kl=args.alpha_kl,
        alpha_epoch_anneal=args.alpha_epoch_anneal,
        weight_decay=args.weight_decay,
        seed=args.seed,
        use_early_stopping=args.early_stopping,
    )

    latent = model.get_latent(mean=args.mean_latent)
    if isinstance(latent, tuple):
        latent = latent[0]
    query_mapped = model.adata
    terms = list(map(str, query_mapped.uns.get("terms", [f"term_{i}" for i in range(latent.shape[1])])))[: latent.shape[1]]

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
        "outputs": {
            "query_model": str(model_dir),
            "query_scores": str(scores_path),
            "mapped_h5ad": str(mapped_h5ad),
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
    parser.add_argument("--alpha-kl", type=float, default=0.35)
    parser.add_argument("--alpha-epoch-anneal", type=int, default=10)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--mean-latent", action="store_true")
    parser.add_argument("--no-freeze", action="store_true")
    parser.add_argument("--no-freeze-expression", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
