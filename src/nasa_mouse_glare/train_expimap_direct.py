"""Train direct expiMap models on prepared OSDR tissue AnnData inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import require_import


def parse_hidden_layers(value: str) -> list[int]:
    return [int(token) for token in str(value).split(",") if token.strip()]


def infer_recon_loss(adata, requested: str) -> str:
    if requested != "auto":
        return requested
    preprocessing = adata.uns.get("expimap_preprocessing", {})
    return str(preprocessing.get("recommended_recon_loss", "nb"))


def term_columns(terms, latent, active_indices=None):
    if active_indices is None:
        return list(map(str, terms))
    return [str(terms[index]) for index in active_indices]


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    sca = require_import("scarches", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")

    adata = ad.read_h5ad(args.input)
    if "I" not in adata.varm:
        raise SystemExit("Input AnnData is missing adata.varm['I'].")
    if adata.varm["I"].shape[0] != adata.n_vars:
        raise SystemExit(
            f"adata.varm['I'] must be gene x term; got {adata.varm['I'].shape} "
            f"for {adata.n_vars} genes."
        )
    if args.condition_key not in adata.obs:
        raise SystemExit(f"condition_key {args.condition_key!r} not in adata.obs.")

    recon_loss = infer_recon_loss(adata, args.recon_loss)
    if recon_loss == "nb":
        if "counts" not in adata.layers:
            raise SystemExit("NB expiMap requires raw counts in adata.layers['counts'].")
        adata.X = adata.layers["counts"].copy()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    model_dir = output_dir / "model"

    scarches_file = str(getattr(sca, "__file__", ""))
    if "src/expiMap_scarches" in scarches_file:
        raise SystemExit(f"Refusing to use shadowed local scarches: {scarches_file}")
    if not hasattr(sca.models, "EXPIMAP"):
        raise SystemExit("Installed scarches does not expose sca.models.EXPIMAP.")

    conditions = sorted(adata.obs[args.condition_key].astype(str).unique().tolist())
    conditions_path = output_dir / "conditions.tsv"
    pd.DataFrame({args.condition_key: conditions}).to_csv(
        conditions_path,
        sep="\t",
        index=False,
    )
    model = sca.models.EXPIMAP(
        adata=adata,
        condition_key=args.condition_key,
        conditions=conditions,
        hidden_layer_sizes=parse_hidden_layers(args.hidden_layers),
        recon_loss=recon_loss,
        mask_key="I",
        use_ln=not args.no_layer_norm,
        use_bn=args.batch_norm,
        dr_rate=args.dropout,
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
        batch_size=args.batch_size,
        train_frac=args.train_frac,
        monitor=not args.no_monitor,
        monitor_only_val=args.monitor_only_val,
        early_stopping_kwargs={
            "early_stopping_metric": args.early_stopping_metric,
            "threshold": args.early_stopping_threshold,
            "patience": args.early_stopping_patience,
            "reduce_lr": not args.no_reduce_lr,
            "lr_patience": args.early_stopping_lr_patience,
            "lr_factor": args.early_stopping_lr_factor,
        },
    )

    latent = model.get_latent(mean=args.mean_latent)
    if isinstance(latent, tuple):
        latent = latent[0]
    active_indices = None
    if args.only_active:
        active_indices = list(map(int, model.nonzero_terms()))
        latent = latent[:, active_indices]

    terms = list(map(str, adata.uns["terms"]))
    columns = term_columns(terms, latent, active_indices)
    adata.obsm["X_expimap"] = latent
    score_frame = adata.obs.reset_index(names="obs_name").copy()
    score_values = pd.DataFrame(latent, columns=columns)
    scores = pd.concat([score_frame, score_values], axis=1)
    scores_path = output_dir / "pathway_scores.tsv"
    scores.to_csv(scores_path, sep="\t", index=False)

    term_table = pd.DataFrame(
        {
            "term": terms,
            "description": list(map(str, adata.uns.get("term_descriptions", terms))),
            "active": [
                bool(index in set(active_indices)) if active_indices is not None else True
                for index in range(len(terms))
            ],
        }
    )
    terms_path = output_dir / "terms.tsv"
    term_table.to_csv(terms_path, sep="\t", index=False)

    trained_h5ad = output_dir / "trained_input_with_scores.h5ad"
    adata.write_h5ad(trained_h5ad)
    model.save(str(model_dir), overwrite=True, save_anndata=False)

    history = {
        key: [float(value) for value in values]
        for key, values in model.trainer.logs.items()
    }
    history_path = output_dir / "training_history.json"
    history_path.write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
    final_metrics = {
        key: values[-1]
        for key, values in history.items()
        if values
    }
    summary = {
        "input": str(args.input),
        "output_dir": str(output_dir),
        "scarches_file": scarches_file,
        "has_EXPIMAP": bool(hasattr(sca.models, "EXPIMAP")),
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
        "condition_key": args.condition_key,
        "n_conditions": int(len(conditions)),
        "recon_loss": recon_loss,
        "epochs": args.epochs,
        "hidden_layers": parse_hidden_layers(args.hidden_layers),
        "n_samples": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_terms": int(len(terms)),
        "n_score_terms": int(latent.shape[1]),
        "posterior_mean_latent": bool(args.mean_latent),
        "training": {
            "epochs_completed": int(len(model.trainer.logs.get("epoch_loss", []))),
            "best_epoch": model.trainer.best_epoch,
            "training_seconds": float(model.trainer.training_time),
            "batch_size": args.batch_size,
            "train_fraction": args.train_frac,
            "early_stopping": bool(args.early_stopping),
            "early_stopping_patience": args.early_stopping_patience,
            "early_stopping_metric": args.early_stopping_metric,
            "early_stopping_threshold": args.early_stopping_threshold,
            "early_stopping_reduce_lr": not args.no_reduce_lr,
            "early_stopping_lr_patience": args.early_stopping_lr_patience,
            "early_stopping_lr_factor": args.early_stopping_lr_factor,
            "monitor_only_val": bool(args.monitor_only_val),
            "final_metrics": final_metrics,
        },
        "outputs": {
            "model": str(model_dir),
            "conditions": str(conditions_path),
            "pathway_scores": str(scores_path),
            "terms": str(terms_path),
            "trained_h5ad": str(trained_h5ad),
            "training_history": str(history_path),
        },
    }
    summary_path = output_dir / "training_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train direct OSDR expiMap model.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--condition-key", default="id.accession")
    parser.add_argument("--recon-loss", choices=["auto", "nb", "mse"], default="auto")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--hidden-layers", "--hidden-layer-sizes", dest="hidden_layers", default="128,128")
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--alpha", type=float, default=0.7)
    parser.add_argument("--alpha-kl", type=float, default=0.35)
    parser.add_argument("--alpha-epoch-anneal", type=int, default=25)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--batch-norm", action="store_true")
    parser.add_argument("--no-layer-norm", action="store_true")
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--early-stopping-metric", default="val_unweighted_loss")
    parser.add_argument("--early-stopping-threshold", type=float, default=0.0)
    parser.add_argument("--early-stopping-patience", type=int, default=30)
    parser.add_argument("--early-stopping-lr-patience", type=int, default=15)
    parser.add_argument("--early-stopping-lr-factor", type=float, default=0.1)
    parser.add_argument("--no-reduce-lr", action="store_true")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--train-frac", type=float, default=0.9)
    parser.add_argument("--no-monitor", action="store_true")
    parser.add_argument("--monitor-only-val", action="store_true")
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
    parser.add_argument("--only-active", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
