"""Command-line interface for the paper-parity mouse RNA diffusion project."""

from __future__ import annotations

import argparse

from .data import prepare
from .evaluate import evaluate
from .train import train
from .conditional_data import prepare_conditional
from .conditional_evaluate import evaluate_conditional
from .conditional_train import train_conditional
from .real_effect_ceiling import run_ceiling
from .factorized_train import train_factorized
from .factorized_evaluate import evaluate_factorized


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare", help="Prepare paper-matched ARCHS4 data")
    prepare_parser.add_argument(
        "--config", default="configs/rna_diffusion/archs4_mouse_paper_parity.yaml"
    )
    prepare_parser.add_argument("--force", action="store_true")
    train_parser = subparsers.add_parser("train", help="Train the exact upstream DDIM")
    train_parser.add_argument(
        "--config", default="configs/rna_diffusion/archs4_mouse_paper_parity.yaml"
    )
    train_parser.add_argument("--restart", action="store_true")
    evaluate_parser = subparsers.add_parser(
        "evaluate", help="Generate the paper-style trajectory and held-out metrics"
    )
    evaluate_parser.add_argument(
        "--config", default="configs/rna_diffusion/archs4_mouse_paper_parity.yaml"
    )
    conditional_prepare = subparsers.add_parser(
        "prepare-osdr", help="Prepare API-derived OSDR data for upstream ModelDDIM"
    )
    conditional_prepare.add_argument(
        "--config",
        default="configs/rna_diffusion/osdr_pooled_flt_gc_paper_architecture.yaml",
    )
    conditional_prepare.add_argument("--force", action="store_true")
    conditional_train = subparsers.add_parser(
        "train-osdr", help="Train upstream ModelDDIM on tissue and FLT/GC classes"
    )
    conditional_train.add_argument(
        "--config",
        default="configs/rna_diffusion/osdr_pooled_flt_gc_paper_architecture.yaml",
    )
    conditional_train.add_argument("--restart", action="store_true")
    conditional_evaluate = subparsers.add_parser(
        "evaluate-osdr", help="Evaluate conditional ModelDDIM on held-out accessions"
    )
    conditional_evaluate.add_argument(
        "--config",
        default="configs/rna_diffusion/osdr_pooled_flt_gc_paper_architecture.yaml",
    )
    conditional_evaluate.add_argument("--unlock-test", action="store_true")
    conditional_evaluate.add_argument(
        "--eta", type=float, default=None, help="Override DDIM sampling stochasticity"
    )
    conditional_evaluate.add_argument(
        "--evaluation-variant",
        default="",
        help="Write evaluation to a named validation/test variant directory",
    )
    ceiling = subparsers.add_parser(
        "real-ceiling", help="Measure real cross-accession FLT/GC reproducibility"
    )
    ceiling.add_argument("--prepared-h5", required=True)
    ceiling.add_argument("--samples-tsv", required=True)
    ceiling.add_argument("--output-dir", required=True)
    ceiling.add_argument("--roles", nargs="+", default=["train"])
    ceiling.add_argument("--transform", choices=("log1p", "none"), default="log1p")
    ceiling.add_argument("--permutation-repeats", type=int, default=100)
    ceiling.add_argument("--seed", type=int, default=1234)
    adapter_train = subparsers.add_parser(
        "train-adapter", help="Train staged factorized residual DDIM adapters"
    )
    adapter_train.add_argument("--config", required=True)
    adapter_train.add_argument("--restart", action="store_true")
    adapter_evaluate = subparsers.add_parser(
        "evaluate-adapter", help="Evaluate factorized adapters on validation only"
    )
    adapter_evaluate.add_argument("--config", required=True)
    adapter_evaluate.add_argument("--guidance-scales", nargs="+", type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare(args.config, force=args.force)
    elif args.command == "train":
        train(args.config, restart=args.restart)
    elif args.command == "evaluate":
        evaluate(args.config)
    elif args.command == "prepare-osdr":
        prepare_conditional(args.config, force=args.force)
    elif args.command == "train-osdr":
        train_conditional(args.config, restart=args.restart)
    elif args.command == "evaluate-osdr":
        evaluate_conditional(
            args.config,
            unlock_test=args.unlock_test,
            eta_override=args.eta,
            evaluation_variant=args.evaluation_variant,
        )
    elif args.command == "real-ceiling":
        run_ceiling(
            args.prepared_h5,
            args.samples_tsv,
            args.output_dir,
            roles=args.roles,
            transform=args.transform,
            permutation_repeats=args.permutation_repeats,
            seed=args.seed,
        )
    elif args.command == "train-adapter":
        train_factorized(args.config, restart=args.restart)
    elif args.command == "evaluate-adapter":
        evaluate_factorized(
            args.config, guidance_scales=args.guidance_scales
        )


if __name__ == "__main__":
    main()
