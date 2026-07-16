"""Command-line interface for the paper-parity mouse RNA diffusion project."""

from __future__ import annotations

import argparse

from .data import prepare
from .evaluate import evaluate
from .train import train


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare(args.config, force=args.force)
    elif args.command == "train":
        train(args.config, restart=args.restart)
    elif args.command == "evaluate":
        evaluate(args.config)


if __name__ == "__main__":
    main()
