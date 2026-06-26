"""Train expiMap reference models on prepared ARCHS4 tissue AnnData inputs."""

from __future__ import annotations

import argparse

from .train_expimap_direct import run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train an ARCHS4 tissue reference expiMap model."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--condition-key", default="archs4_condition")
    parser.add_argument("--recon-loss", choices=["auto", "nb", "mse"], default="auto")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument(
        "--hidden-layers",
        "--hidden-layer-sizes",
        dest="hidden_layers",
        default="128,128",
    )
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
    parser.add_argument("--mean-latent", action="store_true")
    parser.add_argument("--only-active", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
