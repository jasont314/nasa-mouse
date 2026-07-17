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
from .factorized_calibrate import calibrate_factorized
from .factorized_mean_calibrate import calibrate_factorized_means
from .factorized_distribution_calibrate import calibrate_factorized_distribution
from .factorized_subset import subset_factorized_data
from .factorized_final_evaluate import evaluate_factorized_finalist_test


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
    adapter_evaluate.add_argument("--model-artifact", default="model.pt")
    adapter_evaluate.add_argument("--validation-sampling-seed", type=int)
    adapter_evaluate.add_argument("--train-sampling-seed", type=int)
    adapter_evaluate.add_argument("--evaluation-variant", default="")
    adapter_evaluate.add_argument(
        "--sampling-noise",
        choices=["pseudo_random", "stratified_antithetic"],
    )
    adapter_calibrate = subparsers.add_parser(
        "calibrate-adapter", help="Fit train-only covariance calibration"
    )
    adapter_calibrate.add_argument("--config", required=True)
    adapter_calibrate.add_argument("--guidance-scale", type=float, default=1.0)
    adapter_calibrate.add_argument(
        "--ridge-fractions", nargs="+", type=float, default=[0.001, 0.01, 0.1]
    )
    mean_calibrate = subparsers.add_parser(
        "calibrate-mean-adapter",
        help="Fit train-only condition-blind hierarchical mean calibration",
    )
    mean_calibrate.add_argument("--config", required=True)
    mean_calibrate.add_argument("--guidance-scale", type=float, default=1.0)
    mean_calibrate.add_argument(
        "--group-columns", nargs="+", default=["accession", "tissue"]
    )
    mean_calibrate.add_argument(
        "--prior-strengths", nargs="+", type=float, default=[2.0, 5.0, 10.0]
    )
    mean_calibrate.add_argument("--evaluation-variant", default="")
    distribution_calibrate = subparsers.add_parser(
        "calibrate-distribution-adapter",
        help="Fit and evaluate repeated-seed train-only distribution calibration",
    )
    distribution_calibrate.add_argument("--config", required=True)
    distribution_calibrate.add_argument("--guidance-scale", type=float, default=1.0)
    distribution_calibrate.add_argument(
        "--fit-variants", nargs="+", default=["base", "seed3022", "seed3023"]
    )
    distribution_calibrate.add_argument(
        "--evaluation-variants",
        nargs="+",
        default=["base", "seed3021", "seed3022", "seed3023"],
    )
    distribution_calibrate.add_argument(
        "--group-columns", nargs="+", default=["accession", "tissue"]
    )
    distribution_calibrate.add_argument("--prior-strength", type=float, default=5.0)
    distribution_calibrate.add_argument("--residual-scale", type=float, default=0.5)
    distribution_calibrate.add_argument("--residual-seed", type=int, default=9100)
    distribution_calibrate.add_argument(
        "--noise-group-columns",
        nargs="+",
        default=["accession", "tissue", "condition"],
    )
    distribution_calibrate.add_argument(
        "--minimum-repeat-pass-fraction", type=float, default=0.75
    )
    subset = subparsers.add_parser(
        "subset-factorized-data",
        help="Create a train/validation tissue subset without opening locked test",
    )
    subset.add_argument("--source-h5", required=True)
    subset.add_argument("--samples-tsv", required=True)
    subset.add_argument("--output-h5", required=True)
    subset.add_argument("--tissues", nargs="+", required=True)
    final_test = subparsers.add_parser(
        "evaluate-finalist-test",
        help="Run the one-time repeated locked-test evaluation",
    )
    final_test.add_argument("--config", required=True)
    final_test.add_argument("--calibrator-dir", required=True)
    final_test.add_argument("--unlock-test", action="store_true")
    final_test.add_argument(
        "--sampling-seeds", nargs="+", type=int, default=[5020, 5021, 5022, 5023]
    )
    final_test.add_argument("--residual-seed", type=int, default=15020)
    final_test.add_argument(
        "--minimum-repeat-pass-fraction", type=float, default=0.75
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
            args.config,
            guidance_scales=args.guidance_scales,
            model_artifact=args.model_artifact,
            validation_sampling_seed=args.validation_sampling_seed,
            train_sampling_seed=args.train_sampling_seed,
            evaluation_variant=args.evaluation_variant,
            sampling_noise=args.sampling_noise,
        )
    elif args.command == "calibrate-adapter":
        calibrate_factorized(
            args.config,
            guidance_scale=args.guidance_scale,
            ridge_fractions=args.ridge_fractions,
        )
    elif args.command == "calibrate-mean-adapter":
        calibrate_factorized_means(
            args.config,
            guidance_scale=args.guidance_scale,
            group_columns=args.group_columns,
            prior_strengths=args.prior_strengths,
            evaluation_variant=args.evaluation_variant,
        )
    elif args.command == "calibrate-distribution-adapter":
        calibrate_factorized_distribution(
            args.config,
            guidance_scale=args.guidance_scale,
            fit_variants=args.fit_variants,
            evaluation_variants=args.evaluation_variants,
            group_columns=args.group_columns,
            prior_strength=args.prior_strength,
            residual_scale=args.residual_scale,
            residual_seed=args.residual_seed,
            noise_group_columns=args.noise_group_columns,
            minimum_repeat_pass_fraction=args.minimum_repeat_pass_fraction,
        )
    elif args.command == "subset-factorized-data":
        subset_factorized_data(
            args.source_h5,
            args.samples_tsv,
            args.output_h5,
            tissues=args.tissues,
        )
    elif args.command == "evaluate-finalist-test":
        evaluate_factorized_finalist_test(
            args.config,
            args.calibrator_dir,
            unlock_test=args.unlock_test,
            sampling_seeds=args.sampling_seeds,
            residual_seed=args.residual_seed,
            minimum_repeat_pass_fraction=args.minimum_repeat_pass_fraction,
        )


if __name__ == "__main__":
    main()
