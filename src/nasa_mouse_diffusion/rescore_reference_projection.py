"""Rescore OSDR query samples through a frozen ARCHS4-pretrained diffusion model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from nasa_mouse_glare.io import require_import

from . import analyze_features
from .data import (
    prepare_diffusion_data,
    reference_projection_categories,
    reference_projection_obs,
)
from .evaluate import generated_quality
from .model import ConditionalDiffusionMLP
from .reconstruction import train_linear_reconstructor
from .train_diffusion import encode_features, generate_full, write_scores


def checkpoint_args(checkpoint: dict) -> SimpleNamespace:
    values = dict(checkpoint.get("args", {}))
    if not values:
        raise SystemExit("Checkpoint does not contain training args.")
    return SimpleNamespace(**values)


def run(args) -> Path:
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    model_dir = Path(args.model_dir)
    checkpoint_path = model_dir / "reference_pretrained_model.pt"
    if not checkpoint_path.exists():
        raise SystemExit(f"Missing reference pretrained checkpoint: {checkpoint_path}")
    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    train_args = checkpoint_args(checkpoint)
    prepared = prepare_diffusion_data(
        query_h5ad=train_args.query_h5ad,
        query_source=train_args.query_source,
        reference_h5ad=train_args.reference_h5ad,
        reference_source=train_args.reference_source,
        max_genes=train_args.max_genes,
        n_landmarks=train_args.n_landmarks,
        landmark_strategy=train_args.landmark_strategy,
        l1000_map=train_args.l1000_map,
    )
    model = ConditionalDiffusionMLP(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    betas = checkpoint["betas"].to(device)
    projection_categories = reference_projection_categories(prepared)
    projection_obs = reference_projection_obs(prepared)
    projection_path = model_dir / "pretrained_query_projection_profiles.tsv"
    projection_obs[list(prepared.categorical_covariates)].drop_duplicates().to_csv(
        projection_path,
        sep="\t",
        index=False,
    )
    features = encode_features(
        model,
        prepared.query_landmark,
        projection_categories,
        batch_size=args.batch_size,
        device=device,
    )
    score_path = write_scores(model_dir / "pretrained_query_diffusion_feature_scores.tsv", prepared.query_obs, features)
    recon_train = prepared.reference_full if prepared.reference_full is not None else prepared.query_full
    reconstructor = train_linear_reconstructor(
        recon_train,
        prepared.landmark_indices,
        prepared.target_indices,
        alpha=train_args.reconstruction_alpha,
    )
    fake = generate_full(
        model,
        projection_categories,
        reconstructor,
        betas=betas,
        sample_steps=args.sample_steps or train_args.sample_steps,
        eta=train_args.eta,
        batch_size=args.batch_size,
        device=device,
        seed=args.seed,
    )
    quality = generated_quality(prepared.query_full, fake, max_pr_samples=train_args.max_metric_samples)
    pd.DataFrame([quality]).to_csv(model_dir / "pretrained_query_generated_quality.tsv", sep="\t", index=False)
    (model_dir / "pretrained_query_generated_quality.json").write_text(json.dumps(quality, indent=2) + "\n", encoding="utf-8")
    summary_path = model_dir / "training_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary["pretrained_query_quality"] = quality
        summary.setdefault("outputs", {})["pretrained_query_scores"] = str(score_path)
        summary.setdefault("outputs", {})["pretrained_query_projection_profiles"] = str(projection_path)
        summary["pretrained_query_projection_note"] = (
            "Frozen ARCHS4 projection uses query expression and tissue, but replaces condition/accession/source/"
            "platform/material/sex covariates with trained ARCHS4 reference defaults by tissue."
        )
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.analyze:
        analyze_features.run(
            SimpleNamespace(
                scores=str(score_path),
                output_dir=str(model_dir / "pretrained_query_analysis"),
                top_features=args.top_features,
            )
        )
    result = {
        "model_dir": str(model_dir),
        "checkpoint": str(checkpoint_path),
        "device": str(device),
        "scores": str(score_path),
        "projection_profiles": str(projection_path),
        "quality": quality,
    }
    output_path = model_dir / "pretrained_query_rescore_summary.json"
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", default="outputs/diffusion_conditional_generation/archs4_pretrain_osdr_finetune")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--sample-steps", type=int, default=0)
    parser.add_argument("--top-features", type=int, default=30)
    parser.add_argument("--seed", type=int, default=2037)
    parser.add_argument("--skip-analysis", action="store_true")
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()
    args.analyze = not args.skip_analysis
    return args


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
