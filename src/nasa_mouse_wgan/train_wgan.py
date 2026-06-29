"""Train conditional WGAN-GP on API-derived OSDR/ARCHS4 bulk RNA-seq inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .data import CONDITIONAL_GENERATION_COVARIATES, load_prepared_data
from .model import ConditionalWGANGP
from .training import (
    TrainConfig,
    critic_features,
    generate_samples,
    generation_quality,
    train_model,
)


def parse_hidden_dims(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item.strip())


def parse_covariates(value: str) -> tuple[str, ...]:
    if value == "conditional_generation":
        return CONDITIONAL_GENERATION_COVARIATES
    return tuple(item.strip() for item in value.split(",") if item.strip())


def cardinalities(vocabularies: dict[str, list[str]], covariates: tuple[str, ...]) -> list[int]:
    return [len(vocabularies[covariate]) for covariate in covariates]


def device_summary(torch, device) -> dict:
    return {
        "version": str(torch.__version__),
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
        "cuda_device_name": (
            torch.cuda.get_device_name(0) if torch.cuda.is_available() else ""
        ),
        "model_device": str(device),
    }


def write_feature_scores(path: Path, obs, critic_score, features) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    metadata_columns = [
        column
        for column in [
            "profile_id",
            "id.accession",
            "condition_inferred",
            "tissue_final",
            "muscle_group",
            "study.characteristics.material type",
            "wgan_condition",
            "wgan_tissue",
            "wgan_material_type",
            "wgan_muscle_group",
            "wgan_accession",
            "wgan_source",
            "wgan_sex",
            "wgan_assay",
            "wgan_platform",
            "wgan_data_source",
        ]
        if column in obs
    ]
    frame = obs[metadata_columns].copy()
    score_frame = pd.DataFrame({"CRITIC_SCORE": critic_score}, index=frame.index)
    feature_frame = pd.DataFrame(
        features,
        columns=[f"WGAN_FEATURE_{idx:03d}" for idx in range(features.shape[1])],
        index=frame.index,
    )
    frame = pd.concat([frame, score_frame, feature_frame], axis=1)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
    return path


def write_quality(path: Path, quality: dict) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([quality]).to_csv(path.with_suffix(".tsv"), sep="\t", index=False)
    path.write_text(json.dumps(quality, indent=2) + "\n", encoding="utf-8")


def write_readme(output_dir: Path, summary: dict) -> None:
    lines = [
        "# WGAN-GP Run",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Query samples: {summary['counts']['query_samples']}",
        f"- Reference samples: {summary['counts']['reference_samples']}",
        f"- Genes: {summary['counts']['genes']}",
        f"- Critic features: {summary['counts']['critic_features']}",
        f"- Device: `{summary['torch']['model_device']}` {summary['torch']['cuda_device_name']}",
        "",
        "Primary outputs:",
        "",
        "- `critic_feature_scores.tsv`",
        "- `pretrained_query_critic_feature_scores.tsv` when applicable",
        "- `generated_quality.json` / `generated_quality.tsv`",
        "- `training_summary.json`",
        "- `model.pt`",
    ]
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_observed_profiles(output_dir: Path, prepared) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frames = [prepared.query_obs]
    if prepared.reference_obs is not None:
        frames.insert(0, prepared.reference_obs)
    columns = list(prepared.categorical_covariates)
    rows = []
    for label, frame in zip(
        ["reference", "query"] if prepared.reference_obs is not None else ["query"],
        frames,
    ):
        item = frame[columns].drop_duplicates().copy()
        item.insert(0, "training_source", label)
        rows.append(item)
    profiles = pd.concat(rows, axis=0, ignore_index=True).drop_duplicates()
    path = output_dir / "observed_conditioning_profiles.tsv"
    profiles.to_csv(path, sep="\t", index=False)
    return path


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    categorical_covariates = parse_covariates(
        getattr(args, "categorical_covariates", "wgan_condition,wgan_accession")
    )

    prepared = load_prepared_data(
        query_h5ad=args.query_h5ad,
        reference_h5ad=args.reference_h5ad,
        query_source=getattr(args, "query_source", "osdr"),
        reference_source=getattr(args, "reference_source", "archs4"),
        categorical_covariates=categorical_covariates,
        clip=args.clip,
        max_genes=args.max_genes,
    )
    card = cardinalities(prepared.vocabularies, prepared.categorical_covariates)
    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    model = ConditionalWGANGP(
        expression_dim=len(prepared.genes),
        categorical_cardinalities=card,
        noise_dim=args.noise_dim,
        hidden_dims=parse_hidden_dims(args.hidden_dims),
    ).to(device)
    model_config = {
        "expression_dim": int(len(prepared.genes)),
        "categorical_cardinalities": [int(value) for value in card],
        "noise_dim": int(args.noise_dim),
        "hidden_dims": list(parse_hidden_dims(args.hidden_dims)),
    }

    histories = {}
    pretrained_scores_path = None
    pretrained_quality = None
    if prepared.reference_x is not None:
        histories["reference"] = train_model(
            model,
            prepared.reference_x,
            prepared.reference_categories,
            config=TrainConfig(
                epochs=args.reference_epochs,
                batch_size=args.batch_size,
                learning_rate=args.learning_rate,
                critic_steps=args.critic_steps,
                gradient_penalty=args.gradient_penalty,
                seed=args.seed,
            ),
            device=device,
        )
        torch.save(
            {
                "model_state_dict": model.state_dict(),
                "genes": prepared.genes,
                "vocabularies": prepared.vocabularies,
                "categorical_covariates": prepared.categorical_covariates,
                "model_config": model_config,
            },
            output_dir / "reference_pretrained_model.pt",
        )
        critic_score, features = critic_features(
            model,
            prepared.query_x,
            prepared.query_categories,
            batch_size=args.batch_size,
            device=device,
        )
        pretrained_scores_path = output_dir / "pretrained_query_critic_feature_scores.tsv"
        write_feature_scores(pretrained_scores_path, prepared.query_obs, critic_score, features)
        fake = generate_samples(
            model, prepared.query_categories, batch_size=args.batch_size, device=device
        )
        pretrained_quality = generation_quality(prepared.query_x, fake)
        write_quality(output_dir / "pretrained_query_generated_quality.json", pretrained_quality)

        histories["query_finetune"] = train_model(
            model,
            prepared.query_x,
            prepared.query_categories,
            config=TrainConfig(
                epochs=args.query_epochs,
                batch_size=args.batch_size,
                learning_rate=args.finetune_learning_rate or args.learning_rate,
                critic_steps=args.critic_steps,
                gradient_penalty=args.gradient_penalty,
                seed=args.seed + 17,
            ),
            device=device,
        )
    else:
        histories["direct"] = train_model(
            model,
            prepared.query_x,
            prepared.query_categories,
            config=TrainConfig(
                epochs=args.query_epochs,
                batch_size=args.batch_size,
                learning_rate=args.learning_rate,
                critic_steps=args.critic_steps,
                gradient_penalty=args.gradient_penalty,
                seed=args.seed,
            ),
            device=device,
        )

    critic_score, features = critic_features(
        model,
        prepared.query_x,
        prepared.query_categories,
        batch_size=args.batch_size,
        device=device,
    )
    scores_path = output_dir / "critic_feature_scores.tsv"
    write_feature_scores(scores_path, prepared.query_obs, critic_score, features)
    fake = generate_samples(model, prepared.query_categories, batch_size=args.batch_size, device=device)
    quality = generation_quality(prepared.query_x, fake)
    write_quality(output_dir / "generated_quality.json", quality)

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "genes": prepared.genes,
            "vocabularies": prepared.vocabularies,
            "categorical_covariates": prepared.categorical_covariates,
            "model_config": model_config,
        },
        output_dir / "model.pt",
    )
    np.savez_compressed(output_dir / "normalization_stats.npz", mean=prepared.mean, std=prepared.std)
    (output_dir / "genes.tsv").write_text(
        "gene\n" + "\n".join(prepared.genes) + "\n", encoding="utf-8"
    )
    vocab_path = output_dir / "categorical_vocabularies.json"
    vocab_path.write_text(json.dumps(prepared.vocabularies, indent=2) + "\n", encoding="utf-8")
    profiles_path = write_observed_profiles(output_dir, prepared)

    mode = getattr(args, "run_mode", "") or (
        "archs4_pretrain_osdr_finetune" if prepared.reference_x is not None else "direct_osdr"
    )
    summary = {
        "method": "conditional WGAN-GP",
        "paper": {
            "title": "Adversarial generation of gene expression data",
            "pmc": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8756177/",
            "doi": "https://doi.org/10.1093/bioinformatics/btab035",
        },
        "mode": mode,
        "query_h5ad": str(args.query_h5ad),
        "reference_h5ad": str(args.reference_h5ad or ""),
        "query_source": str(getattr(args, "query_source", "osdr")),
        "reference_source": str(getattr(args, "reference_source", "archs4")),
        "output_dir": str(output_dir),
        "normalization": "log1p(CPM) followed by reference-gene z-score for pretrained runs or query-gene z-score for direct runs",
        "categorical_covariates": list(prepared.categorical_covariates),
        "training_design": {
            "reference_epochs": int(args.reference_epochs if prepared.reference_x is not None else 0),
            "query_epochs": int(args.query_epochs),
            "noise_dim": int(args.noise_dim),
            "hidden_dims": list(parse_hidden_dims(args.hidden_dims)),
            "critic_steps": int(args.critic_steps),
            "gradient_penalty": float(args.gradient_penalty),
            "learning_rate": float(args.learning_rate),
            "finetune_learning_rate": float(args.finetune_learning_rate or args.learning_rate),
        },
        "torch": device_summary(torch, device),
        "counts": {
            "query_samples": int(prepared.query_x.shape[0]),
            "reference_samples": int(prepared.reference_x.shape[0]) if prepared.reference_x is not None else 0,
            "genes": int(prepared.query_x.shape[1]),
            "critic_features": int(features.shape[1]),
        },
        "training": histories,
        "quality": quality,
        "pretrained_query_quality": pretrained_quality,
        "outputs": {
            "scores": str(scores_path),
            "pretrained_query_scores": str(pretrained_scores_path or ""),
            "generated_quality": str(output_dir / "generated_quality.json"),
            "model": str(output_dir / "model.pt"),
            "reference_pretrained_model": str(output_dir / "reference_pretrained_model.pt")
            if prepared.reference_x is not None
            else "",
            "normalization_stats": str(output_dir / "normalization_stats.npz"),
            "genes": str(output_dir / "genes.tsv"),
            "vocabularies": str(vocab_path),
            "observed_conditioning_profiles": str(profiles_path),
            "summary": str(output_dir / "training_summary.json"),
        },
    }
    summary_path = output_dir / "training_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_readme(output_dir, summary)
    printed_summary = {
        "method": summary["method"],
        "mode": summary["mode"],
        "output_dir": summary["output_dir"],
        "torch": summary["torch"],
        "counts": summary["counts"],
        "training_design": summary["training_design"],
        "quality": summary["quality"],
        "pretrained_query_quality": summary["pretrained_query_quality"],
        "outputs": summary["outputs"],
    }
    print(json.dumps(printed_summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--reference-h5ad", default="")
    parser.add_argument("--query-source", default="osdr", choices=["osdr", "archs4"])
    parser.add_argument("--reference-source", default="archs4", choices=["osdr", "archs4"])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-mode", default="")
    parser.add_argument(
        "--categorical-covariates",
        default="wgan_condition,wgan_accession",
        help="Comma-separated covariates, or 'conditional_generation'.",
    )
    parser.add_argument("--reference-epochs", type=int, default=100)
    parser.add_argument("--query-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=0.0)
    parser.add_argument("--critic-steps", type=int, default=5)
    parser.add_argument("--gradient-penalty", type=float, default=10.0)
    parser.add_argument("--noise-dim", type=int, default=128)
    parser.add_argument("--hidden-dims", default="256,256")
    parser.add_argument("--clip", type=float, default=10.0)
    parser.add_argument("--max-genes", type=int, default=0)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
