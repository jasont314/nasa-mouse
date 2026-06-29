"""Train conditional DDPM/DDIM on API-derived OSDR/ARCHS4 bulk RNA-seq inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .data import DIFFUSION_COVARIATES, prepare_diffusion_data, write_observed_profiles
from .diffusion import beta_schedule, noise_estimation_loss, sample
from .evaluate import generated_quality, reverse_validation
from .model import ConditionalDiffusionMLP
from .reconstruction import reconstruction_metrics, train_linear_reconstructor


def torch_info(device, torch):
    return {
        "version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device_count": int(torch.cuda.device_count()),
        "cuda_device_name": torch.cuda.get_device_name(device) if device.type == "cuda" else "",
        "model_device": str(device),
    }


def parse_hidden(value: str) -> tuple[int, int]:
    parts = [int(part) for part in str(value).split(",") if part]
    if not parts:
        return 512, 2
    if len(parts) == 1:
        return parts[0], 2
    return parts[0], len(parts)


def train_epochs(model, x, categories, *, betas, epochs: int, batch_size: int, lr: float, device, amp: bool, seed: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    data_utils = require_import("torch.utils.data", "pip install -r requirements-nasa-mouse-glare.txt")
    dataset = data_utils.TensorDataset(
        torch.as_tensor(x, dtype=torch.float32),
        torch.as_tensor(categories, dtype=torch.long),
    )
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    loader = data_utils.DataLoader(dataset, batch_size=int(batch_size), shuffle=True, generator=generator)
    optimizer = torch.optim.Adam(model.parameters(), lr=float(lr), eps=1e-8)
    scaler = torch.cuda.amp.GradScaler(enabled=amp and device.type == "cuda")
    history = []
    for epoch in range(int(epochs)):
        losses = []
        avgs = []
        for xb, cb in loader:
            xb = xb.to(device)
            cb = cb.to(device)
            optimizer.zero_grad(set_to_none=True)
            with torch.cuda.amp.autocast(enabled=amp and device.type == "cuda"):
                loss, avg = noise_estimation_loss(model, xb, cb, betas)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            losses.append(float(loss.detach().cpu()))
            avgs.append(float(avg.detach().cpu()))
        history.append({"epoch": epoch, "loss": float(np.mean(losses)), "noise_abs_error": float(np.mean(avgs))})
    return history


def encode_features(model, matrix, categories, *, batch_size: int, device):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    outputs = []
    model.eval()
    with torch.no_grad():
        for start in range(0, matrix.shape[0], batch_size):
            end = min(start + batch_size, matrix.shape[0])
            xb = torch.as_tensor(matrix[start:end], dtype=torch.float32, device=device)
            cb = torch.as_tensor(categories[start:end], dtype=torch.long, device=device)
            tb = torch.zeros(end - start, dtype=torch.long, device=device)
            outputs.append(model.features(xb, tb, cb).detach().cpu().numpy())
    return np.concatenate(outputs, axis=0).astype("float32")


def write_scores(path: Path, obs, features):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = obs.copy().reset_index(drop=True)
    for idx in range(features.shape[1]):
        frame[f"DIFFUSION_FEATURE_{idx:03d}"] = features[:, idx]
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
    return path


def generate_full(model, categories, reconstructor, *, betas, sample_steps: int, eta: float, batch_size: int, device, seed: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    generated = []
    gen = torch.Generator(device=device)
    gen.manual_seed(int(seed))
    for start in range(0, categories.shape[0], batch_size):
        end = min(start + batch_size, categories.shape[0])
        noise = torch.randn((end - start, model.expression_dim), generator=gen, device=device)
        cats = categories[start:end]
        landmarks = sample(model, cats, betas=betas, sample_steps=sample_steps, eta=eta, noise=noise, device=device)
        generated.append(reconstructor.reconstruct_full(landmarks.detach().cpu().numpy()))
    return np.concatenate(generated, axis=0).astype("float32")


def save_model(path: Path, model, prepared, args, *, model_config: dict, betas, reconstructor) -> None:
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    payload = {
        "model_state_dict": model.state_dict(),
        "model_config": model_config,
        "categorical_covariates": prepared.categorical_covariates,
        "vocabularies": prepared.vocabularies,
        "genes": prepared.genes,
        "landmark_genes": prepared.landmark_genes,
        "target_genes": prepared.target_genes,
        "landmark_indices": prepared.landmark_indices,
        "target_indices": prepared.target_indices,
        "betas": betas.detach().cpu(),
        "args": vars(args),
        "reconstruction": {
            "kind": "ridge_lr",
            "coef": getattr(reconstructor.model, "coef_", None),
            "intercept": getattr(reconstructor.model, "intercept_", None),
            "alpha": getattr(reconstructor.model, "alpha", None),
        },
    }
    torch.save(payload, path)


def write_matrix(path: Path, values, genes: list[str]) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = pd.DataFrame(values, columns=genes)
    frame.insert(0, "synthetic_sample_id", [f"synthetic_{idx:05d}" for idx in range(values.shape[0])])
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def write_training_readme(output_dir: Path, summary: dict) -> None:
    lines = [
        "# Diffusion Run",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Query samples: {summary['counts']['query_samples']}",
        f"- Reference samples: {summary['counts']['reference_samples']}",
        f"- Full genes: {summary['counts']['genes']}",
        f"- Landmark genes: {summary['counts']['landmark_genes']}",
        f"- Landmark source: `{summary['landmarks']['source']}`",
        f"- Device: `{summary['torch']['model_device']}` {summary['torch']['cuda_device_name']}",
        "",
        "Outputs include `model.pt`, feature scores, generated-quality metrics, reconstruction metrics, and observed conditioning profiles.",
    ]
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(int(args.seed))
    np.random.seed(int(args.seed))
    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    prepared = prepare_diffusion_data(
        query_h5ad=args.query_h5ad,
        query_source=args.query_source,
        reference_h5ad=args.reference_h5ad,
        reference_source=args.reference_source,
        max_genes=args.max_genes,
        n_landmarks=args.n_landmarks,
        landmark_strategy=args.landmark_strategy,
        l1000_map=args.l1000_map,
    )
    hidden_dim, n_blocks = parse_hidden(args.hidden_dims)
    model_config = {
        "expression_dim": int(len(prepared.landmark_genes)),
        "categorical_cardinalities": [len(prepared.vocabularies[cov]) for cov in prepared.categorical_covariates],
        "hidden_dim": int(hidden_dim),
        "n_blocks": int(n_blocks),
        "dropout": float(args.dropout),
        "time_embedding_dim": int(args.time_embedding_dim),
        "categorical_embedding_dim": int(args.categorical_embedding_dim),
        "sinusoidal_time": bool(args.sinusoidal_time),
        "num_timesteps": int(args.diffusion_timesteps),
    }
    model = ConditionalDiffusionMLP(**model_config).to(device)
    betas = torch.as_tensor(
        beta_schedule(
            args.beta_schedule,
            beta_start=args.beta_start,
            beta_end=args.beta_end,
            timesteps=args.diffusion_timesteps,
        ),
        dtype=torch.float32,
        device=device,
    )
    recon_train = prepared.reference_full if prepared.reference_full is not None else prepared.query_full
    reconstructor = train_linear_reconstructor(
        recon_train,
        prepared.landmark_indices,
        prepared.target_indices,
        alpha=args.reconstruction_alpha,
    )
    reference_history = []
    pretrained_query_quality = None
    if prepared.reference_landmark is not None and int(args.reference_epochs) > 0:
        reference_history = train_epochs(
            model,
            prepared.reference_landmark,
            prepared.reference_categories,
            betas=betas,
            epochs=args.reference_epochs,
            batch_size=args.batch_size,
            lr=args.learning_rate,
            device=device,
            amp=not args.no_amp,
            seed=args.seed,
        )
        pretrained_features = encode_features(
            model,
            prepared.query_landmark,
            prepared.query_categories,
            batch_size=args.batch_size,
            device=device,
        )
        write_scores(output_dir / "pretrained_query_diffusion_feature_scores.tsv", prepared.query_obs, pretrained_features)
        pretrained_fake = generate_full(
            model,
            prepared.query_categories,
            reconstructor,
            betas=betas,
            sample_steps=args.sample_steps,
            eta=args.eta,
            batch_size=args.batch_size,
            device=device,
            seed=args.seed + 17,
        )
        pretrained_query_quality = generated_quality(prepared.query_full, pretrained_fake, max_pr_samples=args.max_metric_samples)
        save_model(output_dir / "reference_pretrained_model.pt", model, prepared, args, model_config=model_config, betas=betas, reconstructor=reconstructor)

    query_lr = args.finetune_learning_rate if args.finetune_learning_rate > 0 else args.learning_rate
    query_history = train_epochs(
        model,
        prepared.query_landmark,
        prepared.query_categories,
        betas=betas,
        epochs=args.query_epochs,
        batch_size=args.batch_size,
        lr=query_lr,
        device=device,
        amp=not args.no_amp,
        seed=args.seed + 1,
    )
    features = encode_features(model, prepared.query_landmark, prepared.query_categories, batch_size=args.batch_size, device=device)
    write_scores(output_dir / "diffusion_feature_scores.tsv", prepared.query_obs, features)
    fake_full = generate_full(
        model,
        prepared.query_categories,
        reconstructor,
        betas=betas,
        sample_steps=args.sample_steps,
        eta=args.eta,
        batch_size=args.batch_size,
        device=device,
        seed=args.seed + 29,
    )
    quality = generated_quality(prepared.query_full, fake_full, max_pr_samples=args.max_metric_samples)
    labels = prepared.query_obs["condition_inferred"].astype(str).to_numpy()
    rv = reverse_validation(prepared.query_full, fake_full, labels)
    recon_metrics = reconstruction_metrics(reconstructor, recon_train)
    np.savez(
        output_dir / "normalization_stats.npz",
        center=prepared.full_center,
        scale=prepared.full_scale,
        landmark_indices=prepared.landmark_indices,
        target_indices=prepared.target_indices,
    )
    pd.DataFrame({"gene_id": prepared.genes}).to_csv(output_dir / "genes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": prepared.landmark_genes}).to_csv(output_dir / "landmark_genes.tsv", sep="\t", index=False)
    (output_dir / "categorical_vocabularies.json").write_text(json.dumps(prepared.vocabularies, indent=2) + "\n", encoding="utf-8")
    write_observed_profiles(output_dir / "observed_conditioning_profiles.tsv", prepared)
    save_model(output_dir / "model.pt", model, prepared, args, model_config=model_config, betas=betas, reconstructor=reconstructor)
    pd.DataFrame(reference_history).to_csv(output_dir / "reference_training_history.tsv", sep="\t", index=False)
    pd.DataFrame(query_history).to_csv(output_dir / "query_training_history.tsv", sep="\t", index=False)
    pd.DataFrame([quality]).to_csv(output_dir / "generated_quality.tsv", sep="\t", index=False)
    (output_dir / "generated_quality.json").write_text(json.dumps(quality, indent=2) + "\n", encoding="utf-8")
    if pretrained_query_quality is not None:
        pd.DataFrame([pretrained_query_quality]).to_csv(output_dir / "pretrained_query_generated_quality.tsv", sep="\t", index=False)
        (output_dir / "pretrained_query_generated_quality.json").write_text(
            json.dumps(pretrained_query_quality, indent=2) + "\n", encoding="utf-8"
        )
    summary = {
        "method": "conditional DDPM/DDIM diffusion",
        "mode": args.run_mode,
        "output_dir": str(output_dir),
        "torch": torch_info(device, torch),
        "counts": {
            "query_samples": int(prepared.query_full.shape[0]),
            "reference_samples": int(prepared.reference_full.shape[0]) if prepared.reference_full is not None else 0,
            "genes": int(len(prepared.genes)),
            "landmark_genes": int(len(prepared.landmark_genes)),
            "target_genes": int(len(prepared.target_genes)),
            "diffusion_features": int(features.shape[1]),
        },
        "landmarks": {
            "strategy": prepared.landmark_strategy,
            "source": prepared.landmark_source,
            "l1000_map": str(args.l1000_map),
        },
        "training_design": {
            "reference_epochs": int(args.reference_epochs),
            "query_epochs": int(args.query_epochs),
            "diffusion_timesteps": int(args.diffusion_timesteps),
            "sample_steps": int(args.sample_steps),
            "eta": float(args.eta),
            "beta_schedule": args.beta_schedule,
            "learning_rate": float(args.learning_rate),
            "finetune_learning_rate": float(query_lr),
            "batch_size": int(args.batch_size),
            "hidden_dim": int(hidden_dim),
            "n_blocks": int(n_blocks),
        },
        "quality": quality,
        "pretrained_query_quality": pretrained_query_quality,
        "reconstruction": recon_metrics,
        "reverse_validation": rv,
        "outputs": {
            "scores": str(output_dir / "diffusion_feature_scores.tsv"),
            "pretrained_query_scores": str(output_dir / "pretrained_query_diffusion_feature_scores.tsv")
            if (output_dir / "pretrained_query_diffusion_feature_scores.tsv").exists()
            else "",
            "model": str(output_dir / "model.pt"),
            "reference_pretrained_model": str(output_dir / "reference_pretrained_model.pt")
            if (output_dir / "reference_pretrained_model.pt").exists()
            else "",
            "observed_conditioning_profiles": str(output_dir / "observed_conditioning_profiles.tsv"),
        },
    }
    summary_path = output_dir / "training_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_training_readme(output_dir, summary)
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--reference-h5ad", default="")
    parser.add_argument("--query-source", choices=["osdr", "archs4"], default="osdr")
    parser.add_argument("--reference-source", choices=["osdr", "archs4"], default="archs4")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-mode", default="diffusion")
    parser.add_argument("--reference-epochs", type=int, default=100)
    parser.add_argument("--query-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--learning-rate", type=float, default=5e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=0.0)
    parser.add_argument("--hidden-dims", default="512,512")
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--time-embedding-dim", type=int, default=64)
    parser.add_argument("--categorical-embedding-dim", type=int, default=8)
    parser.add_argument("--sinusoidal-time", action="store_true")
    parser.add_argument("--diffusion-timesteps", type=int, default=1000)
    parser.add_argument("--sample-steps", type=int, default=50)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--beta-schedule", default="quad")
    parser.add_argument("--beta-start", type=float, default=0.0001)
    parser.add_argument("--beta-end", type=float, default=0.02)
    parser.add_argument("--n-landmarks", type=int, default=512)
    parser.add_argument("--landmark-strategy", choices=["l1000", "l1000_or_hvg", "hvg"], default="l1000_or_hvg")
    parser.add_argument("--l1000-map", default="data/diffusion/l1000_human_to_mouse_ensembl.tsv")
    parser.add_argument("--reconstruction-alpha", type=float, default=1.0)
    parser.add_argument("--max-genes", type=int, default=0)
    parser.add_argument("--max-metric-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--no-amp", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
