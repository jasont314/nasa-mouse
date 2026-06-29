"""Generate synthetic expression from a trained conditional diffusion model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .diffusion import sample
from .model import ConditionalDiffusionMLP


def parse_set(values: list[str]) -> dict[str, str]:
    profile = {}
    for item in values:
        if "=" not in item:
            raise SystemExit(f"Expected --set key=value, got {item!r}")
        key, value = item.split("=", 1)
        profile[key] = value
    return profile


def load_checkpoint(model_dir: Path, *, device):
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    checkpoint = torch.load(model_dir / "model.pt", map_location=device, weights_only=False)
    model = ConditionalDiffusionMLP(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return checkpoint, model


def default_profile(model_dir: Path, covariates: tuple[str, ...]) -> dict[str, str]:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    path = model_dir / "observed_conditioning_profiles.tsv"
    if not path.exists():
        return {}
    profiles = pd.read_csv(path, sep="\t")
    query = profiles.loc[profiles.get("training_source", "").eq("query")]
    row = (query if not query.empty else profiles).iloc[0]
    return {covariate: str(row[covariate]) for covariate in covariates if covariate in row}


def encode_profile(profile: dict[str, str], vocabularies: dict[str, list[str]], covariates: tuple[str, ...], *, n: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    encoded = []
    for covariate in covariates:
        value = profile.get(covariate)
        if value is None:
            raise SystemExit(f"Missing covariate {covariate!r}; pass --set {covariate}=...")
        vocabulary = vocabularies[covariate]
        if value not in vocabulary:
            preview = ", ".join(vocabulary[:12])
            raise SystemExit(
                f"Value {value!r} not in vocabulary for {covariate}. First allowed values: {preview}"
            )
        encoded.append(vocabulary.index(value))
    return np.tile(np.asarray(encoded, dtype="int64"), (int(n), 1))


def reconstruct_full(landmarks, checkpoint):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    landmark_indices = np.asarray(checkpoint["landmark_indices"], dtype=int)
    target_indices = np.asarray(checkpoint["target_indices"], dtype=int)
    full = np.zeros((landmarks.shape[0], len(checkpoint["genes"])), dtype="float32")
    full[:, landmark_indices] = landmarks
    if len(target_indices):
        rec = checkpoint.get("reconstruction", {})
        coef = np.asarray(rec.get("coef"), dtype="float32")
        intercept = np.asarray(rec.get("intercept"), dtype="float32")
        full[:, target_indices] = landmarks.dot(coef.T) + intercept.reshape(1, -1)
    return full


def generate_landmarks(model, categories, betas, *, n: int, sample_steps: int, eta: float, seed: int, batch_size: int, device):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    categories = torch.as_tensor(categories, dtype=torch.long)
    rng = torch.Generator(device=device)
    rng.manual_seed(int(seed))
    generated = []
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        noise = torch.randn((end - start, model.expression_dim), generator=rng, device=device)
        batch = categories[start:end]
        landmarks = sample(model, batch, betas=betas, sample_steps=sample_steps, eta=eta, noise=noise, device=device)
        generated.append(landmarks.detach().cpu().numpy())
    return np.concatenate(generated, axis=0).astype("float32")


def write_matrix(path: Path, matrix, genes: list[str]) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = pd.DataFrame(matrix, columns=genes)
    frame.insert(0, "synthetic_sample_id", [f"synthetic_{idx:05d}" for idx in range(matrix.shape[0])])
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def scaled_to_log1p_cpm(full_scaled, center, scale):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    raw = full_scaled * scale.reshape(1, -1) + center.reshape(1, -1)
    max_log1p_cpm = float(np.log1p(1_000_000.0))
    finite = np.isfinite(raw)
    finite_values = raw[finite]
    report = {
        "n_values": int(raw.size),
        "n_nonfinite": int((~finite).sum()),
        "n_clipped_low": int((raw < 0.0).sum()),
        "n_clipped_high": int((raw > max_log1p_cpm).sum()),
        "raw_log1p_cpm_min": float(finite_values.min()) if finite_values.size else float("nan"),
        "raw_log1p_cpm_max": float(finite_values.max()) if finite_values.size else float("nan"),
        "max_valid_log1p_cpm": max_log1p_cpm,
    }
    clipped = np.nan_to_num(raw, nan=0.0, posinf=max_log1p_cpm, neginf=0.0)
    return np.clip(clipped, 0.0, max_log1p_cpm).astype("float32"), report


def write_one(output_dir: Path, stem: str, full_scaled, center, scale, genes: list[str], profile: dict) -> dict:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    log1p_cpm, clip_report = scaled_to_log1p_cpm(full_scaled, center, scale)
    cpm = np.expm1(log1p_cpm).astype("float32")
    paths = {
        "scaled": output_dir / f"{stem}_scaled.tsv.gz",
        "log1p_cpm": output_dir / f"{stem}_log1p_cpm.tsv.gz",
        "cpm": output_dir / f"{stem}_cpm.tsv.gz",
        "profile": output_dir / f"{stem}_profile.json",
        "clip_report": output_dir / f"{stem}_clip_report.json",
    }
    write_matrix(paths["scaled"], full_scaled, genes)
    write_matrix(paths["log1p_cpm"], log1p_cpm, genes)
    write_matrix(paths["cpm"], cpm, genes)
    paths["profile"].write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
    paths["clip_report"].write_text(json.dumps(clip_report, indent=2) + "\n", encoding="utf-8")
    return {key: str(path) for key, path in paths.items()}


def run(args) -> Path:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    model_dir = Path(args.model_dir)
    output_dir = Path(args.output_dir)
    checkpoint, model = load_checkpoint(model_dir, device=device)
    covariates = tuple(checkpoint["categorical_covariates"])
    vocabularies = checkpoint["vocabularies"]
    genes = list(checkpoint["genes"])
    stats = np.load(model_dir / "normalization_stats.npz")
    center = stats["center"].astype("float32")
    scale = stats["scale"].astype("float32")
    betas = checkpoint["betas"].to(device)

    base_profile = default_profile(model_dir, covariates)
    base_profile.update(parse_set(args.set))
    if args.condition:
        base_profile["wgan_condition"] = args.condition
    n = int(args.n)
    outputs = {}
    if args.counterfactual:
        first_condition, second_condition = args.counterfactual
        first_profile = dict(base_profile, wgan_condition=first_condition)
        second_profile = dict(base_profile, wgan_condition=second_condition)
        first_categories = encode_profile(first_profile, vocabularies, covariates, n=n)
        second_categories = encode_profile(second_profile, vocabularies, covariates, n=n)
        first_landmarks = generate_landmarks(
            model,
            first_categories,
            betas,
            n=n,
            sample_steps=args.sample_steps,
            eta=args.eta,
            seed=args.seed,
            batch_size=args.batch_size,
            device=device,
        )
        second_landmarks = generate_landmarks(
            model,
            second_categories,
            betas,
            n=n,
            sample_steps=args.sample_steps,
            eta=args.eta,
            seed=args.seed,
            batch_size=args.batch_size,
            device=device,
        )
        first_full = reconstruct_full(first_landmarks, checkpoint)
        second_full = reconstruct_full(second_landmarks, checkpoint)
        outputs[first_condition] = write_one(output_dir, first_condition, first_full, center, scale, genes, first_profile)
        outputs[second_condition] = write_one(output_dir, second_condition, second_full, center, scale, genes, second_profile)
        first_log1p, _ = scaled_to_log1p_cpm(first_full, center, scale)
        second_log1p, _ = scaled_to_log1p_cpm(second_full, center, scale)
        delta = second_log1p - first_log1p
        delta_path = output_dir / f"{second_condition}_minus_{first_condition}_mean_log1p_cpm_delta.tsv"
        pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
        pd.DataFrame({"gene": genes, "mean_log1p_cpm_delta": delta.mean(axis=0)}).to_csv(delta_path, sep="\t", index=False)
        outputs["delta"] = str(delta_path)
    else:
        categories = encode_profile(base_profile, vocabularies, covariates, n=n)
        landmarks = generate_landmarks(
            model,
            categories,
            betas,
            n=n,
            sample_steps=args.sample_steps,
            eta=args.eta,
            seed=args.seed,
            batch_size=args.batch_size,
            device=device,
        )
        full = reconstruct_full(landmarks, checkpoint)
        stem = base_profile.get("wgan_condition", "synthetic")
        outputs[stem] = write_one(output_dir, stem, full, center, scale, genes, base_profile)
    summary = {
        "model_dir": str(model_dir),
        "output_dir": str(output_dir),
        "device": str(device),
        "n": n,
        "sample_steps": int(args.sample_steps),
        "eta": float(args.eta),
        "covariates": list(covariates),
        "outputs": outputs,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "synthetic_generation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n", type=int, default=32)
    parser.add_argument("--condition", default="")
    parser.add_argument("--counterfactual", nargs=2, metavar=("FROM", "TO"))
    parser.add_argument("--set", action="append", default=[], help="Covariate override as key=value.")
    parser.add_argument("--sample-steps", type=int, default=50)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
