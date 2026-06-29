"""Generate synthetic expression from a trained conditional WGAN-GP."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .model import ConditionalWGANGP


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
    config = checkpoint.get("model_config")
    if config is None:
        raise SystemExit(f"{model_dir / 'model.pt'} does not contain model_config.")
    model = ConditionalWGANGP(
        expression_dim=int(config["expression_dim"]),
        categorical_cardinalities=[int(value) for value in config["categorical_cardinalities"]],
        noise_dim=int(config["noise_dim"]),
        hidden_dims=tuple(int(value) for value in config["hidden_dims"]),
    ).to(device)
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


def encode_profile(profile: dict[str, str], vocabularies: dict[str, list[str]],
                   covariates: tuple[str, ...], *, n: int):
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
                f"Value {value!r} not in vocabulary for {covariate}. "
                f"First allowed values: {preview}"
            )
        encoded.append(vocabulary.index(value))
    return np.tile(np.asarray(encoded, dtype="int64"), (int(n), 1))


def generate_with_noise(model, categories, noise, *, batch_size: int, device):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    categories = torch.as_tensor(categories, dtype=torch.long)
    generated = []
    with torch.no_grad():
        for start in range(0, categories.shape[0], batch_size):
            end = min(start + batch_size, categories.shape[0])
            cats = categories[start:end].to(device)
            batch_noise = noise[start:end].to(device)
            fake = model.generator(batch_noise, cats)
            generated.append(fake.detach().cpu().numpy())
    return np.concatenate(generated, axis=0).astype("float32")


def write_matrix(path: Path, matrix, genes: list[str]) -> None:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = pd.DataFrame(matrix, columns=genes)
    frame.insert(0, "synthetic_sample_id", [f"synthetic_{idx:04d}" for idx in range(matrix.shape[0])])
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def write_one(output_dir: Path, stem: str, zscore, mean, std, genes: list[str], profile: dict) -> dict:
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    log1p_cpm = zscore * std.reshape(1, -1) + mean.reshape(1, -1)
    cpm = np.expm1(log1p_cpm)
    cpm = np.maximum(cpm, 0.0)
    z_path = output_dir / f"{stem}_zscore.tsv.gz"
    log_path = output_dir / f"{stem}_log1p_cpm.tsv.gz"
    cpm_path = output_dir / f"{stem}_cpm.tsv.gz"
    write_matrix(z_path, zscore, genes)
    write_matrix(log_path, log1p_cpm, genes)
    write_matrix(cpm_path, cpm, genes)
    profile_path = output_dir / f"{stem}_profile.json"
    profile_path.write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
    return {
        "zscore": str(z_path),
        "log1p_cpm": str(log_path),
        "cpm": str(cpm_path),
        "profile": str(profile_path),
    }


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
    mean = stats["mean"].astype("float32")
    std = stats["std"].astype("float32")

    base_profile = default_profile(model_dir, covariates)
    base_profile.update(parse_set(args.set))
    if args.condition:
        base_profile["wgan_condition"] = args.condition
    n = int(args.n)
    rng = torch.Generator(device=device)
    rng.manual_seed(int(args.seed))
    noise = torch.randn((n, int(checkpoint["model_config"]["noise_dim"])), generator=rng, device=device)

    outputs = {}
    if args.counterfactual:
        conditions = args.counterfactual
        if len(conditions) != 2:
            raise SystemExit("--counterfactual requires exactly two condition values.")
        first_profile = dict(base_profile, wgan_condition=conditions[0])
        second_profile = dict(base_profile, wgan_condition=conditions[1])
        first_categories = encode_profile(first_profile, vocabularies, covariates, n=n)
        second_categories = encode_profile(second_profile, vocabularies, covariates, n=n)
        first = generate_with_noise(model, first_categories, noise, batch_size=args.batch_size, device=device)
        second = generate_with_noise(model, second_categories, noise, batch_size=args.batch_size, device=device)
        outputs[conditions[0]] = write_one(output_dir, conditions[0], first, mean, std, genes, first_profile)
        outputs[conditions[1]] = write_one(output_dir, conditions[1], second, mean, std, genes, second_profile)
        delta = (second * std.reshape(1, -1) + mean.reshape(1, -1)) - (
            first * std.reshape(1, -1) + mean.reshape(1, -1)
        )
        delta_path = output_dir / f"{conditions[1]}_minus_{conditions[0]}_mean_log1p_cpm_delta.tsv"
        frame = __import__("pandas").DataFrame(
            {"gene": genes, "mean_log1p_cpm_delta": delta.mean(axis=0)}
        )
        frame.to_csv(delta_path, sep="\t", index=False)
        outputs["delta"] = str(delta_path)
    else:
        categories = encode_profile(base_profile, vocabularies, covariates, n=n)
        generated = generate_with_noise(model, categories, noise, batch_size=args.batch_size, device=device)
        stem = base_profile.get("wgan_condition", "synthetic")
        outputs[stem] = write_one(output_dir, stem, generated, mean, std, genes, base_profile)

    summary = {
        "model_dir": str(model_dir),
        "output_dir": str(output_dir),
        "device": str(device),
        "n": n,
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
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
