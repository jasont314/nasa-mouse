"""Refresh reverse-validation metrics for saved diffusion production models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nasa_mouse_glare.io import require_import

from .data import prepare_diffusion_data
from .evaluate import reverse_validation
from .model import ConditionalDiffusionMLP
from .paths import OUTPUT_ROOT, SUMMARY_DIR
from .reconstruction import train_linear_reconstructor
from .train_diffusion import generate_full
from .summarize_results import run as summarize_run


DEFAULT_MODEL_DIRS = (
    "osdr_only",
    "archs4_pretrain_osdr_finetune",
)


def refresh_one(model_dir: Path, *, batch_size: int, sample_steps: int, seed: int, cpu: bool) -> dict:
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    model_path = model_dir / "model.pt"
    summary_path = model_dir / "training_summary.json"
    if not model_path.exists():
        return {"model_dir": str(model_dir), "status": "missing_model"}
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    args = checkpoint.get("args", {})
    prepared = prepare_diffusion_data(
        query_h5ad=args["query_h5ad"],
        query_source=args["query_source"],
        reference_h5ad=args.get("reference_h5ad", ""),
        reference_source=args.get("reference_source", "archs4"),
        max_genes=args.get("max_genes", 0),
        n_landmarks=args.get("n_landmarks", 512),
        landmark_strategy=args.get("landmark_strategy", "l1000_or_hvg"),
        l1000_map=args.get("l1000_map", "data/diffusion/l1000_human_to_mouse_ensembl.tsv"),
    )
    labels = prepared.query_obs["condition_inferred"].astype(str).to_numpy()
    if len(set(labels.tolist())) < 2:
        return {"model_dir": str(model_dir), "status": "skipped_single_label"}
    device = torch.device("cuda" if torch.cuda.is_available() and not cpu else "cpu")
    model = ConditionalDiffusionMLP(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    recon_train = prepared.reference_full if prepared.reference_full is not None else prepared.query_full
    reconstructor = train_linear_reconstructor(
        recon_train,
        prepared.landmark_indices,
        prepared.target_indices,
        alpha=args.get("reconstruction_alpha", 1.0),
    )
    fake = generate_full(
        model,
        prepared.query_categories,
        reconstructor,
        betas=checkpoint["betas"].to(device),
        sample_steps=sample_steps or args.get("sample_steps", 50),
        eta=args.get("eta", 0.0),
        batch_size=batch_size,
        device=device,
        seed=seed,
    )
    metrics = reverse_validation(prepared.query_full, fake, labels)
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        summary["reverse_validation"] = metrics
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return {"model_dir": str(model_dir), "status": "completed", "device": str(device), **metrics}


def run(args) -> Path:
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    root = Path(args.root)
    rows = [
        refresh_one(root / model_dir, batch_size=args.batch_size, sample_steps=args.sample_steps, seed=args.seed, cpu=args.cpu)
        for model_dir in args.model_dirs
    ]
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    path = SUMMARY_DIR / "diffusion_reverse_validation_refresh.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    (SUMMARY_DIR / "diffusion_reverse_validation_refresh.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    summarize_run(argparse.Namespace(root=str(root), summary_dir=str(SUMMARY_DIR)))
    print(json.dumps({"summary": str(path), "rows": rows}, indent=2))
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(OUTPUT_ROOT))
    parser.add_argument("--model-dirs", nargs="+", default=list(DEFAULT_MODEL_DIRS))
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--sample-steps", type=int, default=50)
    parser.add_argument("--seed", type=int, default=2051)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
