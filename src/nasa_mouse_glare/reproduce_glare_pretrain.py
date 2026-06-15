"""Reproduce GLARE single-cell pretraining with the released config."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.io import mmread
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader


class SparseAutoEncoder(nn.Module):
    """Architecture used by GLARE's released `representation_learning.py`."""

    def __init__(self, input_dim: int):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.LayerNorm(128),
            nn.ELU(),
            nn.Linear(128, 64),
            nn.LayerNorm(64),
            nn.ELU(),
            nn.Linear(64, 32),
            nn.LayerNorm(32),
            nn.ELU(),
            nn.Linear(32, 16),
            nn.ELU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(16, 32),
            nn.LayerNorm(32),
            nn.ELU(),
            nn.Linear(32, 64),
            nn.LayerNorm(64),
            nn.ELU(),
            nn.Linear(64, 128),
            nn.LayerNorm(128),
            nn.ELU(),
            nn.Linear(128, input_dim),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.decoder(self.encoder(x))


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def matrixmarket_to_dense(path: Path) -> np.ndarray:
    matrix = mmread(path)
    if hasattr(matrix, "toarray"):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)


def write_epoch_logs(output_dir: Path, records: list[dict]) -> None:
    csv_path = output_dir / "epoch_losses.csv"
    json_path = output_dir / "epoch_losses.json"
    fields = [
        "epoch",
        "loss",
        "best_loss",
        "elapsed_seconds",
        "elapsed",
        "epoch_seconds",
        "epoch_elapsed",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    json_path.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GLARE's released single-cell pretraining config."
    )
    parser.add_argument(
        "--input",
        default=(
            "assets/glare_original/"
            "E-CURD-5.aggregated_filtered_normalised_counts.mtx"
        ),
        help="Path to E-CURD-5 normalized MatrixMarket file.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/glare_original_pretrain_config5",
        help="Directory for weights and epoch logs.",
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1996)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_start = time.perf_counter()

    torch.manual_seed(args.seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    config = {
        "hidden_layers": [128, 64, 32, 16],
        "layer_norm": True,
        "activation_fn": "ELU",
        "sparsity_penalty": 1e-5,
        "learning_rate": 1e-3,
        "weight_decay": 1e-4,
        "batch_size": args.batch_size,
        "epochs": args.epochs,
        "seed": args.seed,
        "device": str(device),
        "input": str(input_path),
    }
    log(f"Loading MatrixMarket: {input_path}")
    data = matrixmarket_to_dense(input_path)
    log(f"Loaded data shape={data.shape}")

    scaler = StandardScaler()
    x = scaler.fit_transform(data)
    x = torch.tensor(x, dtype=torch.float32)
    loader = DataLoader(
        x,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=4,
        pin_memory=torch.cuda.is_available(),
    )

    model = SparseAutoEncoder(x.shape[1]).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(
        model.parameters(),
        lr=config["learning_rate"],
        weight_decay=config["weight_decay"],
    )

    records = []
    best_loss = float("inf")
    log(f"Starting GLARE original pretraining config on device={device}")
    for epoch in range(args.epochs):
        epoch_start = time.perf_counter()
        total_loss = 0.0
        model.train()
        for batch_data in loader:
            batch_data = batch_data.to(device)
            optimizer.zero_grad()
            outputs = model(batch_data)
            loss = criterion(outputs, batch_data)
            encoded = model.encoder[-1](batch_data)
            l1_regularization = torch.mean(torch.abs(encoded))
            loss += config["sparsity_penalty"] * l1_regularization
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
            optimizer.step()
            total_loss += loss.item()

        avg_loss = total_loss / len(loader)
        best_loss = min(best_loss, avg_loss)
        epoch_seconds = time.perf_counter() - epoch_start
        elapsed_seconds = time.perf_counter() - run_start
        record = {
            "epoch": epoch + 1,
            "loss": round(float(avg_loss), 8),
            "best_loss": round(float(best_loss), 8),
            "elapsed_seconds": round(elapsed_seconds, 3),
            "elapsed": format_elapsed(elapsed_seconds),
            "epoch_seconds": round(epoch_seconds, 3),
            "epoch_elapsed": format_elapsed(epoch_seconds),
        }
        records.append(record)
        write_epoch_logs(output_dir, records)
        log(
            f"epoch {epoch + 1}/{args.epochs} "
            f"loss={avg_loss:.8f} best={best_loss:.8f} "
            f"epoch_time={format_elapsed(epoch_seconds)}"
        )

    weights_path = output_dir / "sc_shulse_pretrained_reproduced.pth"
    torch.save(model.state_dict(), weights_path)
    summary = {
        "config": config,
        "data_shape": list(data.shape),
        "best_loss": round(float(best_loss), 8),
        "final_loss": records[-1]["loss"],
        "weights": str(weights_path),
        "epoch_logs_csv": str(output_dir / "epoch_losses.csv"),
        "epoch_logs_json": str(output_dir / "epoch_losses.json"),
        "elapsed_seconds": round(time.perf_counter() - run_start, 3),
        "elapsed": format_elapsed(time.perf_counter() - run_start),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    log(f"Saved weights: {weights_path}")
    log(f"Saved summary: {output_dir / 'summary.json'}")


if __name__ == "__main__":
    main()
