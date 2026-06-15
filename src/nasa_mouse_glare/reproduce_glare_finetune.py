"""Reproduce GLARE OSD-120 fine-tuning with the released config."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader


KMEANS_KWARGS = {
    "init": "k-means++",
    "n_init": 10,
    "random_state": 1,
}
OUTLIER_GENE_IDS = ["AT3G41768", "ATMG00020", "AT1G07590"]


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


class Adapter(nn.Module):
    def __init__(self, input_dim: int, output_dim: int):
        super().__init__()
        self.linear = nn.Linear(input_dim, output_dim)
        self.activation = nn.ELU()

    def forward(self, x):
        return self.activation(self.linear(x))


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


def restructure_data(
    df: pd.DataFrame,
    condition_keywords: list[str],
    first_keywords: list[str],
    secondary_keywords: list[str] | None = None,
    replicate_identifier: str = "Rep",
    id_col: str = "gene_id",
) -> pd.DataFrame:
    """Match GLARE's CARA wide-to-location restructuring."""
    expr_cols = [col for col in df.columns if col != id_col]
    metadata = []
    for col in expr_cols:
        col_str = col if isinstance(col, str) else col[0]
        components = re.split(r"[,_]", col_str)
        cond = next(
            (
                c
                for c in condition_keywords
                if any(c.lower() in comp.lower() for comp in components)
            ),
            None,
        )
        genotype = next(
            (
                g
                for g in first_keywords
                if any(g.lower() in comp.lower() for comp in components)
            ),
            None,
        )
        secondary = None
        if secondary_keywords:
            secondary = next(
                (
                    s
                    for s in secondary_keywords
                    if any(s.lower() in comp.lower() for comp in components)
                ),
                None,
            )
        rep = None
        for comp in components:
            if replicate_identifier.lower() in comp.lower():
                match = re.search(r"(\d+)", comp)
                if match:
                    rep = f"{replicate_identifier}{match.group(1)}"
                    break
        if not all([cond, genotype, rep]):
            raise ValueError(f"Could not parse all metadata from column: '{col_str}'")
        label_parts = [rep, genotype]
        if secondary:
            label_parts.append(secondary)
        metadata.append(
            {
                "original_col": col,
                "condition": cond,
                "label": "_".join(label_parts),
            }
        )

    meta_df = pd.DataFrame(metadata)
    long_df = df.melt(
        id_vars=[id_col],
        value_vars=meta_df["original_col"].tolist(),
        var_name="original_col",
        value_name="expression",
    )
    long_df = long_df.merge(meta_df, on="original_col")
    wide_df = long_df.pivot_table(
        index=[id_col, "condition"],
        columns="label",
        values="expression",
    ).reset_index()
    wide_df.rename(columns={"condition": "Location"}, inplace=True)
    labels = sorted([col for col in wide_df.columns if col not in [id_col, "Location"]])
    return wide_df[[id_col] + labels + ["Location"]]


def preprocess_location(df: pd.DataFrame, location: str, output_dir: Path) -> pd.DataFrame:
    """Run GLARE's PCA/k-means outlier step and remove fixed CARA outliers."""
    loc_value = 1 if location == "FLT" else 0
    loc_df = df[df["Location"] == loc_value].reset_index(drop=True)
    pca_values = PCA(n_components=3).fit_transform(loc_df.iloc[:, 1:-1])
    num_cluster = 5 if location == "FLT" else 4
    clusters = KMeans(n_clusters=num_cluster, **KMEANS_KWARGS).fit_predict(pca_values)
    pca_df = pd.DataFrame(pca_values, columns=["x", "y", "z"])
    pca_df["cluster"] = clusters.astype(str)
    pca_df["gene_id"] = loc_df["gene_id"]
    pca_df.to_csv(output_dir / f"outlier_detection_df_{location}.csv", index=False)
    clean_df = loc_df.drop(
        loc_df[loc_df["gene_id"].isin(OUTLIER_GENE_IDS)].index
    ).reset_index(drop=True)
    clean_df.to_csv(output_dir / f"clean_glds120_{location}.csv", index=False)
    return clean_df


def load_state_dict(path: Path):
    checkpoint = torch.load(path, map_location="cpu")
    if isinstance(checkpoint, dict) and "model_state_dict" in checkpoint:
        return checkpoint["model_state_dict"]
    return checkpoint


def write_epoch_logs(output_dir: Path, location: str, records: list[dict]) -> None:
    csv_path = output_dir / f"{location}_epoch_losses.csv"
    json_path = output_dir / f"{location}_epoch_losses.json"
    fields = [
        "location",
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


def finetune_location(
    clean_df: pd.DataFrame,
    location: str,
    pretrained_weights: Path,
    output_dir: Path,
    device: torch.device,
    pi_dim: int,
    epochs: int,
    batch_size: int,
) -> dict:
    run_start = time.perf_counter()
    x = StandardScaler().fit_transform(clean_df.iloc[:, 1:-1])
    x = torch.tensor(x, dtype=torch.float32)
    adapter = Adapter(x.shape[1], pi_dim)
    x = adapter(x).clone().detach()
    loader = DataLoader(
        x,
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        pin_memory=torch.cuda.is_available(),
    )

    model = SparseAutoEncoder(pi_dim).to(device)
    model.load_state_dict(load_state_dict(pretrained_weights))
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    sparsity_penalty = 1e-5

    records = []
    best_loss = float("inf")
    best_epoch = 0
    log(
        f"{location}: starting fine-tune rows={clean_df.shape[0]} "
        f"features={clean_df.shape[1] - 2} adapter_dim={pi_dim}"
    )
    for epoch in range(epochs):
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
            loss += sparsity_penalty * l1_regularization
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        avg_loss = total_loss / len(loader)
        if avg_loss < best_loss:
            best_loss = avg_loss
            best_epoch = epoch + 1
        epoch_seconds = time.perf_counter() - epoch_start
        elapsed_seconds = time.perf_counter() - run_start
        record = {
            "location": location,
            "epoch": epoch + 1,
            "loss": round(float(avg_loss), 8),
            "best_loss": round(float(best_loss), 8),
            "elapsed_seconds": round(elapsed_seconds, 3),
            "elapsed": format_elapsed(elapsed_seconds),
            "epoch_seconds": round(epoch_seconds, 3),
            "epoch_elapsed": format_elapsed(epoch_seconds),
        }
        records.append(record)
        write_epoch_logs(output_dir, location, records)
        if (epoch + 1) % 5 == 0:
            torch.save(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                },
                output_dir / f"{location}_finetune_checkpoint_epoch_{epoch + 1}.pth",
            )
        log(
            f"{location}: epoch {epoch + 1}/{epochs} "
            f"loss={avg_loss:.8f} best={best_loss:.8f} "
            f"epoch_time={format_elapsed(epoch_seconds)}"
        )

    final_weights = output_dir / f"OSD120_{location}_finetuned_reproduced.pth"
    torch.save(
        {
            "epoch": epochs - 1,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
        },
        final_weights,
    )
    model.eval()
    with torch.no_grad():
        representation = model.encoder(x.to(device)).detach().cpu().numpy()
    representation_path = output_dir / f"FTSAE_{location}_GLDS120_reproduced.npy"
    np.save(representation_path, representation)

    return {
        "location": location,
        "rows": int(clean_df.shape[0]),
        "features": int(clean_df.shape[1] - 2),
        "best_loss": round(float(best_loss), 8),
        "best_epoch": best_epoch,
        "final_loss": records[-1]["loss"],
        "final_weights": str(final_weights),
        "representation": str(representation_path),
        "epoch_logs_csv": str(output_dir / f"{location}_epoch_losses.csv"),
        "epoch_logs_json": str(output_dir / f"{location}_epoch_losses.json"),
        "elapsed_seconds": round(time.perf_counter() - run_start, 3),
        "elapsed": format_elapsed(time.perf_counter() - run_start),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GLARE's released OSD-120 fine-tuning config."
    )
    parser.add_argument(
        "--data",
        default="assets/glare_original/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
        help="Path to GLDS-120 normalized counts CSV.",
    )
    parser.add_argument(
        "--pretrained-weights",
        default="outputs/glare_original_pretrain_config5/sc_shulse_pretrained_reproduced.pth",
        help="Path to GLARE pretraining weights.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/glare_original_finetune_osd120",
        help="Directory for weights, representations, and logs.",
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--pretrained-input-dim", type=int, default=3552)
    parser.add_argument("--seed", type=int, default=1996)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_start = time.perf_counter()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    data_path = Path(args.data)
    pretrained_weights = Path(args.pretrained_weights)
    log(f"Loading GLDS-120 normalized counts: {data_path}")
    df = pd.read_csv(data_path)
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "gene_id"})
    log(f"Loaded GLDS-120 shape={df.shape}")

    glds120 = restructure_data(
        df=df,
        condition_keywords=["FLT", "GC"],
        first_keywords=["Col-0-PhyD", "Col-0", "Ws"],
        secondary_keywords=["Alight", "dark"],
        replicate_identifier="Rep",
        id_col="gene_id",
    )
    glds120["Location"] = glds120["Location"].map({"FLT": 1, "GC": 0})
    glds120 = glds120.infer_objects(copy=False)
    glds120.to_csv(output_dir / "glds120_restructured.csv", index=False)
    log(f"Restructured GLDS-120 shape={glds120.shape}")

    clean_flt = preprocess_location(glds120, "FLT", output_dir)
    clean_gc = preprocess_location(glds120, "GC", output_dir)

    config = {
        "pretrained_weights": str(pretrained_weights),
        "pretrained_input_dim": args.pretrained_input_dim,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "device": str(device),
        "learning_rate": 1e-3,
        "weight_decay": 0,
        "sparsity_penalty": 1e-5,
        "gradient_clipping": False,
        "adapter": "Linear(target_features -> pretrained_input_dim) + ELU",
    }
    summaries = [
        finetune_location(
            clean_flt,
            "FLT",
            pretrained_weights,
            output_dir,
            device,
            args.pretrained_input_dim,
            args.epochs,
            args.batch_size,
        ),
        finetune_location(
            clean_gc,
            "GC",
            pretrained_weights,
            output_dir,
            device,
            args.pretrained_input_dim,
            args.epochs,
            args.batch_size,
        ),
    ]
    summary = {
        "config": config,
        "data": str(data_path),
        "raw_shape": list(df.shape),
        "restructured_shape": list(glds120.shape),
        "locations": summaries,
        "elapsed_seconds": round(time.perf_counter() - run_start, 3),
        "elapsed": format_elapsed(time.perf_counter() - run_start),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    log(f"Saved summary: {summary_path}")


if __name__ == "__main__":
    main()
