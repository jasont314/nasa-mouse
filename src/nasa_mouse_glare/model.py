"""GLARE-style sparse autoencoder pretraining and fine-tuning for mouse data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import dense_matrix, load_matrix_bundle, require_import


def _imports():
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
    optim = require_import("torch.optim", "pip install -r requirements-nasa-mouse-glare.txt")
    sklearn = require_import("sklearn.preprocessing", "pip install -r requirements-nasa-mouse-glare.txt")
    return np, torch, nn, optim, sklearn.StandardScaler


def build_classes():
    _, torch, nn, _, _ = _imports()

    class SparseAutoEncoder(nn.Module):
        def __init__(self, input_dim: int, hidden_layers: tuple[int, ...] = (128, 64, 32, 16)):
            super().__init__()
            encoder_layers = []
            in_dim = input_dim
            for out_dim in hidden_layers:
                encoder_layers.extend([nn.Linear(in_dim, out_dim), nn.LayerNorm(out_dim), nn.ELU()])
                in_dim = out_dim
            decoder_layers = []
            reversed_layers = list(hidden_layers[::-1][1:]) + [input_dim]
            for out_dim in reversed_layers[:-1]:
                decoder_layers.extend([nn.Linear(in_dim, out_dim), nn.LayerNorm(out_dim), nn.ELU()])
                in_dim = out_dim
            decoder_layers.extend([nn.Linear(in_dim, reversed_layers[-1]), nn.Sigmoid()])
            self.encoder = nn.Sequential(*encoder_layers)
            self.decoder = nn.Sequential(*decoder_layers)

        def forward(self, x):
            encoded = self.encoder(x)
            return self.decoder(encoded)

    class Adapter(nn.Module):
        def __init__(self, input_dim: int, output_dim: int):
            super().__init__()
            self.linear = nn.Linear(input_dim, output_dim)
            self.activation = nn.ELU()

        def forward(self, x):
            return self.activation(self.linear(x))

    return SparseAutoEncoder, Adapter


def standardize_matrix(matrix, max_dense_gb: float):
    np, _, _, _, StandardScaler = _imports()
    X = dense_matrix(matrix, max_dense_gb=max_dense_gb)
    X = StandardScaler().fit_transform(X)
    return np.asarray(X, dtype="float32")


def train_autoencoder(
    X,
    output_path: str | Path,
    epochs: int = 30,
    batch_size: int = 16,
    learning_rate: float = 1e-3,
    weight_decay: float = 1e-4,
    sparsity_penalty: float = 1e-5,
    device_name: str = "auto",
    seed: int = 2026,
    pretrained_state: str | Path | None = None,
    adapter_input_dim: int | None = None,
):
    np, torch, nn, optim, _ = _imports()
    SparseAutoEncoder, Adapter = build_classes()

    torch.manual_seed(seed)
    device = torch.device(
        "cuda:0" if device_name == "auto" and torch.cuda.is_available() else "cpu"
    )
    if device_name != "auto":
        device = torch.device(device_name)

    X = torch.tensor(np.asarray(X, dtype="float32"))
    adapter = None
    if adapter_input_dim is not None:
        adapter = Adapter(X.shape[1], adapter_input_dim)
        X = adapter(X).detach()

    loader = torch.utils.data.DataLoader(X, batch_size=batch_size, shuffle=True)
    model = SparseAutoEncoder(X.shape[1]).to(device)
    if pretrained_state is not None:
        model.load_state_dict(torch.load(pretrained_state, map_location=device))

    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.MSELoss()
    history = []
    for epoch in range(epochs):
        model.train()
        total_loss = 0.0
        for batch in loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            output = model(batch)
            loss = criterion(output, batch)
            encoded = model.encoder(batch)
            loss = loss + sparsity_penalty * torch.mean(torch.abs(encoded))
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += float(loss.item())
        average_loss = total_loss / max(1, len(loader))
        history.append({"epoch": epoch + 1, "loss": average_loss})
        print(f"epoch {epoch + 1}/{epochs}\tloss={average_loss:.6f}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), output_path)
    history_path = output_path.with_suffix(".history.json")
    history_path.write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
    return output_path


def encode_bundle(
    bundle_manifest: str | Path,
    model_state: str | Path,
    output_npy: str | Path,
    max_dense_gb: float = 8.0,
    device_name: str = "auto",
):
    np, torch, _, _, _ = _imports()
    SparseAutoEncoder, _ = build_classes()

    bundle = load_matrix_bundle(bundle_manifest)
    X = standardize_matrix(bundle.matrix, max_dense_gb=max_dense_gb)
    device = torch.device(
        "cuda:0" if device_name == "auto" and torch.cuda.is_available() else "cpu"
    )
    if device_name != "auto":
        device = torch.device(device_name)

    model = SparseAutoEncoder(X.shape[1]).to(device)
    model.load_state_dict(torch.load(model_state, map_location=device))
    model.eval()
    with torch.no_grad():
        encoded = model.encoder(torch.tensor(X, dtype=torch.float32).to(device))
    output_npy = Path(output_npy)
    output_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_npy, encoded.cpu().numpy())
    output_npy.with_suffix(".genes.tsv").write_text(
        "\n".join(bundle.genes) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Train GLARE-style SAE models.")
    sub = parser.add_subparsers(dest="command", required=True)

    pre = sub.add_parser("pretrain")
    pre.add_argument("--matrix", required=True, help="Pretraining matrix manifest")
    pre.add_argument("--output", required=True)
    pre.add_argument("--epochs", type=int, default=30)
    pre.add_argument("--batch-size", type=int, default=16)
    pre.add_argument("--max-dense-gb", type=float, default=8.0)
    pre.add_argument("--device", default="auto")

    ft = sub.add_parser("finetune")
    ft.add_argument("--target-matrix", required=True)
    ft.add_argument("--pretrained", required=True)
    ft.add_argument("--pretrain-input-dim", type=int, required=True)
    ft.add_argument("--output", required=True)
    ft.add_argument("--epochs", type=int, default=30)
    ft.add_argument("--batch-size", type=int, default=16)
    ft.add_argument("--max-dense-gb", type=float, default=8.0)
    ft.add_argument("--device", default="auto")

    enc = sub.add_parser("encode")
    enc.add_argument("--matrix", required=True)
    enc.add_argument("--model", required=True)
    enc.add_argument("--output", required=True)
    enc.add_argument("--max-dense-gb", type=float, default=8.0)
    enc.add_argument("--device", default="auto")

    args = parser.parse_args()
    if args.command == "pretrain":
        bundle = load_matrix_bundle(args.matrix)
        X = standardize_matrix(bundle.matrix, max_dense_gb=args.max_dense_gb)
        train_autoencoder(
            X,
            args.output,
            epochs=args.epochs,
            batch_size=args.batch_size,
            device_name=args.device,
        )
        print(f"pretrain_input_dim={X.shape[1]}")
    elif args.command == "finetune":
        bundle = load_matrix_bundle(args.target_matrix)
        X = standardize_matrix(bundle.matrix, max_dense_gb=args.max_dense_gb)
        train_autoencoder(
            X,
            args.output,
            epochs=args.epochs,
            batch_size=args.batch_size,
            device_name=args.device,
            pretrained_state=args.pretrained,
            adapter_input_dim=args.pretrain_input_dim,
            weight_decay=0.0,
        )
    elif args.command == "encode":
        encode_bundle(args.matrix, args.model, args.output, args.max_dense_gb, args.device)


if __name__ == "__main__":
    main()
