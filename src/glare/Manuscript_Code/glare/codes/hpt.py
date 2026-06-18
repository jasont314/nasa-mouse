# Utils
import pandas as pd
import numpy as np
import argparse
import csv
import itertools
import json
import time
from pathlib import Path

# Representation learning via SAE
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.preprocessing import StandardScaler

from scipy.io import mmread

# Set environment
# os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:64'
torch.cuda.empty_cache()
torch.manual_seed(2023)


def log(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def format_elapsed(seconds):
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def activation_name(activation_fn):
    return getattr(activation_fn, "__name__", str(activation_fn))


ACTIVATION_FNS = {
    "ELU": nn.ELU,
    "ReLU": nn.ReLU,
    "LeakyReLU": nn.LeakyReLU,
}


def format_params(params):
    if params is None:
        return None
    formatted = dict(params)
    formatted["activation_fn"] = activation_name(formatted["activation_fn"])
    return formatted


def restore_params(params):
    restored = dict(params)
    activation_fn = restored["activation_fn"]
    if isinstance(activation_fn, str):
        if activation_fn not in ACTIVATION_FNS:
            raise ValueError(f"Unsupported activation function in results: {activation_fn}")
        restored["activation_fn"] = ACTIVATION_FNS[activation_fn]
    return restored


class ConfigResultLogger:
    csv_fields = [
        "timestamp",
        "stage",
        "config_index",
        "total_configs",
        "data_shape",
        "device",
        "hidden_layers",
        "layer_norm",
        "activation_fn",
        "sparsity_penalty",
        "learning_rate",
        "weight_decay",
        "best_loss",
        "epochs_run",
        "elapsed_seconds",
        "elapsed",
        "early_stopped",
        "is_stage_best",
    ]

    def __init__(self, output_dir, prefix="hpt_config_results", resume=False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.csv_path = self.output_dir / f"{prefix}.csv"
        self.json_path = self.output_dir / f"{prefix}.json"
        self.records = []
        if resume:
            if self.json_path.exists():
                self.records = json.loads(self.json_path.read_text(encoding="utf-8"))
            if self.records and not self.csv_path.exists():
                self._rewrite_csv()
        else:
            self.csv_path.unlink(missing_ok=True)
            self.json_path.unlink(missing_ok=True)

    def write(self, record):
        self.records.append(record)
        write_header = not self.csv_path.exists()
        with self.csv_path.open("a", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.csv_fields)
            if write_header:
                writer.writeheader()
            writer.writerow(self._csv_row(record))
        self.json_path.write_text(
            json.dumps(self.records, indent=2) + "\n", encoding="utf-8"
        )

    def _csv_row(self, record):
        params = record["params"]
        row = {
            "timestamp": record["timestamp"],
            "stage": record["stage"],
            "config_index": record["config_index"],
            "total_configs": record["total_configs"],
            "data_shape": json.dumps(record["data_shape"]),
            "device": record["device"],
            "hidden_layers": json.dumps(params["hidden_layers"]),
            "layer_norm": params["layer_norm"],
            "activation_fn": params["activation_fn"],
            "sparsity_penalty": params["sparsity_penalty"],
            "learning_rate": params["learning_rate"],
            "weight_decay": params["weight_decay"],
            "best_loss": record["best_loss"],
            "epochs_run": record["epochs_run"],
            "elapsed_seconds": record["elapsed_seconds"],
            "elapsed": record["elapsed"],
            "early_stopped": record["early_stopped"],
            "is_stage_best": record["is_stage_best"],
        }
        return row

    def _rewrite_csv(self):
        with self.csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.csv_fields)
            writer.writeheader()
            for record in self.records:
                writer.writerow(self._csv_row(record))

    def stage_records(self, stage_name):
        return [record for record in self.records if record["stage"] == stage_name]

    def paths(self):
        return self.csv_path, self.json_path

# Build Sparse AutoEncoder with chosen hyperparmeter
class BuildSAE(nn.Module):
    def __init__(self, input_dim, hidden_layers, layer_norm, activation_fn):
        super(BuildSAE, self).__init__()
        self.encoder = self._build_encoder(input_dim, hidden_layers, layer_norm, activation_fn)
        self.decoder = self._build_decoder(hidden_layers[-1], hidden_layers[::-1][1:] + [input_dim], layer_norm, activation_fn)

    def _build_encoder(self, input_dim, layers, layer_norm, activation_fn):
        seq = []
        in_dim = input_dim
        for out_dim in layers:
            seq.append(nn.Linear(in_dim, out_dim))
            if layer_norm:
                seq.append(nn.LayerNorm(out_dim))
            seq.append(activation_fn())
            in_dim = out_dim
        return nn.Sequential(*seq)

    def _build_decoder(self, input_dim, layers, layer_norm, activation_fn):
        seq = []
        in_dim = input_dim
        for out_dim in layers[:-1]:
            seq.append(nn.Linear(in_dim, out_dim))
            if layer_norm:
                seq.append(nn.LayerNorm(out_dim))
            seq.append(activation_fn())
            in_dim = out_dim
        # Add final layer with Sigmoid activation
        seq.append(nn.Linear(in_dim, layers[-1]))
        seq.append(nn.Sigmoid())
        return nn.Sequential(*seq)

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded


class Adapter(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(Adapter, self).__init__()
        self.linear = nn.Linear(input_dim, output_dim)
        self.activation = nn.ELU()

    def forward(self, x):
        x = self.linear(x)
        x = self.activation(x)
        return x


# Grid Search for Hyperparameter Tuning
def grid_search(
    data,
    device,
    pretrained_model=None,
    pi_dim=None,
    fixed_architecture=None,
    stage_name="grid_search",
    log_every_epochs=1,
    result_logger=None,
    resume=False,
    start_config=1,
    num_workers=0,
):
    param_grid = {
        "hidden_layers": [[128, 64, 32, 16], [256, 128, 64, 32], [512, 256, 128, 64]],
        "layer_norm": [True, False],
        "activation_fn": [nn.ELU, nn.ReLU, nn.LeakyReLU],
        "sparsity_penalty": [1e-5, 1e-4, 1e-3],
        "learning_rate": [1e-4, 1e-3, 1e-2],
        "weight_decay": [1e-5, 1e-4, 1e-3],
    }
    if fixed_architecture is not None:
        param_grid["hidden_layers"] = [fixed_architecture["hidden_layers"]]
        param_grid["layer_norm"] = [fixed_architecture["layer_norm"]]
        param_grid["activation_fn"] = [fixed_architecture["activation_fn"]]

    param_combinations = list(itertools.product(*param_grid.values()))
    param_keys = list(param_grid.keys())
    total_configs = len(param_combinations)
    if start_config < 1 or start_config > total_configs:
        raise ValueError(
            f"{stage_name}: start_config must be between 1 and {total_configs}; "
            f"got {start_config}"
        )

    best_params = None
    best_loss = float("inf")
    completed_indices = set()
    prior_records = result_logger.stage_records(stage_name) if result_logger else []
    if resume and prior_records:
        completed_indices = {record["config_index"] for record in prior_records}
        best_record = min(prior_records, key=lambda record: record["best_loss"])
        best_loss = float(best_record["best_loss"])
        best_params = restore_params(best_record["params"])
        log(
            f"{stage_name}: loaded {len(prior_records)} prior records from "
            f"{result_logger.json_path}; prior best_loss={best_loss:.6f} "
            f"params={format_params(best_params)}"
        )
    if resume and start_config > 1 and not prior_records:
        raise SystemExit(
            f"{stage_name}: --start-config/--osdr-start-config was set to "
            f"{start_config}, but no prior records for this stage were found in "
            f"{result_logger.json_path}. Copy the prior hpt_config_results.json "
            "into --output-dir or start from config 1."
        )

    stage_start = time.perf_counter()
    skipped_before_start = start_config - 1
    skipped_completed = len(
        [idx for idx in completed_indices if idx >= start_config]
    ) if resume else 0
    log(
        f"{stage_name}: starting grid search with {total_configs} configs, "
        f"data_shape={data.shape}, device={device}, start_config={start_config}, "
        f"resume={resume}, skipped_before_start={skipped_before_start}, "
        f"completed_to_skip={skipped_completed}"
    )

    for config_idx, params in enumerate(param_combinations, start=1):
        if config_idx < start_config:
            continue
        if resume and config_idx in completed_indices:
            continue

        param_dict = dict(zip(param_keys, params))
        config_start = time.perf_counter()
        log(
            f"{stage_name}: config {config_idx}/{total_configs} start "
            f"{format_params(param_dict)}"
        )

        # Prepare Data
        scaler = StandardScaler()
        X = scaler.fit_transform(data)
        X = torch.tensor(X, dtype=torch.float32)
        # For fine-tuning step when the input dimension needs to be changed
        if pi_dim is not None:
            adapter = Adapter(X.shape[1], pi_dim)
            X = adapter(X).clone().detach()
        # Set batch size and make the data to tensor # batch size can change depending on your device
        data_loader = DataLoader(
            X,
            batch_size=16,
            shuffle=True,
            num_workers=num_workers,
            pin_memory=device.type == "cuda",
        )

        # Model Initialization
        input_dim = X.shape[1]
        model = BuildSAE(input_dim, param_dict["hidden_layers"], param_dict["layer_norm"], param_dict["activation_fn"]).to(device)
        if pretrained_model is not None:
            model.load_state_dict(pretrained_model.state_dict())
        optimizer = optim.Adam(model.parameters(), lr=param_dict["learning_rate"], weight_decay=param_dict["weight_decay"])
        criterion = nn.MSELoss()

        # Training Loop with Early Stopping
        best_epoch_loss = float("inf")
        patience = 5
        patience_counter = 0
        early_stopped = False
        epochs_run = 0

        for epoch in range(50):  # Maximum epochs set for 50
            epochs_run = epoch + 1
            epoch_start = time.perf_counter()
            total_loss = 0.0
            model.train()   
            for batch_data in data_loader:
                optimizer.zero_grad()
                outputs = model(batch_data.to(device))
                loss = criterion(outputs, batch_data.to(device))
                # Get the activations (outputs) from the bottleneck layer
                encoded = model.encoder[-1](batch_data.to(device))
                # Add sparsity-inducing regularizer (L1 regularization) to the loss
                l1_regularization = torch.mean(torch.abs(encoded))
                loss += param_dict["sparsity_penalty"] * l1_regularization
                loss.backward()
                # Apply gradient clipping here
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
                optimizer.step()
                total_loss += loss.item()

            avg_loss = total_loss / len(data_loader)

            if avg_loss < best_epoch_loss:
                best_epoch_loss = avg_loss
                patience_counter = 0
            else:
                patience_counter += 1

            if log_every_epochs and (
                epoch == 0 or (epoch + 1) % log_every_epochs == 0
            ):
                log(
                    f"{stage_name}: config {config_idx}/{total_configs} "
                    f"epoch {epoch + 1}/50 loss={avg_loss:.6f} "
                    f"best={best_epoch_loss:.6f} patience={patience_counter}/5 "
                    f"epoch_time={format_elapsed(time.perf_counter() - epoch_start)}"
                )

            if patience_counter >= patience:
                early_stopped = True
                log(
                    f"{stage_name}: config {config_idx}/{total_configs} "
                    f"early stopped at epoch {epoch + 1}"
                )
                break

        elapsed_seconds = time.perf_counter() - config_start
        is_stage_best = best_epoch_loss < best_loss
        if result_logger is not None:
            result_logger.write(
                {
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "stage": stage_name,
                    "config_index": config_idx,
                    "total_configs": total_configs,
                    "data_shape": list(data.shape),
                    "device": str(device),
                    "params": format_params(param_dict),
                    "best_loss": float(best_epoch_loss),
                    "epochs_run": epochs_run,
                    "elapsed_seconds": round(elapsed_seconds, 3),
                    "elapsed": format_elapsed(elapsed_seconds),
                    "early_stopped": early_stopped,
                    "is_stage_best": is_stage_best,
                }
            )

        log(
            f"{stage_name}: config {config_idx}/{total_configs} complete "
            f"best_loss={best_epoch_loss:.6f} "
            f"epochs_run={epochs_run} early_stopped={early_stopped} "
            f"elapsed={format_elapsed(elapsed_seconds)}"
        )

        # Update Best Parameters
        if is_stage_best:
            best_loss = best_epoch_loss
            best_params = param_dict
            log(
                f"{stage_name}: new best config {config_idx}/{total_configs} "
                f"loss={best_loss:.6f} params={format_params(best_params)}"
            )

    log(
        f"{stage_name}: complete in {format_elapsed(time.perf_counter() - stage_start)} "
        f"best_loss={best_loss:.6f} best_params={format_params(best_params)}"
    )
    if best_params is None:
        raise SystemExit(f"{stage_name}: no configs were run and no prior best was loaded")
    return best_params

# Pre-Training and Fine-Tuning Workflow
def pretrain_sae(data, device, best_params, log_every_epochs=1, num_workers=0):
    # Extract best parameters
    hidden_layers = best_params["hidden_layers"]
    layer_norm = best_params["layer_norm"]
    activation_fn = best_params["activation_fn"]
    learning_rate = best_params["learning_rate"]
    weight_decay = best_params["weight_decay"]
    sparsity_penalty = best_params["sparsity_penalty"]
    run_start = time.perf_counter()
    log(
        f"pretrain: starting data_shape={data.shape}, device={device}, "
        f"params={format_params(best_params)}"
    )

    scaler = StandardScaler()
    X = scaler.fit_transform(data)
    X = torch.tensor(X, dtype=torch.float32)
    # Set batch size and make the data to tensor # batch size can change depending on your device
    data_loader = DataLoader(
        X,
        batch_size=16,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
    )

    # Model Initialization
    input_dim = X.shape[1]
    model = BuildSAE(input_dim, hidden_layers, layer_norm, activation_fn).to(device)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.MSELoss()

    best_epoch_loss = float("inf")
    patience = 5
    patience_counter = 0
    early_stopped = False
    epochs_run = 0

    for epoch in range(50):  # Maximum epochs set for 50
        epochs_run = epoch + 1
        epoch_start = time.perf_counter()
        total_loss = 0.0
        model.train()    
        for batch_data in data_loader:
            optimizer.zero_grad()
            outputs = model(batch_data.to(device))
            loss = criterion(outputs, batch_data.to(device))
            # Get the activations (outputs) from the bottleneck layer
            encoded = model.encoder[-1](batch_data.to(device))
            # Add sparsity-inducing regularizer (L1 regularization) to the loss
            l1_regularization = torch.mean(torch.abs(encoded))
            loss += sparsity_penalty * l1_regularization
            loss.backward()
            # Apply gradient clipping here
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
            optimizer.step()
            total_loss += loss.item()

        avg_loss = total_loss / len(data_loader)

        if avg_loss < best_epoch_loss:
            best_epoch_loss = avg_loss
            patience_counter = 0
        else:
            patience_counter += 1

        if log_every_epochs and (
            epoch == 0 or (epoch + 1) % log_every_epochs == 0
        ):
            log(
                f"pretrain: epoch {epoch + 1}/50 loss={avg_loss:.6f} "
                f"best={best_epoch_loss:.6f} patience={patience_counter}/5 "
                f"epoch_time={format_elapsed(time.perf_counter() - epoch_start)}"
            )

        if patience_counter >= patience:
            early_stopped = True
            log(f"pretrain: early stopped at epoch {epoch + 1}")
            break

    elapsed_seconds = time.perf_counter() - run_start
    log(
        f"pretrain: complete best_loss={best_epoch_loss:.6f} "
        f"elapsed={format_elapsed(elapsed_seconds)}"
    )
    return model, {
        "best_loss": float(best_epoch_loss),
        "epochs_run": epochs_run,
        "early_stopped": early_stopped,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "elapsed": format_elapsed(elapsed_seconds),
    }

def finetune_sae(
    data,
    pretrained_model,
    device,
    best_params,
    pi_dim,
    log_every_epochs=1,
    num_workers=0,
):
    # Extract best parameters
    hidden_layers = best_params["hidden_layers"]
    layer_norm = best_params["layer_norm"]
    activation_fn = best_params["activation_fn"]
    learning_rate = best_params["learning_rate"]
    weight_decay = best_params["weight_decay"]
    sparsity_penalty = best_params["sparsity_penalty"]
    run_start = time.perf_counter()
    log(
        f"finetune: starting data_shape={data.shape}, adapter_output_dim={pi_dim}, "
        f"device={device}, params={format_params(best_params)}"
    )

    scaler = StandardScaler()
    X = scaler.fit_transform(data)
    X = torch.tensor(X, dtype=torch.float32)

    adapter = Adapter(X.shape[1], pi_dim)
    X = adapter(X).clone().detach()
    # Set batch size and make the data to tensor # batch size can change depending on your device
    data_loader = DataLoader(
        X,
        batch_size=16,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
    )

    # Model Initialization
    input_dim = X.shape[1]
    model = BuildSAE(input_dim, hidden_layers, layer_norm, activation_fn).to(device)
    model.load_state_dict(pretrained_model.state_dict())  # Load pretrained weights
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.MSELoss()

    best_epoch_loss = float("inf")
    patience = 5
    patience_counter = 0
    early_stopped = False
    epochs_run = 0

    for epoch in range(50):  # Maximum epochs set for 50
        epochs_run = epoch + 1
        epoch_start = time.perf_counter()
        model.train()
        total_loss = 0.0
        for batch_data in data_loader:
            optimizer.zero_grad()
            outputs = model(batch_data.to(device))
            loss = criterion(outputs, batch_data.to(device))
            # Get the activations (outputs) from the bottleneck layer
            encoded = model.encoder[-1](batch_data.to(device))
            # Add sparsity-inducing regularizer (L1 regularization) to the loss
            l1_regularization = torch.mean(torch.abs(encoded))
            loss += sparsity_penalty * l1_regularization
            loss.backward()
            # Apply gradient clipping here
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
            optimizer.step()
            total_loss += loss.item()

        avg_loss = total_loss / len(data_loader)

        if avg_loss < best_epoch_loss:
            best_epoch_loss = avg_loss
            patience_counter = 0
        else:
            patience_counter += 1

        if log_every_epochs and (
            epoch == 0 or (epoch + 1) % log_every_epochs == 0
        ):
            log(
                f"finetune: epoch {epoch + 1}/50 loss={avg_loss:.6f} "
                f"best={best_epoch_loss:.6f} patience={patience_counter}/5 "
                f"epoch_time={format_elapsed(time.perf_counter() - epoch_start)}"
            )

        if patience_counter >= patience:
            early_stopped = True
            log(f"finetune: early stopped at epoch {epoch + 1}")
            break

    elapsed_seconds = time.perf_counter() - run_start
    log(
        f"finetune: complete best_loss={best_epoch_loss:.6f} "
        f"elapsed={format_elapsed(elapsed_seconds)}"
    )
    return model, adapter, {
        "best_loss": float(best_epoch_loss),
        "epochs_run": epochs_run,
        "early_stopped": early_stopped,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "elapsed": format_elapsed(elapsed_seconds),
    }

def get_representation(SAE_model, data, device, adapter=None):
    run_start = time.perf_counter()
    log(f"representation: starting data_shape={data.shape}, device={device}")
    # Train SAE
    SAE_model.eval()
    # Create a StandardScaler instance
    scaler = StandardScaler()
    # Fit the scaler to your data and transform it
    X = scaler.fit_transform(np.array(data))
    # To tensor
    X = torch.tensor(X, dtype=torch.float32)
    if adapter is not None:
        adapter.eval()
        with torch.no_grad():
            X = adapter(X).clone().detach()
    X = X.to(device)
    # Retrieve data representation from bottleneck layer
    with torch.no_grad():
        encoded_data = SAE_model.encoder(X)
    # Convert the encoded data tensor to NumPy array
    SAE_representation = encoded_data.detach().cpu().numpy()

    log(
        f"representation: complete shape={SAE_representation.shape} "
        f"elapsed={format_elapsed(time.perf_counter() - run_start)}"
    )
    return SAE_representation

def parse_arguments():
    # Command-Line Argument Parsing
    parser = argparse.ArgumentParser(description="Hyperparameter Tuning for Sparse AutoEncoder")
    parser.add_argument("--data1", required=True, type=str, help="Path to the first dataset (e.g., GeneLab)")
    parser.add_argument("--data2", required=True, type=str, help="Path to the second dataset (e.g., Single-Cell)")
    parser.add_argument(
        "--log-every-epochs",
        type=int,
        default=1,
        help="Print epoch progress every N epochs. Use 0 to disable epoch logs.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="",
        help="Directory for HPT outputs. Defaults to outputs/glare_hpt_<timestamp>.",
    )
    parser.add_argument(
        "--results-prefix",
        type=str,
        default="hpt_config_results",
        help="Filename prefix for per-config CSV/JSON logs.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Append to existing result logs in --output-dir and restore prior "
            "best configs from the JSON log."
        ),
    )
    parser.add_argument(
        "--start-config",
        type=int,
        default=None,
        help=(
            "First single-cell HPT config to run, 1-indexed. Alias for "
            "--single-cell-start-config; use 250 after configs 1-249 completed."
        ),
    )
    parser.add_argument(
        "--single-cell-start-config",
        type=int,
        default=None,
        help="First single-cell HPT config to run, 1-indexed.",
    )
    parser.add_argument(
        "--osdr-start-config",
        type=int,
        default=1,
        help="First OSDR HPT config to run, 1-indexed.",
    )
    parser.add_argument(
        "--reuse-best-configs-from",
        type=str,
        default="",
        help=(
            "Load single_cell_best_params and osdr_best_params from a prior "
            "hpt_summary.json and skip both hyperparameter sweeps."
        ),
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help=(
            "PyTorch DataLoader worker processes. Use 0 on macOS or restricted "
            "environments that do not permit torch shared-memory workers."
        ),
    )
    
    return parser.parse_args()

def matrixmarket_to_dense(matrix):
    if hasattr(matrix, "toarray"):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)

# Usage Example
if __name__ == "__main__":
    pipeline_start = time.perf_counter()
    # Set seed
    torch.manual_seed(1996)
    # Set gpu # cuda required, if you are using different gpu change.
    gpu = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # Get arguments
    args = parse_arguments()
    single_cell_start_config = (
        args.single_cell_start_config
        if args.single_cell_start_config is not None
        else args.start_config
    )
    if single_cell_start_config is None:
        single_cell_start_config = 1
    resume_run = (
        args.resume
        or single_cell_start_config > 1
        or args.osdr_start_config > 1
    )
    run_id = time.strftime("glare_hpt_%Y%m%d_%H%M%S")
    output_dir = Path(args.output_dir) if args.output_dir else Path("outputs") / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    fixed_config_source = Path(args.reuse_best_configs_from) if args.reuse_best_configs_from else None
    if fixed_config_source:
        fixed_config_summary = json.loads(
            fixed_config_source.read_text(encoding="utf-8")
        )
        best_params_sc = restore_params(
            fixed_config_summary["single_cell_best_params"]
        )
        best_params_gl = restore_params(
            fixed_config_summary["osdr_best_params"]
        )
        result_logger = None
        config_csv_output = None
        config_json_output = None
    else:
        result_logger = ConfigResultLogger(
            output_dir,
            prefix=args.results_prefix,
            resume=resume_run,
        )
        config_csv_output, config_json_output = result_logger.paths()
    log(f"Starting GLARE HPT pipeline on device={gpu}")
    log(f"GeneLab/OSDR input: {args.data1}")
    log(f"single-cell input: {args.data2}")
    log(f"Output directory: {output_dir}")
    log(
        f"Resume: {resume_run}; single_cell_start_config={single_cell_start_config}; "
        f"osdr_start_config={args.osdr_start_config}"
    )
    if fixed_config_source:
        log(f"Fixed-config mode: reusing best configs from {fixed_config_source}")
        log(f"Single-cell config: {format_params(best_params_sc)}")
        log(f"OSDR config: {format_params(best_params_gl)}")
    else:
        log(f"Config result CSV: {config_csv_output}")
        log(f"Config result JSON: {config_json_output}")
    # Load data 
    # # RECOMMENDED:Process the csv to include appropriate matrix for model training
    load_start = time.perf_counter()
    gl_csv = pd.read_csv(args.data1)
    log(
        f"Loaded GeneLab/OSDR CSV shape={gl_csv.shape} "
        f"elapsed={format_elapsed(time.perf_counter() - load_start)}"
    )
    # # RECOMMENDED: Convert your single-cell data to .mtx data so it could be read appropriately
    load_start = time.perf_counter()
    sc_matrix = mmread(args.data2)
    sc_dense = matrixmarket_to_dense(sc_matrix)
    log(
        f"Loaded single-cell MatrixMarket shape={sc_dense.shape} "
        f"elapsed={format_elapsed(time.perf_counter() - load_start)}"
    )
    # # If sparse matrix change to sparse tensor
    # sparse_tensor = torch.sparse.FloatTensor(
    #     torch.LongTensor([sc_matrix.row, sc_matrix.col]),
    #     torch.FloatTensor(sc_matrix.data),
    #     torch.Size(sc_matrix.shape)
    # )

    # Hyperparameter Tuning and Pretraining for chosen single-cell Dataset
    if fixed_config_source:
        log("Skipping single-cell hyperparameter sweep")
    else:
        log("Tuning hyperparameters for single-cell data")
        best_params_sc = grid_search(
            sc_dense,
            gpu,
            stage_name="single-cell",
            log_every_epochs=args.log_every_epochs,
            result_logger=result_logger,
            resume=resume_run,
            start_config=single_cell_start_config,
            num_workers=args.num_workers,
        )
    log(f"Best Hyperparameters for single-cell data: {format_params(best_params_sc)}")

    log("Pretraining SAE with single-cell data")
    pretrained_model, pretrain_metrics = pretrain_sae(
        sc_dense,
        gpu,
        best_params_sc,
        log_every_epochs=args.log_every_epochs,
        num_workers=args.num_workers,
    )
    # Save pre-trained weights
    pretrained_output = output_dir / "pretrained_sae_sc.pth"
    torch.save(pretrained_model.state_dict(), pretrained_output) # rename it with your project/data name
    log(f"Saved pretrained weights: {pretrained_output}")

    # Load pretrained weights and tune hyperparameters for GeneLab Dataset
    pretrained_input_dim = sc_matrix.shape[1] 
    if fixed_config_source:
        log("Skipping OSDR hyperparameter sweep")
    else:
        log("Loading pretrained weights and tuning hyperparameters for GeneLab")
        best_params_gl = grid_search(
            gl_csv,
            gpu,
            pretrained_model,
            pi_dim=pretrained_input_dim,
            fixed_architecture=best_params_sc,
            stage_name="OSDR",
            log_every_epochs=args.log_every_epochs,
            result_logger=result_logger,
            resume=resume_run,
            start_config=args.osdr_start_config,
            num_workers=args.num_workers,
        )
    log(f"Best Hyperparameters for GeneLab data: {format_params(best_params_gl)}")

    log("Fine-tuning SAE with GeneLab data")
    finetuned_model, finetune_adapter, finetune_metrics = finetune_sae(
        gl_csv,
        pretrained_model,
        gpu,
        best_params_gl,
        pi_dim=pretrained_input_dim,
        log_every_epochs=args.log_every_epochs,
        num_workers=args.num_workers,
    )
    # Save Fine-tuned weights
    finetuned_output = output_dir / "finetuned_sae_gl.pth"
    torch.save(finetuned_model.state_dict(), finetuned_output) # rename it with your project/data name
    log(f"Saved fine-tuned weights: {finetuned_output}")
    adapter_output = output_dir / "finetune_adapter.pth"
    torch.save(finetune_adapter.state_dict(), adapter_output)
    log(f"Saved fine-tuning adapter: {adapter_output}")
    # Get Final data representation
    FTSAE_representation = get_representation(
        finetuned_model,
        gl_csv,
        gpu,
        adapter=finetune_adapter,
    )
    representation_output = output_dir / "FTSAE_representation.npy"
    np.save(representation_output, FTSAE_representation)
    log(f"Saved final representation: {representation_output}")
    total_elapsed_seconds = time.perf_counter() - pipeline_start
    summary_output = output_dir / "hpt_summary.json"
    summary = {
        "data1": args.data1,
        "data2": args.data2,
        "device": str(gpu),
        "num_workers": args.num_workers,
        "mode": "fixed_best_configs" if fixed_config_source else "hyperparameter_sweep",
        "best_config_source": str(fixed_config_source) if fixed_config_source else "",
        "resume": resume_run,
        "single_cell_start_config": single_cell_start_config,
        "osdr_start_config": args.osdr_start_config,
        "single_cell_best_params": format_params(best_params_sc),
        "osdr_best_params": format_params(best_params_gl),
        "pretrained_weights": str(pretrained_output),
        "finetuned_weights": str(finetuned_output),
        "finetune_adapter": str(adapter_output),
        "representation": str(representation_output),
        "config_results_csv": str(config_csv_output) if config_csv_output else "",
        "config_results_json": str(config_json_output) if config_json_output else "",
        "training_metrics": {
            "pretrain": pretrain_metrics,
            "finetune": finetune_metrics,
        },
        "total_elapsed_seconds": round(total_elapsed_seconds, 3),
        "total_elapsed": format_elapsed(total_elapsed_seconds),
    }
    summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    log(f"Saved HPT summary: {summary_output}")
    log(f"GLARE HPT pipeline complete in {format_elapsed(total_elapsed_seconds)}")
     
    # # The gl_csv might have different data structure to CARA data which is `tsne4viz` is based on.
    # # Edit the function accordingly on utils.py
    # viz_df = tsne4viz(gl_csv, FTSAE_representation, dim=2)
    # # Clustering
    # # The gl_csv might have different data structure to CARA data which is `GLARECluseter` is based on.
    # # Edit the class accordingly on clustering.py
    # glare_cst = GLARECluster(viz_df, FTSAE_representation)
    # # Run ensemble clustering
    # labeled_df = glare_cst.eac()
    # # Export
    # labeled_df.to_csv('final_clustering_df.csv', index=False)
