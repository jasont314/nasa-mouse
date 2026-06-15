# Utils
import pandas as pd
import numpy as np
import argparse
import itertools
import time

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


def format_params(params):
    if params is None:
        return None
    formatted = dict(params)
    formatted["activation_fn"] = activation_name(formatted["activation_fn"])
    return formatted

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

    best_params = None
    best_loss = float("inf")

    param_combinations = list(itertools.product(*param_grid.values()))
    param_keys = list(param_grid.keys())
    total_configs = len(param_combinations)
    stage_start = time.perf_counter()
    log(
        f"{stage_name}: starting grid search with {total_configs} configs, "
        f"data_shape={data.shape}, device={device}"
    )

    for config_idx, params in enumerate(param_combinations, start=1):
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
        data_loader = DataLoader(X, batch_size=16, shuffle=True, num_workers=4, pin_memory=True)

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

        for epoch in range(50):  # Maximum epochs set for 50
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
                log(
                    f"{stage_name}: config {config_idx}/{total_configs} "
                    f"early stopped at epoch {epoch + 1}"
                )
                break

        log(
            f"{stage_name}: config {config_idx}/{total_configs} complete "
            f"best_loss={best_epoch_loss:.6f} "
            f"elapsed={format_elapsed(time.perf_counter() - config_start)}"
        )

        # Update Best Parameters
        if best_epoch_loss < best_loss:
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
    return best_params

# Pre-Training and Fine-Tuning Workflow
def pretrain_sae(data, device, best_params, log_every_epochs=1):
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
    data_loader = DataLoader(X, batch_size=16, shuffle=True, num_workers=4, pin_memory=True)

    # Model Initialization
    input_dim = X.shape[1]
    model = BuildSAE(input_dim, hidden_layers, layer_norm, activation_fn).to(device)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.MSELoss()

    best_epoch_loss = float("inf")
    patience = 5
    patience_counter = 0

    for epoch in range(50):  # Maximum epochs set for 50
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
            log(f"pretrain: early stopped at epoch {epoch + 1}")
            break

    log(
        f"pretrain: complete best_loss={best_epoch_loss:.6f} "
        f"elapsed={format_elapsed(time.perf_counter() - run_start)}"
    )
    return model

def finetune_sae(data, pretrained_model, device, best_params, pi_dim, log_every_epochs=1):
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
    data_loader = DataLoader(X, batch_size=16, shuffle=True, num_workers=4, pin_memory=True)

    # Model Initialization
    input_dim = X.shape[1]
    model = BuildSAE(input_dim, hidden_layers, layer_norm, activation_fn).to(device)
    model.load_state_dict(pretrained_model.state_dict())  # Load pretrained weights
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    criterion = nn.MSELoss()

    best_epoch_loss = float("inf")
    patience = 5
    patience_counter = 0

    for epoch in range(50):  # Maximum epochs set for 50
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
            log(f"finetune: early stopped at epoch {epoch + 1}")
            break

    log(
        f"finetune: complete best_loss={best_epoch_loss:.6f} "
        f"elapsed={format_elapsed(time.perf_counter() - run_start)}"
    )
    return model, adapter

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
    log(f"Starting GLARE HPT pipeline on device={gpu}")
    log(f"GeneLab/OSDR input: {args.data1}")
    log(f"single-cell input: {args.data2}")
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
    log("Tuning hyperparameters for single-cell data")
    best_params_sc = grid_search(
        sc_dense,
        gpu,
        stage_name="single_cell_hpt",
        log_every_epochs=args.log_every_epochs,
    )
    log(f"Best Hyperparameters for single-cell data: {format_params(best_params_sc)}")

    log("Pretraining SAE with single-cell data")
    pretrained_model = pretrain_sae(
        sc_dense,
        gpu,
        best_params_sc,
        log_every_epochs=args.log_every_epochs,
    )
    # Save pre-trained weights
    pretrained_output = "pretrained_sae_sc.pth"
    torch.save(pretrained_model.state_dict(), pretrained_output) # rename it with your project/data name
    log(f"Saved pretrained weights: {pretrained_output}")

    # Load pretrained weights and tune hyperparameters for GeneLab Dataset
    pretrained_input_dim = sc_matrix.shape[1] 
    log("Loading pretrained weights and tuning hyperparameters for GeneLab")
    best_params_gl = grid_search(
        gl_csv,
        gpu,
        pretrained_model,
        pi_dim=pretrained_input_dim,
        fixed_architecture=best_params_sc,
        stage_name="genelab_hpt",
        log_every_epochs=args.log_every_epochs,
    )
    log(f"Best Hyperparameters for GeneLab data: {format_params(best_params_gl)}")

    log("Fine-tuning SAE with GeneLab data")
    finetuned_model, finetune_adapter = finetune_sae(
        gl_csv,
        pretrained_model,
        gpu,
        best_params_gl,
        pi_dim=pretrained_input_dim,
        log_every_epochs=args.log_every_epochs,
    )
    # Save Fine-tuned weights
    finetuned_output = "finetuned_sae_gl.pth"
    torch.save(finetuned_model.state_dict(), finetuned_output) # rename it with your project/data name
    log(f"Saved fine-tuned weights: {finetuned_output}")
    # Get Final data representation
    FTSAE_representation = get_representation(
        finetuned_model,
        gl_csv,
        gpu,
        adapter=finetune_adapter,
    )
    representation_output = "FTSAE_representation.npy"
    np.save(representation_output, FTSAE_representation)
    log(f"Saved final representation: {representation_output}")
    log(f"GLARE HPT pipeline complete in {format_elapsed(time.perf_counter() - pipeline_start)}")
     
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
