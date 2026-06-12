# Utils
import pandas as pd
import numpy as np
import argparse
import itertools

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
def grid_search(data, device, pretrained_model=None, pi_dim=None, fixed_architecture=None):
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

    for params in param_combinations:
        param_dict = dict(zip(param_keys, params))

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

            if patience_counter >= patience:
                break

        # Update Best Parameters
        if best_epoch_loss < best_loss:
            best_loss = best_epoch_loss
            best_params = param_dict

    return best_params

# Pre-Training and Fine-Tuning Workflow
def pretrain_sae(data, device, best_params):
    # Extract best parameters
    hidden_layers = best_params["hidden_layers"]
    layer_norm = best_params["layer_norm"]
    activation_fn = best_params["activation_fn"]
    learning_rate = best_params["learning_rate"]
    weight_decay = best_params["weight_decay"]
    sparsity_penalty = best_params["sparsity_penalty"]

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

        if patience_counter >= patience:
            break

    return model

def finetune_sae(data, pretrained_model, device, best_params, pi_dim):
    # Extract best parameters
    hidden_layers = best_params["hidden_layers"]
    layer_norm = best_params["layer_norm"]
    activation_fn = best_params["activation_fn"]
    learning_rate = best_params["learning_rate"]
    weight_decay = best_params["weight_decay"]
    sparsity_penalty = best_params["sparsity_penalty"]

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

        if patience_counter >= patience:
            break

    return model, adapter

def get_representation(SAE_model, data, device, adapter=None):
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

    return SAE_representation

def parse_arguments():
    # Command-Line Argument Parsing
    parser = argparse.ArgumentParser(description="Hyperparameter Tuning for Sparse AutoEncoder")
    parser.add_argument("--data1", required=True, type=str, help="Path to the first dataset (e.g., GeneLab)")
    parser.add_argument("--data2", required=True, type=str, help="Path to the second dataset (e.g., Single-Cell)")
    
    return parser.parse_args()

def matrixmarket_to_dense(matrix):
    if hasattr(matrix, "toarray"):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)

# Usage Example
if __name__ == "__main__":
    # Set seed
    torch.manual_seed(1996)
    # Set gpu # cuda required, if you are using different gpu change.
    gpu = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # Get arguments
    args = parse_arguments()
    # Load data 
    # # RECOMMENDED:Process the csv to include appropriate matrix for model training
    gl_csv = pd.read_csv(args.data1)
    # # RECOMMENDED: Convert your single-cell data to .mtx data so it could be read appropriately
    sc_matrix = mmread(args.data2)
    sc_dense = matrixmarket_to_dense(sc_matrix)
    # # If sparse matrix change to sparse tensor
    # sparse_tensor = torch.sparse.FloatTensor(
    #     torch.LongTensor([sc_matrix.row, sc_matrix.col]),
    #     torch.FloatTensor(sc_matrix.data),
    #     torch.Size(sc_matrix.shape)
    # )

    # Hyperparameter Tuning and Pretraining for chosen single-cell Dataset
    print("Tuning hyperparameters for single-cell data")
    best_params_sc = grid_search(sc_dense, gpu)
    print("Best Hyperparameters for single-cell data:", best_params_sc)

    print("Pretraining SAE with single-cell data")
    pretrained_model = pretrain_sae(sc_dense, gpu, best_params_sc)
    # Save pre-trained weights
    torch.save(pretrained_model.state_dict(), "pretrained_sae_sc.pth") # rename it with your project/data name

    # Load pretrained weights and tune hyperparameters for GeneLab Dataset
    pretrained_input_dim = sc_matrix.shape[1] 
    print("Loading pretrained weights and tuning hyperparameters for GeneLab")
    best_params_gl = grid_search(
        gl_csv,
        gpu,
        pretrained_model,
        pi_dim=pretrained_input_dim,
        fixed_architecture=best_params_sc,
    )
    print("Best Hyperparameters for GeneLab data:", best_params_gl)

    print("Fine-tuning SAE with GeneLab data")
    finetuned_model, finetune_adapter = finetune_sae(
        gl_csv,
        pretrained_model,
        gpu,
        best_params_gl,
        pi_dim=pretrained_input_dim,
    )
    # Save Fine-tuned weights
    torch.save(finetuned_model.state_dict(), "finetuned_sae_gl.pth") # rename it with your project/data name
    # Get Final data representation
    FTSAE_representation = get_representation(
        finetuned_model,
        gl_csv,
        gpu,
        adapter=finetune_adapter,
    )
    np.save('FTSAE_representation.npy', FTSAE_representation)
     
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
