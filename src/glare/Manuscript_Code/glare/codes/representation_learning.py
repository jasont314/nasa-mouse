# Utils
import pandas as pd
from .utils import concat_df
# Representation learning via Dimensionality Reduction
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
# Representation learning via SAE
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.preprocessing import StandardScaler

# Set environment
# os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:64'
torch.cuda.empty_cache()
torch.manual_seed(2023)


def get_pca(nc_df):
    # PCA, 2-dimension
    glds120_pca = PCA(n_components=2).fit_transform(nc_df)
    # Get full pca_df
    pca_df = concat_df(glds120_pca, nc_df, dimension=2)

    return pca_df


def get_tsne(nc_df):
    # t-SNE, 2-dimension
    glds120_tsne = TSNE(n_components=2, random_state=1996, n_jobs=-1,
                        learning_rate='auto').fit_transform(nc_df)
    # Get full tsne_df
    tsne_df = concat_df(glds120_tsne, nc_df, dimension=2)

    return tsne_df


def get_umap(nc_df):
    # UMAP, 2-dimension
    glds120_umap = UMAP(n_neighbors=5, min_dist=0.0, spread=2, n_components=2,
                        random_state=1996).fit_transform(nc_df)
    # Get full umap_df
    umap_df = concat_df(glds120_umap, nc_df, dimension=2)

    return umap_df


# Create sparse autoencoder model architecture
class SparseAutoEncoder(nn.Module):
    def __init__(self, input_dim):
        super(SparseAutoEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.LayerNorm(128),
            nn.ELU(),  # LeakyReLU(0.2)
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
            nn.Sigmoid()
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded


def train_SAE(X, device, exp_type):
    # Create a StandardScaler instance
    scaler = StandardScaler()
    # Fit the scaler to your data and transform it
    X = scaler.fit_transform(X)
    # Set batch size and make the data to tensor # batch size can change depending on your device
    batch_size = 16
    X = torch.tensor(X, dtype=torch.float32)
    data_loader = DataLoader(X, batch_size=batch_size, shuffle=True,
                             num_workers=4, pin_memory=True)

    # Initialize the sparse autoencoder model, loss function, and optimizer
    input_dim = X.shape[1]
    sparse_autoencoder = SparseAutoEncoder(input_dim).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(sparse_autoencoder.parameters(), lr=0.001, weight_decay=0.0001)
    # LR scheduler if needed
    # scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.5)

    # tuned early stopping #28 for FLT # 24 for GC # 34 for single-cell pretraining
    num_epochs = 30 # 28 if exp_type == 'FLT' else (24 if exp_type == 'GC' else 34)
    sparsity_penalty = 1e-5  # Adjust the sparsity penalty coefficient as needed
    for epoch in range(num_epochs):
        total_loss = 0.0
        sparse_autoencoder.train()
        for batch_data in data_loader:
            optimizer.zero_grad()
            outputs = sparse_autoencoder(batch_data.to(device))
            loss = criterion(outputs, batch_data.to(device))
            # Get the activations (outputs) from the bottleneck layer
            encoded = sparse_autoencoder.encoder[-1](batch_data.to(device))
            # Add sparsity-inducing regularizer (L1 regularization) to the loss
            l1_regularization = torch.mean(torch.abs(encoded))  # Applying L1 to the bottleneck layer
            loss += sparsity_penalty * l1_regularization
            loss.backward()
            # Apply gradient clipping here
            torch.nn.utils.clip_grad_norm_(sparse_autoencoder.parameters(), max_norm=1)
            optimizer.step()
            # scheduler.step()
            total_loss += loss.item()
        # average loss
        average_loss = total_loss / len(data_loader)
        print(f'Epoch [{epoch + 1}/{num_epochs}], Training Loss: {average_loss:.6f}')

    torch.save(sparse_autoencoder.state_dict(),
               './weights/sc_shulse_pretrained.pth') if exp_type == 'sc_pretrain' else None

    return sparse_autoencoder


# Define adapter layers to adjust the dimensions
class Adapter(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(Adapter, self).__init__()
        self.linear = nn.Linear(input_dim, output_dim)
        self.activation = nn.ELU()

    def forward(self, x):
        x = self.linear(x)
        x = self.activation(x)
        return x


def finetune_SAE_sc(X, pi_dim, weights, device, ft_epoch):
    # Create a StandardScaler instance
    scaler = StandardScaler()
    # Fit the scaler to your data and transform it
    X = scaler.fit_transform(X)
    # Set batch size and make the data to tensor # batch size can change depending on your device
    batch_size = 16
    X = torch.tensor(X, dtype=torch.float32)
    # Use adapter layer to set dimension
    adapter = Adapter(X.shape[1], pi_dim)
    X = adapter(X).clone().detach()
    # data loader
    data_loader = DataLoader(X, batch_size=batch_size, shuffle=True,
                             num_workers=4, pin_memory=True)

    # Initialize the sparse autoencoder model, loss function, and optimizer
    input_dim = X.shape[1]
    sparse_autoencoder = SparseAutoEncoder(input_dim).to(device)
    # Load weights
    sparse_autoencoder.load_state_dict(torch.load(weights))
    criterion = nn.MSELoss()
    optimizer = optim.Adam(sparse_autoencoder.parameters(), lr=0.001) # Turn off weight decay #, weight_decay=0.0001
    # LR scheduler if needed
    # scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.5)

    # tuned early stopping
    num_epochs = ft_epoch #12
    sparsity_penalty = 1e-5  # Adjust the sparsity penalty coefficient as needed
    for epoch in range(num_epochs):
        total_loss = 0.0
        sparse_autoencoder.train()
        for batch_data in data_loader:
            optimizer.zero_grad()
            outputs = sparse_autoencoder(batch_data.to(device))
            loss = criterion(outputs, batch_data.to(device))
            # Get the activations (outputs) from the bottleneck layer
            encoded = sparse_autoencoder.encoder[-1](batch_data.to(device))
            # Add sparsity-inducing regularizer (L1 regularization) to the loss
            l1_regularization = torch.mean(torch.abs(encoded))  # Applying L1 to the bottleneck layer
            loss += sparsity_penalty * l1_regularization
            loss.backward()
            # Apply gradient clipping here # Turn off gradient clipping for fine-tuning
            # torch.nn.utils.clip_grad_norm_(sparse_autoencoder.parameters(), max_norm=1)
            optimizer.step()
            # scheduler.step()
            total_loss += loss.item()
        # average loss
        average_loss = total_loss / len(data_loader)
        print(f'Epoch [{epoch + 1}/{num_epochs}], Training Loss: {average_loss:.7f}')

        # Get checkpoint for every 5 epochs # Take the best performing checkpoint weights based on the loss
        if (epoch + 1) % 5 == 0:
            torch.save({
                'epoch': epoch,
                'model_state_dict': sparse_autoencoder.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                # Add other information as needed
            }, f'./weights/finetune_checkpoint_epoch_{epoch + 1}.pth')

    return sparse_autoencoder
