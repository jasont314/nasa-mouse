"""Conditional residual MLP denoiser for tabular DDPM/DDIM."""

from __future__ import annotations

import math

from nasa_mouse_glare.io import require_import

torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")


def sinusoidal_embedding(timesteps, dim: int):
    half = dim // 2
    scale = math.log(10000.0) / max(half - 1, 1)
    freqs = torch.exp(torch.arange(half, device=timesteps.device, dtype=torch.float32) * -scale)
    values = timesteps.float().unsqueeze(1) * freqs.unsqueeze(0)
    emb = torch.cat([torch.sin(values), torch.cos(values)], dim=1)
    if dim % 2:
        emb = torch.nn.functional.pad(emb, (0, 1))
    return emb


class ResidualBlock(nn.Module):
    def __init__(self, dim: int, context_dim: int, *, dropout: float):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim + context_dim)
        self.fc1 = nn.Linear(dim + context_dim, dim)
        self.norm2 = nn.LayerNorm(dim + context_dim)
        self.fc2 = nn.Linear(dim + context_dim, dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, context):
        y = torch.cat([x, context], dim=1)
        y = self.norm1(y)
        y = torch.nn.functional.silu(self.fc1(y))
        y = self.dropout(y)
        y = torch.cat([y, context], dim=1)
        y = self.norm2(y)
        y = torch.nn.functional.silu(self.fc2(y))
        y = self.dropout(y)
        return x + y


class ConditionalDiffusionMLP(nn.Module):
    def __init__(
        self,
        *,
        expression_dim: int,
        categorical_cardinalities: list[int],
        hidden_dim: int = 512,
        n_blocks: int = 2,
        dropout: float = 0.1,
        time_embedding_dim: int = 64,
        categorical_embedding_dim: int = 8,
        sinusoidal_time: bool = False,
        num_timesteps: int = 1000,
    ):
        super().__init__()
        self.expression_dim = int(expression_dim)
        self.hidden_dim = int(hidden_dim)
        self.num_timesteps = int(num_timesteps)
        self.sinusoidal_time = bool(sinusoidal_time)
        self.time_embedding_dim = int(time_embedding_dim if sinusoidal_time else 1)
        self.category_embeddings = nn.ModuleList(
            [
                nn.Embedding(max(1, int(cardinality)), int(categorical_embedding_dim))
                for cardinality in categorical_cardinalities
            ]
        )
        context_dim = self.time_embedding_dim + len(self.category_embeddings) * int(categorical_embedding_dim)
        self.input = nn.Linear(self.expression_dim, self.hidden_dim)
        self.blocks = nn.ModuleList(
            [ResidualBlock(self.hidden_dim, context_dim, dropout=dropout) for _ in range(int(n_blocks))]
        )
        self.norm = nn.LayerNorm(self.hidden_dim)
        self.output = nn.Linear(self.hidden_dim, self.expression_dim)

    def context(self, timesteps, categories):
        if self.sinusoidal_time:
            pieces = [sinusoidal_embedding(timesteps, self.time_embedding_dim)]
        else:
            pieces = [(timesteps.float() / max(self.num_timesteps - 1, 1)).unsqueeze(1)]
        if categories is not None and len(self.category_embeddings):
            for idx, embedding in enumerate(self.category_embeddings):
                pieces.append(embedding(categories[:, idx]))
        return torch.cat(pieces, dim=1)

    def features(self, x, timesteps, categories):
        context = self.context(timesteps, categories)
        h = torch.nn.functional.silu(self.input(x))
        for block in self.blocks:
            h = block(h, context)
        return torch.nn.functional.silu(self.norm(h))

    def forward(self, x, timesteps, categories=None):
        return self.output(self.features(x, timesteps, categories))
