"""Conditional WGAN-GP model components."""

from __future__ import annotations

import math

from nasa_mouse_glare.io import require_import


def embedding_dim(cardinality: int) -> int:
    """Return the exact Viñas source rule ``int(sqrt(vocab)) + 1``."""

    return int(math.sqrt(max(0, int(cardinality)))) + 1


class CovariateEmbeddings(require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt").Module):
    """Embedding block for categorical covariates."""

    def __init__(self, cardinalities: list[int]):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.cardinalities = [int(cardinality) for cardinality in cardinalities]
        self.embedding_dims = [embedding_dim(cardinality) for cardinality in cardinalities]
        self.embeddings = nn.ModuleList(
            [
                nn.Embedding(num_embeddings=max(1, cardinality), embedding_dim=dim)
                for cardinality, dim in zip(self.cardinalities, self.embedding_dims)
            ]
        )

    @property
    def output_dim(self) -> int:
        return int(sum(self.embedding_dims))

    def forward(self, categories):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        if categories.ndim == 1:
            categories = categories[:, None]
        pieces = []
        for idx, embedding in enumerate(self.embeddings):
            pieces.append(embedding(categories[:, idx].long()))
        if not pieces:
            return torch.empty((categories.shape[0], 0), device=categories.device)
        return torch.cat(pieces, dim=1)


class Generator(require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt").Module):
    """Conditional MLP generator."""

    def __init__(
        self,
        *,
        noise_dim: int,
        output_dim: int,
        categorical_cardinalities: list[int],
        numeric_dim: int = 0,
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.noise_dim = int(noise_dim)
        self.output_dim = int(output_dim)
        self.numeric_dim = int(numeric_dim)
        self.covariates = CovariateEmbeddings(categorical_cardinalities)
        layers = []
        in_dim = self.noise_dim + self.numeric_dim + self.covariates.output_dim
        for hidden_dim in hidden_dims:
            layers.extend([nn.Linear(in_dim, hidden_dim), nn.ReLU()])
            in_dim = hidden_dim
        layers.append(nn.Linear(in_dim, self.output_dim))
        self.network = nn.Sequential(*layers)

    def forward(self, noise, categories, numeric=None):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        pieces = [noise]
        if self.numeric_dim:
            if numeric is None:
                numeric = torch.zeros(
                    (len(noise), self.numeric_dim),
                    dtype=noise.dtype,
                    device=noise.device,
                )
            pieces.append(numeric)
        pieces.append(self.covariates(categories))
        x = torch.cat(pieces, dim=1)
        return self.network(x)


class Critic(require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt").Module):
    """Conditional MLP critic returning scalar scores and hidden features."""

    def __init__(
        self,
        *,
        input_dim: int,
        categorical_cardinalities: list[int],
        numeric_dim: int = 0,
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.input_dim = int(input_dim)
        self.numeric_dim = int(numeric_dim)
        self.covariates = CovariateEmbeddings(categorical_cardinalities)
        layers = []
        in_dim = self.input_dim + self.numeric_dim + self.covariates.output_dim
        for hidden_dim in hidden_dims:
            layers.extend([nn.Linear(in_dim, hidden_dim), nn.ReLU()])
            in_dim = hidden_dim
        self.features = nn.Sequential(*layers)
        self.output = nn.Linear(in_dim, 1)
        self.feature_dim = int(in_dim)

    def forward(
        self, expression, categories, numeric=None, *, return_features: bool = False
    ):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        pieces = [expression]
        if self.numeric_dim:
            if numeric is None:
                numeric = torch.zeros(
                    (len(expression), self.numeric_dim),
                    dtype=expression.dtype,
                    device=expression.device,
                )
            pieces.append(numeric)
        pieces.append(self.covariates(categories))
        x = torch.cat(pieces, dim=1)
        features = self.features(x)
        score = self.output(features).view(-1)
        if return_features:
            return score, features
        return score


class ConditionalWGANGP(require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt").Module):
    """Container for generator and critic."""

    def __init__(
        self,
        *,
        expression_dim: int,
        categorical_cardinalities: list[int],
        noise_dim: int = 128,
        numeric_dim: int = 0,
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.expression_dim = int(expression_dim)
        self.noise_dim = int(noise_dim)
        self.numeric_dim = int(numeric_dim)
        self.generator = Generator(
            noise_dim=noise_dim,
            output_dim=expression_dim,
            categorical_cardinalities=categorical_cardinalities,
            numeric_dim=numeric_dim,
            hidden_dims=hidden_dims,
        )
        self.critic = Critic(
            input_dim=expression_dim,
            categorical_cardinalities=categorical_cardinalities,
            numeric_dim=numeric_dim,
            hidden_dims=hidden_dims,
        )

    def sample_noise(self, n: int, device):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        return torch.randn((int(n), self.noise_dim), device=device)
