"""Conditional WGAN-GP model components."""

from __future__ import annotations

import math

from nasa_mouse_glare.io import require_import


def embedding_dim(cardinality: int) -> int:
    """Bounded categorical embedding size inspired by the paper's rule of thumb."""

    if cardinality <= 1:
        return 1
    return min(50, max(2, int(math.sqrt(cardinality)) + 1))


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
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.noise_dim = int(noise_dim)
        self.output_dim = int(output_dim)
        self.covariates = CovariateEmbeddings(categorical_cardinalities)
        layers = []
        in_dim = self.noise_dim + self.covariates.output_dim
        for hidden_dim in hidden_dims:
            layers.extend([nn.Linear(in_dim, hidden_dim), nn.ReLU()])
            in_dim = hidden_dim
        layers.append(nn.Linear(in_dim, self.output_dim))
        self.network = nn.Sequential(*layers)

    def forward(self, noise, categories):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        x = torch.cat([noise, self.covariates(categories)], dim=1)
        return self.network(x)


class Critic(require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt").Module):
    """Conditional MLP critic returning scalar scores and hidden features."""

    def __init__(
        self,
        *,
        input_dim: int,
        categorical_cardinalities: list[int],
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.input_dim = int(input_dim)
        self.covariates = CovariateEmbeddings(categorical_cardinalities)
        layers = []
        in_dim = self.input_dim + self.covariates.output_dim
        for hidden_dim in hidden_dims:
            layers.extend([nn.Linear(in_dim, hidden_dim), nn.ReLU()])
            in_dim = hidden_dim
        self.features = nn.Sequential(*layers)
        self.output = nn.Linear(in_dim, 1)
        self.feature_dim = int(in_dim)

    def forward(self, expression, categories, *, return_features: bool = False):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        x = torch.cat([expression, self.covariates(categories)], dim=1)
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
        hidden_dims: tuple[int, ...] = (256, 256),
    ):
        nn = require_import("torch.nn", "pip install -r requirements-nasa-mouse-glare.txt")
        super().__init__()
        self.expression_dim = int(expression_dim)
        self.noise_dim = int(noise_dim)
        self.generator = Generator(
            noise_dim=noise_dim,
            output_dim=expression_dim,
            categorical_cardinalities=categorical_cardinalities,
            hidden_dims=hidden_dims,
        )
        self.critic = Critic(
            input_dim=expression_dim,
            categorical_cardinalities=categorical_cardinalities,
            hidden_dims=hidden_dims,
        )

    def sample_noise(self, n: int, device):
        torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
        return torch.randn((int(n), self.noise_dim), device=device)
