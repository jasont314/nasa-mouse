"""Training utilities for conditional WGAN-GP."""

from __future__ import annotations

from dataclasses import dataclass
import time

from nasa_mouse_glare.io import require_import


@dataclass
class TrainConfig:
    epochs: int = 100
    batch_size: int = 128
    learning_rate: float = 1e-4
    critic_steps: int = 5
    gradient_penalty: float = 10.0
    augmentation_probability: float = 0.0
    augmentation_noise_scale: float = 0.5
    seed: int = 2020


class ExpressionDataset(require_import("torch.utils.data", "pip install -r requirements.txt").Dataset):
    """Torch dataset with expression vectors and categorical codes."""

    def __init__(self, expression, categories):
        torch = require_import("torch", "pip install -r requirements.txt")
        self.expression = torch.as_tensor(expression, dtype=torch.float32)
        self.categories = torch.as_tensor(categories, dtype=torch.long)

    def __len__(self):
        return int(self.expression.shape[0])

    def __getitem__(self, index):
        return self.expression[index], self.categories[index]


def make_loader(expression, categories, *, batch_size: int, seed: int, shuffle: bool = True):
    torch = require_import("torch", "pip install -r requirements.txt")
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    dataset = ExpressionDataset(expression, categories)
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=min(int(batch_size), len(dataset)),
        shuffle=shuffle,
        drop_last=False,
        generator=generator,
    )


def gradient_penalty(model, real, fake, categories, device):
    torch = require_import("torch", "pip install -r requirements.txt")
    alpha = torch.rand((real.shape[0], 1), device=device)
    interpolated = (alpha * real + (1.0 - alpha) * fake).requires_grad_(True)
    score = model.critic(interpolated, categories)
    gradients = torch.autograd.grad(
        outputs=score,
        inputs=interpolated,
        grad_outputs=torch.ones_like(score),
        create_graph=True,
        retain_graph=True,
        only_inputs=True,
    )[0]
    gradients = gradients.view(gradients.shape[0], -1)
    return ((gradients.norm(2, dim=1) - 1.0) ** 2).mean()


def augment_expression(expression, *, probability: float, noise_scale: float):
    """Apply the released WGAN's per-profile Gaussian augmentation."""

    torch = require_import("torch", "pip install -r requirements.txt")
    probability = float(probability)
    noise_scale = float(noise_scale)
    if not 0.0 <= probability <= 1.0:
        raise ValueError("augmentation_probability must be in [0, 1]")
    if noise_scale < 0.0:
        raise ValueError("augmentation_noise_scale cannot be negative")
    if probability == 0.0 or noise_scale == 0.0:
        return expression
    selected = (
        torch.rand((len(expression), 1), device=expression.device) < probability
    ).to(expression.dtype)
    return expression + selected * torch.randn_like(expression) * noise_scale


def _gradient_norm(parameters):
    torch = require_import("torch", "pip install -r requirements.txt")
    squared = torch.zeros((), dtype=torch.float32)
    for parameter in parameters:
        if parameter.grad is not None:
            squared = squared.to(parameter.grad.device)
            squared = squared + parameter.grad.detach().float().square().sum()
    return squared.sqrt()


def train_epoch(model, loader, *, config: TrainConfig, optim_g, optim_d, device):
    torch = require_import("torch", "pip install -r requirements.txt")
    model.train()
    totals = {
        "critic_loss": 0.0,
        "generator_loss": 0.0,
        "wasserstein_estimate": 0.0,
        "gradient_penalty": 0.0,
        "critic_gradient_norm": 0.0,
        "generator_gradient_norm": 0.0,
        "batches": 0,
    }
    for real, categories in loader:
        real = real.to(device)
        categories = categories.to(device)
        batch_size = real.shape[0]
        critic_loss = None
        gp_value = None
        wasserstein = None
        for _ in range(int(config.critic_steps)):
            optim_d.zero_grad(set_to_none=True)
            noise = model.sample_noise(batch_size, device)
            fake = model.generator(noise, categories).detach()
            augmented_real = augment_expression(
                real,
                probability=config.augmentation_probability,
                noise_scale=config.augmentation_noise_scale,
            )
            augmented_fake = augment_expression(
                fake,
                probability=config.augmentation_probability,
                noise_scale=config.augmentation_noise_scale,
            )
            real_score = model.critic(augmented_real, categories)
            fake_score = model.critic(augmented_fake, categories)
            gp_value = gradient_penalty(
                model,
                augmented_real,
                augmented_fake,
                categories,
                device,
            )
            critic_loss = (
                fake_score.mean()
                - real_score.mean()
                + float(config.gradient_penalty) * gp_value
            )
            critic_loss.backward()
            critic_gradient_norm = _gradient_norm(model.critic.parameters())
            optim_d.step()
            wasserstein = real_score.mean() - fake_score.mean()

        optim_g.zero_grad(set_to_none=True)
        noise = model.sample_noise(batch_size, device)
        fake = model.generator(noise, categories)
        augmented_fake = augment_expression(
            fake,
            probability=config.augmentation_probability,
            noise_scale=config.augmentation_noise_scale,
        )
        generator_loss = -model.critic(augmented_fake, categories).mean()
        generator_loss.backward()
        generator_gradient_norm = _gradient_norm(model.generator.parameters())
        optim_g.step()

        totals["critic_loss"] += float(critic_loss.detach().cpu())
        totals["generator_loss"] += float(generator_loss.detach().cpu())
        totals["wasserstein_estimate"] += float(wasserstein.detach().cpu())
        totals["gradient_penalty"] += float(gp_value.detach().cpu())
        totals["critic_gradient_norm"] += float(
            critic_gradient_norm.detach().cpu()
        )
        totals["generator_gradient_norm"] += float(
            generator_gradient_norm.detach().cpu()
        )
        totals["batches"] += 1

    batches = max(1, totals.pop("batches"))
    return {key: value / batches for key, value in totals.items()}


def train_model(model, expression, categories, *, config: TrainConfig, device):
    torch = require_import("torch", "pip install -r requirements.txt")
    torch.manual_seed(int(config.seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(config.seed))
    loader = make_loader(
        expression,
        categories,
        batch_size=config.batch_size,
        seed=config.seed,
        shuffle=True,
    )
    optim_g = torch.optim.Adam(
        model.generator.parameters(), lr=float(config.learning_rate), betas=(0.5, 0.9)
    )
    optim_d = torch.optim.Adam(
        model.critic.parameters(), lr=float(config.learning_rate), betas=(0.5, 0.9)
    )
    started = time.time()
    history = []
    for epoch in range(1, int(config.epochs) + 1):
        row = train_epoch(
            model,
            loader,
            config=config,
            optim_g=optim_g,
            optim_d=optim_d,
            device=device,
        )
        row["epoch"] = epoch
        history.append(row)
    return {
        "epochs_completed": int(config.epochs),
        "training_seconds": float(time.time() - started),
        "history": history,
    }


def critic_features(model, expression, categories, *, batch_size: int, device):
    np = require_import("numpy", "pip install -r requirements.txt")
    torch = require_import("torch", "pip install -r requirements.txt")
    loader = make_loader(
        expression,
        categories,
        batch_size=batch_size,
        seed=0,
        shuffle=False,
    )
    model.eval()
    scores = []
    features = []
    with torch.no_grad():
        for real, cats in loader:
            real = real.to(device)
            cats = cats.to(device)
            score, feature = model.critic(real, cats, return_features=True)
            scores.append(score.detach().cpu().numpy())
            features.append(feature.detach().cpu().numpy())
    return np.concatenate(scores), np.concatenate(features, axis=0)


def generate_samples(model, categories, *, batch_size: int, device):
    np = require_import("numpy", "pip install -r requirements.txt")
    torch = require_import("torch", "pip install -r requirements.txt")
    categories = torch.as_tensor(categories, dtype=torch.long)
    dataset = torch.utils.data.TensorDataset(categories)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    model.eval()
    generated = []
    with torch.no_grad():
        for (cats,) in loader:
            cats = cats.to(device)
            noise = model.sample_noise(cats.shape[0], device)
            fake = model.generator(noise, cats)
            generated.append(fake.detach().cpu().numpy())
    return np.concatenate(generated, axis=0).astype("float32")


def generation_quality(real, fake):
    np = require_import("numpy", "pip install -r requirements.txt")

    def corr(a, b):
        if np.std(a) == 0 or np.std(b) == 0:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    real_mean = real.mean(axis=0)
    fake_mean = fake.mean(axis=0)
    real_std = real.std(axis=0)
    fake_std = fake.std(axis=0)
    return {
        "n_real": int(real.shape[0]),
        "n_fake": int(fake.shape[0]),
        "genes": int(real.shape[1]),
        "gene_mean_correlation": corr(real_mean, fake_mean),
        "gene_std_correlation": corr(real_std, fake_std),
        "mean_rmse": float(np.sqrt(np.mean((real_mean - fake_mean) ** 2))),
        "std_rmse": float(np.sqrt(np.mean((real_std - fake_std) ** 2))),
        "real_global_mean": float(real.mean()),
        "fake_global_mean": float(fake.mean()),
        "real_global_std": float(real.std()),
        "fake_global_std": float(fake.std()),
    }
