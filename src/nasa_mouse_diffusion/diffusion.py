"""DDPM training loss and DDIM/DDPM sampling utilities."""

from __future__ import annotations

from nasa_mouse_glare.io import require_import


def beta_schedule(schedule: str, *, beta_start: float, beta_end: float, timesteps: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    if schedule in {"quad", "quadratic"}:
        betas = np.linspace(beta_start ** 0.5, beta_end ** 0.5, timesteps, dtype=np.float64) ** 2
    elif schedule == "linear":
        betas = np.linspace(beta_start, beta_end, timesteps, dtype=np.float64)
    elif schedule in {"const", "constant"}:
        betas = np.full(timesteps, beta_end, dtype=np.float64)
    elif schedule == "sigmoid":
        x = np.linspace(-6, 6, timesteps)
        betas = 1 / (np.exp(-x) + 1) * (beta_end - beta_start) + beta_start
    else:
        raise ValueError(f"Unknown beta schedule: {schedule}")
    return betas.astype("float32")


def compute_alpha(betas, timesteps):
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    padded = torch.cat([torch.zeros(1, device=betas.device), betas], dim=0)
    return (1 - padded).cumprod(dim=0).index_select(0, timesteps + 1).view(-1, 1)


def noise_estimation_loss(model, x0, categories, betas):
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    n = x0.shape[0]
    noise = torch.randn_like(x0)
    t = torch.randint(low=0, high=betas.shape[0], size=(n // 2 + 1,), device=x0.device)
    t = torch.cat([t, betas.shape[0] - t - 1], dim=0)[:n]
    alpha = (1 - betas).cumprod(dim=0).index_select(0, t).view(-1, 1)
    xt = x0 * alpha.sqrt() + noise * (1.0 - alpha).sqrt()
    predicted = model(xt, t, categories)
    return (noise - predicted).square().sum(dim=1).mean(), (noise - predicted).abs().mean()


def sample(model, categories, *, betas, sample_steps: int, eta: float, noise=None, device=None):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    if device is None:
        device = next(model.parameters()).device
    categories = torch.as_tensor(categories, dtype=torch.long, device=device)
    n = categories.shape[0]
    if noise is None:
        x = torch.randn(n, model.expression_dim, device=device)
    else:
        x = noise.to(device)
    timesteps = betas.shape[0]
    if sample_steps >= timesteps:
        seq = list(range(timesteps))
    else:
        seq = [int(value) for value in np.linspace(0, timesteps - 1, int(sample_steps))]
    seq_next = [-1] + seq[:-1]
    with torch.no_grad():
        for i, j in zip(reversed(seq), reversed(seq_next)):
            t = torch.full((n,), i, device=device, dtype=torch.long)
            next_t = torch.full((n,), j, device=device, dtype=torch.long)
            at = compute_alpha(betas, t)
            at_next = compute_alpha(betas, next_t)
            et = model(x, t, categories)
            x0 = (x - et * (1 - at).sqrt()) / at.sqrt()
            c1 = eta * ((1 - at / at_next) * (1 - at_next) / (1 - at)).clamp_min(0).sqrt()
            c2 = ((1 - at_next) - c1 ** 2).clamp_min(0).sqrt()
            x = at_next.sqrt() * x0 + c1 * torch.randn_like(x) + c2 * et
    return x
