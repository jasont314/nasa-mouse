"""Pinned access to the Lacan et al. DDIM architecture and equations."""

from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib.util
from pathlib import Path
import subprocess
from types import ModuleType, SimpleNamespace
from typing import Iterable

import numpy as np
import torch


SOURCE_URL = "https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git"
SOURCE_COMMIT = "cde890154698fcea96c924804aaff04af3351b48"
SOURCE_FILES = {
    "src/generation/ddim/models/diffusion_ddim.py": (
        "7f06b79b89dc08efab5e2f0de319ae7420d97d51a5de11a4eb8ce6e264987f91"
    ),
    "src/generation/ddim/functions/losses.py": (
        "e87b2597b73c522e0e45570555cad2316f681c95fff000ce815106cd4a681968"
    ),
    "src/generation/ddim/functions/denoising.py": (
        "680b5eab75c51bcc0bffb63b6892b4860060eb50ad981732c7b2d32f9f6d418c"
    ),
}


def repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def default_source_root() -> Path:
    return repository_root() / "assets/model_sources/rna-diffusion"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_source(source_root: str | Path | None = None) -> dict[str, object]:
    root = Path(source_root or default_source_root()).resolve()
    observed_commit = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if observed_commit != SOURCE_COMMIT:
        raise RuntimeError(
            f"Expected RNA diffusion source {SOURCE_COMMIT}, observed {observed_commit}"
        )
    hashes: dict[str, str] = {}
    for relative, expected in SOURCE_FILES.items():
        path = root / relative
        observed = _sha256(path)
        if observed != expected:
            raise RuntimeError(
                f"Pinned source file changed: {relative} expected {expected}, "
                f"observed {observed}"
            )
        hashes[relative] = observed
    return {
        "source_url": SOURCE_URL,
        "source_commit": observed_commit,
        "source_root": str(root),
        "source_file_sha256": hashes,
    }


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Python module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=1)
def upstream_model_class():
    root = default_source_root()
    verify_source(root)
    module = _load_module(
        "nasa_mouse_pinned_lacan_model",
        root / "src/generation/ddim/models/diffusion_ddim.py",
    )
    return module.ModelDDIM


@lru_cache(maxsize=1)
def upstream_loss_module() -> ModuleType:
    root = default_source_root()
    verify_source(root)
    return _load_module(
        "nasa_mouse_pinned_lacan_losses",
        root / "src/generation/ddim/functions/losses.py",
    )


@lru_cache(maxsize=1)
def upstream_denoising_module() -> ModuleType:
    root = default_source_root()
    verify_source(root)
    return _load_module(
        "nasa_mouse_pinned_lacan_denoising",
        root / "src/generation/ddim/functions/denoising.py",
    )


def model_config(*, expression_dim: int, num_classes: int, model: dict) -> SimpleNamespace:
    """Build the namespace consumed by the unmodified upstream ModelDDIM."""

    return SimpleNamespace(
        data=SimpleNamespace(
            image_size=int(expression_dim),
            num_classes=int(num_classes),
        ),
        diffusion=SimpleNamespace(
            beta_schedule=str(model["beta_schedule"]),
            beta_start=float(model["beta_start"]),
            beta_end=float(model["beta_end"]),
            num_diffusion_timesteps=int(model["diffusion_timesteps"]),
        ),
        model=SimpleNamespace(
            dropout=float(model["dropout"]),
            dim_t=int(model["time_embedding_dim"]),
            is_time_embed=bool(model["sinusoidal_time"]),
            is_y_cond=True,
            num_res_blocks=0,
            attn_resolutions=[],
            d_layers=[int(value) for value in model["hidden_dims"]],
            with_attn=False,
            use_y_emb=True,
            dim_y_emb=int(model["tissue_embedding_dim"]),
            precision="single",
            type="simple",
            var_type="fixedlarge",
            model="ddim",
            ema=True,
            ema_rate=float(model["ema_decay"]),
            parallel=False,
        ),
    )


def quadratic_beta_schedule(
    *, beta_start: float, beta_end: float, timesteps: int
) -> torch.Tensor:
    values = np.linspace(
        float(beta_start) ** 0.5,
        float(beta_end) ** 0.5,
        int(timesteps),
        dtype=np.float64,
    ) ** 2
    return torch.from_numpy(values).float()


def antithetic_timesteps(batch_size: int, timesteps: int, device: torch.device) -> torch.Tensor:
    first = torch.randint(
        low=0,
        high=int(timesteps),
        size=(int(batch_size) // 2 + 1,),
        device=device,
    )
    return torch.cat([first, int(timesteps) - first - 1], dim=0)[:batch_size]


def noise_estimation_loss(
    model: torch.nn.Module,
    clean: torch.Tensor,
    timesteps: torch.Tensor,
    noise: torch.Tensor,
    betas: torch.Tensor,
    labels: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    alpha = (1 - betas).cumprod(dim=0).index_select(0, timesteps).view(-1, 1)
    noisy = clean * alpha.sqrt() + noise * (1.0 - alpha).sqrt()
    prediction = model(noisy, timesteps, labels)
    absolute_error = (noise - prediction).abs().mean(dim=1).mean(dim=0)
    loss = (noise - prediction).square().sum(dim=1).mean(dim=0)
    return loss, absolute_error


def compute_alpha(betas: torch.Tensor, timesteps: torch.Tensor) -> torch.Tensor:
    padded = torch.cat([torch.zeros(1, device=betas.device), betas], dim=0)
    return (1 - padded).cumprod(dim=0).index_select(0, timesteps + 1).view(-1, 1)


def ddim_trajectory(
    initial_noise: torch.Tensor,
    labels: torch.Tensor,
    model: torch.nn.Module,
    betas: torch.Tensor,
    *,
    sequence: Iterable[int],
    snapshot_timesteps: Iterable[int] = (1000, 200, 0),
    eta: float = 0.0,
    generator: torch.Generator | None = None,
) -> dict[int, torch.Tensor]:
    """Run upstream-equivalent generalized steps while retaining requested states."""

    requested = {int(value) for value in snapshot_timesteps}
    seq = list(map(int, sequence))
    seq_next = [-1] + seq[:-1]
    total = len(seq)
    current = initial_noise
    snapshots: dict[int, torch.Tensor] = {}
    if total in requested:
        snapshots[total] = current.detach().cpu()
    with torch.no_grad():
        for current_t, next_t in zip(reversed(seq), reversed(seq_next)):
            t = torch.full(
                (len(current),), current_t, dtype=torch.long, device=current.device
            )
            next_values = torch.full(
                (len(current),), next_t, dtype=torch.long, device=current.device
            )
            alpha_t = compute_alpha(betas, t)
            alpha_next = compute_alpha(betas, next_values)
            with torch.autocast(
                current.device.type,
                dtype=torch.float16,
                enabled=current.device.type == "cuda",
            ):
                predicted_noise = model(current, t.float(), labels)
            predicted_clean = (
                current - predicted_noise * (1 - alpha_t).sqrt()
            ) / alpha_t.sqrt()
            c1 = float(eta) * (
                (1 - alpha_t / alpha_next)
                * (1 - alpha_next)
                / (1 - alpha_t)
            ).sqrt()
            c2 = ((1 - alpha_next) - c1**2).sqrt()
            current = (
                alpha_next.sqrt() * predicted_clean
                + c1
                * torch.randn(
                    initial_noise.shape,
                    dtype=initial_noise.dtype,
                    device=initial_noise.device,
                    generator=generator,
                )
                + c2 * predicted_noise
            )
            display_timestep = next_t + 1
            if display_timestep in requested:
                snapshots[display_timestep] = current.detach().cpu()
    missing = requested.difference(snapshots)
    if missing:
        raise ValueError(f"Trajectory did not visit requested timesteps: {sorted(missing)}")
    return snapshots


class EMA:
    """Parameter-only EMA matching the upstream EMAHelper implementation."""

    def __init__(self, model: torch.nn.Module, decay: float) -> None:
        self.decay = float(decay)
        self.shadow = {
            name: parameter.data.clone()
            for name, parameter in model.named_parameters()
            if parameter.requires_grad
        }

    def update(self, model: torch.nn.Module) -> None:
        for name, parameter in model.named_parameters():
            if parameter.requires_grad:
                self.shadow[name].data = (
                    (1.0 - self.decay) * parameter.data
                    + self.decay * self.shadow[name].data
                )

    def copy_to(self, model: torch.nn.Module) -> None:
        for name, parameter in model.named_parameters():
            if parameter.requires_grad:
                parameter.data.copy_(self.shadow[name].data)

    def state_dict(self) -> dict[str, torch.Tensor]:
        return self.shadow

    def load_state_dict(self, state: dict[str, torch.Tensor], device: torch.device) -> None:
        self.shadow = {name: value.to(device) for name, value in state.items()}
