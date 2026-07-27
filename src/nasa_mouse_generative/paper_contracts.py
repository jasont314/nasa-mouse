"""Pinned source identities and immutable paper-native hyperparameter contracts."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
from typing import Any


PAPER_SOURCES: dict[str, dict[str, Any]] = {
    "vinas_wgan_gp": {
        "url": "https://github.com/rvinas/adversarial-gene-expression.git",
        "commit": "94fa44dd1bd52d924efd3af0fcd8eeb18bd141a8",
        "files": {
            "gtex_tcga_gan.py": "beeb3bc4eb8de47d78ad263d356814de9e68a3eec8308025535ec922082dd7aa",
            "tf_utils.py": "950fb137c3e42e7e04b73080c68456869f7a4a253eea2da79e7be9715460a954",
            "rnaseqdb.py": "e5bed742f0cde9d2bfdb2ce76d28ae4f3acb878627608d63dda8f88c6fb0dcf1",
        },
    },
    "lacan_diffusion": {
        "url": "https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git",
        "commit": "cde890154698fcea96c924804aaff04af3351b48",
        "files": {
            "src/generation/ddim/models/diffusion_ddim.py": "7f06b79b89dc08efab5e2f0de319ae7420d97d51a5de11a4eb8ce6e264987f91",
            "src/generation/ddim/functions/losses.py": "e87b2597b73c522e0e45570555cad2316f681c95fff000ce815106cd4a681968",
            "src/generation/ddim/functions/denoising.py": "680b5eab75c51bcc0bffb63b6892b4860060eb50ad981732c7b2d32f9f6d418c",
        },
    },
    "genejepa": {
        "url": "https://github.com/BiostateAI/GeneJEPA.git",
        "commit": "a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265",
        "files": {
            "genejepa/models.py": "2538c2cede1eabbfd3e2e726254d9fbcabff620377b0d208dbabc939a817385d",
            "genejepa/tokenizer.py": "8708c2c6c20d8258badf5fb4aff7ee50689cce556d49075e4f765d19f5a0abb3",
            "genejepa/train.py": "f0225a013aefc138eb98c2200af4a84068330d574a7bf87e15a1eca853080175",
            "genejepa/configs.py": "e4c995d06028cdd46165a00bbc2738689a050a45f6a737737542f17bdf524d55",
        },
    },
    "mbatch": {
        "url": "https://github.com/MD-Anderson-Bioinformatics/MBatch.git",
        "commit": "93cddd2ba18ed8781b9865ba0259fafa057bcc17",
        "files": {
            "apps/MBatch/R/BEA_CorrectionsMP.R": "706084efe5c2bad7bad2ec24b7013c9105cfab7e1f865d72946c0cf3b411dc0f",
            "apps/MBatch/R/BEA_CorrectionsEB.R": "897c03e2a375e0e3602330099a5d2ca3f02a593ea20fb792bb928f0a8ad52544",
            "apps/MBatch/R/BEA_CorrectionsAN.R": "9f7afddf224d8fc5d46e74c5a55b2a98457a6c7aabf7a3e4328f6a66b2f8ff5b",
        },
    },
}


PAPER_NATIVE_LOCKS: dict[str, dict[str, Any]] = {
    "vinas_wgan_gp": {
        "noise_dim": 64,
        "numeric_dim": 1,
        "hidden_dims": [256, 256],
        "batch_size": 32,
        "critic_steps": 5,
        "gradient_penalty": 10.0,
        "optimizer": "rmsprop",
        "rmsprop_alpha": 0.9,
        "rmsprop_epsilon": 1e-7,
        "learning_rate": 5e-4,
        "epochs": 2000,
        "reference_epochs": 2000,
        "weighted_sampling": False,
        "early_stopping": True,
        "early_stopping_first_epoch": 1,
        "early_stopping_evaluate_every_epochs": 5,
    },
    "lacan_diffusion": {
        "batch_size": 2048,
        "learning_rate": 0.0004783833151836702,
        "epochs": 15000,
        "hidden_dim": 8192,
        "n_blocks": 2,
        "dropout": 0.1,
        "time_embedding_dim": 1,
        "categorical_embedding_dim": 2,
        "sinusoidal_time": False,
        "diffusion_timesteps": 1000,
        "sample_steps": 1000,
        "beta_schedule": "quad",
        "beta_start": 0.0001,
        "beta_end": 0.02,
        "optimizer": "adam",
        "weight_decay": 0.0,
        "scheduler": "one_cycle",
        "use_amp": True,
        "use_ema": True,
        "ema_decay": 0.999,
        "landmark_strategy": "l1000",
        "n_landmarks": 974,
        "gradient_clipping": False,
        "weighted_sampling": False,
        "loss_reduction": "sum_genes_mean_batch",
        "antithetic_timesteps": True,
    },
    "genejepa": {
        "batch_size": 92,
        "learning_rate": 1e-4,
        "weight_decay": 2e-4,
        "accumulate_grad_batches": 2,
        "warmup_ratio": 0.05,
        "scheduler": "cosine",
        "epochs": 50,
        "reference_epochs": 50,
        "samples_per_epoch": 1000000,
        "num_workers": 8,
        "weighted_sampling": False,
        "use_amp": True,
        "amp_dtype": "bfloat16",
        "d": 768,
        "latents_L": 512,
        "blocks_D": 24,
        "heads_h": 12,
        "mask_ratio": 0.45,
        "ema_start_decay": 0.992,
        "ema_end_decay": 0.9995,
        "ema_warmup_steps": 2000,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_pinned_source(model: str, source_root: str | Path) -> dict[str, Any]:
    """Verify both the Git identity and selected source files for one paper."""

    contract = PAPER_SOURCES[model]
    root = Path(source_root).resolve()
    if not (root / ".git").exists():
        raise FileNotFoundError(f"Pinned {model} source checkout is missing: {root}")
    observed_commit = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if observed_commit != contract["commit"]:
        raise RuntimeError(
            f"Expected {model} source {contract['commit']}, observed {observed_commit}"
        )
    hashes: dict[str, str] = {}
    for relative, expected in contract["files"].items():
        path = root / relative
        if not path.exists():
            raise FileNotFoundError(f"Pinned {model} source file is missing: {path}")
        observed = _sha256(path)
        if observed != expected:
            raise RuntimeError(
                f"Pinned {model} source file changed: {relative}; "
                f"expected {expected}, observed {observed}"
            )
        hashes[relative] = observed
    return {
        "model": model,
        "source_url": contract["url"],
        "source_root": str(root),
        "source_commit": observed_commit,
        "source_file_sha256": hashes,
    }


def validate_paper_native_parameters(
    model: str, parameters: dict[str, Any], profile: str = "paper_native"
) -> None:
    """Reject a profile labeled paper-native when a locked value was changed."""

    expected = dict(PAPER_NATIVE_LOCKS[model])
    if model == "vinas_wgan_gp":
        if profile == "paper_native":
            expected.update(
                {
                    "early_stopping_variant": "released_code",
                    "early_stopping_patience_checks": 10,
                }
            )
        elif profile == "paper_native_paper_text":
            expected.update(
                {
                    "early_stopping_variant": "paper_text",
                    "early_stopping_patience_epochs": 30,
                }
            )
    mismatches = {
        key: {"expected": value, "observed": parameters.get(key)}
        for key, value in expected.items()
        if parameters.get(key) != value
    }
    if mismatches:
        details = ", ".join(
            f"{key}={item['observed']!r} (expected {item['expected']!r})"
            for key, item in mismatches.items()
        )
        raise ValueError(
            f"{model} paper-native contract was modified: {details}. "
            "Use a practical/tunable profile for architecture or training changes."
        )
