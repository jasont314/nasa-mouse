"""Capabilities and provenance for the three requested model families."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelCapabilities:
    model_id: str
    display_name: str
    paper_url: str
    code_url: str
    inspected_commit: str
    native_input: str
    supports_expression_generation: bool
    supports_categorical_conditioning: bool
    supports_pretrain_finetune: bool
    supports_counterfactual_pairing: bool
    notes: str


MODEL_REGISTRY = {
    "vinas_wgan_gp": ModelCapabilities(
        model_id="vinas_wgan_gp",
        display_name="Vinas et al. conditional WGAN-GP",
        paper_url="https://doi.org/10.1093/bioinformatics/btab035",
        code_url="https://github.com/rvinas/adversarial-gene-expression",
        inspected_commit="94fa44dd1bd52d924efd3af0fcd8eeb18bd141a8",
        native_input="log1p expression followed by train-fitted gene-wise z-score",
        supports_expression_generation=True,
        supports_categorical_conditioning=True,
        supports_pretrain_finetune=False,
        supports_counterfactual_pairing=True,
        notes=(
            "The paper conditions on tissue and dataset with learned embeddings. "
            "ARCHS4 pretraining and OSDR fine-tuning require a local extension."
        ),
    ),
    "lacan_diffusion": ModelCapabilities(
        model_id="lacan_diffusion",
        display_name="Lacan et al. landmark-space DDIM",
        paper_url="https://doi.org/10.1186/s12859-026-06470-8",
        code_url="https://forge.ibisc.univ-evry.fr/alacan/rna-diffusion.git",
        inspected_commit="cde890154698fcea96c924804aaff04af3351b48",
        native_input="TPM expression, L1000 landmark selection, and MaxAbs scaling",
        supports_expression_generation=True,
        supports_categorical_conditioning=True,
        supports_pretrain_finetune=False,
        supports_counterfactual_pairing=False,
        notes=(
            "The paper generates landmarks with tissue-conditioned DDIM and "
            "reconstructs target genes. Pretrain/fine-tune is a local extension."
        ),
    ),
    "genejepa": ModelCapabilities(
        model_id="genejepa",
        display_name="GeneJEPA",
        paper_url="https://doi.org/10.1101/2025.10.14.682378",
        code_url="https://github.com/BiostateAI/GeneJEPA",
        inspected_commit="a2f4d7218b17f2f52cc5f1cc94420c8ef1ae3265",
        native_input="ragged nonzero scRNA-seq genes with global log1p z-score",
        supports_expression_generation=False,
        supports_categorical_conditioning=False,
        supports_pretrain_finetune=True,
        supports_counterfactual_pairing=False,
        notes=(
            "GeneJEPA predicts masked-set embeddings and explicitly does not decode "
            "counts. It is a representation benchmark unless a separately validated "
            "decoder is added, which would no longer be the paper model."
        ),
    ),
}


def require_generation_capability(model_id: str) -> None:
    model = MODEL_REGISTRY[model_id]
    if not model.supports_expression_generation:
        raise ValueError(
            f"{model.display_name} does not generate expression matrices. "
            "Use task=representation or choose a native generator."
        )
