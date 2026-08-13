"""Gene annotation helpers shared by the GLARE analyses."""

from __future__ import annotations

from pathlib import Path

from .io import require_import


DEFAULT_TMS_H5AD = "assets/tms/be2af593-fb71-4c76-85a8-3c8400783c2a.h5ad"


def load_tms_gene_symbols(
    tms_h5ad: str | Path = DEFAULT_TMS_H5AD,
) -> dict[str, str]:
    """Map mouse Ensembl IDs to symbols from the retained TMS reference."""
    anndata = require_import(
        "anndata", "pip install -r requirements.txt"
    )
    adata = anndata.read_h5ad(tms_h5ad, backed="r")
    if "feature_name" not in adata.var.columns:
        raise ValueError(f"TMS reference has no feature_name annotation: {tms_h5ad}")
    symbols = {
        str(gene).split(".", 1)[0]: str(symbol)
        for gene, symbol in zip(adata.var_names, adata.var["feature_name"])
        if str(symbol) and str(symbol).lower() != "nan"
    }
    try:
        adata.file.close()
    except AttributeError:
        pass
    return symbols
