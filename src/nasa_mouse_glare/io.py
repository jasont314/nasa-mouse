"""Shared matrix-bundle IO for the mouse GLARE workflow.

The bundle format keeps large matrices out of CSV while preserving row and
column labels in small TSV sidecars. Matrices are always gene x profile.
"""

from __future__ import annotations

from dataclasses import dataclass
import importlib
import json
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class MatrixBundle:
    matrix: object
    genes: list[str]
    profiles: list[str]
    manifest_path: Path
    profile_metadata: object | None = None


def require_import(module_name: str, install_hint: str | None = None):
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        hint = install_hint or f"pip install {module_name}"
        raise SystemExit(
            f"Missing dependency '{module_name}'. Install with: {hint}"
        ) from exc


def normalize_prefix(output_prefix: str | Path) -> Path:
    prefix = Path(output_prefix)
    if prefix.suffix == ".json":
        return prefix.with_suffix("")
    return prefix


def append_suffix(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def manifest_for(output_prefix: str | Path) -> Path:
    prefix = normalize_prefix(output_prefix)
    if str(output_prefix).endswith(".manifest.json"):
        return Path(output_prefix)
    return append_suffix(prefix, ".manifest.json")


def write_lines(path: Path, values: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def read_lines(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def write_matrix_bundle(
    output_prefix: str | Path,
    matrix,
    genes: Iterable[str],
    profiles: Iterable[str],
    profile_metadata=None,
    description: str = "",
) -> Path:
    """Write a gene x profile matrix bundle and return its manifest path."""
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    prefix = normalize_prefix(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    genes = [str(g) for g in genes]
    profiles = [str(p) for p in profiles]

    matrix_path = append_suffix(prefix, ".matrix.npz")
    genes_path = append_suffix(prefix, ".genes.tsv")
    profiles_path = append_suffix(prefix, ".profiles.tsv")
    metadata_path = append_suffix(prefix, ".profile_metadata.tsv")
    manifest_path = append_suffix(prefix, ".manifest.json")

    if scipy_sparse.issparse(matrix):
        scipy_sparse.save_npz(matrix_path, matrix.tocsr())
        matrix_format = "scipy_sparse_npz"
    else:
        np.savez_compressed(matrix_path, matrix=np.asarray(matrix, dtype="float32"))
        matrix_format = "numpy_dense_npz"

    write_lines(genes_path, genes)
    write_lines(profiles_path, profiles)

    has_metadata = profile_metadata is not None
    if has_metadata:
        profile_metadata.to_csv(metadata_path, sep="\t", index=False)

    manifest = {
        "description": description,
        "matrix_path": str(matrix_path),
        "matrix_format": matrix_format,
        "genes_path": str(genes_path),
        "profiles_path": str(profiles_path),
        "profile_metadata_path": str(metadata_path) if has_metadata else "",
        "shape": [len(genes), len(profiles)],
        "orientation": "genes_x_profiles",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def load_matrix_bundle(manifest_path: str | Path) -> MatrixBundle:
    """Load a matrix bundle from its manifest."""
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    matrix_path = Path(manifest["matrix_path"])
    if manifest["matrix_format"] == "scipy_sparse_npz":
        matrix = scipy_sparse.load_npz(matrix_path)
    elif manifest["matrix_format"] == "numpy_dense_npz":
        matrix = np.load(matrix_path)["matrix"]
    else:
        raise ValueError(f"Unsupported matrix format: {manifest['matrix_format']}")

    genes = read_lines(Path(manifest["genes_path"]))
    profiles = read_lines(Path(manifest["profiles_path"]))

    metadata_path = manifest.get("profile_metadata_path") or ""
    metadata = pd.read_csv(metadata_path, sep="\t") if metadata_path else None
    return MatrixBundle(matrix, genes, profiles, manifest_path, metadata)


def dense_matrix(matrix, max_dense_gb: float = 8.0):
    """Return a dense float32 matrix with a guardrail around accidental blowups."""
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )

    n_bytes = matrix.shape[0] * matrix.shape[1] * 4
    n_gb = n_bytes / (1024**3)
    if n_gb > max_dense_gb:
        raise SystemExit(
            f"Refusing to densify matrix {matrix.shape}; estimated float32 size "
            f"is {n_gb:.2f} GB. Use --max-cells, or raise "
            f"--max-dense-gb if this is intentional."
        )
    if scipy_sparse.issparse(matrix):
        return matrix.toarray().astype("float32", copy=False)
    return np.asarray(matrix, dtype="float32")
