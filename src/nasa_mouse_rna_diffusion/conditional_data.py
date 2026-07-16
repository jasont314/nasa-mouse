"""Leakage-safe NASA OSDR input for the exact Lacan ModelDDIM architecture."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from nasa_mouse_generative.config import load_config as load_generative_config
from nasa_mouse_generative.config import PreprocessingConfig
from nasa_mouse_generative.preprocessing import FittedPreprocessor
from nasa_mouse_generative.training_data import _osdr_rows

from .conditional_config import load_conditional_config


ROLES = ("train", "validation", "test")
FORMAT = "nasa_mouse_lacan_conditional_osdr_v1"


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode(values: Iterable[object]) -> list[str]:
    return [
        value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _joint_class_labels(
    rows: pd.DataFrame, covariates: Iterable[str]
) -> pd.Series:
    columns = tuple(map(str, covariates))
    missing = [column for column in columns if column not in rows]
    if missing:
        raise ValueError(f"Conditioning columns are absent from OSDR metadata: {missing}")
    normalized = rows.loc[:, columns].fillna("__missing__").astype(str)
    return normalized.apply(
        lambda row: "||".join(
            f"{column}={row[column]}" for column in columns
        ),
        axis=1,
    )


def _within_study_roles(
    rows: pd.DataFrame,
    *,
    seed: int,
    validation_fraction: float,
    test_fraction: float,
) -> pd.Series:
    """Split samples within each study/tissue/condition stratum."""

    roles = pd.Series("train", index=rows.index, dtype="object")
    strata = ["accession", "tissue", "condition"]
    for keys, group in rows.groupby(strata, dropna=False, sort=True):
        ordered = sorted(
            group.index,
            key=lambda index: int.from_bytes(
                hashlib.sha256(
                    f"{seed}|{keys}|{rows.loc[index, 'profile_id']}".encode(
                        "utf-8", "replace"
                    )
                ).digest()[:8],
                "big",
            ),
        )
        count = len(ordered)
        if count < 3:
            continue
        test_count = min(max(1, round(count * test_fraction)), count - 2)
        validation_count = min(
            max(1, round(count * validation_fraction)), count - test_count - 1
        )
        roles.loc[ordered[:test_count]] = "test"
        roles.loc[
            ordered[test_count : test_count + validation_count]
        ] = "validation"
    return roles


def _explicit_accession_roles(
    rows: pd.DataFrame,
    *,
    validation_accessions: Iterable[str],
    test_accessions: Iterable[str],
) -> tuple[pd.Series, dict[str, object]]:
    validation = set(map(str, validation_accessions))
    test = set(map(str, test_accessions))
    available = set(rows["accession"].astype(str))
    missing = sorted((validation | test) - available)
    if missing:
        raise ValueError(f"Explicit split accessions are absent from the cohort: {missing}")
    roles = pd.Series("train", index=rows.index, dtype="object")
    roles.loc[rows["accession"].astype(str).isin(validation)] = "validation"
    roles.loc[rows["accession"].astype(str).isin(test)] = "test"
    accessions_by_role = {
        role: sorted(rows.loc[roles.eq(role), "accession"].astype(str).unique())
        for role in ROLES
    }
    if any(not values for values in accessions_by_role.values()):
        raise ValueError("Explicit accession split produced an empty role")
    return roles, {
        "kind": "explicit_accession_holdout",
        "split_unit": "NASA OSDR accession",
        "accessions_by_role": accessions_by_role,
    }


def _dense(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=np.float64)


def _full_transcriptome_tpm(
    adata: ad.AnnData,
    row_indices: np.ndarray,
    landmark_indices: np.ndarray,
    gene_lengths: np.ndarray,
    *,
    chunk_size: int = 64,
) -> np.ndarray:
    lengths = np.asarray(gene_lengths, dtype=np.float64)
    valid = np.isfinite(lengths) & (lengths > 0)
    if not valid.any():
        raise ValueError("No OSDR genes have a positive gene length")
    if not valid[np.asarray(landmark_indices, dtype=int)].all():
        raise ValueError("One or more landmark genes lack a positive gene length")
    valid_indices = np.flatnonzero(valid)
    valid_lengths_kb = lengths[valid_indices] / 1000.0
    landmark_lengths_kb = lengths[np.asarray(landmark_indices)] / 1000.0
    output = np.empty(
        (len(row_indices), len(landmark_indices)), dtype=np.float32
    )
    for start in range(0, len(row_indices), int(chunk_size)):
        end = min(start + int(chunk_size), len(row_indices))
        source_rows = row_indices[start:end]
        full = _dense(adata.X[source_rows][:, valid_indices])
        if (full < 0).any():
            raise ValueError("OSDR API matrix contains negative raw expression")
        denominator = np.sum(full / valid_lengths_kb.reshape(1, -1), axis=1)
        if not np.isfinite(denominator).all() or (denominator <= 0).any():
            raise ValueError("OSDR API matrix contains an empty or invalid library")
        landmark = _dense(adata.X[source_rows][:, landmark_indices])
        output[start:end] = (
            landmark
            / landmark_lengths_kb.reshape(1, -1)
            / denominator.reshape(-1, 1)
            * 1_000_000.0
        ).astype(np.float32)
    return output


def _source_identity(path: str | Path) -> dict[str, object]:
    source = Path(path)
    stat = source.stat()
    return {
        "path": str(source.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": _sha256(source),
    }


def prepare_conditional(config_path: str | Path, *, force: bool = False) -> Path:
    config = load_conditional_config(config_path)
    options = config["data"]
    output = Path(options["prepared_h5"])
    manifest_path = output.with_suffix(".manifest.json")
    config_sha256 = _sha256(config_path)
    data_contract_sha256 = hashlib.sha256(
        json.dumps(
            {"seed": config["run"]["seed"], "data": options},
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    source_identity = _source_identity(options["osdr_h5ad"])
    custom_preprocessing = options.get("preprocessing")
    scale_source_option = str(options.get("maxabs_scale_source", "osdr_train"))
    pretrained_classes_option = str(options.get("pretrained_classes_h5", ""))
    optional_source_identities: dict[str, dict[str, object]] = {}
    if scale_source_option != "osdr_train" and not custom_preprocessing:
        optional_source_identities["maxabs_scale_source"] = _source_identity(
            scale_source_option
        )
    if pretrained_classes_option:
        optional_source_identities["pretrained_classes_h5"] = (
            optional_source_identities.get("maxabs_scale_source")
            if pretrained_classes_option == scale_source_option
            else _source_identity(pretrained_classes_option)
        )
    if output.exists() and manifest_path.exists() and not force:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("data_contract_sha256") == data_contract_sha256
            and manifest.get("source", {}).get("osdr_api_h5ad") == source_identity
            and all(
                manifest.get("source", {}).get(name) == identity
                for name, identity in optional_source_identities.items()
            )
        ):
            print(f"[conditional-ddim:data] using existing {output}", flush=True)
            return output

    base = load_generative_config(options["generative_config"])
    requested_tissues = tuple(map(str, options.get("tissues", [])))
    base = replace(
        base,
        data=replace(
            base.data,
            osdr_h5ad=str(options["osdr_h5ad"]),
            osdr_accession_scope=str(options["accession_scope"]),
            osdr_include_accessions=tuple(map(str, options.get("include_accessions", []))),
            osdr_exclude_accessions=tuple(map(str, options.get("exclude_accessions", []))),
            osdr_tissues=requested_tissues,
        ),
        training=replace(
            base.training,
            tissue_mode=str(options["tissue_mode"]),
            condition_on_flight=True,
            seed=int(config["run"]["seed"]),
        ),
    )
    tissue_override = None
    if str(options["tissue_mode"]) == "per_tissue":
        if len(requested_tissues) != 1:
            raise ValueError("Per-tissue conditional DDIM requires exactly one tissue")
        tissue_override = requested_tissues[0]
    adata, rows, split_metadata = _osdr_rows(base, tissue_override)
    split_strategy = str(options.get("split_strategy", "accession_holdout"))
    if options.get("validation_accessions") or options.get("test_accessions"):
        rows = rows.copy()
        rows["role"], split_metadata = _explicit_accession_roles(
            rows,
            validation_accessions=options.get("validation_accessions", []),
            test_accessions=options.get("test_accessions", []),
        )
    elif split_strategy == "within_study_stratified":
        rows = rows.copy()
        rows["role"] = _within_study_roles(
            rows,
            seed=int(config["run"]["seed"]),
            validation_fraction=float(base.validation.pooled_validation_fraction),
            test_fraction=float(base.validation.pooled_test_fraction),
        )
        split_metadata = {
            "kind": "within_study_stratified",
            "split_unit": "sample within accession/tissue/condition",
            "accessions": int(rows["accession"].nunique()),
            "limitation": (
                "Evaluates generation for studies represented during training; "
                "it does not measure unseen-study generalization."
            ),
        }
    genes = [str(value).split(".", 1)[0] for value in adata.var_names]
    if len(set(genes)) != len(genes):
        raise ValueError("OSDR API matrix has duplicate versionless Ensembl genes")
    gene_map = {gene: index for index, gene in enumerate(genes)}
    panel = pd.read_csv(options["mouse_landmark_panel"], sep="\t")
    landmark_genes = panel["mouse_ensembl_gene"].astype(str).tolist()
    if len(landmark_genes) != int(options["landmark_dimensions"]):
        raise ValueError("Mouse landmark panel does not contain exactly 974 genes")
    missing_landmarks = [gene for gene in landmark_genes if gene not in gene_map]
    if missing_landmarks:
        raise ValueError(
            f"OSDR API matrix lacks {len(missing_landmarks)} mouse landmarks"
        )
    landmark_indices = np.asarray(
        [gene_map[gene] for gene in landmark_genes], dtype=np.int64
    )
    length_table = pd.read_csv(options["gene_lengths"], sep="\t").set_index(
        "gene_id"
    )["length_bp"]
    lengths = np.asarray(
        [length_table.get(gene, np.nan) for gene in genes], dtype=np.float64
    )
    source_rows = rows["_row_index"].to_numpy(dtype=np.int64)
    tpm = _full_transcriptome_tpm(
        adata,
        source_rows,
        landmark_indices,
        lengths,
    )

    covariates = tuple(map(str, options["conditioning_covariates"]))
    samples = rows.copy()
    samples["class_label"] = _joint_class_labels(samples, covariates)
    train_labels = set(
        samples.loc[samples["role"].eq("train"), "class_label"].astype(str)
    )
    heldout_labels = set(
        samples.loc[~samples["role"].eq("train"), "class_label"].astype(str)
    )
    unseen = sorted(heldout_labels - train_labels)
    if unseen:
        raise ValueError(
            "Upstream ModelDDIM uses a joint class vocabulary, but held-out OSDR "
            f"profiles contain {len(unseen)} classes absent from training: {unseen[:5]}"
        )
    additional_classes: list[str] = []
    pretrained_classes_h5 = str(options.get("pretrained_classes_h5", ""))
    if pretrained_classes_h5:
        with h5py.File(pretrained_classes_h5, "r") as handle:
            pretrained_tissues = _decode(handle["classes"][:])
        additional_classes = [
            f"tissue={tissue}||condition=reference"
            for tissue in pretrained_tissues
        ]
    classes = sorted(train_labels.union(additional_classes))
    class_map = {label: index for index, label in enumerate(classes)}
    samples["class_index"] = samples["class_label"].map(class_map).astype(int)
    role_indices = {
        role: np.flatnonzero(samples["role"].astype(str).eq(role).to_numpy())
        for role in ROLES
    }
    if any(len(indices) == 0 for indices in role_indices.values()):
        raise ValueError("Conditional DDIM requires nonempty train/validation/test roles")
    preprocessing_audit: dict[str, object] | None = None
    preprocessing_dir: Path | None = None
    scale_source = scale_source_option
    if custom_preprocessing:
        spec_options = dict(custom_preprocessing)
        if "harmonization_covariates" in spec_options:
            spec_options["harmonization_covariates"] = tuple(
                spec_options["harmonization_covariates"]
            )
        processor = FittedPreprocessor(
            PreprocessingConfig(**spec_options), seed=int(config["run"]["seed"])
        )
        studies = samples["study"].astype(str).to_numpy()
        scaled = np.empty_like(tpm, dtype=np.float32)
        train_indices = role_indices["train"]
        scaled[train_indices] = processor.fit_transform(
            tpm[train_indices], studies[train_indices], metadata=samples.iloc[train_indices]
        )
        allow_transductive = bool(
            spec_options.get("unseen_study_policy") == "transductive_unlabeled"
        )
        if allow_transductive:
            unseen_mask = np.asarray(
                [study not in processor.study_stats for study in studies],
                dtype=bool,
            )
            unseen_indices = np.flatnonzero(unseen_mask)
            processor.fit_additional_study_stats(
                tpm[unseen_indices], studies[unseen_indices]
            )
        for role in ("validation", "test"):
            indices = role_indices[role]
            scaled[indices] = processor.transform(
                tpm[indices],
                studies[indices],
                allow_transductive=allow_transductive,
                metadata=samples.iloc[indices],
            )
        preprocessing_dir = output.with_suffix(".preprocessing")
        processor.save(preprocessing_dir)
        preprocessing_audit = processor.audit()
        preprocessing_audit["allow_transductive"] = allow_transductive
        scale = np.ones(len(landmark_genes), dtype=np.float32)
        scale_source = "fitted_preprocessor"
    else:
        if scale_source == "osdr_train":
            scale = np.max(np.abs(tpm[role_indices["train"]]), axis=0).astype(
                np.float32
            )
        else:
            with h5py.File(scale_source, "r") as handle:
                source_genes = _decode(handle["genes"][:])
                if source_genes != landmark_genes:
                    raise ValueError("Pretraining MaxAbs gene order differs from OSDR")
                scale = np.asarray(handle["maxabs_scale"][:], dtype=np.float32)
        scale[scale == 0] = 1.0
        scaled = (tpm / scale.reshape(1, -1)).astype(np.float32)
    if not np.isfinite(scaled).all():
        raise FloatingPointError("Conditional DDIM prepared matrix is not finite")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.h5")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["format"] = FORMAT
        handle.attrs["normalization"] = str(options["normalization"])
        handle.attrs["preprocessing_dir"] = (
            str(preprocessing_dir) if preprocessing_dir is not None else ""
        )
        handle.create_dataset(
            "genes",
            data=np.asarray(landmark_genes, dtype=h5py.string_dtype("utf-8")),
        )
        handle.create_dataset(
            "classes", data=np.asarray(classes, dtype=h5py.string_dtype("utf-8"))
        )
        handle.create_dataset(
            "conditioning_covariates",
            data=np.asarray(covariates, dtype=h5py.string_dtype("utf-8")),
        )
        handle.create_dataset("maxabs_scale", data=scale)
        for role, indices in role_indices.items():
            group = handle.create_group(role)
            group.create_dataset(
                "expression",
                data=scaled[indices],
                chunks=(min(256, len(indices)), len(landmark_genes)),
                compression="lzf",
            )
            group.create_dataset(
                "tpm",
                data=tpm[indices],
                chunks=(min(256, len(indices)), len(landmark_genes)),
                compression="lzf",
            )
            group.create_dataset(
                "class_index", data=samples.iloc[indices]["class_index"].to_numpy()
            )
            group.create_dataset(
                "source_row", data=samples.iloc[indices]["_row_index"].to_numpy()
            )
    temporary.replace(output)
    samples_path = output.with_suffix(".samples.tsv.gz")
    samples.to_csv(samples_path, sep="\t", index=False, compression="gzip")
    coverage = (
        samples.groupby(["role", "tissue", "condition"], dropna=False)
        .size()
        .rename("profiles")
        .reset_index()
    )
    coverage_path = output.with_suffix(".coverage.tsv")
    coverage.to_csv(coverage_path, sep="\t", index=False)
    sources = {
        "osdr_api_h5ad": source_identity,
        "mouse_landmark_panel": _source_identity(options["mouse_landmark_panel"]),
        "gene_lengths": _source_identity(options["gene_lengths"]),
    }
    if not custom_preprocessing and scale_source != "osdr_train":
        sources["maxabs_scale_source"] = optional_source_identities[
            "maxabs_scale_source"
        ]
    if pretrained_classes_h5:
        sources["pretrained_classes_h5"] = optional_source_identities[
            "pretrained_classes_h5"
        ]
    manifest = {
        "format": FORMAT,
        "prepared_h5": str(output),
        "prepared_h5_sha256": _sha256(output),
        "samples": str(samples_path),
        "coverage": str(coverage_path),
        "config": str(Path(config_path).resolve()),
        "config_sha256": config_sha256,
        "data_contract_sha256": data_contract_sha256,
        "raw_integrated_osdr_h5_used": False,
        "source": sources,
        "split": split_metadata,
        "profiles": {role: int(len(indices)) for role, indices in role_indices.items()},
        "genes": len(landmark_genes),
        "full_transcriptome_genes": len(genes),
        "full_transcriptome_genes_with_lengths": int(np.isfinite(lengths).sum()),
        "classes": classes,
        "additional_pretraining_classes": additional_classes,
        "conditioning_covariates": list(covariates),
        "maxabs_scale_source": scale_source,
        "normalization": (
            "TPM denominator uses every OSDR API gene with a positive GENCODE M39 "
            "length; the 974 landmarks are selected afterward; the serialized "
            "fold-aware preprocessor supplies transformation and harmonization."
            if custom_preprocessing
            else (
                "TPM denominator uses every OSDR API gene with a positive GENCODE "
                "M39 length; the 974 landmarks are selected afterward; MaxAbs "
                "comes from "
                + (
                    "OSDR training accessions only."
                    if scale_source == "osdr_train"
                    else "the declared ARCHS4 pretraining partition."
                )
            )
        ),
        "preprocessing": preprocessing_audit,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return output


def load_conditional_prepared(path: str | Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("format") != FORMAT:
            raise ValueError(f"Unsupported conditional DDIM data: {path}")
        result["genes"] = _decode(handle["genes"][:])
        result["classes"] = _decode(handle["classes"][:])
        result["conditioning_covariates"] = _decode(
            handle["conditioning_covariates"][:]
        )
        result["maxabs_scale"] = np.asarray(
            handle["maxabs_scale"][:], dtype=np.float32
        )
        preprocessing_dir = str(handle.attrs.get("preprocessing_dir", ""))
        result["preprocessing"] = (
            FittedPreprocessor.load(preprocessing_dir)
            if preprocessing_dir
            else None
        )
        for role in ROLES:
            result[role] = {
                "expression": np.asarray(
                    handle[f"{role}/expression"][:], dtype=np.float32
                ),
                "tpm": np.asarray(handle[f"{role}/tpm"][:], dtype=np.float32),
                "class_index": np.asarray(
                    handle[f"{role}/class_index"][:], dtype=np.int64
                ),
                "source_row": np.asarray(
                    handle[f"{role}/source_row"][:], dtype=np.int64
                ),
            }
    if len(result["genes"]) != 974:
        raise ValueError("Conditional DDIM prepared data must have 974 genes")
    return result
