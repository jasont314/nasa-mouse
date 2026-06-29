"""Train and score OntoVAE on API-derived mouse OSDR/ARCHS4 AnnData inputs."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

from .io import require_import


ONTO_VAE_SOURCE = Path(__file__).resolve().parents[1] / "onto-vae"
KEY = "0_0"


def ensure_onto_vae_on_path() -> None:
    source = str(ONTO_VAE_SOURCE)
    if source not in sys.path:
        sys.path.insert(0, source)


def matrix_to_numpy(values):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    scipy_sparse = require_import(
        "scipy.sparse", "pip install -r requirements-nasa-mouse-glare.txt"
    )
    if scipy_sparse.issparse(values):
        return values.toarray().astype("float32", copy=False)
    return np.asarray(values, dtype="float32")


def counts_matrix(adata):
    if "counts" in adata.layers:
        return matrix_to_numpy(adata.layers["counts"])
    return matrix_to_numpy(adata.X)


def log1p_cpm(counts):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    counts = np.asarray(counts, dtype="float32")
    library = counts.sum(axis=1, keepdims=True)
    library[library <= 0] = 1.0
    return np.log1p(counts / library * 1_000_000.0).astype("float32", copy=False)


def standardize(train, *others, clip: float = 10.0):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    mean = train.mean(axis=0, keepdims=True)
    std = train.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    arrays = []
    for array in (train, *others):
        transformed = (array - mean) / std
        if clip:
            transformed = np.clip(transformed, -clip, clip)
        arrays.append(transformed.astype("float32", copy=False))
    return arrays, mean.ravel().astype("float32"), std.ravel().astype("float32")


def align_inputs(query, reference=None, min_term_genes: int = 5):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    query_genes = query.var_names.astype(str).tolist()
    if reference is None:
        common_genes = query_genes
        query_positions = list(range(len(query_genes)))
        reference_positions = None
    else:
        reference_genes = reference.var_names.astype(str).tolist()
        reference_set = set(reference_genes)
        common_genes = [gene for gene in query_genes if gene in reference_set]
        reference_lookup = {gene: index for index, gene in enumerate(reference_genes)}
        query_positions = [query_genes.index(gene) for gene in common_genes]
        reference_positions = [reference_lookup[gene] for gene in common_genes]

    mask = np.asarray(query.varm["I"], dtype="float32")[query_positions, :]
    terms = list(map(str, query.uns["terms"]))
    descriptions = list(map(str, query.uns.get("term_descriptions", terms)))
    term_gene_counts = mask.sum(axis=0)
    term_keep = term_gene_counts >= min_term_genes
    if not bool(term_keep.any()):
        raise SystemExit("No OntoVAE terms retained after term gene-count filtering.")
    mask = mask[:, term_keep]
    terms = [term for term, keep in zip(terms, term_keep) if keep]
    descriptions = [
        description for description, keep in zip(descriptions, term_keep) if keep
    ]
    gene_keep = mask.sum(axis=1) > 0
    if not bool(gene_keep.any()):
        raise SystemExit("No OntoVAE genes retained after dropping unconnected genes.")

    query_aligned = query[:, [query_positions[i] for i, keep in enumerate(gene_keep) if keep]].copy()
    if reference is None:
        reference_aligned = None
    else:
        reference_aligned = reference[
            :,
            [reference_positions[i] for i, keep in enumerate(gene_keep) if keep],
        ].copy()
    mask = mask[gene_keep, :].astype("float32", copy=False)
    genes = [gene for gene, keep in zip(common_genes, gene_keep) if keep]
    return query_aligned, reference_aligned, genes, terms, descriptions, mask


def build_flat_ontobj(genes, terms, descriptions, mask, datasets):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")

    gene_counts = np.asarray(mask.sum(axis=0), dtype=int)
    annot = pd.DataFrame(
        {
            "ID": terms,
            "Name": descriptions,
            "depth": 0,
            "children": gene_counts,
            "parents": 0,
            "descendants": 0,
            "desc_genes": gene_counts,
            "genes": gene_counts,
        }
    )
    graph = {
        gene: [terms[index] for index in np.flatnonzero(mask[row_index] > 0)]
        for row_index, gene in enumerate(genes)
    }
    desc_genes = {
        term: [genes[index] for index in np.flatnonzero(mask[:, term_index] > 0)]
        for term_index, term in enumerate(terms)
    }
    return SimpleNamespace(
        description="flat_mouse_reactome_ensembl",
        genes={KEY: list(genes)},
        annot={KEY: annot},
        graph={KEY: graph},
        desc_genes={KEY: desc_genes},
        masks={KEY: {"decoder": [mask.astype("float32", copy=False)]}},
        data={KEY: {name: value for name, value in datasets.items()}},
    )


def split_indices(n_samples: int, train_frac: float, seed: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    indices = np.random.RandomState(seed=seed).permutation(n_samples)
    if n_samples < 3:
        return indices, indices
    split = max(1, min(n_samples - 1, round(n_samples * train_frac)))
    return indices[:split], indices[split:]


def train_loop(
    model,
    data,
    *,
    epochs: int,
    batch_size: int,
    lr: float,
    kl_coeff: float,
    train_frac: float,
    seed: int,
):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    optim = require_import("torch.optim", "pip install -r requirements-nasa-mouse-glare.txt")
    functional = require_import("torch.nn.functional", "pip install -r requirements-nasa-mouse-glare.txt")

    model.to(model.device)
    train_indices, val_indices = split_indices(data.shape[0], train_frac, seed)
    train = torch.tensor(data[train_indices], dtype=torch.float32, device=model.device)
    val = torch.tensor(data[val_indices], dtype=torch.float32, device=model.device)
    optimizer = optim.AdamW(model.parameters(), lr=lr)
    history = []
    best_state = None
    best_val = float("inf")
    start = time.time()
    rng = np.random.RandomState(seed)

    for epoch in range(epochs):
        model.train()
        order = rng.permutation(train.shape[0])
        train_loss = 0.0
        n_train = 0
        for start_index in range(0, len(order), batch_size):
            batch_index = order[start_index : start_index + batch_size]
            batch = train[batch_index]
            optimizer.zero_grad()
            reconstruction, mu, log_var = model.forward(batch)
            kl_loss = -0.5 * torch.sum(1.0 + log_var - mu.pow(2) - log_var.exp())
            rec_loss = functional.mse_loss(reconstruction, batch, reduction="sum")
            loss = rec_loss + kl_coeff * kl_loss
            loss.backward()
            for layer_index in range(len(model.decoder.decoder)):
                layer = model.decoder.decoder[layer_index][0]
                if layer.weight.grad is not None:
                    layer.weight.grad = torch.mul(
                        layer.weight.grad,
                        model.decoder.masks[layer_index],
                    )
            optimizer.step()
            for layer_index in range(len(model.decoder.decoder)):
                model.decoder.decoder[layer_index][0].weight.data = (
                    model.decoder.decoder[layer_index][0].weight.data.clamp(0)
                )
            train_loss += float(loss.item())
            n_train += int(batch.shape[0])

        model.eval()
        with torch.no_grad():
            reconstruction, mu, log_var = model.forward(val)
            kl_loss = -0.5 * torch.sum(1.0 + log_var - mu.pow(2) - log_var.exp())
            rec_loss = functional.mse_loss(reconstruction, val, reduction="sum")
            val_loss = float((rec_loss + kl_coeff * kl_loss).item())
        train_loss_per_sample = train_loss / max(n_train, 1)
        val_loss_per_sample = val_loss / max(int(val.shape[0]), 1)
        history.append(
            {
                "epoch": epoch + 1,
                "train_loss_per_sample": train_loss_per_sample,
                "val_loss_per_sample": val_loss_per_sample,
            }
        )
        if val_loss_per_sample < best_val:
            best_val = val_loss_per_sample
            best_state = copy.deepcopy(model.state_dict())

    if best_state is not None:
        model.load_state_dict(best_state)
    return {
        "epochs_completed": int(epochs),
        "best_val_loss_per_sample": float(best_val),
        "training_seconds": float(time.time() - start),
        "history": history,
    }


def posterior_scores(model, data, terms, batch_size: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")

    model.eval()
    chunks = []
    with torch.no_grad():
        for start in range(0, data.shape[0], batch_size):
            batch = torch.tensor(
                data[start : start + batch_size],
                dtype=torch.float32,
                device=model.device,
            )
            mu, _ = model.encoder(batch)
            array = mu.detach().cpu().numpy()
            array = array.reshape(array.shape[0], len(terms), model.neuronnum).mean(axis=2)
            chunks.append(array)
    return np.vstack(chunks).astype("float32", copy=False)


def write_scores(path: Path, obs, scores, terms):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    frame = obs.reset_index(names="obs_name").copy()
    frame = pd.concat([frame, pd.DataFrame(scores, columns=terms)], axis=1)
    frame.to_csv(path, sep="\t", index=False)
    return frame


def write_terms(path: Path, terms, descriptions, mask):
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")
    pd.DataFrame(
        {
            "term": terms,
            "description": descriptions,
            "n_genes": mask.sum(axis=0).astype(int),
        }
    ).to_csv(path, sep="\t", index=False)


def write_top_gene_weights(path: Path, model, genes, terms, mask, top_n: int):
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    pd = require_import("pandas", "pip install -r requirements-nasa-mouse-glare.txt")

    weight = model.decoder.decoder[-1][0].weight.detach().cpu().numpy()
    weight = weight.reshape(len(genes), len(terms), model.neuronnum).mean(axis=2)
    rows = []
    for term_index, term in enumerate(terms):
        connected = np.flatnonzero(mask[:, term_index] > 0)
        if len(connected) == 0:
            continue
        order = connected[np.argsort(-weight[connected, term_index])[:top_n]]
        for rank, gene_index in enumerate(order, start=1):
            rows.append(
                {
                    "term": term,
                    "rank": rank,
                    "gene": genes[gene_index],
                    "decoder_weight": float(weight[gene_index, term_index]),
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def run(args) -> Path:
    ensure_onto_vae_on_path()
    ad = require_import("anndata", "pip install -r requirements-nasa-mouse-glare.txt")
    np = require_import("numpy", "pip install -r requirements-nasa-mouse-glare.txt")
    torch = require_import("torch", "pip install -r requirements-nasa-mouse-glare.txt")
    from onto_vae.vae_model import OntoVAE

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    query = ad.read_h5ad(args.query_h5ad)
    reference = ad.read_h5ad(args.reference_h5ad) if args.reference_h5ad else None
    query, reference, genes, terms, descriptions, mask = align_inputs(
        query,
        reference=reference,
        min_term_genes=args.min_term_genes,
    )
    query_log = log1p_cpm(counts_matrix(query))
    if reference is not None:
        reference_log = log1p_cpm(counts_matrix(reference))
        (reference_x, query_x), mean, std = standardize(
            reference_log,
            query_log,
            clip=args.clip,
        )
        datasets = {"reference": reference_x, "query": query_x}
        train_dataset = "reference"
    else:
        (query_x,), mean, std = standardize(query_log, clip=args.clip)
        reference_x = None
        datasets = {"query": query_x}
        train_dataset = "query"

    ontobj = build_flat_ontobj(genes, terms, descriptions, mask, datasets)
    model = OntoVAE(
        ontobj=ontobj,
        dataset=train_dataset,
        top_thresh=0,
        bottom_thresh=0,
        neuronnum=args.neuronnum,
        drop=args.drop,
        z_drop=args.z_drop,
    )
    model.to(model.device)

    histories = {}
    pretrained_scores_path = None
    if reference_x is not None:
        histories["reference"] = train_loop(
            model,
            reference_x,
            epochs=args.reference_epochs,
            batch_size=args.batch_size,
            lr=args.learning_rate,
            kl_coeff=args.kl_coeff,
            train_frac=args.train_frac,
            seed=args.seed,
        )
        torch.save(
            {"model_state_dict": model.state_dict(), "terms": terms, "genes": genes},
            output_dir / "reference_pretrained_model.pt",
        )
        pretrained_scores = posterior_scores(model, query_x, terms, args.batch_size)
        pretrained_scores_path = output_dir / "pretrained_query_pathway_scores.tsv"
        write_scores(pretrained_scores_path, query.obs, pretrained_scores, terms)
        histories["query_finetune"] = train_loop(
            model,
            query_x,
            epochs=args.query_epochs,
            batch_size=args.batch_size,
            lr=args.finetune_learning_rate or args.learning_rate,
            kl_coeff=args.kl_coeff,
            train_frac=args.train_frac,
            seed=args.seed + 17,
        )
    else:
        histories["direct"] = train_loop(
            model,
            query_x,
            epochs=args.query_epochs,
            batch_size=args.batch_size,
            lr=args.learning_rate,
            kl_coeff=args.kl_coeff,
            train_frac=args.train_frac,
            seed=args.seed,
        )

    torch.save(
        {"model_state_dict": model.state_dict(), "terms": terms, "genes": genes},
        output_dir / "model.pt",
    )
    scores = posterior_scores(model, query_x, terms, args.batch_size)
    scores_path = output_dir / "pathway_scores.tsv"
    write_scores(scores_path, query.obs, scores, terms)
    terms_path = output_dir / "terms.tsv"
    write_terms(terms_path, terms, descriptions, mask)
    gene_weights_path = output_dir / "term_gene_weights_top.tsv"
    write_top_gene_weights(gene_weights_path, model, genes, terms, mask, args.top_genes)
    np.savez_compressed(output_dir / "normalization_stats.npz", mean=mean, std=std)

    mode = getattr(args, "run_mode", "") or (
        "archs4_pretrain_osdr_finetune" if reference_x is not None else "direct_osdr"
    )
    summary = {
        "method": "OntoVAE",
        "mode": mode,
        "query_h5ad": str(args.query_h5ad),
        "reference_h5ad": str(args.reference_h5ad or ""),
        "output_dir": str(output_dir),
        "normalization": "log1p(CPM) followed by reference-gene z-score for pretrained runs or query-gene z-score for direct runs",
        "ontology": {
            "type": "flat Reactome/GMT decoder mask",
            "limitations": "OntoVAE's native GO/HPO DAG hierarchy is not used; each retained pathway/module is a root-level latent program.",
            "min_term_genes": int(args.min_term_genes),
        },
        "training_design": {
            "query_mapping_mode": "not native in OntoVAE",
            "implemented_equivalent": (
                "ARCHS4 pretraining followed by OSDR fine-tuning from pretrained weights"
                if reference_x is not None
                else "OSDR-only direct training"
            ),
            "reference_epochs": int(args.reference_epochs if reference_x is not None else 0),
            "query_epochs": int(args.query_epochs),
            "neuronnum": int(args.neuronnum),
            "kl_coeff": float(args.kl_coeff),
            "learning_rate": float(args.learning_rate),
        },
        "torch": {
            "version": str(torch.__version__),
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_count": int(torch.cuda.device_count()),
            "cuda_device_name": (
                torch.cuda.get_device_name(0) if torch.cuda.is_available() else ""
            ),
            "model_device": str(model.device),
        },
        "counts": {
            "query_samples": int(query.n_obs),
            "reference_samples": int(reference.n_obs) if reference is not None else 0,
            "genes": int(len(genes)),
            "terms": int(len(terms)),
        },
        "training": histories,
        "outputs": {
            "scores": str(scores_path),
            "pretrained_query_scores": str(pretrained_scores_path or ""),
            "terms": str(terms_path),
            "top_gene_weights": str(gene_weights_path),
            "model": str(output_dir / "model.pt"),
            "reference_pretrained_model": str(output_dir / "reference_pretrained_model.pt")
            if reference_x is not None
            else "",
            "normalization_stats": str(output_dir / "normalization_stats.npz"),
        },
    }
    summary_path = output_dir / "training_summary.json"
    summary["outputs"]["summary"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    readme = [
        "# OntoVAE Run",
        "",
        f"- Mode: {summary['mode']}",
        f"- Query samples: {summary['counts']['query_samples']}",
        f"- Reference samples: {summary['counts']['reference_samples']}",
        f"- Genes: {summary['counts']['genes']}",
        f"- Terms: {summary['counts']['terms']}",
        f"- Device: {summary['torch']['model_device']} {summary['torch']['cuda_device_name']}",
        "",
        "This run uses a flat Reactome/GMT OntoVAE decoder mask and MSE reconstruction on log1p(CPM) z-scored expression.",
    ]
    (output_dir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train OntoVAE on mouse OSDR/ARCHS4 inputs.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--reference-h5ad", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--reference-epochs", type=int, default=60)
    parser.add_argument("--query-epochs", type=int, default=60)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--finetune-learning-rate", type=float, default=None)
    parser.add_argument("--kl-coeff", type=float, default=1e-4)
    parser.add_argument("--train-frac", type=float, default=0.9)
    parser.add_argument("--seed", type=int, default=2020)
    parser.add_argument("--neuronnum", type=int, default=1)
    parser.add_argument("--drop", type=float, default=0.1)
    parser.add_argument("--z-drop", type=float, default=0.1)
    parser.add_argument("--clip", type=float, default=10.0)
    parser.add_argument("--min-term-genes", type=int, default=5)
    parser.add_argument("--top-genes", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
