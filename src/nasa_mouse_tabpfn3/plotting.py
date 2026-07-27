"""Plot helpers for TabPFN3 OSDR classification outputs."""

from __future__ import annotations

from pathlib import Path

from nasa_mouse_glare.io import require_import


def _setup():
    matplotlib = require_import("matplotlib", "pip install matplotlib")

    matplotlib.use("Agg")
    plt = require_import("matplotlib.pyplot", "pip install matplotlib")
    sns = require_import("seaborn", "pip install seaborn")
    return plt, sns


def plot_confusion(predictions, output_dir: Path, *, label: str) -> Path | None:
    metrics = require_import("sklearn.metrics", "pip install scikit-learn")
    plt, sns = _setup()
    if predictions.empty:
        return None
    cm = metrics.confusion_matrix(
        predictions["y_true"], predictions["y_pred"], labels=[0, 1]
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{label}_confusion_matrix.png"
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=["GC", "FLT"],
        yticklabels=["GC", "FLT"],
        cbar=False,
        ax=ax,
    )
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Observed")
    ax.set_title(label.replace("_", " "))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_roc_pr(predictions, output_dir: Path, *, label: str) -> list[Path]:
    metrics = require_import("sklearn.metrics", "pip install scikit-learn")
    plt, _ = _setup()
    if predictions.empty or predictions["y_true"].nunique() < 2:
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    fpr, tpr, _ = metrics.roc_curve(predictions["y_true"], predictions["p_flight"])
    auroc = metrics.auc(fpr, tpr)
    path = output_dir / f"{label}_roc.png"
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    ax.plot(fpr, tpr, color="#1f77b4", lw=2, label=f"AUROC={auroc:.3f}")
    ax.plot([0, 1], [0, 1], color="#777777", lw=1, linestyle="--")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(frameon=False)
    ax.set_title(label.replace("_", " "))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    precision, recall, _ = metrics.precision_recall_curve(
        predictions["y_true"], predictions["p_flight"]
    )
    auprc = metrics.average_precision_score(
        predictions["y_true"], predictions["p_flight"]
    )
    path = output_dir / f"{label}_precision_recall.png"
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    ax.plot(recall, precision, color="#2ca02c", lw=2, label=f"AUPRC={auprc:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(frameon=False)
    ax.set_title(label.replace("_", " "))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_importance(importance, output_dir: Path, *, label: str, top_n: int = 25) -> Path | None:
    plt, sns = _setup()
    if importance.empty:
        return None
    summary = (
        importance.groupby("gene_id", as_index=False)["mean_decrease_balanced_accuracy"]
        .mean()
        .sort_values("mean_decrease_balanced_accuracy", ascending=False)
        .head(int(top_n))
    )
    if summary.empty:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{label}_top_gene_importance.png"
    fig_height = max(4.0, 0.24 * len(summary) + 1.0)
    fig, ax = plt.subplots(figsize=(7.0, fig_height))
    sns.barplot(
        data=summary,
        x="mean_decrease_balanced_accuracy",
        y="gene_id",
        color="#4c78a8",
        ax=ax,
    )
    ax.set_xlabel("Mean decrease in balanced accuracy")
    ax.set_ylabel("Gene")
    ax.set_title(label.replace("_", " "))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_plots(predictions, importance, output_dir: Path, *, label: str) -> list[str]:
    paths: list[Path] = []
    path = plot_confusion(predictions, output_dir, label=label)
    if path is not None:
        paths.append(path)
    paths.extend(plot_roc_pr(predictions, output_dir, label=label))
    path = plot_importance(importance, output_dir, label=label)
    if path is not None:
        paths.append(path)
    return [str(path) for path in paths]

