"""Evaluation and visualization functions for pseudo-space results."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def evaluate_pseudo_space(
    sc_adata: ad.AnnData,
    ground_truth_key: str = "scSpace",
) -> dict:
    """Compute clustering metrics (ARI, NMI) for pseudo-space results.

    Requires ``sc_adata.obs`` to contain a ground-truth label column
    (``ground_truth_key``) and the ``"scSpace"`` column produced by
    :func:`scSpace.scspace.spatial_cluster`.

    Parameters
    ----------
    sc_adata
        AnnData object with clustering results.
    ground_truth_key
        Column name in ``sc_adata.obs`` holding ground-truth labels.

    Returns
    -------
    Dict with keys ``"ari"``, ``"nmi"``.
    """
    if ground_truth_key not in sc_adata.obs:
        raise KeyError(f"Ground-truth column '{ground_truth_key}' not found in sc_adata.obs")
    if "scSpace" not in sc_adata.obs:
        raise KeyError("'scSpace' column not found in sc_adata.obs; run spatial_cluster first")

    true_labels = sc_adata.obs[ground_truth_key].astype(str)
    pred_labels = sc_adata.obs["scSpace"].astype(str)

    ari = adjusted_rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)

    return {"ari": round(ari, 4), "nmi": round(nmi, 4)}


def plot_pseudo_space(
    sc_adata: ad.AnnData,
    output_dir: str | Path,
) -> Path:
    """Generate pseudo-space visualisation plots.

    Saves a scatter plot of the 2-D pseudo-space coordinates and a
    clustering comparison figure.

    Parameters
    ----------
    sc_adata
        AnnData object with ``obsm["pseudo_space"]``.
    output_dir
        Directory to save figures into.

    Returns
    -------
    Path to the output directory.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pseudo = sc_adata.obsm["pseudo_space"]
    if pseudo.shape[1] < 2:
        raise ValueError(f"pseudo_space has {pseudo.shape[1]} dimensions; need at least 2")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Pseudo-space scatter (first 2 dims)
    ax = axes[0]
    scatter = ax.scatter(pseudo[:, 0], pseudo[:, 1], c=pseudo[:, 0], cmap="viridis", s=5, alpha=0.7)
    ax.set_title("Pseudo-space (first 2 dims)")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    fig.colorbar(scatter, ax=ax)

    # Clustering view
    ax = axes[1]
    if "scSpace" in sc_adata.obs:
        clusters = sc_adata.obs["scSpace"].astype(str).values
        unique = np.unique(clusters)
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique)))
        for i, cl in enumerate(unique):
            mask = clusters == cl
            ax.scatter(pseudo[mask, 0], pseudo[mask, 1], c=[colors[i]], label=cl, s=5, alpha=0.7)
        ax.legend(markerscale=3, fontsize="small", ncol=2)
    ax.set_title("Pseudo-space coloured by cluster")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")

    fig.tight_layout()
    fig_path = output_dir / "pseudo_space_plot.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    return fig_path


def benchmark_compare(
    results_dir: str | Path,
    baseline_json: str | Path,
) -> dict:
    """Compare pipeline results against a baseline JSON.

    Parameters
    ----------
    results_dir
        Directory containing pipeline output (``clustering.csv``,
        ``pseudo_space.csv``).
    baseline_json
        Path to baseline ``evaluation_metrics.json``.

    Returns
    -------
    Dict of metric deviations (current - baseline).
    """
    results_dir = Path(results_dir)
    with open(baseline_json) as f:
        baseline = json.load(f)

    current: dict = {}
    clustering_path = results_dir / "clustering.csv"
    if clustering_path.exists():
        df = pd.read_csv(clustering_path)
        n_clusters = df["scSpace"].nunique()
        current["n_clusters"] = n_clusters

    deviations = {}
    for key in baseline:
        if key in current:
            deviations[f"{key}_delta"] = round(current[key] - baseline[key], 4)

    return deviations