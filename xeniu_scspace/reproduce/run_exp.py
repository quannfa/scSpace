"""Experiment harness: run scSpace pipeline and evaluate pseudo-space quality.

Usage:
    uv run python xeniu_scspace/reproduce/run_exp.py --help
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to sys.path
PROJ = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJ))


def evaluate_pseudo_space(
    pseudo_path: Path,
    clust_path: Path,
    sc_meta_path: Path,
) -> dict:
    """Compute quality metrics for pseudo-space results."""
    metrics: dict[str, float] = {}

    pseudo = pd.read_csv(pseudo_path, index_col=0)
    clust = pd.read_csv(clust_path, index_col=0)
    sc_meta = pd.read_csv(sc_meta_path, index_col=0)

    n_cells = len(pseudo)
    n_clusters = clust["scSpace"].nunique()
    metrics["n_cells"] = n_cells
    metrics["n_clusters"] = n_clusters

    # 1. Pseudo-space dimension correlation
    #    If r ≈ 1, the two dims are degenerate (diagonal line)
    r = np.corrcoef(pseudo["pseudo_dim_0"], pseudo["pseudo_dim_1"])[0, 1]
    metrics["pseudo_dim_corr"] = round(r, 4)

    # 2. Spatial vs pseudo distance correlation (PCC, scSpace paper metric)
    if "cell_position_x" in sc_meta.columns and "cell_position_y" in sc_meta.columns:
        from scipy.stats import pearsonr
        from scipy.spatial.distance import pdist

        real_coords = sc_meta.loc[pseudo.index, ["cell_position_x", "cell_position_y"]].values
        pseudo_coords = pseudo.values

        # Subsample for speed if > 5000 cells
        if n_cells > 5000:
            idx = np.random.RandomState(42).choice(n_cells, 5000, replace=False)
            real_coords = real_coords[idx]
            pseudo_coords = pseudo_coords[idx]

        real_dist = pdist(real_coords)
        pseudo_dist = pdist(pseudo_coords)
        pcc, _ = pearsonr(real_dist, pseudo_dist)
        metrics["pcc_dist"] = round(pcc, 4)

    # 3. Cluster separation in pseudo-space
    # Handle clustering CSV: might have 'cell_id' col or be indexed by cell_id
    if "cell_id" in clust.columns:
        clust = clust.set_index("cell_id")
    df = pseudo.join(clust)
    cluster_centers = df.groupby("scSpace")[["pseudo_dim_0", "pseudo_dim_1"]].mean()
    intra_dists = []
    for cl in df["scSpace"].unique():
        sub = df[df["scSpace"] == cl][["pseudo_dim_0", "pseudo_dim_1"]].values
        center = cluster_centers.loc[cl].values
        intra_dists.append(np.mean(np.sqrt(((sub - center) ** 2).sum(axis=1))))
    metrics["cluster_intra_dist_mean"] = round(np.mean(intra_dists), 4)
    metrics["cluster_intra_dist_std"] = round(np.std(intra_dists), 4)

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Run scSpace experiment")
    parser.add_argument("--bundle-dir", type=Path, default="_l_result/processed")
    parser.add_argument("--output-dir", type=Path, default="_l_result/xenium_output")
    parser.add_argument("--epoch-num", type=int, default=2000)
    parser.add_argument("--log-epoch", type=int, default=200)
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--hidden-size", type=int, default=128)
    parser.add_argument("--activation", default="sigmoid")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--common-size", type=int, default=2)
    parser.add_argument("--no-normalize-coords", action="store_true")
    parser.add_argument("--label", default=None, help="Experiment label for logging")
    args = parser.parse_args()

    from xenium_scspace.pipeline import run_scspace_pipeline

    label = args.label or f"epoch{args.epoch_num}_lr{args.lr}_hs{args.hidden_size}_act{args.activation}_norm{not args.no_normalize_coords}"

    print(f"\n{'='*60}")
    print(f"Experiment: {label}")
    print(f"{'='*60}")
    print(f"  bundle_dir = {args.bundle_dir}")
    print(f"  epoch_num  = {args.epoch_num}")
    print(f"  lr         = {args.lr}")
    print(f"  hidden_size= {args.hidden_size}")
    print(f"  activation = {args.activation}")
    print(f"  batch_size = {args.batch_size}")
    print(f"  common_size= {args.common_size}")
    print(f"  normalize_coords = {not args.no_normalize_coords}")

    # Run pipeline
    t0 = time.time()
    result = run_scspace_pipeline(
        bundle_dir=args.bundle_dir,
        output_dir=args.output_dir,
        st_type="spot",
        n_features=2000,
        normalize=True,
        normalize_coords=not args.no_normalize_coords,
        epoch_num=args.epoch_num,
        log_epoch=args.log_epoch,
        lr=args.lr,
        hidden_size=args.hidden_size,
        activation=args.activation,
        batch_size=args.batch_size,
        common_size=args.common_size,
        dim=50,
        kernel_type="primal",
        res=0.5,
        random_seed=123,
    )
    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")

    # Evaluate
    sc_meta_path = args.bundle_dir / "xenium_sc_meta.csv"
    metrics = evaluate_pseudo_space(
        pseudo_path=args.output_dir / "pseudo_space.csv",
        clust_path=args.output_dir / "clustering.csv",
        sc_meta_path=sc_meta_path,
    )
    metrics["elapsed_s"] = round(elapsed, 1)
    metrics["label"] = label

    print(f"\n{'='*60}")
    print(f"Results for: {label}")
    print(f"{'='*60}")
    for k, v in metrics.items():
        print(f"  {k}: {v}")

    # Save metrics
    exp_dir = PROJ / "_l_result" / "experiments"
    exp_dir.mkdir(parents=True, exist_ok=True)
    report_path = exp_dir / f"{label}.json"
    with open(report_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nReport saved to {report_path}")

    # Print pass/fail
    corr = metrics["pseudo_dim_corr"]
    if abs(corr) > 0.99:
        print(f"\n❌ FAIL: pseudo_dim_corr = {corr} (degenerate diagonal)")
    else:
        print(f"\n✅ PASS: pseudo_dim_corr = {corr} (good separation)")
    if metrics.get("pcc_dist", 0) > 0.3:
        print(f"✅ PASS: pcc_dist = {metrics.get('pcc_dist', 0):.4f}")
    else:
        print(f"⚠️  pcc_dist = {metrics.get('pcc_dist', 0):.4f} (low spatial correlation)")


if __name__ == "__main__":
    main()