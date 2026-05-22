"""Pipeline orchestration: Xenium → bundle → scSpace → results."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .bundle import ScSpaceBundlePaths, bundle_files


def prepare_from_xenium(
    xenium_sample_dir: str | Path,
    output_dir: str | Path,
    st_type: str = "spot",
    **kwargs,
) -> ScSpaceBundlePaths:
    """Convert Xenium sample output to scSpace-ready bundle."""
    import xenium.trans_sc as _trans_sc
    import xenium.trans_visium as _trans_visium

    # Patch esypro‑derived OUTPUT_ROOT with a real Path so mkdir() works
    _scratch = output_dir / ".scratch"
    _trans_sc.OUTPUT_ROOT = _scratch / "scrna"
    _trans_visium.OUTPUT_ROOT = _scratch / "visium"

    sample_dir = Path(xenium_sample_dir)

    # Build sc-like h5ad
    sc_output = _trans_sc.convert_xenium_to_scrna(sample_dir, **kwargs)
    sc_h5ad = sc_output / "xenium_scrna.h5ad"
    sc_adata = ad.read_h5ad(sc_h5ad)

    # Build st-like h5ad
    st_output = _trans_visium.convert_xenium_to_visium(sample_dir, **kwargs)
    st_h5ad = st_output / "xenium_visium.h5ad"
    st_adata = ad.read_h5ad(st_h5ad)

    # Export as bundle CSV files
    return bundle_files(sc_adata, st_adata, output_dir=output_dir, prefix="xenium")


def run_scspace_pipeline(
    bundle_dir: str | Path,
    output_dir: str | Path,
    st_type: str = "spot",
    n_features: int = 2000,
    normalize: bool = True,
    normalize_coords: bool = True,
    kernel_type: str = "primal",
    dim: int = 50,
    lamb: int = 1,
    gamma: int = 1,
    batch_size: int = 16,
    hidden_size: int = 128,
    common_size: int = 2,
    activation: str = "sigmoid",
    lr: float = 0.001,
    epoch_num: int = 1000,
    log_epoch: int = 100,
    use_neighbors: bool = True,
    Ks: int = 10,
    Kg: int = 20,
    n_comps: int = 50,
    alpha: int = 0,
    beta: int = 0,
    res: float = 0.5,
    target_num: int | None = None,
    random_seed: int = 123,
) -> dict:
    """Run the full scSpace pipeline on a bundle directory."""
    from scSpace.utils import load_data, preporcess
    from scSpace.scspace import construct_pseudo_space, spatial_cluster

    bundle_dir = Path(bundle_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = "xenium"

    # Locate bundle files
    bundle = ScSpaceBundlePaths(
        sc_data_path=bundle_dir / f"{prefix}_sc_data.csv",
        sc_meta_path=bundle_dir / f"{prefix}_sc_meta.csv",
        st_data_path=bundle_dir / f"{prefix}_st_data.csv",
        st_meta_path=bundle_dir / f"{prefix}_st_meta.csv",
        output_dir=output_dir,
    )

    # Load data via scSpace's loader
    sc_adata, st_adata = load_data(
        sc_data_path=str(bundle.sc_data_path),
        sc_meta_path=str(bundle.sc_meta_path),
        st_data_path=str(bundle.st_data_path),
        st_meta_path=str(bundle.st_meta_path),
    )

    # ── Coordinate normalisation ──────────────────────────────────────
    # scSpace MLP trains to predict raw pixel coordinates (0-5500) from
    # normalised TCA features.  The scale mismatch makes optimisation hard
    # and frequently yields degenerate solutions (both pseudo-dims equal).
    # We z-score the st spatial coordinates so that the target has unit
    # variance; the pseudo-space can later be interpreted in the normalised
    # frame.
    coords_mean = coords_std = None
    if normalize_coords:
        spatial = st_adata.obsm["spatial"]
        coords_mean = spatial.mean(axis=0)
        coords_std = spatial.std(axis=0)
        st_adata.obsm["spatial"] = (spatial - coords_mean) / coords_std
        print(f"Normalised st spatial coords: mean=({coords_mean[0]:.1f}, {coords_mean[1]:.1f}), "
              f"std=({coords_std[0]:.1f}, {coords_std[1]:.1f})")

    # Preprocess — auto‑clamp n_features to available gene count
    n_genes = sc_adata.shape[1]
    if n_features > n_genes:
        print(f"n_features clamped from {n_features} → {n_genes} (only {n_genes} genes available)")
        n_features = n_genes
    if normalize and n_genes < 200:
        print(f"Disabling normalization for small gene set ({n_genes} genes)")
        normalize = False
    sc_adata, st_adata = preporcess(
        sc_adata, st_adata, st_type=st_type,
        n_features=n_features, normalize=normalize,
    )

    # Construct pseudo-space
    sc_adata, st_adata = construct_pseudo_space(
        sc_adata, st_adata,
        kernel_type=kernel_type, dim=dim, lamb=lamb, gamma=gamma,
        batch_size=batch_size, hidden_size=hidden_size,
        common_size=common_size, activation=activation,
        lr=lr, epoch_num=epoch_num, log_epoch=log_epoch,
    )

    # Spatial clustering
    sc_adata = spatial_cluster(
        sc_adata,
        use_neighbors=use_neighbors, Ks=Ks, Kg=Kg,
        n_comps=n_comps, alpha=alpha, beta=beta,
        res=res, target_num=target_num, random_seed=random_seed,
    )

    # Save results
    pseudo_space = pd.DataFrame(
        sc_adata.obsm["pseudo_space"],
        index=sc_adata.obs_names,
        columns=[f"pseudo_dim_{i}" for i in range(sc_adata.obsm["pseudo_space"].shape[1])],
    )
    pseudo_space.to_csv(output_dir / "pseudo_space.csv")

    clustering = pd.DataFrame({"cell_id": sc_adata.obs_names, "scSpace": sc_adata.obs["scSpace"]})
    clustering.to_csv(output_dir / "clustering.csv", index=False)

    return {"sc_adata": sc_adata, "st_adata": st_adata, "output_dir": str(output_dir)}


def prepare(
    sc_counts: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_counts: pd.DataFrame,
    st_meta: pd.DataFrame,
    output_dir: str | Path,
) -> ScSpaceBundlePaths:
    """Create a bundle from pre-built DataFrames."""
    sc_adata = ad.AnnData(X=sc_counts, obs=sc_meta)
    st_adata = ad.AnnData(X=st_counts, obs=st_meta)
    return bundle_files(sc_adata, st_adata, output_dir=output_dir)