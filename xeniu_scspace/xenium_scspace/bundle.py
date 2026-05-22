"""Bundle data structures and I/O for scSpace pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


@dataclass
class ScSpaceBundlePaths:
    """Paths for the four CSV files required by scSpace."""

    sc_data_path: Path
    sc_meta_path: Path
    st_data_path: Path
    st_meta_path: Path
    output_dir: Path


def bundle_files(
    sc_adata: ad.AnnData,
    st_adata: ad.AnnData,
    output_dir: str | Path,
    prefix: str = "xenium",
) -> ScSpaceBundlePaths:
    """Write sc- and st-like AnnData objects to four CSV files.

    Parameters
    ----------
    sc_adata
        Single-cell-style AnnData (cells × genes).
    st_adata
        Spatial-transcriptomics-style AnnData (spots × genes).
    output_dir
        Directory to write the CSV files into.
    prefix
        Filename prefix for the four output files.

    Returns
    -------
    ScSpaceBundlePaths with paths to the written files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # sc data: genes × cells (scSpace convention)
    sc_data = pd.DataFrame(
        sc_adata.X.T.toarray() if hasattr(sc_adata.X, "toarray") else np.asarray(sc_adata.X.T),
        index=sc_adata.var_names,
        columns=sc_adata.obs_names,
    )
    sc_meta = sc_adata.obs.copy()

    # st data: genes × spots (scSpace convention)
    st_data = pd.DataFrame(
        st_adata.X.T.toarray() if hasattr(st_adata.X, "toarray") else np.asarray(st_adata.X.T),
        index=st_adata.var_names,
        columns=st_adata.obs_names,
    )
    st_meta = st_adata.obs.copy()

    paths = ScSpaceBundlePaths(
        sc_data_path=output_dir / f"{prefix}_sc_data.csv",
        sc_meta_path=output_dir / f"{prefix}_sc_meta.csv",
        st_data_path=output_dir / f"{prefix}_st_data.csv",
        st_meta_path=output_dir / f"{prefix}_st_meta.csv",
        output_dir=output_dir,
    )

    sc_data.to_csv(paths.sc_data_path)
    sc_meta.to_csv(paths.sc_meta_path)
    st_data.to_csv(paths.st_data_path)
    st_meta.to_csv(paths.st_meta_path)

    return paths


def validate_bundle(bundle: ScSpaceBundlePaths) -> list[str]:
    """Validate a bundle: check file existence, shape alignment, index matching.

    Parameters
    ----------
    bundle
        Bundle paths to validate.

    Returns
    -------
    List of error messages (empty if valid).
    """
    errors: list[str] = []

    required = [
        ("sc_data", bundle.sc_data_path),
        ("sc_meta", bundle.sc_meta_path),
        ("st_data", bundle.st_data_path),
        ("st_meta", bundle.st_meta_path),
    ]

    for name, path in required:
        if not path.exists():
            errors.append(f"{name} file not found: {path}")

    if errors:
        return errors

    sc_data = pd.read_csv(bundle.sc_data_path, index_col=0)
    sc_meta = pd.read_csv(bundle.sc_meta_path, index_col=0)
    st_data = pd.read_csv(bundle.st_data_path, index_col=0)
    st_meta = pd.read_csv(bundle.st_meta_path, index_col=0)

    # sc: data columns (cells) should match meta index
    if set(sc_data.columns) != set(sc_meta.index):
        errors.append(
            f"sc_data columns ({len(sc_data.columns)}) do not match "
            f"sc_meta index ({len(sc_meta.index)})"
        )

    # st: data columns (spots) should match meta index
    if set(st_data.columns) != set(st_meta.index):
        errors.append(
            f"st_data columns ({len(st_data.columns)}) do not match "
            f"st_meta index ({len(st_meta.index)})"
        )

    # sc data rows (genes) should match sc_data.columns length in meta
    sc_genes_in_meta = "n_genes" in sc_meta.columns
    if sc_genes_in_meta:
        expected = sc_meta["n_genes"].iloc[0]
        if len(sc_data.index) != expected:
            errors.append(f"sc_data n_genes ({len(sc_data.index)}) != meta ({expected})")

    return errors