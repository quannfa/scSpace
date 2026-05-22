"""Tests for the bundle module."""

from __future__ import annotations

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_scspace.bundle import ScSpaceBundlePaths, bundle_files, validate_bundle


@pytest.fixture
def sample_adata():
    """Create small synthetic AnnData objects for testing."""
    n_cells, n_spots, n_genes = 20, 5, 10
    np.random.seed(42)

    sc_adata = ad.AnnData(
        X=np.random.poisson(1, size=(n_cells, n_genes)).astype(np.float32),
        obs=pd.DataFrame({"cell_type": ["type_A"] * 10 + ["type_B"] * 10}, index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
    )

    st_adata = ad.AnnData(
        X=np.random.poisson(2, size=(n_spots, n_genes)).astype(np.float32),
        obs=pd.DataFrame({"xcoord": np.random.rand(n_spots), "ycoord": np.random.rand(n_spots)}, index=[f"spot_{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
    )

    return sc_adata, st_adata


class TestBundleFiles:
    def test_bundle_files_creates_four_csvs(self, sample_adata):
        sc_adata, st_adata = sample_adata
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = bundle_files(sc_adata, st_adata, output_dir=tmpdir)

            assert bundle.sc_data_path.exists()
            assert bundle.sc_meta_path.exists()
            assert bundle.st_data_path.exists()
            assert bundle.st_meta_path.exists()

    def test_bundle_files_correct_shapes(self, sample_adata):
        sc_adata, st_adata = sample_adata
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = bundle_files(sc_adata, st_adata, output_dir=tmpdir)

            sc_data = pd.read_csv(bundle.sc_data_path, index_col=0)
            assert sc_data.shape == (sc_adata.n_vars, sc_adata.n_obs)

            sc_meta = pd.read_csv(bundle.sc_meta_path, index_col=0)
            assert sc_meta.shape == (sc_adata.n_obs, len(sc_adata.obs.columns))

    def test_bundle_custom_prefix(self, sample_adata):
        sc_adata, st_adata = sample_adata
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = bundle_files(sc_adata, st_adata, output_dir=tmpdir, prefix="test")
            assert bundle.sc_data_path.name == "test_sc_data.csv"


class TestValidateBundle:
    def test_validate_valid_bundle(self, sample_adata):
        sc_adata, st_adata = sample_adata
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = bundle_files(sc_adata, st_adata, output_dir=tmpdir)
            errors = validate_bundle(bundle)
            assert errors == []

    def test_validate_missing_file(self):
        bundle = ScSpaceBundlePaths(
            sc_data_path=Path("/nonexistent/sc_data.csv"),
            sc_meta_path=Path("/nonexistent/sc_meta.csv"),
            st_data_path=Path("/nonexistent/st_data.csv"),
            st_meta_path=Path("/nonexistent/st_meta.csv"),
            output_dir=Path("/nonexistent"),
        )
        errors = validate_bundle(bundle)
        assert len(errors) == 4  # all four files missing

    def test_validate_index_mismatch(self, sample_adata):
        sc_adata, st_adata = sample_adata
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = bundle_files(sc_adata, st_adata, output_dir=tmpdir)

            # Tamper with sc_meta to cause mismatch
            sc_meta = pd.read_csv(bundle.sc_meta_path, index_col=0)
            sc_meta = sc_meta.iloc[:-1]  # drop one row
            sc_meta.to_csv(bundle.sc_meta_path)

            errors = validate_bundle(bundle)
            assert len(errors) > 0