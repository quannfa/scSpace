"""Tests for the pipeline module."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from xenium_scspace.bundle import bundle_files
from xenium_scspace.pipeline import prepare


@pytest.fixture
def sample_dataframes():
    """Create small synthetic DataFrames for testing prepare()."""
    np.random.seed(42)
    n_cells, n_spots, n_genes = 20, 5, 10

    sc_counts = pd.DataFrame(
        np.random.poisson(1, size=(n_cells, n_genes)).astype(np.float32),
        index=[f"cell_{i}" for i in range(n_cells)],
        columns=[f"gene_{i}" for i in range(n_genes)],
    )
    sc_meta = pd.DataFrame(
        {"cell_type": ["type_A"] * 10 + ["type_B"] * 10},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    st_counts = pd.DataFrame(
        np.random.poisson(2, size=(n_spots, n_genes)).astype(np.float32),
        index=[f"spot_{i}" for i in range(n_spots)],
        columns=[f"gene_{i}" for i in range(n_genes)],
    )
    st_meta = pd.DataFrame(
        {"xcoord": np.random.rand(n_spots), "ycoord": np.random.rand(n_spots)},
        index=[f"spot_{i}" for i in range(n_spots)],
    )

    return sc_counts, sc_meta, st_counts, st_meta


class TestPrepare:
    def test_prepare_creates_bundle(self, sample_dataframes):
        sc_counts, sc_meta, st_counts, st_meta = sample_dataframes
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = prepare(sc_counts, sc_meta, st_counts, st_meta, output_dir=tmpdir)

            assert bundle.sc_data_path.exists()
            assert bundle.sc_meta_path.exists()
            assert bundle.st_data_path.exists()
            assert bundle.st_meta_path.exists()

    def test_prepare_correct_shapes(self, sample_dataframes):
        sc_counts, sc_meta, st_counts, st_meta = sample_dataframes
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = prepare(sc_counts, sc_meta, st_counts, st_meta, output_dir=tmpdir)

            sc_data = pd.read_csv(bundle.sc_data_path, index_col=0)
            # scSpace convention: data is genes × cells
            assert sc_data.shape == (sc_counts.shape[1], sc_counts.shape[0])

            st_data = pd.read_csv(bundle.st_data_path, index_col=0)
            assert st_data.shape == (st_counts.shape[1], st_counts.shape[0])