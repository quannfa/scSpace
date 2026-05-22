# xeniu-scspace

Xenium → scSpace data processing and analysis pipeline.

Convert Xenium raw output (`cells.csv.gz`, `transcripts.csv.gz`) into scRNA-seq and Visium-style representations, then run the scSpace pseudo-space construction and spatial clustering pipeline.

## Quick Start

### Setup

```bash
# From project root
uv sync
```

### Prepare bundle from Xenium data

```bash
uv run xenium-scspace prepare --from-xenium _l_data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs --output-dir _l_result/xenium_bundle
```

Or prepare from pre-built CSV files:

```bash
uv run xenium-scspace prepare \
  --sc-data sc_data.csv --sc-meta sc_meta.csv \
  --st-data st_data.csv --st-meta st_meta.csv \
  --output-dir _l_result/bundle
```

### Run scSpace pipeline

```bash
uv run xenium-scspace run --bundle-dir _l_result/xenium_bundle --output-dir _l_result/xenium_output
```

### Run tests

```bash
uv run python -m pytest xeniu_scspace/tests/ -v
```

## Package Structure

```
xeniu_scspace/
├── pyproject.toml
├── README.md
├── design.md
├── xenium_scspace/
│   ├── __init__.py
│   ├── bundle.py         # Bundle data structures & I/O
│   ├── pipeline.py       # Pipeline orchestration
│   ├── cli.py            # Command-line interface
│   └── advanced_analysis.py  # Evaluation & visualization
└── tests/
    ├── __init__.py
    ├── test_bundle.py
    └── test_pipeline.py
```

## Dependencies

- [scSpace](https://github.com/ZJUFanLab/scSpace) — pseudo-space reconstruction
- [t20260507-xenium-dataset](../t20260507_xenium_dataset/) — Xenium data conversion tools

## Parameters

See `uv run xenium-scspace prepare --help` and `uv run xenium-scspace run --help` for full parameter descriptions.