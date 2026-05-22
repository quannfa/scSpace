"""Reproduce Figure 1a-style visualization: original space vs pseudo-space.

Reads existing scSpace pipeline outputs and generates a two-panel scatter plot
colored by scSpace cluster labels, analogous to scSpace paper Figure 1a.
Uses the **improved** pipeline run (coordinate-normalised MLP training).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────
SC_META = Path("_l_result/processed/xenium_sc_meta.csv")
PSEUDO = Path("_l_result/xenium_output/pseudo_space.csv")
CLUSTERING = Path("_l_result/xenium_output/clustering.csv")
OUTPUT = Path("_l_result/figure1/fig1_original_vs_pseudo_space.png")
DPI = 150
POINT_SIZE = 3
ALPHA = 0.7

# ── Load and merge ─────────────────────────────────────────────────────────
print("Loading data...")
sc_meta = pd.read_csv(SC_META, index_col=0)
pseudo = pd.read_csv(PSEUDO, index_col=0)
clustering = pd.read_csv(CLUSTERING, index_col=0)

assert sc_meta.index.name == "cell_id"
assert pseudo.index.name == "cell_id"

df = (
    sc_meta[["cell_position_x", "cell_position_y"]]
    .join(pseudo)
    .join(clustering)
)
df["scSpace"] = df["scSpace"].astype(str)

n_total = len(df)
n_clusters = df["scSpace"].nunique()
print(f"Loaded {n_total} cells, {n_clusters} clusters.")

assert n_total == len(sc_meta), "Row count mismatch after merge"
assert n_clusters >= 2, f"Expected >=2 clusters, got {n_clusters}"

# ── Plotting ────────────────────────────────────────────────────────────────
clusters = sorted(df["scSpace"].unique(), key=int)
colors = plt.cm.tab20(np.linspace(0, 1, len(clusters)))
color_map = dict(zip(clusters, colors))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: original Xenium space
for cl in clusters:
    mask = df["scSpace"] == cl
    ax1.scatter(
        df.loc[mask, "cell_position_x"],
        df.loc[mask, "cell_position_y"],
        c=[color_map[cl]],
        s=POINT_SIZE,
        alpha=ALPHA,
        label=cl,
    )
ax1.set_title("Original space (Xenium TgCRND8)")
ax1.set_xlabel("cell_position_x")
ax1.set_ylabel("cell_position_y")
ax1.legend(markerscale=5, fontsize="small", ncol=2, loc="upper right")

# Right: pseudo space
for cl in clusters:
    mask = df["scSpace"] == cl
    ax2.scatter(
        df.loc[mask, "pseudo_dim_0"],
        df.loc[mask, "pseudo_dim_1"],
        c=[color_map[cl]],
        s=POINT_SIZE,
        alpha=ALPHA,
        label=cl,
    )
ax2.set_title("Pseudo space (scSpace)")
ax2.set_xlabel("pseudo_dim_0")
ax2.set_ylabel("pseudo_dim_1")
ax2.legend(markerscale=5, fontsize="small", ncol=2, loc="upper right")

fig.tight_layout()
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT, dpi=DPI)
plt.close(fig)

print(f"Figure saved to {OUTPUT}")