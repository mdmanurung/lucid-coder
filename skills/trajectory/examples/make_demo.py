"""Purpose: build the synthetic demo.h5ad for trajectory.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   One linear trajectory: 600 cells with a known time t in [0, 1]
         (obs["true_time"]), 4 stages cut at quartiles of t, 40 genes that rise,
         fall or peak along t. Log-normalised values in X; PCA, neighbours and
         UMAP are computed here as sc-embed would, so the demo is "embedded".
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
N_CELLS = 600
N_GENES = 40
STAGES = ["progenitor", "early", "mid", "late"]

rng = np.random.default_rng(SEED)
true_time = np.sort(rng.uniform(0, 1, N_CELLS))
# Each gene follows one of three shapes along t, with its own centre.
centres = rng.uniform(0.1, 0.9, N_GENES)
shape = np.arange(N_GENES) % 3
logistic = 1 / (1 + np.exp(-12 * (true_time[:, None] - centres)))
profiles = np.select(
    [shape == 0, shape == 1],
    [logistic, 1 - logistic],
    np.exp(-((true_time[:, None] - centres) ** 2) / 0.02),
)
counts = rng.poisson(10 * profiles + 0.5)

# Shuffle so obs order carries no time information; the root rule takes the
# first root-group cell by index.
order = rng.permutation(N_CELLS)
adata = ad.AnnData(
    X=counts[order].astype(np.float32),
    obs=pd.DataFrame(
        {
            "true_time": true_time[order],
            "stage": pd.Categorical(
                pd.qcut(true_time[order], 4, labels=STAGES), categories=STAGES
            ),
        },
        index=[f"cell{i:04d}" for i in range(N_CELLS)],
    ),
    var=pd.DataFrame(index=[f"G{i:02d}" for i in range(N_GENES)]),
)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=10, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
adata.X = adata.X.astype(np.float32)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
