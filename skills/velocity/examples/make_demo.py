"""Purpose: build the synthetic demo.h5ad for velocity.py --demo.
Inputs:  none (seeded scvelo simulation, offline)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   scv.datasets.simulation: 500 cells on one linear time axis
         (obs["true_t"]), 50 genes with splicing kinetics. 4 stages cut at
         quartiles of true_t (s1..s4). X and the spliced/unspliced layers hold
         the simulated abundances; X_umap comes from a log-normalised copy.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

SEED = 0
N_CELLS = 500
N_GENES = 50
STAGES = ["s1", "s2", "s3", "s4"]

adata = scv.datasets.simulation(
    n_obs=N_CELLS,
    n_vars=N_GENES,
    t_max=25,
    alpha=5,
    beta=0.3,
    gamma=0.5,
    noise_level=1,
    random_seed=SEED,
)
adata.obs_names = [f"cell{i:04d}" for i in range(N_CELLS)]
adata.var_names = [f"G{i:02d}" for i in range(N_GENES)]
adata.obs["stage"] = pd.Categorical(
    pd.qcut(adata.obs["true_t"], 4, labels=STAGES), categories=STAGES
)
adata.var = adata.var[[]]

# UMAP for plotting only; velocity.py recomputes PCA and neighbours itself.
embed = adata.copy()
sc.pp.normalize_total(embed)
sc.pp.log1p(embed)
sc.pp.pca(embed, n_comps=10, random_state=SEED)
sc.pp.neighbors(embed, n_neighbors=30, random_state=SEED)
sc.tl.umap(embed, random_state=SEED)
adata.obsm["X_umap"] = embed.obsm["X_umap"].astype(np.float32)

for key in ("spliced", "unspliced"):
    adata.layers[key] = adata.layers[key].astype(np.float32)
adata.X = adata.layers["spliced"].copy()
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
