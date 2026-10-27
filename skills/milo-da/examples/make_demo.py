"""Purpose: build the synthetic demo.h5ad for milo_da.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   6 samples (A1-A3 in condition A, B1-B3 in B), 200 cells each, 4 cell
         types as Gaussian blobs in a 10-D latent (obsm["X_pca"]). Mono is
         planted at 3x its condition-A share in B. Per-sample proportions jitter
         by Dirichlet(200 * p). X holds no genes: Milo reads only the graph.
         Each sample has its own donor (D1-D6), so --donor-col passes the check.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
N_PER_SAMPLE = 200
N_DIMS = 10
N_NEIGHBORS = 30
CELL_TYPES = ["T", "B", "NK", "Mono"]
SHARE_A = np.array([0.40, 0.25, 0.20, 0.15])
MONO_FOLD = 3.0
SAMPLES = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "B3": "B"}

rng = np.random.default_rng(SEED)
share_b = SHARE_A * np.where(np.array(CELL_TYPES) == "Mono", MONO_FOLD, 1.0)
share_b /= share_b.sum()
centres = rng.normal(0, 4, (len(CELL_TYPES), N_DIMS))
latent, obs = [], []
for i, (sample, condition) in enumerate(SAMPLES.items(), start=1):
    share = rng.dirichlet(200 * (SHARE_A if condition == "A" else share_b))
    types = rng.choice(len(CELL_TYPES), N_PER_SAMPLE, p=share)
    latent.append(centres[types] + rng.normal(0, 1, (N_PER_SAMPLE, N_DIMS)))
    obs.append(
        pd.DataFrame(
            {
                "sample": sample,
                "donor": f"D{i}",
                "condition": condition,
                "cell_type": np.array(CELL_TYPES)[types],
            }
        )
    )

obs = pd.concat(obs, ignore_index=True).astype("category")
obs.index = [f"cell{i:04d}" for i in range(len(obs))]
adata = ad.AnnData(obs=obs, obsm={"X_pca": np.vstack(latent).astype(np.float32)})
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca", random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
