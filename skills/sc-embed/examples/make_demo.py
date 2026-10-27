"""Purpose: build the synthetic demo.h5ad for sc_embed.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   2 batches x 3 true groups (obs["true_group"]), 100 cells each, 300 genes.
         Genes G000-G059 mark the groups (20 each, 4x up). Batch b2 multiplies
         every gene by its own lognormal(0, 0.8) factor, so the shift survives a
         per-batch HVG choice and, without integration, clusters split by batch.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
N_GENES = 300
N_PER_CELL_GROUP = 100
GROUPS = ["g1", "g2", "g3"]
BATCHES = ["b1", "b2"]
MARKERS_PER_GROUP = 20
BATCH_SIGMA = 0.8
FOLD = 4.0

rng = np.random.default_rng(SEED)
base_mean = rng.lognormal(mean=0.0, sigma=0.8, size=N_GENES)
batch_shift = rng.lognormal(mean=0.0, sigma=BATCH_SIGMA, size=N_GENES)
blocks, obs = [], []
for batch in BATCHES:
    for i, group in enumerate(GROUPS):
        mean = base_mean.copy()
        mean[i * MARKERS_PER_GROUP : (i + 1) * MARKERS_PER_GROUP] *= FOLD
        if batch == "b2":
            mean = mean * batch_shift
        size = rng.lognormal(0, 0.3, (N_PER_CELL_GROUP, 1))
        blocks.append(rng.poisson(rng.gamma(2.0, mean * size / 2.0)))
        obs.append(
            pd.DataFrame({"batch": [batch] * N_PER_CELL_GROUP, "true_group": group})
        )

obs = pd.concat(obs, ignore_index=True)
obs.index = [f"cell{i:04d}" for i in range(len(obs))]
adata = ad.AnnData(
    X=sparse.csr_matrix(np.vstack(blocks).astype(np.int32)),
    obs=obs.astype("category"),
    var=pd.DataFrame(index=[f"G{i:03d}" for i in range(N_GENES)]),
)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
