"""Purpose: build the synthetic demo inputs for sc_normalize.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad (RNA), demo_adt.h5ad (ADT of the same cells),
         demo_empty.h5ad (ADT of empty droplets)
Run:     python examples/make_demo.py
Notes:   3 true groups (obs["true_group"]) x 200 cells. RNA: 300 genes, G000-G059
         mark the groups (20 each, 4x up); per-cell depth is lognormal(0, 0.7), so
         log1p keeps a depth signal that Pearson residuals remove. ADT: CD3 on g1,
         CD19 on g2, CD14 on g3, CD16 on none, three isotypes. Every ADT count is
         ambient background times a per-droplet technical factor, plus the
         specific signal on positive cells. Empty barcodes differ from cell ones.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
N_PER_GROUP = 200
GROUPS = ["g1", "g2", "g3"]
N_GENES = 300
MARKERS_PER_GROUP = 20
FOLD = 4.0
DEPTH_SIGMA = 0.7
POSITIVE = {"CD3": "g1", "CD19": "g2", "CD14": "g3", "CD16": None}
ISOTYPES = ["IgG1", "IgG2a", "IgG2b"]
N_EMPTY = 2000
SIGNAL = 400.0

here = Path(__file__).resolve().parent
rng = np.random.default_rng(SEED)
groups = np.repeat(GROUPS, N_PER_GROUP)
obs = pd.DataFrame(
    {"true_group": pd.Categorical(groups)},
    index=[f"cell{i:04d}" for i in range(len(groups))],
)

# RNA: gamma-Poisson counts, group markers, wide depth spread.
base_mean = rng.lognormal(mean=-1.0, sigma=1.0, size=N_GENES)
mean = np.tile(base_mean, (len(groups), 1))
for i, group in enumerate(GROUPS):
    mean[groups == group, i * MARKERS_PER_GROUP : (i + 1) * MARKERS_PER_GROUP] *= FOLD
depth = rng.lognormal(0, DEPTH_SIGMA, (len(groups), 1))
rna = rng.poisson(rng.gamma(5.0, mean * depth / 5.0))
ad.AnnData(
    X=sparse.csr_matrix(rna.astype(np.int32)),
    obs=obs,
    var=pd.DataFrame(index=[f"G{i:03d}" for i in range(N_GENES)]),
).write_h5ad(here / "demo.h5ad", compression="gzip")

# ADT: ambient background per protein, scaled per droplet, plus specific signal.
proteins = [*POSITIVE, *ISOTYPES]
ambient = rng.lognormal(np.log(20), 0.5, len(proteins))
var = pd.DataFrame(index=proteins)


def droplets(n):
    technical = rng.lognormal(0, 0.3, (n, 1))
    return ambient * technical


cell_mean = droplets(len(groups))
for j, group in enumerate(POSITIVE.values()):
    if group:
        cell_mean[groups == group, j] += SIGNAL * rng.lognormal(0, 0.3, N_PER_GROUP)
adt = rng.poisson(cell_mean).astype(np.int32)
ad.AnnData(X=sparse.csr_matrix(adt), obs=obs, var=var).write_h5ad(
    here / "demo_adt.h5ad", compression="gzip"
)
empty = rng.poisson(droplets(N_EMPTY)).astype(np.int32)
ad.AnnData(
    X=sparse.csr_matrix(empty),
    obs=pd.DataFrame(index=[f"empty{i:04d}" for i in range(N_EMPTY)]),
    var=var,
).write_h5ad(here / "demo_empty.h5ad", compression="gzip")
