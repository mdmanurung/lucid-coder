"""Purpose: build the synthetic demo.h5ad for sc_markers.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   3 clusters x 50 cells, 200 genes. Genes G000-G014 are planted markers,
         5 per cluster (var["marker_of"]), 8-fold up in their cluster. Background
         genes sit at a moderate level in every cell so no null gene has a near-zero
         rest mean. Raw counts in layers["counts"]; X = normalize_total + log1p.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SEED = 0
N_GENES = 200
N_CELLS = 50
N_MARKERS = 5
MARKER_FOLD = 8
CLUSTERS = ["c1", "c2", "c3"]

rng = np.random.default_rng(SEED)
genes = [f"G{i:03d}" for i in range(N_GENES)]
base_mean = rng.uniform(1, 4, N_GENES)
marker_of = np.full(N_GENES, "", dtype=object)
for k, cluster in enumerate(CLUSTERS):
    marker_of[k * N_MARKERS : (k + 1) * N_MARKERS] = cluster

blocks = []
for cluster in CLUSTERS:
    mean = base_mean * np.where(marker_of == cluster, MARKER_FOLD, 1)
    size = rng.lognormal(0, 0.2, (N_CELLS, 1))
    blocks.append(rng.poisson(rng.gamma(4.0, mean * size / 4.0)))

counts = sparse.csr_matrix(np.vstack(blocks).astype(np.float32))
adata = ad.AnnData(
    X=counts.copy(),
    obs=pd.DataFrame(
        {"cluster": pd.Categorical(np.repeat(CLUSTERS, N_CELLS), categories=CLUSTERS)},
        index=[f"cell{i:03d}" for i in range(len(CLUSTERS) * N_CELLS)],
    ),
    var=pd.DataFrame({"marker_of": marker_of.astype(str)}, index=genes),
)
adata.layers["counts"] = counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
