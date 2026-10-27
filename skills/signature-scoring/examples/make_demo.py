"""Purpose: build the synthetic demo.h5ad and demo_gene_sets.tsv for --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad, examples/demo_gene_sets.tsv
Run:     python examples/make_demo.py
Notes:   2 groups (healthy, disease) x 4 samples x 50 cells, 300 genes. sig_up
         (G000-G009) is 4-fold up in disease cells only; sig_null (G010-G019) is
         flat and lists GMISSING, a gene absent from the data. X = normalize_total
         + log1p of raw counts in layers["counts"]; obsm["X_umap"] is two blobs.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SEED = 0
N_GENES = 300
N_CELLS = 50
FOLD = 4
GROUPS = ["healthy", "disease"]
SAMPLES = {"healthy": ["H1", "H2", "H3", "H4"], "disease": ["D1", "D2", "D3", "D4"]}
HERE = Path(__file__).resolve().parent

rng = np.random.default_rng(SEED)
genes = [f"G{i:03d}" for i in range(N_GENES)]
base_mean = rng.uniform(0.5, 8, N_GENES)
up = np.arange(N_GENES) < 10

blocks, obs = [], []
for group in GROUPS:
    mean = base_mean * np.where(up & (group == "disease"), FOLD, 1)
    for sample in SAMPLES[group]:
        size = rng.lognormal(0, 0.2, (N_CELLS, 1))
        blocks.append(rng.poisson(rng.gamma(4.0, mean * size / 4.0)))
        obs += [(sample, group)] * N_CELLS

counts = sparse.csr_matrix(np.vstack(blocks).astype(np.float32))
obs = pd.DataFrame(obs, columns=["sample", "group"])
obs.index = [f"cell{i:03d}" for i in range(len(obs))]
obs["group"] = pd.Categorical(obs["group"], categories=GROUPS)
adata = ad.AnnData(X=counts.copy(), obs=obs, var=pd.DataFrame(index=genes))
adata.layers["counts"] = counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
shift = np.where(obs["group"] == "disease", 4.0, 0.0)
adata.obsm["X_umap"] = rng.normal(0, 1, (len(obs), 2)) + shift[:, None]
adata.write_h5ad(HERE / "demo.h5ad", compression="gzip")

sets = [("sig_up", g) for g in genes[:10]]
sets += [("sig_null", g) for g in [*genes[10:20], "GMISSING"]]
pd.DataFrame(sets, columns=["set", "gene"]).to_csv(
    HERE / "demo_gene_sets.tsv", sep="\t", index=False
)
