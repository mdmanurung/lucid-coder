"""Reproduce the box test answer; writes answer.md in the cwd."""

from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

CELL_TYPE, GENE = "NK cells", "CCL5"

adata = ad.read_h5ad("kang2018.h5ad")
adata = adata[adata.obs["cell_type"] == CELL_TYPE].copy()
adata.X = adata.X.astype("float64")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
df = adata.obs[["replicate", "label"]].copy()
df["expr"] = np.asarray(adata[:, GENE].X.todense()).ravel()
per_sample = df.groupby(["replicate", "label"], observed=True)["expr"].mean()
medians = per_sample.groupby(level="label", observed=True).median()
Path("answer.md").write_text(medians.idxmax() + "\n")
