# Naive answer: concatenates X from both files and treats it as counts.
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

a = ad.read_h5ad("data/lab_a.h5ad")
b = ad.read_h5ad("data/lab_b.h5ad")
for x in (a, b):
    x.obs["cell_id"] = x.obs_names
adata = ad.concat({"lab_a": a, "lab_b": b}, label="lab", index_unique="-")
adata.layers["counts"] = adata.X.copy()
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(out / "combined.h5ad")

res = adata.obs[["lab", "cell_id", "total_counts", "n_genes_by_counts"]]
res.columns = ["lab", "cell_id", "library_size", "n_genes"]
res.to_csv(out / "cells.tsv", sep="\t", index=False)
