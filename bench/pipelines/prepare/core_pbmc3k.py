"""Write core_pbmc3k.h5ad: 10x PBMC 3k raw counts with obs `sample`, for the qc task.

Kang 2018 ships no mitochondrial genes, so QC needs another library. Copies the plots
bench's pbmc3k.h5ad when present, else downloads scanpy.datasets.pbmc3k().
"""

import shutil
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

out = Path(sys.argv[1])
local = Path(__file__).resolve().parents[2] / "plots" / "cache" / "pbmc3k.h5ad"
if local.is_file():
    adata = ad.read_h5ad(local)
else:
    # scanpy downloads into datasetdir; keep it out of the cache, which acumen freezes.
    sc.settings.datasetdir = out / ".download"
    adata = sc.datasets.pbmc3k()
    shutil.rmtree(out / ".download")
adata.var_names_make_unique()
adata = ad.AnnData(X=adata.X, obs=adata.obs[[]].copy(), var=adata.var[[]].copy())
adata.obs["sample"] = "pbmc3k"
print(adata)
adata.write_h5ad(out / "core_pbmc3k.h5ad", compression="gzip")
