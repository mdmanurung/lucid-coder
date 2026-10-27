"""Write pbmc3k.h5ad: 10x PBMC 3k raw counts from scanpy.datasets.pbmc3k()."""

import shutil
import sys
from pathlib import Path

import scanpy as sc

out = Path(sys.argv[1])
# scanpy downloads into datasetdir; keep it out of the cache, which acumen freezes.
sc.settings.datasetdir = out / ".download"
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
print(adata)
adata.write_h5ad(out / "pbmc3k.h5ad", compression="gzip")
shutil.rmtree(out / ".download")
