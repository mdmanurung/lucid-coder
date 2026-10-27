"""Write kang2018.h5ad: Kang 2018 PBMC, IFN-beta stimulated vs control, raw counts.

decoupler 2.x ships no Kang loader and figshare sits behind a bot challenge, so this
copies a local h5ad (the pertpy `kang_2018` file) and falls back to downloading it.
Kept obs: label (ctrl/stim), replicate (8 patients), cell_type (8 types).
"""

import os
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc
import scipy.sparse as sp

# ponytail: host path baked in; set LUCID_KANG_H5AD on another host.
SOURCE = os.environ.get(
    "LUCID_KANG_H5AD",
    "/exports/para-lipg-hpc/mdmanurung/.cache/cellimo-docs-data/kang_2018.h5ad",
)
URL = "https://figshare.com/ndownloader/files/34464122"

out = Path(sys.argv[1])
if Path(SOURCE).is_file():
    adata = ad.read_h5ad(SOURCE)
else:
    sc.settings.datasetdir = out / ".download"
    adata = sc.read(out / ".download" / "kang.h5ad", backup_url=URL)

cols = ["label", "replicate", "cell_type"]
adata = ad.AnnData(
    X=sp.csr_matrix(adata.X).astype("int32"),
    obs=adata.obs[cols].copy(),
    var=adata.var[[]].copy(),
)
for col in cols:
    print(col, list(adata.obs[col].cat.categories))
print(adata)
adata.write_h5ad(out / "kang2018.h5ad", compression="gzip")
