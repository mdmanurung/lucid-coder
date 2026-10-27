"""Write kang2018_pb.h5ad: kang2018.h5ad plus a `sample` column (replicate x label).

Runs after kang2018.py (prepare scripts run in filename order) in the cache dir.
"""

import sys
from pathlib import Path

import anndata as ad

cache = Path(sys.argv[1])
adata = ad.read_h5ad(cache / "kang2018.h5ad")
adata.obs["sample"] = (
    adata.obs["replicate"].astype(str) + "_" + adata.obs["label"].astype(str)
).astype("category")
print(adata.obs["sample"].nunique(), "samples")
adata.write_h5ad(cache / "kang2018_pb.h5ad", compression="gzip")
