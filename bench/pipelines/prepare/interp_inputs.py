"""Write the interp-group inputs into the cache (argv[1]); runs after kang2018*.py.

- kang2018_lognorm.h5ad: kang2018_pb.h5ad with X = log1p(normalize_total(target_sum=1e4)).
- collectri_human.tsv: CollecTRI source, target, weight. Copied from a local decoupler
  export when present, else downloaded once with decoupler.
- liana_consensus.tsv: liana's bundled consensus ligand-receptor resource (offline).
"""

import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

cache = Path(sys.argv[1])

adata = ad.read_h5ad(cache / "kang2018_pb.h5ad")
adata.X = adata.X.astype("float32")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.uns.pop("log1p", None)
adata.write_h5ad(cache / "kang2018_lognorm.h5ad", compression="gzip")
print(adata)

# ponytail: host path baked in; set LUCID_COLLECTRI_CSV on another host.
local = Path(
    os.environ.get(
        "LUCID_COLLECTRI_CSV",
        "/home/mdmanurung/.cache/decoupler/collectri_human.csv",
    )
)
if local.is_file():
    net = pd.read_csv(local)
else:
    import decoupler as dc

    net = dc.op.collectri(organism="human", remove_complexes=False, license="academic")
net = net[["source", "target", "weight"]].drop_duplicates(["source", "target"])
net.to_csv(cache / "collectri_human.tsv", sep="\t", index=False)
print(len(net), "CollecTRI edges,", net["source"].nunique(), "sources")

import liana as li  # noqa: E402

res = li.rs.select_resource("consensus")
res.to_csv(cache / "liana_consensus.tsv", sep="\t", index=False)
print(len(res), "consensus pairs")
