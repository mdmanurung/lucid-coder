"""Write pbmc3k_umap.tsv and pbmc3k_markers.tsv from scanpy's processed PBMC 3k.

scanpy.datasets.pbmc3k_processed() holds the tutorial's UMAP, its louvain cell-type
labels and log1p(CP10k) expression in .raw.
  pbmc3k_umap.tsv     cell, umap_1, umap_2, cell_type (umap)
  pbmc3k_markers.tsv  cell, umap_1, umap_2, one log1p(CP10k) column per marker (marker_feature)
"""

import shutil
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

MARKERS = ["CD14", "MS4A1", "NKG7", "LYZ", "CD3E", "FCGR3A", "PPBP"]

out = Path(sys.argv[1])
# scanpy downloads into datasetdir; keep it out of the cache, which acumen freezes.
sc.settings.datasetdir = out / ".download"
adata = sc.datasets.pbmc3k_processed()
shutil.rmtree(out / ".download")
umap = pd.DataFrame(adata.obsm["X_umap"], columns=["umap_1", "umap_2"])
umap.insert(0, "cell", adata.obs_names)
umap["cell_type"] = adata.obs["louvain"].astype(str).to_numpy()
umap.to_csv(out / "pbmc3k_umap.tsv", sep="\t", index=False, float_format="%.4f")
expr = adata.raw[:, MARKERS].X
expr = pd.DataFrame(expr.toarray() if hasattr(expr, "toarray") else expr, columns=MARKERS)
markers = pd.concat([umap[["cell", "umap_1", "umap_2"]], expr], axis=1)
markers.to_csv(out / "pbmc3k_markers.tsv", sep="\t", index=False, float_format="%.4f")
print(umap.shape, markers.shape)
