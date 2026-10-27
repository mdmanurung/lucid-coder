"""Write the annotate task inputs from kang2018.h5ad.

core_kang_stim_ref.h5ad: stim cells, raw counts, obs cell_type (the reference).
core_kang_ctrl_query.h5ad: ctrl cells, raw counts, obs replicate and cell_type, plus a
UMAP in obsm["X_umap"] (the label-transfer tool draws on it; no answer depends on it).
core_pbmc_markers.tsv: two canonical markers per Kang cell type (label, gene).
`core_` sorts before `kang2018`, so this builds kang2018.h5ad first when it is missing.
"""

import runpy
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

out = Path(sys.argv[1])
if not (out / "kang2018.h5ad").is_file():
    runpy.run_path(str(Path(__file__).with_name("kang2018.py")), run_name="__main__")
kang = ad.read_h5ad(out / "kang2018.h5ad")

ref = kang[kang.obs["label"] == "stim"].copy()
ref.obs = ref.obs[["cell_type"]]
ref.write_h5ad(out / "core_kang_stim_ref.h5ad", compression="gzip")

query = kang[kang.obs["label"] == "ctrl"].copy()
query.obs = query.obs[["replicate", "cell_type"]]
lognorm = query.copy()
sc.pp.normalize_total(lognorm, target_sum=1e4)
sc.pp.log1p(lognorm)
sc.pp.highly_variable_genes(lognorm, n_top_genes=2000, subset=True)
sc.pp.scale(lognorm)
sc.tl.pca(lognorm, n_comps=30, random_state=0)
sc.pp.neighbors(lognorm, random_state=0)
sc.tl.umap(lognorm, random_state=0)
query.obsm["X_umap"] = lognorm.obsm["X_umap"]
query.write_h5ad(out / "core_kang_ctrl_query.h5ad", compression="gzip")

markers = {
    "B cells": ["MS4A1", "CD79A"],
    "CD14+ Monocytes": ["CD14", "LYZ"],
    "CD4 T cells": ["CD3E", "IL7R"],
    "CD8 T cells": ["CD8A", "CD8B"],
    "Dendritic cells": ["FCER1A", "CLEC10A"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
    "Megakaryocytes": ["PPBP", "PF4"],
    "NK cells": ["GNLY", "NKG7"],
}
rows = [(label, gene) for label, genes in markers.items() for gene in genes]
pd.DataFrame(rows, columns=["label", "gene"]).to_csv(
    out / "core_pbmc_markers.tsv", sep="\t", index=False
)
print(ref, query, sep="\n")
