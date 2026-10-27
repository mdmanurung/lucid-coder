"""Write struct_mofa/: per-cell-type pseudobulk views of Kang 2018 for MOFA tasks.

One TSV per cell type, samples x features: counts summed per patient x condition
(16 samples), log2(CPM + 1), then the 1000 genes with the largest variance across
samples among genes with mean log2(CPM + 1) >= 3. Feature names carry the view
prefix (`b_ISG15`), so names are unique across views. meta.tsv holds `sample`,
`label` and `replicate`. Reads kang2018_pb.h5ad from the cache dir argv[1].
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

VIEWS = {
    "b": "B cells",
    "cd4t": "CD4 T cells",
    "cd8t": "CD8 T cells",
    "mono14": "CD14+ Monocytes",
    "mono16": "FCGR3A+ Monocytes",
    "nk": "NK cells",
}
N_FEATURES, MIN_MEAN = 1000, 3.0

cache = Path(sys.argv[1])
out = cache / "struct_mofa"
out.mkdir(exist_ok=True)
adata = ad.read_h5ad(cache / "kang2018_pb.h5ad")
obs = adata.obs[["sample", "label", "replicate"]].astype(str)
obs.drop_duplicates("sample").sort_values("sample").to_csv(
    out / "meta.tsv", sep="\t", index=False
)
for view, cell_type in VIEWS.items():
    in_type = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
    sample = obs["sample"].to_numpy()[in_type]
    x = adata.X[in_type]
    summed = pd.DataFrame(
        {s: np.asarray(x[sample == s].sum(axis=0)).ravel() for s in sorted(set(sample))},
        index=adata.var_names,
    ).T
    logcpm = np.log2(summed.div(summed.sum(axis=1), axis=0) * 1e6 + 1)
    logcpm = logcpm.loc[:, logcpm.mean() >= MIN_MEAN]
    top = logcpm.var().sort_values(ascending=False, kind="stable").index[:N_FEATURES]
    table = logcpm[sorted(top)].add_prefix(f"{view}_").round(6)
    table.rename_axis("sample").to_csv(out / f"{view}.tsv", sep="\t")
    print(view, table.shape)
