"""Reproduce the heatmap test answer; writes answer.md in the cwd."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ANCHOR, N_GENES, LIMIT = "CD4 T cells", 25, 2
TYPES = [
    "B cells",
    "CD14+ Monocytes",
    "CD4 T cells",
    "CD8 T cells",
    "FCGR3A+ Monocytes",
    "NK cells",
]

adata = ad.read_h5ad("kang2018.h5ad")
adata = adata[adata.obs["cell_type"].isin(TYPES)]
key = adata.obs["cell_type"].astype(str) + "|" + adata.obs["label"].astype(str)
pb = pd.DataFrame(adata.X.toarray(), columns=adata.var_names, dtype="float64")
pb = pb.groupby(key.values).sum()
log_cpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)
lfc = pd.DataFrame(
    {t: log_cpm.loc[f"{t}|stim"] - log_cpm.loc[f"{t}|ctrl"] for t in TYPES}
)
total = pb.loc[f"{ANCHOR}|stim"] + pb.loc[f"{ANCHOR}|ctrl"]
top = lfc.loc[total >= 100, ANCHOR].abs().nlargest(N_GENES).index
centred = lfc.loc[top].sub(lfc.loc[top].mean(axis=1), axis=0)
n_sat = int((centred.abs() > LIMIT).to_numpy().sum())
Path("answer.md").write_text(f"{n_sat}\n")
