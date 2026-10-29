# Purpose: sum Kang 2018 raw counts per patient x condition x cell type (CD14+ Monocytes,
#          CD4 T cells) for task 14.
# Run:     uv run --python 3.12 --with anndata --with pandas \
#            python -I 01_pseudobulk.py <kang_2018.h5ad> <work_dir>
# Inputs:  kang_2018.h5ad   [https://exampledata.scverse.org/pertpy/kang_2018.h5ad,
#                            38,356,412 B, GEO GSE96583, sha256
#                            e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830]
# Outputs: <work_dir>/pseudobulk_counts.tsv.gz   genes x 32 pseudobulk samples, summed counts
#          <work_dir>/pseudobulk_samples.tsv     one row per sample: patient, label, cell type, cells
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# The two types the interaction compares.
CELL_TYPES = ["CD14+ Monocytes", "CD4 T cells"]

adata = ad.read_h5ad(sys.argv[1])
work_dir = Path(sys.argv[2])
work_dir.mkdir(parents=True, exist_ok=True)

# Sum counts per sample --------------------------------------------------------
adata = adata[adata.obs["cell_type"].isin(CELL_TYPES)]
counts = adata.X.tocsr()
if not np.allclose(counts.data, np.round(counts.data)):
    raise ValueError("X does not hold raw counts")
obs = adata.obs[["replicate", "label", "cell_type"]].astype(str)
sample_key = obs.agg("|".join, axis=1)
samples = sorted(sample_key.unique())
pseudobulk = np.vstack(
    [np.asarray(counts[(sample_key == s).to_numpy()].sum(0)).ravel() for s in samples]
)

# Write tables -----------------------------------------------------------------
pd.DataFrame(
    pseudobulk.T.astype(np.int64),
    index=pd.Index(adata.var_names, name="gene"),
    columns=samples,
).to_csv(work_dir / "pseudobulk_counts.tsv.gz", sep="\t")
sample_sheet = (
    obs.assign(sample=sample_key).value_counts().rename("n_cells").reset_index()
)
sample_sheet.sort_values("sample").to_csv(
    work_dir / "pseudobulk_samples.tsv", sep="\t", index=False
)
