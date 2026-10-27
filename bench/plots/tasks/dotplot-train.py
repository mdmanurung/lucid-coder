"""Reproduce the dotplot train answer; writes answer.md in the cwd."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

CELL_TYPE, GENES = "B cells", ["MS4A1", "CD79A", "CD79B", "CD19", "TCL1A"]

adata = ad.read_h5ad("kang2018.h5ad")
sub = adata[adata.obs["cell_type"] == CELL_TYPE, GENES]
frac = pd.Series(np.asarray((sub.X > 0).mean(axis=0)).ravel(), index=GENES)
Path("answer.md").write_text(frac.idxmax() + "\n")
