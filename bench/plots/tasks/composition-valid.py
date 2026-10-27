"""Reproduce the composition valid answer; writes answer.md in the cwd."""

from pathlib import Path

import anndata as ad
import pandas as pd

REPLICATE, LABEL = "patient_1016", "stim"

obs = ad.read_h5ad("kang2018.h5ad").obs
share = pd.crosstab(
    [obs["replicate"], obs["label"]], obs["cell_type"], normalize="index"
)
Path("answer.md").write_text(share.loc[(REPLICATE, LABEL)].idxmax() + "\n")
