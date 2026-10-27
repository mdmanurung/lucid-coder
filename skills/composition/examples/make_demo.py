"""Purpose: write examples/demo.h5ad, synthetic cells with a planted composition shift.
Inputs:  none
Outputs: examples/demo.h5ad (obs only: donor, condition, sample, cell_type; no genes)
Run:     python examples/make_demo.py
Notes:   8 donors (D1-D8), each sampled ctrl and stim, 300-500 cells per sample.
         Mono expands 3-fold and B shrinks 3-fold in stim (planted truth); pDC is
         rare (about 2 cells per sample), so --min-count drops it. Each donor has
         its own baseline per cell type (SD 0.2 on the log scale).
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 0
N_DONORS = 8
BASE = {"CD4 T": 0.35, "CD8 T": 0.2, "NK": 0.12, "B": 0.12, "Mono": 0.15, "DC": 0.05}
BASE["pDC"] = 0.005
STIM_FOLD = {"Mono": 3.0, "B": 1 / 3}

rng = np.random.default_rng(SEED)
labels = list(BASE)
donor_effect = rng.normal(0, 0.2, (N_DONORS, len(labels)))
cells = []
for d in range(N_DONORS):
    for condition in ["ctrl", "stim"]:
        fold = [STIM_FOLD.get(lv, 1.0) if condition == "stim" else 1.0 for lv in labels]
        weight = np.array(list(BASE.values())) * fold * np.exp(donor_effect[d])
        n = rng.multinomial(rng.integers(300, 501), weight / weight.sum())
        cells.append(
            pd.DataFrame(
                {
                    "donor": f"D{d + 1}",
                    "condition": condition,
                    "sample": f"D{d + 1}_{condition}",
                    "cell_type": np.repeat(labels, n),
                }
            )
        )
obs = pd.concat(cells, ignore_index=True).astype("category")
obs.index = [f"cell{i}" for i in range(len(obs))]
ad.AnnData(obs=obs).write_h5ad(
    Path(__file__).resolve().parent / "demo.h5ad", compression="gzip"
)
