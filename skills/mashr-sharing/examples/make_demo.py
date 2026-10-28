"""Purpose: write four synthetic per-cell-type DE tables with planted effects.
Inputs:  none
Outputs: examples/demo_{B,CD4T,CD8T,Mono}.tsv (gene_id, gene_symbol, log2_fc, lfc_se)
Run:     python examples/make_demo.py
Notes:   800 genes. G0000-G0079 are shared: one effect, +-1.2, in all four cell
         types (each type adds noise of sd 0.1). G0080-G0159 are specific: 20 per
         cell type, in table order, +-1.5 in that type only. The rest are null.
         Observed log2_fc = true + N(0, lfc_se), lfc_se ~ U(0.2, 0.35).
         gene_id is a made-up Ensembl-like id (ENSG9..., not a real gene).
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N_GENES, N_SHARED, N_SPECIFIC = 800, 80, 20
CELL_TYPES = ["B", "CD4T", "CD8T", "Mono"]
SHARED_EFFECT, SPECIFIC_EFFECT = 1.2, 1.5

rng = np.random.default_rng(SEED)
genes = [f"G{i:04d}" for i in range(N_GENES)]
gene_ids = [f"ENSG9{i:010d}.1" for i in range(N_GENES)]
shared = SHARED_EFFECT * rng.choice([-1, 1], N_SHARED)
for k, cell_type in enumerate(CELL_TYPES):
    truth = np.zeros(N_GENES)
    truth[:N_SHARED] = shared + rng.normal(0, 0.1, N_SHARED)
    start = N_SHARED + k * N_SPECIFIC
    truth[start : start + N_SPECIFIC] = SPECIFIC_EFFECT * rng.choice(
        [-1, 1], N_SPECIFIC
    )
    se = rng.uniform(0.2, 0.35, N_GENES)
    pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_symbol": genes,
            "log2_fc": (truth + rng.normal(0, se)).round(4),
            "lfc_se": se.round(4),
        }
    ).to_csv(
        Path(__file__).resolve().parent / f"demo_{cell_type}.tsv", sep="\t", index=False
    )
