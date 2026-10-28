"""Purpose: build the synthetic demo.h5ad for pseudobulk_de.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   7 patients; P1-P6 have ctrl and stim samples, P7 has one ctrl sample
         with 5 T cells so --min-cells 10 drops it. 30 of 300 genes are truly DE
         in T cells (var["true_de"]); 40 genes are near-silent so --min-frac drops them.
         var["gene_ids"] holds Ensembl-like ids (ENSG90000000000.1, ...).
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
N_GENES = 300
N_TRUE_DE = 30
N_SILENT = 40
LOG2_FC = 1.5
CELLS = {"T": 40, "B": 30}

rng = np.random.default_rng(SEED)
genes = [f"G{i:03d}" for i in range(N_GENES)]
base_mean = rng.lognormal(mean=-0.5, sigma=1.0, size=N_GENES)
base_mean[-N_SILENT:] = 0.002
fold = np.ones(N_GENES)
fold[:N_TRUE_DE] = 2.0 ** np.where(np.arange(N_TRUE_DE) % 2 == 0, LOG2_FC, -LOG2_FC)

samples = [(f"P{p}", cond) for p in range(1, 7) for cond in ("ctrl", "stim")]
samples.append(("P7", "ctrl"))
blocks, obs = [], []
for patient, condition in samples:
    # Patient effect per gene, shared by both of a patient's samples (pairing).
    patient_rng = np.random.default_rng([SEED, int(patient[1:])])
    patient_effect = patient_rng.lognormal(0, 0.3, N_GENES)
    for cell_type, n_cells in CELLS.items():
        if patient == "P7" and cell_type == "T":
            n_cells = 5
        mean = base_mean * patient_effect
        if condition == "stim" and cell_type == "T":
            mean = mean * fold
        size = rng.lognormal(0, 0.3, (n_cells, 1))
        blocks.append(rng.poisson(rng.gamma(2.0, mean * size / 2.0)))
        obs.append(
            pd.DataFrame(
                {
                    "sample": f"{patient}_{condition}",
                    "patient": patient,
                    "condition": condition,
                    "cell_type": cell_type,
                },
                index=range(n_cells),
            )
        )

obs = pd.concat(obs, ignore_index=True)
obs.index = [f"cell{i:04d}" for i in range(len(obs))]
adata = ad.AnnData(
    X=sparse.csr_matrix(np.vstack(blocks).astype(np.int32)),
    obs=obs.astype("category"),
    # Ensembl-like ids for --gene-id-col; var_names are the symbols.
    var=pd.DataFrame(
        {
            "true_de": np.arange(N_GENES) < N_TRUE_DE,
            "gene_ids": [f"ENSG9{i:010d}.1" for i in range(N_GENES)],
        },
        index=genes,
    ),
)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
