"""Purpose: write the two seeded demo inputs for sc_de.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad (mast mode), examples/demo_npx.tsv (lm mode)
Run:     python examples/make_demo.py
Notes:   demo.h5ad: 6 donors x ctrl/stim x cell types T and Mono, 80 cells each,
         200 genes. In T cells G000-G009 are 3-fold up and G010-G014 3-fold down
         in stim (var["true_de"]); Mono has no DE. X = normalize_total(1e4) + log1p
         (natural log), sparse; var["gene_ids"] holds Ensembl-like ids
         (ENSG90000000000.1, ...). demo_npx.tsv: 40 samples (20 healthy, 20 disease),
         50 proteins P01-P50, NPX (log2). P01-P08 are 1 NPX higher in disease;
         P09-P12 rise 0.05 NPX per year of age; the rest are null.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(SEED)

# Single-cell demo.
N_GENES, N_CELLS, FOLD = 200, 80, 3.0
genes = [f"G{i:03d}" for i in range(N_GENES)]
true_de = np.r_[np.full(10, "up"), np.full(5, "down"), np.full(N_GENES - 15, "")]
base_mean = rng.uniform(0.2, 3, N_GENES)
blocks, obs = [], []
for cell_type in ["Mono", "T"]:
    for donor in [f"D{i}" for i in range(1, 7)]:
        donor_shift = rng.lognormal(0, 0.2, N_GENES)
        for condition in ["ctrl", "stim"]:
            mean = base_mean * donor_shift
            if cell_type == "T" and condition == "stim":
                mean = mean * np.select(
                    [true_de == "up", true_de == "down"], [FOLD, 1 / FOLD], 1
                )
            size = rng.lognormal(0, 0.3, (N_CELLS, 1))
            blocks.append(rng.poisson(rng.gamma(2.0, mean * size / 2.0)))
            obs += [(cell_type, donor, condition)] * N_CELLS

counts = np.vstack(blocks).astype(np.float32)
lib = counts.sum(axis=1, keepdims=True)
x = np.log1p(counts / np.maximum(lib, 1) * 1e4).astype(np.float32)
obs = pd.DataFrame(obs, columns=["cell_type", "donor", "condition"])
obs.index = [f"cell{i:04d}" for i in range(len(obs))]
adata = ad.AnnData(
    X=sparse.csr_matrix(x),
    obs=obs,
    # Ensembl-like ids for --gene-id-col; var_names are the symbols.
    var=pd.DataFrame(
        {"true_de": true_de, "gene_ids": [f"ENSG9{i:010d}.1" for i in range(N_GENES)]},
        index=genes,
    ),
)
adata.write_h5ad(HERE / "demo.h5ad", compression="gzip")

# Olink NPX demo.
samples = pd.DataFrame(
    {
        "sample": [f"S{i:02d}" for i in range(1, 41)],
        "status": ["healthy"] * 20 + ["disease"] * 20,
        "age": rng.integers(25, 76, 40),
        "sex": rng.choice(["F", "M"], 40),
    }
)
rows = samples.merge(
    pd.DataFrame({"feature": [f"P{i:02d}" for i in range(1, 51)]}), how="cross"
)
number = rows["feature"].str[1:].astype(int)
rows["value"] = (
    rng.normal(5, 1, 50)[number - 1]
    + rng.normal(0, 0.5, len(rows))
    + 1.0 * ((number <= 8) & (rows["status"] == "disease"))
    + 0.05 * (rows["age"] - 50) * number.between(9, 12)
).round(3)
cols = ["sample", "feature", "value", "status", "age", "sex"]
rows[cols].sort_values(["feature", "sample"]).to_csv(
    HERE / "demo_npx.tsv", sep="\t", index=False
)
