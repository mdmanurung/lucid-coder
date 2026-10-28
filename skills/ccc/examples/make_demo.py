"""Purpose: build the synthetic demo.h5ad and demo_resource.tsv for ccc.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad, examples/demo_resource.tsv
Run:     python examples/make_demo.py
Notes:   3 cell types x 40 cells (sender, receiver, bystander), 2 samples
         alternating, 60 genes. Planted pairs: LIG1-REC1 and LIG2-REC2, ligands
         8-fold up in sender, receptors 8-fold up in receiver. The other 18
         resource pairs use background genes at the same level in every type.
         Raw counts in layers["counts"]; X = normalize_total + log1p.
         var["gene_ids"] holds made-up Ensembl-like ids.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SEED = 0
N_CELLS = 40
FOLD = 8
CELL_TYPES = ["sender", "receiver", "bystander"]
PLANTED = {"sender": ["LIG1", "LIG2"], "receiver": ["REC1", "REC2"]}
N_BACKGROUND = 18

examples = Path(__file__).resolve().parent
rng = np.random.default_rng(SEED)
background = [(f"BL{i:02d}", f"BR{i:02d}") for i in range(1, N_BACKGROUND + 1)]
genes = ["LIG1", "LIG2", "REC1", "REC2", *[g for pair in background for g in pair]]
genes += [f"G{i:02d}" for i in range(60 - len(genes))]
base_mean = rng.uniform(1, 4, len(genes))

blocks = []
for cell_type in CELL_TYPES:
    up = np.isin(genes, PLANTED.get(cell_type, []))
    mean = base_mean * np.where(up, FOLD, 1)
    size = rng.lognormal(0, 0.2, (N_CELLS, 1))
    blocks.append(rng.poisson(rng.gamma(4.0, mean * size / 4.0)))

counts = sparse.csr_matrix(np.vstack(blocks).astype(np.float32))
n_cells = len(CELL_TYPES) * N_CELLS
adata = ad.AnnData(
    X=counts.copy(),
    obs=pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                np.repeat(CELL_TYPES, N_CELLS), categories=CELL_TYPES
            ),
            "sample": pd.Categorical(np.tile(["s1", "s2"], n_cells // 2)),
        },
        index=[f"cell{i:03d}" for i in range(n_cells)],
    ),
    # Ensembl-like ids so --gene-id-col gene_ids has a column to read.
    var=pd.DataFrame(
        {"gene_ids": [f"ENSG9{i:010d}.1" for i in range(len(genes))]}, index=genes
    ),
)
adata.layers["counts"] = counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(examples / "demo.h5ad", compression="gzip")

pairs = [("LIG1", "REC1"), ("LIG2", "REC2"), *background]
pd.DataFrame(pairs, columns=["ligand", "receptor"]).to_csv(
    examples / "demo_resource.tsv", sep="\t", index=False
)
