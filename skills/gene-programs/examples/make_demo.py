"""Purpose: build the synthetic demo.h5ad for gene_programs.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   3 groups x 60 cells, 300 genes, raw counts in X. Three planted programs:
         genes G000-G009 belong to P1, G010-G019 to P2, G020-G029 to P3
         (var["program_of"]). Each cell mixes the programs with Dirichlet weights
         biased to its group's program, so group g1 mostly uses P1, and so on.
         var["gene_ids"] holds made-up Ensembl-like ids.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
N_GENES = 300
N_CELLS = 60
N_PLANTED = 10
PLANTED_FOLD = 20
GROUPS = ["g1", "g2", "g3"]
# Dirichlet concentration: own program 8, the others 1, so cells are mixtures.
OWN_ALPHA, OTHER_ALPHA = 8.0, 1.0

rng = np.random.default_rng(SEED)
genes = [f"G{i:03d}" for i in range(N_GENES)]
background = rng.uniform(0.5, 2, N_GENES)
program_of = np.full(N_GENES, "", dtype=object)
programs = np.tile(background, (len(GROUPS), 1))
for k in range(len(GROUPS)):
    planted = slice(k * N_PLANTED, (k + 1) * N_PLANTED)
    program_of[planted] = f"P{k + 1}"
    programs[k, planted] *= PLANTED_FOLD
programs /= programs.sum(axis=1, keepdims=True)

blocks = []
for k in range(len(GROUPS)):
    alpha = np.full(len(GROUPS), OTHER_ALPHA)
    alpha[k] = OWN_ALPHA
    usage = rng.dirichlet(alpha, N_CELLS)
    library = rng.lognormal(np.log(2000), 0.3, (N_CELLS, 1))
    blocks.append(rng.poisson(library * usage @ programs))

counts = sparse.csr_matrix(np.vstack(blocks).astype(np.float32))
adata = ad.AnnData(
    X=counts,
    obs=pd.DataFrame(
        {"group": pd.Categorical(np.repeat(GROUPS, N_CELLS), categories=GROUPS)},
        index=[f"cell{i:03d}" for i in range(len(GROUPS) * N_CELLS)],
    ),
    var=pd.DataFrame(
        {
            "program_of": program_of.astype(str),
            # Made-up Ensembl-like ids, so --gene-id-col gene_ids has a column.
            "gene_ids": [f"ENSG9{i:010d}.1" for i in range(N_GENES)],
        },
        index=genes,
    ),
)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
