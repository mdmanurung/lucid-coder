"""Purpose: build the synthetic demo files for tf_activity.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad, examples/demo_network.tsv, examples/demo_de.tsv
Run:     python examples/make_demo.py
Notes:   3 groups x 60 cells, 300 genes. TF1-TF6 have 15 targets each (12 up,
         3 down). TF1 is planted active in g1, TF2 in g2: their up targets are
         4-fold higher there, their down targets 4-fold lower. TF7 has 3 targets
         in the data and TF8 none, so both fall under --tmin 5.
         demo_de.tsv is a Welch t statistic per gene, g1 cells against g2 cells.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

SEED = 0
N_GENES = 300
N_CELLS = 60
N_TARGETS = 15
N_DOWN = 3
FOLD = 4
GROUPS = ["g1", "g2", "g3"]
PLANTED = {"TF1": "g1", "TF2": "g2"}
OUT = Path(__file__).resolve().parent

rng = np.random.default_rng(SEED)
genes = np.array([f"G{i:03d}" for i in range(N_GENES)])

edges = []
for k in range(6):
    targets = genes[k * N_TARGETS : (k + 1) * N_TARGETS]
    weights = np.r_[np.ones(N_TARGETS - N_DOWN), -np.ones(N_DOWN)]
    edges.append(
        pd.DataFrame({"source": f"TF{k + 1}", "target": targets, "weight": weights})
    )
edges.append(pd.DataFrame({"source": "TF7", "target": genes[200:203], "weight": 1.0}))
edges.append(
    pd.DataFrame(
        {"source": "TF8", "target": [f"X{i}" for i in range(6)], "weight": 1.0}
    )
)
net = pd.concat(edges, ignore_index=True)

base_mean = rng.uniform(2, 6, N_GENES)
blocks = []
for group in GROUPS:
    mean = base_mean.copy()
    for tf, home in PLANTED.items():
        if home != group:
            continue
        sub = net[net["source"] == tf]
        idx = np.searchsorted(genes, sub["target"])
        mean[idx] *= np.where(sub["weight"] > 0, FOLD, 1 / FOLD)
    size = rng.lognormal(0, 0.2, (N_CELLS, 1))
    blocks.append(rng.poisson(rng.gamma(4.0, mean * size / 4.0)))
counts = np.vstack(blocks).astype(np.float32)
lognorm = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)

adata = ad.AnnData(
    X=lognorm.astype(np.float32),
    obs=pd.DataFrame(
        {"group": pd.Categorical(np.repeat(GROUPS, N_CELLS), categories=GROUPS)},
        index=[f"cell{i:03d}" for i in range(len(GROUPS) * N_CELLS)],
    ),
    var=pd.DataFrame(index=genes),
)
adata.write_h5ad(OUT / "demo.h5ad", compression="gzip")
net.to_csv(OUT / "demo_network.tsv", sep="\t", index=False)

g1, g2 = lognorm[:N_CELLS], lognorm[N_CELLS : 2 * N_CELLS]
de = pd.DataFrame(
    {
        "gene": genes,
        "log2_fc": np.log2(np.expm1(g1).mean(0) + 1)
        - np.log2(np.expm1(g2).mean(0) + 1),
        "stat": stats.ttest_ind(g1, g2, equal_var=False).statistic,
    }
)
de.round(4).to_csv(OUT / "demo_de.tsv", sep="\t", index=False)
