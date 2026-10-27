"""Purpose: build the synthetic demo for sc_annotate.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo_{reference,query}.h5ad, examples/demo_markers.tsv
Run:     python examples/make_demo.py
Notes:   4 cell types, 2 marker genes each plus 150 background genes; raw counts in X.
         The query has obs true_label (the truth for the golden accuracy) and
         cluster (each type split in two, 5% of cells in another type's cluster),
         and a seeded 2D layout in obsm["X_umap"]. The query has half the depth.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
MARKERS = {
    "B cell": ["MS4A1", "CD79A"],
    "Monocyte": ["CD14", "LYZ"],
    "NK cell": ["NKG7", "GNLY"],
    "T cell": ["CD3E", "CD3D"],
}
N_BACKGROUND = 150
N_REF, N_QUERY = 75, 100  # cells per type
SWAP_FRAC = 0.05

here = Path(__file__).resolve().parent
rng = np.random.default_rng(SEED)
types = list(MARKERS)
markers = [g for genes in MARKERS.values() for g in genes]
genes = markers + [f"G{i:03d}" for i in range(N_BACKGROUND)]

base = rng.lognormal(-1.0, 1.0, len(genes))
base[: len(markers)] = 0.05
profile = {}
for i, cell_type in enumerate(types):
    mean = base * rng.lognormal(0, 0.6, len(genes))
    mean[2 * i : 2 * i + 2] = 3.0
    profile[cell_type] = mean
# NK and T cells share half their programme, so a few cells are ambiguous.
profile["NK cell"][len(markers) :] = profile["T cell"][len(markers) :] * rng.lognormal(
    0, 0.2, N_BACKGROUND
)


def simulate(n_per_type, depth):
    labels = np.repeat(types, n_per_type)
    size = rng.lognormal(np.log(depth), 0.3, (len(labels), 1))
    mean = np.vstack([profile[lab] for lab in labels]) * size
    counts = rng.poisson(rng.gamma(2.0, mean / 2.0)).astype(np.int32)
    return sparse.csr_matrix(counts), labels


x_ref, y_ref = simulate(N_REF, 1.0)
reference = ad.AnnData(
    X=x_ref,
    obs=pd.DataFrame(
        {"cell_type": pd.Categorical(y_ref)},
        index=[f"ref{i:04d}" for i in range(len(y_ref))],
    ),
    var=pd.DataFrame(index=genes),
)

x_q, y_q = simulate(N_QUERY, 0.5)
type_index = np.searchsorted(types, y_q)
cluster = 2 * type_index + rng.integers(0, 2, len(y_q))
swap = rng.random(len(y_q)) < SWAP_FRAC
cluster[swap] = (cluster[swap] + 2 * rng.integers(1, len(types), swap.sum())) % (
    2 * len(types)
)
centres = rng.normal(0, 6, (2 * len(types), 2))
umap = centres[cluster] + rng.normal(0, 1, (len(y_q), 2))
query = ad.AnnData(
    X=x_q,
    obs=pd.DataFrame(
        {
            "true_label": pd.Categorical(y_q),
            "cluster": pd.Categorical(cluster.astype(str)),
        },
        index=[f"cell{i:04d}" for i in range(len(y_q))],
    ),
    var=pd.DataFrame(index=genes),
    obsm={"X_umap": umap},
)

reference.write_h5ad(here / "demo_reference.h5ad", compression="gzip")
query.write_h5ad(here / "demo_query.h5ad", compression="gzip")
pd.DataFrame(
    [(lab, g) for lab, gs in MARKERS.items() for g in gs], columns=["label", "gene"]
).to_csv(here / "demo_markers.tsv", sep="\t", index=False)
