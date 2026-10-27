"""Write struct_traj.h5ad: a seeded synthetic branched trajectory, embedded.

A stem splits into two branches of unequal length. Cells sit on the tree at a known
`true_time`; 60 genes follow smooth curves of the latent position, with Poisson
counts. Seven groups along the tree carry neutral labels (`g1`..`g7`, shuffled), so
the topology cannot be guessed from the names. X is log-normalised; PCA, neighbours
and UMAP are stored as an embedding step would leave them. Obs order is shuffled.
Takes the cache dir as argv[1]; reads nothing from it.
"""

import os

for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[var] = "1"

import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

SEED = 7
rng = np.random.default_rng(SEED)

# Segments: (name, start time, end time, n cells). Stem 0-1, short branch 1-2,
# long branch 1-3. Groups cut each segment into equal-time pieces.
SEGMENTS = [("stem", 0.0, 1.0, 500), ("short", 1.0, 2.0, 400), ("long", 1.0, 3.0, 700)]
PIECES = {"stem": 2, "short": 2, "long": 3}
N_GENES = 60

seg, time = [], []
for name, t0, t1, n in SEGMENTS:
    seg += [name] * n
    time.append(rng.uniform(t0, t1, n))
seg, time = np.array(seg), np.concatenate(time)

# Latent position: stem runs along axis 0, each branch turns into its own axis.
latent = np.zeros((len(time), 3))
latent[:, 0] = np.minimum(time, 1.0)
latent[seg == "short", 1] = time[seg == "short"] - 1
latent[seg == "long", 2] = time[seg == "long"] - 1

# Each gene is a Gaussian bump centred somewhere on the tree.
centres = latent[rng.choice(len(time), N_GENES, replace=False)]
d2 = ((latent[:, None, :] - centres[None]) ** 2).sum(axis=2)
counts = rng.poisson(12 * np.exp(-d2 / 0.15) + 0.3)

group = np.empty(len(time), dtype=object)
for name, t0, t1, _ in SEGMENTS:
    on = seg == name
    k = PIECES[name]
    piece = np.minimum(((time[on] - t0) / (t1 - t0) * k).astype(int), k - 1)
    group[on] = [f"{name}{i}" for i in piece]
names = sorted(set(group))
label = dict(zip(names, [f"g{i + 1}" for i in rng.permutation(len(names))], strict=True))

order = rng.permutation(len(time))
adata = ad.AnnData(
    X=counts[order].astype(np.float32),
    obs=pd.DataFrame(
        {
            "true_time": time[order],
            "cluster": pd.Categorical([label[g] for g in group[order]]),
        },
        index=[f"cell{i:04d}" for i in range(len(time))],
    ),
    var=pd.DataFrame(index=[f"G{i:02d}" for i in range(N_GENES)]),
)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=15, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
print({v: k for k, v in label.items()})
adata.write_h5ad(Path(sys.argv[1]) / "struct_traj.h5ad", compression="gzip")
