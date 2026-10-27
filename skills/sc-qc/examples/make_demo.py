"""Purpose: build the synthetic demo.h5ad for sc_qc.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo.h5ad
Run:     python examples/make_demo.py
Notes:   3 samples with different library sizes, so per-sample MAD thresholds differ.
         Each sample has 300 good cells from 5 populations with disjoint marker
         blocks, 20 low-quality cells (few counts, ~40% MT-) and 25 doublets (the
         sum of two good-like cells from different populations). obs["planted"]
         holds the truth: good, low_quality or doublet.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 0
N_MITO, N_RIBO, N_MARKER, N_BACKGROUND = 13, 20, 40, 300
N_POPULATIONS = 5
CELLS = {"good": 300, "low_quality": 20, "doublet": 25}
LIBRARY = {"S1": 2500, "S2": 1500, "S3": 4000}  # median UMIs of a good cell

rng = np.random.default_rng(SEED)
genes = (
    [f"MT-G{i:02d}" for i in range(N_MITO)]
    + [f"RP{'SL'[i % 2]}{i:02d}" for i in range(N_RIBO)]
    + [f"MK{p}_{i:02d}" for p in range(N_POPULATIONS) for i in range(N_MARKER)]
    + [f"G{i:03d}" for i in range(N_BACKGROUND)]
)
is_mito = np.char.startswith(genes, "MT-")
marker_start = N_MITO + N_RIBO
base = rng.lognormal(0, 1, len(genes))
profiles = []
for p in range(N_POPULATIONS):
    profile = base.copy()
    block = slice(marker_start + p * N_MARKER, marker_start + (p + 1) * N_MARKER)
    profile[block] *= 10
    profiles.append(profile)


def draw(n, population, library, mito_frac):
    """Poisson counts for n cells: library size per cell, fixed mito fraction."""
    profile = profiles[population]
    nuclear = np.where(is_mito, 0, profile)
    mito = np.where(is_mito, profile, 0)
    probs = (1 - mito_frac) * nuclear / nuclear.sum() + mito_frac * mito / mito.sum()
    size = library * rng.lognormal(0, 0.3, (n, 1))
    return rng.poisson(size * probs)


blocks, obs = [], []
for sample, library in LIBRARY.items():
    population = rng.integers(0, N_POPULATIONS, CELLS["good"])
    good = [draw(1, p, library, 0.03) for p in population]
    low = [
        draw(1, p, library * 0.12, 0.4)
        for p in rng.integers(0, N_POPULATIONS, CELLS["low_quality"])
    ]
    # Heterotypic doublets: two parents from different populations, counts summed.
    first = rng.integers(0, N_POPULATIONS, CELLS["doublet"])
    second = (first + rng.integers(1, N_POPULATIONS, CELLS["doublet"])) % N_POPULATIONS
    doublets = [
        draw(1, a, library, 0.03) + draw(1, b, library, 0.03)
        for a, b in zip(first, second, strict=True)
    ]
    cells = good + low + doublets
    blocks.append(np.vstack(cells))
    planted = [kind for kind, n in CELLS.items() for _ in range(n)]
    obs.append(pd.DataFrame({"sample": sample, "planted": planted}))

order = rng.permutation(sum(CELLS.values()) * len(LIBRARY))
obs = pd.concat(obs, ignore_index=True).iloc[order].reset_index(drop=True)
counts = np.vstack(blocks)[order]
obs.index = [f"cell{i:04d}" for i in range(len(obs))]
adata = ad.AnnData(
    X=sparse.csr_matrix(counts.astype(np.int32)),
    obs=obs.astype("category"),
    var=pd.DataFrame(index=genes),
)
adata.write_h5ad(Path(__file__).resolve().parent / "demo.h5ad", compression="gzip")
