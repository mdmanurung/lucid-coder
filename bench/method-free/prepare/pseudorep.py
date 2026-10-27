"""Write pseudorep_<split>.h5ad: a seeded synthetic pseudoreplication trap, one per split.

6 donors (3 ctrl, 3 stim), one sample per donor, ~300 cells each, one cell type, raw
counts. The first 10 genes are named candidates; the rest are background (`BG0001`...).
Candidates come in three kinds, shuffled:
- true DE: a large condition effect, small donor noise -> significant at any unit;
- donor noise: no condition effect, a big donor effect whose pattern puts the stim
  donors above the ctrl donors on average -> significant only if cells are the unit;
- null: small donor noise only, with equal ctrl and stim donor means.
Takes the cache dir as argv[1]; reads nothing from it.
"""

import os

for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[var] = "1"

import string  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# split: (seed, n true DE, n donor-noise); the other candidates are nulls.
SPLITS = {"train": (11, 3, 4), "valid": (22, 2, 5), "test": (33, 4, 3)}
N_CAND, N_BG, CELLS = 10, 200, 300
TRUE_LFC = 1.6  # natural-log condition effect of a true DE gene (~2.3 log2)
# Donor-noise pattern (natural log): within-group spread ~0.7, stim shifted up by 0.35.
NOISE_CTRL = np.array([-0.8, -0.05, 0.6])
NOISE_SHIFT = 0.35
PHI = 0.3  # cell-level gamma-Poisson overdispersion


def simulate(seed, n_true, n_noise):
    rng = np.random.default_rng(seed)
    donors = [f"D{i}" for i in rng.permutation(6) + 1]  # shuffled ids, 3 ctrl then 3 stim
    cond = np.array(["ctrl"] * 3 + ["stim"] * 3)
    stim = cond == "stim"

    letters = np.array(list(string.ascii_uppercase))
    names = set()
    while len(names) < N_CAND:
        names.add("".join(rng.choice(letters, 3)) + str(rng.integers(1, 10)))
    cand = sorted(names)
    kind = rng.permutation(["true"] * n_true + ["noise"] * n_noise + ["null"] * (N_CAND - n_true - n_noise))

    # Donor x gene log effects.
    eff = np.zeros((6, N_CAND + N_BG))
    for g, k in enumerate(kind):
        if k == "true":
            eff[:, g] = rng.normal(0, 0.1, 6) + TRUE_LFC * stim * rng.choice([-1, 1])
        elif k == "noise":
            eff[:3, g] = rng.permutation(NOISE_CTRL) + rng.normal(0, 0.05, 3)
            eff[3:, g] = rng.permutation(NOISE_CTRL) + NOISE_SHIFT + rng.normal(0, 0.05, 3)
        else:  # equal realised group means, so no sound test calls a null by chance
            eff[:, g] = rng.normal(0, 0.15, 6)
            eff[3:, g] += eff[:3, g].mean() - eff[3:, g].mean()
    eff[:, N_CAND:] = rng.normal(0, rng.uniform(0.05, 0.5, N_BG), (6, N_BG))

    mu = np.concatenate([rng.uniform(2, 6, N_CAND), rng.lognormal(0, 1, N_BG)])
    n = rng.integers(CELLS - 30, CELLS + 31, 6)
    donor_idx = np.repeat(np.arange(6), n)
    size = rng.lognormal(0, 0.25, donor_idx.size)[:, None]
    lam = mu * np.exp(eff[donor_idx]) * size
    counts = rng.poisson(rng.gamma(1 / PHI, lam * PHI))

    order = rng.permutation(donor_idx.size)
    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(np.array(donors)[donor_idx][order]),
            "condition": pd.Categorical(cond[donor_idx][order]),
            "cell_type": pd.Categorical(["Monocyte"] * donor_idx.size),
        },
        index=[f"cell{i:04d}" for i in range(donor_idx.size)],
    )
    var = pd.DataFrame(index=cand + [f"BG{i:04d}" for i in range(1, N_BG + 1)])
    print(dict(zip(cand, kind, strict=True)))
    return ad.AnnData(X=counts[order].astype(np.int64), obs=obs, var=var)


for split, (seed, n_true, n_noise) in SPLITS.items():
    simulate(seed, n_true, n_noise).write_h5ad(Path(sys.argv[1]) / f"pseudorep_{split}.h5ad")
