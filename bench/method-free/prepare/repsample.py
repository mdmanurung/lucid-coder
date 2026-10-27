"""Write repsample_<split>.h5ad: a seeded synthetic repeated-sample trap, one per split.

Paired design: each donor gives one `ctrl` and one `stim` sample, ~200 cells each, one
cell type, raw counts. The first 10 genes are named candidates; the rest are background
(`BG0001`...). Every gene has a large donor effect shared by the donor's two samples.
Candidates come in three kinds, shuffled:
- shift: a small condition effect, the same in every donor -> significant only if the
  donor is modelled (paired); the donor spread hides it from `~ condition`;
- decoy: no condition effect, but a large stim-only bump in one donor -> significant
  only if cells are the unit;
- null: no condition effect.
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

# split: (seed, n donors, n shift, n decoy); the other candidates are nulls.
SPLITS = {"train": (101, 7, 5, 2), "valid": (202, 8, 3, 2), "test": (303, 6, 4, 2)}
N_CAND, N_BG, CELLS = 10, 200, 200
DONOR_SD = 1.0  # natural-log between-donor spread of a candidate
SHIFT = 0.5  # natural-log within-donor condition effect of a shift gene
SAMPLE_SD = 0.08  # natural-log sample-level noise around donor + condition
DECOY = 1.5  # natural-log stim-only bump in one donor of a decoy gene
PHI = 0.3  # cell-level gamma-Poisson overdispersion


def simulate(seed, n_donors, n_shift, n_decoy):
    rng = np.random.default_rng(seed)
    donors = np.array([f"D{i}" for i in rng.permutation(n_donors) + 1])
    n_samp = 2 * n_donors
    s_donor = np.repeat(np.arange(n_donors), 2)  # sample -> donor index
    s_stim = np.tile([0, 1], n_donors)

    letters = np.array(list(string.ascii_uppercase))
    names = set()
    while len(names) < N_CAND:
        names.add("".join(rng.choice(letters, 3)) + str(rng.integers(1, 10)))
    cand = sorted(names)
    kind = rng.permutation(["shift"] * n_shift + ["decoy"] * n_decoy + ["null"] * (N_CAND - n_shift - n_decoy))

    # Sample x gene log effects: donor + condition + sample noise.
    sd = np.concatenate([np.full(N_CAND, DONOR_SD), rng.uniform(0.3, 1.0, N_BG)])
    eff = rng.normal(0, sd, (n_donors, sd.size))[s_donor] + rng.normal(0, SAMPLE_SD, (n_samp, sd.size))
    for g, k in enumerate(kind):
        if k == "shift":
            eff[:, g] += SHIFT * s_stim * rng.choice([-1, 1])
        elif k == "decoy":
            eff[2 * rng.integers(n_donors) + 1, g] += DECOY

    mu = np.concatenate([rng.uniform(2, 6, N_CAND), rng.lognormal(0, 1, N_BG)])
    n = rng.integers(CELLS - 20, CELLS + 21, n_samp)
    samp_idx = np.repeat(np.arange(n_samp), n)
    size = rng.lognormal(0, 0.25, samp_idx.size)[:, None]
    lam = mu * np.exp(eff[samp_idx]) * size
    counts = rng.poisson(rng.gamma(1 / PHI, lam * PHI))

    cond = np.array(["ctrl", "stim"])[s_stim]
    sample = np.array([f"{donors[d]}_{c}" for d, c in zip(s_donor, cond, strict=True)])
    order = rng.permutation(samp_idx.size)
    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(donors[s_donor][samp_idx][order]),
            "sample": pd.Categorical(sample[samp_idx][order]),
            "condition": pd.Categorical(cond[samp_idx][order]),
        },
        index=[f"cell{i:04d}" for i in range(samp_idx.size)],
    )
    var = pd.DataFrame(index=cand + [f"BG{i:04d}" for i in range(1, N_BG + 1)])
    print(dict(zip(cand, kind, strict=True)))
    return ad.AnnData(X=counts[order].astype(np.int64), obs=obs, var=var)


for split, (seed, n_donors, n_shift, n_decoy) in SPLITS.items():
    simulate(seed, n_donors, n_shift, n_decoy).write_h5ad(Path(sys.argv[1]) / f"repsample_{split}.h5ad")
