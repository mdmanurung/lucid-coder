"""Purpose: build the synthetic pseudobulk demo for variance_partition.py --demo.
Inputs:  none (seeded simulation)
Outputs: examples/demo_counts.tsv (genes x samples), examples/demo_meta.tsv
Run:     python examples/make_demo.py
Notes:   12 donors x ctrl/stim = 24 samples, 120 genes, negative binomial counts.
         G000-G019 are donor-driven (donor log-SD 1.5), G020-G039 condition-driven
         (log2FC +-2), G040-G059 pure noise. Background genes carry a donor log-SD
         of 0.4 and a small condition effect, so donor has the largest median
         fraction. rin is a per-sample covariate with no effect on expression.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N_DONORS = 12
N_GENES = 120
N_PLANTED = 20
DISPERSION = 0.05

rng = np.random.default_rng(SEED)
out = Path(__file__).resolve().parent
donors = [f"D{i:02d}" for i in range(1, N_DONORS + 1)]
meta = pd.DataFrame(
    [(f"{d}_{c}", d, c) for d in donors for c in ("ctrl", "stim")],
    columns=["sample", "donor", "condition"],
)
meta["rin"] = rng.uniform(6, 9, len(meta)).round(2)

genes = np.array([f"G{i:03d}" for i in range(N_GENES)])
block = np.arange(N_GENES) // N_PLANTED  # 0 donor, 1 condition, 2 noise, 3+ background
donor_sd = np.select([block == 0, block == 1, block == 2], [1.5, 0.0, 0.0], 0.4)
stim_effect = np.select(
    [block == 1, block >= 3],
    [rng.choice([-2.0, 2.0], N_GENES) * np.log(2), rng.normal(0, 0.15, N_GENES)],
    0.0,
)

donor_effect = rng.normal(0, 1, (N_DONORS, N_GENES)) * donor_sd
donor_index = meta["donor"].map({d: i for i, d in enumerate(donors)}).to_numpy()
is_stim = (meta["condition"] == "stim").to_numpy()[:, None]
log_mean = (
    rng.uniform(np.log(100), np.log(1000), N_GENES)
    + donor_effect[donor_index]
    + is_stim * stim_effect
    + np.log(rng.uniform(0.8, 1.2, (len(meta), 1)))
)
mu = np.exp(log_mean)
counts = rng.negative_binomial(1 / DISPERSION, 1 / (1 + mu * DISPERSION))

pd.DataFrame(
    counts.T, index=pd.Index(genes, name="gene"), columns=meta["sample"]
).to_csv(out / "demo_counts.tsv", sep="\t")
meta.to_csv(out / "demo_meta.tsv", sep="\t", index=False)
