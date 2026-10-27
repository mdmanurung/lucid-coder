"""Purpose: write the synthetic MOFA demo: three views and sample metadata.
Inputs:  none
Outputs: examples/demo_{rna,adt,tf}.tsv (samples x features), examples/demo_meta.tsv
Run:     python examples/make_demo.py
Notes:   60 samples, 30 healthy and 30 disease. Three planted factors:
         A is shared by all three views and is 1.5 higher in disease;
         B is shared by rna and adt; C is specific to tf. Each factor loads on a
         random half of a view's features; noise is N(0, 1). age is null.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N = 60
VIEWS = {"rna": 200, "adt": 40, "tf": 80}
# Factor -> views it loads on, with the loading scale.
PLANTED = {"A": {"rna": 1.0, "adt": 1.0, "tf": 1.0}, "B": {"rna": 0.8, "adt": 0.8}}
PLANTED["C"] = {"tf": 1.2}
LOAD_FRAC, GROUP_SHIFT = 0.5, 1.5

ROOT = Path(__file__).resolve().parent
rng = np.random.default_rng(SEED)
samples = [f"S{i:02d}" for i in range(1, N + 1)]
meta = pd.DataFrame(
    {
        "sample": samples,
        "group": ["healthy"] * 30 + ["disease"] * 30,
        "age": rng.integers(30, 81, N),
    }
)
z = {k: rng.normal(0, 1, N) for k in PLANTED}
z["A"] += GROUP_SHIFT * (meta["group"] == "disease").to_numpy()
z["A"] -= z["A"].mean()
for view, n_feat in VIEWS.items():
    y = rng.normal(0, 1, (N, n_feat))
    for k, scales in PLANTED.items():
        if view in scales:
            w = rng.normal(0, scales[view], n_feat) * (rng.random(n_feat) < LOAD_FRAC)
            y += np.outer(z[k], w)
    names = [f"{view}_{j:03d}" for j in range(1, n_feat + 1)]
    pd.DataFrame(
        y.round(4), index=pd.Index(samples, name="sample"), columns=names
    ).to_csv(ROOT / f"demo_{view}.tsv", sep="\t")
meta.to_csv(ROOT / "demo_meta.tsv", sep="\t", index=False)
