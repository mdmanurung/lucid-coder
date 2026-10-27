"""Purpose: write examples/demo.tsv, synthetic sample-level values with planted effects.
Inputs:  none
Outputs: examples/demo.tsv (long: feature, subject, donor, visit, cohort, age, value)
Run:     python examples/make_demo.py
Notes:   20 donors (D01-D10 healthy, D11-D20 disease), each sampled pre and post.
         F01-F06 shift by +1 from pre to post (paired truth); F07-F12 are 2.5
         higher in disease (cohort truth); F13-F18 rise 0.1 per year of age
         (age truth); F19-F30 are null. Each donor has its own baseline per feature.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N_DONORS, N_FEATURES = 20, 30
VISIT_SHIFT, COHORT_SHIFT, AGE_SLOPE = 1.0, 2.5, 0.1

rng = np.random.default_rng(SEED)
donors = pd.DataFrame(
    {
        "donor": [f"D{i:02d}" for i in range(1, N_DONORS + 1)],
        "cohort": ["healthy"] * 10 + ["disease"] * 10,
        "age": rng.integers(30, 81, N_DONORS),
    }
)
rows = donors.merge(pd.DataFrame({"visit": ["pre", "post"]}), how="cross").merge(
    pd.DataFrame({"feature": [f"F{i:02d}" for i in range(1, N_FEATURES + 1)]}),
    how="cross",
)
number = rows["feature"].str[1:].astype(int)
baseline = rng.normal(0, 1, (N_FEATURES, N_DONORS))
rows["value"] = (
    baseline[number - 1, rows["donor"].str[1:].astype(int) - 1]
    + rng.normal(0, 0.3, len(rows))
    + VISIT_SHIFT * ((number <= 6) & (rows["visit"] == "post"))
    + COHORT_SHIFT * (number.between(7, 12) & (rows["cohort"] == "disease"))
    + AGE_SLOPE * (rows["age"] - 55) * number.between(13, 18)
).round(4)
rows["subject"] = rows["donor"] + "_" + rows["visit"]
columns = ["feature", "subject", "donor", "visit", "cohort", "age", "value"]
rows[columns].sort_values(["feature", "subject"]).to_csv(
    Path(__file__).resolve().parent / "demo.tsv", sep="\t", index=False
)
