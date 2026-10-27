"""Purpose: write examples/demo.tsv, synthetic repeated measures with time effects.
Inputs:  none
Outputs: examples/demo.tsv (long: feature, subject, time, age, sex, value)
Run:     python examples/make_demo.py
Notes:   24 subjects (S01-S24) measured at visits 1, 2 and 3; 30 features.
         F01-F06 rise 0.75 per visit (time truth); F07-F29 are null; every feature
         but F30 has a subject baseline (sd 1). F30 has none and its noise is
         centred within each subject, so the subject variance fits at 0 (singular).
         S24 misses visit 2 (every feature), so its figure line must break there.
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N_SUBJECTS, N_FEATURES, VISITS = 24, 30, [1, 2, 3]
SLOPE, SUBJECT_SD, NOISE_SD = 0.75, 1.0, 0.5
GAP_SUBJECT = "S24"

rng = np.random.default_rng(SEED)
subjects = pd.DataFrame(
    {
        "subject": [f"S{i:02d}" for i in range(1, N_SUBJECTS + 1)],
        "age": rng.integers(30, 81, N_SUBJECTS),
        "sex": rng.choice(["female", "male"], N_SUBJECTS),
    }
)
rows = subjects.merge(pd.DataFrame({"time": VISITS}), how="cross").merge(
    pd.DataFrame({"feature": [f"F{i:02d}" for i in range(1, N_FEATURES + 1)]}),
    how="cross",
)
number = rows["feature"].str[1:].astype(int)
baseline = rng.normal(0, SUBJECT_SD, (N_FEATURES, N_SUBJECTS))
baseline[N_FEATURES - 1] = 0
noise = pd.Series(rng.normal(0, NOISE_SD, len(rows)), index=rows.index)
singular = number == N_FEATURES
noise[singular] -= noise[singular].groupby(rows["subject"]).transform("mean")
rows["value"] = (
    baseline[number - 1, rows["subject"].str[1:].astype(int) - 1]
    + noise
    + SLOPE * (rows["time"] - 1) * (number <= 6)
).round(4)
rows = rows[~((rows["subject"] == GAP_SUBJECT) & (rows["time"] == 2))]
columns = ["feature", "subject", "time", "age", "sex", "value"]
rows[columns].sort_values(["feature", "subject", "time"]).to_csv(
    Path(__file__).resolve().parent / "demo.tsv", sep="\t", index=False
)
