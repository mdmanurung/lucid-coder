# Naive answer: treats every cell as a replicate (t-test on log-normalised cells).
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv("data/counts.tsv", sep="\t", index_col="gene")
cells = pd.read_csv("data/cells.tsv", sep="\t").set_index("cell_id").loc[counts.columns]

norm = np.log1p(counts / counts.sum() * 1e4)
disease = (cells["condition"] == "disease").to_numpy()

rows = []
for gene, values in norm.iterrows():
    a = values[disease]
    b = values[~disease]
    rows.append(
        {
            "gene": gene,
            "log2_fc": np.log2(np.expm1(a).mean() + 1) - np.log2(np.expm1(b).mean() + 1),
            "p_value": ttest_ind(a, b, equal_var=False).pvalue,
        }
    )
res = pd.DataFrame(rows)
res["q_value"] = false_discovery_control(res["p_value"], method="bh")
res.to_csv(out / "de.tsv", sep="\t", index=False)
