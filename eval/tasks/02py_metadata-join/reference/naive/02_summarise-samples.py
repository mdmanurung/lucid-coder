# Naive answer: binds the metadata rows onto the count columns by position.
import sys
from pathlib import Path

import pandas as pd

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv("data/counts.tsv", sep="\t", index_col="gene")
meta = pd.read_csv("data/metadata.tsv", sep="\t")

res = pd.DataFrame(
    {
        "sample_id": counts.columns,
        "condition": meta["condition"].to_numpy(),
        "batch": meta["batch"].to_numpy(),
        "library_size": counts.sum().to_numpy(),
    }
)
res.to_csv(out / "sample_summary.tsv", sep="\t", index=False)
