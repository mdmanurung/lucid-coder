# Naive answer: scores the table as given, so only the rows that carry symbols match CollecTRI.
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

de = pd.read_csv("data/de.tsv.gz", sep="\t", index_col="gene")
collectri = pd.read_csv("data/collectri.tsv.gz", sep="\t")

score, padj = dc.mt.ulm(data=de[["stat"]].T, net=collectri)
res = pd.DataFrame(
    {"source": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)
res.to_csv(out / "tf_activity.tsv", sep="\t", index=False)
