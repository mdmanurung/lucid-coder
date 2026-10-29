# Naive answer: scores pathways on log2FoldChange instead of the Wald stat.
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

de = pd.read_csv("data/de.tsv", sep="\t", index_col="gene")
progeny = pd.read_csv("data/progeny.tsv", sep="\t")

score, padj = dc.mt.ulm(data=de[["log2FoldChange"]].T, net=progeny)
res = pd.DataFrame(
    {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)
res.to_csv(out / "pathways.tsv", sep="\t", index=False)
