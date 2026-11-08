# Naive answer: runs PyDESeq2 correctly, then scores pathways on log2FoldChange.
import sys
from pathlib import Path

import decoupler as dc
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv("data/counts.tsv", sep="\t", index_col="gene").T
samples = pd.read_csv("data/samples.tsv", sep="\t", index_col="sample")
progeny = pd.read_csv("data/progeny.tsv", sep="\t")

dds = DeseqDataSet(counts=counts, metadata=samples, design="~ label", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["label", "stim", "ctrl"], quiet=True)
stats.summary()

score, padj = dc.mt.ulm(data=stats.results_df[["log2FoldChange"]].T, net=progeny)
res = pd.DataFrame(
    {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)
res.to_csv(out / "pathways.tsv", sep="\t", index=False)
