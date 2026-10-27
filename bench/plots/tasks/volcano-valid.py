"""Reproduce the volcano valid answer; writes answer.md in the cwd."""

from pathlib import Path

import anndata as ad
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

CELL_TYPE, PADJ_MAX, LFC_MIN = "NK cells", 0.05, 1


def pseudobulk(adata, cell_type):
    """Sum raw counts per replicate x label; keep genes with >= 10 counts in total."""
    sub = adata[adata.obs["cell_type"] == cell_type]
    key = sub.obs["replicate"].astype(str) + "_" + sub.obs["label"].astype(str)
    counts = (
        pd.DataFrame(sub.X.toarray(), columns=sub.var_names).groupby(key.values).sum()
    )
    counts = counts.loc[:, counts.sum() >= 10]
    meta = counts.index.to_series().str.rsplit("_", n=1, expand=True)
    meta.columns = ["replicate", "label"]
    return counts, meta


counts, meta = pseudobulk(ad.read_h5ad("kang2018.h5ad"), CELL_TYPE)
inference = DefaultInference(n_cpus=2)
dds = DeseqDataSet(
    counts=counts, metadata=meta, design="~replicate + label", inference=inference
)
dds.deseq2()
stats = DeseqStats(
    dds, contrast=["label", "stim", "ctrl"], alpha=PADJ_MAX, inference=inference
)
stats.summary()
res = stats.results_df.dropna(subset=["padj"])
n_sig = int(((res["padj"] < PADJ_MAX) & (res["log2FoldChange"].abs() > LFC_MIN)).sum())
Path("answer.md").write_text(f"{n_sig}\n")
