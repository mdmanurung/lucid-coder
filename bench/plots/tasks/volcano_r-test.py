"""Reproduce the volcano_r test answer; writes answer.md in the cwd."""

import subprocess
from pathlib import Path

import anndata as ad
import pandas as pd

CELL_TYPE, PADJ_MAX, LFC_MIN = "CD14+ Monocytes", 0.05, 2


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
counts.T.to_csv("counts.csv")
meta.to_csv("meta.csv")
r_code = f"""
suppressPackageStartupMessages(library(DESeq2))
counts <- as.matrix(read.csv("counts.csv", row.names = 1, check.names = FALSE))
meta <- read.csv("meta.csv", row.names = 1)
meta$label <- factor(meta$label, levels = c("ctrl", "stim"))
design <- ~ replicate + label
dds <- DESeqDataSetFromMatrix(round(counts), meta[colnames(counts), ], design)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("label", "stim", "ctrl"), alpha = {PADJ_MAX})
res <- res[!is.na(res$padj), ]
cat(sum(res$padj < {PADJ_MAX} & abs(res$log2FoldChange) > {LFC_MIN}))
"""
# Plain `Rscript` is the target venv's shim: check prepends the venv bin to PATH.
out = subprocess.run(
    ["Rscript", "-e", r_code], capture_output=True, text=True, check=True
)
Path("answer.md").write_text(out.stdout.strip() + "\n")
