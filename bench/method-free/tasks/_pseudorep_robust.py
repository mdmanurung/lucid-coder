"""Check the pseudorep key is method-robust and that the cell-level trap differs.

Per split: three sound donor-level methods must give the key; a naive cell-level
Wilcoxon must give a larger count. Run with the target venv python:
`<venv>/bin/python tasks/_pseudorep_robust.py cache`.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _pseudorep as pr  # noqa: E402  (sets the BLAS thread env before numpy)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy.stats import ttest_ind  # noqa: E402

EDGER = """
suppressMessages(library(edgeR))
a <- commandArgs(TRUE)
counts <- t(as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)))
cond <- factor(read.delim(a[2], row.names = 1)[colnames(counts), "condition"], c("ctrl", "stim"))
y <- normLibSizes(DGEList(counts))
d <- model.matrix(~ cond)
fit <- glmQLFit(estimateDisp(y, d), d)
res <- glmQLFTest(fit, coef = 2)$table
write.table(data.frame(gene = rownames(res), p = res$PValue), a[3], sep = "\t", row.names = FALSE, quote = FALSE)
"""


def welch(adata):
    counts, meta = pr.pseudobulk(adata)
    logcpm = np.log2(counts.div(counts.sum(axis=1), axis=0) * 1e6 + 1)
    stim = (meta["condition"] == "stim").to_numpy()
    return ttest_ind(logcpm[stim], logcpm[~stim], equal_var=False).pvalue


def edger(adata):
    counts, meta = pr.pseudobulk(adata)
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        counts.to_csv(t / "c.tsv", sep="\t")
        meta.to_csv(t / "m.tsv", sep="\t")
        (t / "e.R").write_text(EDGER)
        subprocess.run(["Rscript", t / "e.R", t / "c.tsv", t / "m.tsv", t / "o.tsv"], check=True)
        return pd.read_csv(t / "o.tsv", sep="\t", index_col=0)["p"]


def wilcoxon_cells(adata):
    import scanpy as sc

    x = adata.copy()
    x.X = x.X.astype(np.float32)
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    sc.tl.rank_genes_groups(x, "condition", groups=["stim"], reference="ctrl", method="wilcoxon")
    df = sc.get.rank_genes_groups_df(x, "stim").set_index("names")
    return df


def main(cache):
    keys = {t["id"]: t for t in yaml.safe_load((Path(__file__).parents[1] / "tasks.yaml").read_text())["tasks"]}
    rows = {}
    for split in pr.SPLITS:
        adata = ad.read_h5ad(Path(cache) / pr.h5ad(split))
        cand = adata.var_names[: pr.N_CAND]
        des = pr.deseq_pvals(adata)
        cells = wilcoxon_cells(adata)
        rows[split] = {
            "deseq2_bh10": pr.n_sig(des.loc[cand, "pvalue"]),
            "deseq2_padj_all": int((des.loc[cand, "padj"] < pr.ALPHA).sum()),
            "welch_logcpm": pr.n_sig(welch(adata)[: pr.N_CAND]),
            "edger_ql": pr.n_sig(edger(adata).loc[cand]),
            "naive_cell_wilcoxon": pr.n_sig(cells.loc[cand, "pvals"]),
            "naive_cell_wilcoxon_padj_all": int((cells.loc[cand, "pvals_adj"] < pr.ALPHA).sum()),
        }
        pv = pd.DataFrame({"deseq2": des.loc[cand, "pvalue"], "welch": welch(adata)[: pr.N_CAND], "cells": cells.loc[cand, "pvals"]})
        print(split, "candidate raw p:\n", pv.map(lambda v: f"{v:.1e}").to_string(), flush=True)
    table = pd.DataFrame(rows)
    table.loc["key"] = [int(keys["pseudorep"][s]["answer"]) for s in pr.SPLITS]
    print(table.to_string())
    for split in pr.SPLITS:
        col = table[split]
        sound = col[["deseq2_bh10", "deseq2_padj_all", "welch_logcpm", "edger_ql"]]
        assert (sound == col["key"]).all(), f"{split}: a sound method misses the key"
        naive = col[["naive_cell_wilcoxon", "naive_cell_wilcoxon_padj_all"]]
        assert (naive > col["key"]).all(), f"{split}: the naive method does not overcall"
    print("robust: ok")


if __name__ == "__main__":
    os.environ.setdefault("ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS", "1")
    main(sys.argv[1])
