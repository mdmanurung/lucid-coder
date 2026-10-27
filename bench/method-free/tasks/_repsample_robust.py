"""Check the repsample key is method-robust and that both naive methods miss it.

Per split: four sound methods that model the donor must give the key; per-sample
PyDESeq2 `~ condition` and a cell-level Wilcoxon must not. Run with the target venv
python, its bin on PATH for Rscript: `<venv>/bin/python tasks/_repsample_robust.py cache`.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _repsample as rs  # noqa: E402  (sets the BLAS thread env before numpy)

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy.stats import ttest_rel  # noqa: E402

EDGER = """
suppressMessages(library(edgeR))
a <- commandArgs(TRUE)
counts <- t(as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)))
m <- read.delim(a[2], row.names = 1)[colnames(counts), ]
donor <- factor(m$donor); cond <- factor(m$condition, c("ctrl", "stim"))
y <- normLibSizes(DGEList(counts))
d <- model.matrix(~ donor + cond)
fit <- glmQLFit(estimateDisp(y, d), d)
res <- glmQLFTest(fit, coef = "condstim")$table
write.table(data.frame(gene = rownames(res), p = res$PValue), a[3], sep = "\t", row.names = FALSE, quote = FALSE)
"""

DREAMLET = """
suppressMessages({library(dreamlet); library(SingleCellExperiment)})
a <- commandArgs(TRUE)
counts <- t(as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)))
obs <- read.delim(a[2], row.names = 1)[colnames(counts), ]
obs$cluster <- "one"
sce <- SingleCellExperiment(list(counts = counts), colData = obs)
pb <- aggregateToPseudoBulk(sce, assay = "counts", cluster_id = "cluster", sample_id = "sample", verbose = FALSE)
pb$condition <- factor(pb$condition, c("ctrl", "stim"))
res <- dreamlet(processAssays(pb, ~ (1 | donor), verbose = FALSE), ~ condition + (1 | donor), verbose = FALSE)
tab <- topTable(res, coef = "conditionstim", number = Inf)
write.table(data.frame(gene = tab$ID, p = tab$P.Value), a[3], sep = "\t", row.names = FALSE, quote = FALSE)
"""


def rscript(code, counts, meta):
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        counts.to_csv(t / "c.tsv", sep="\t")
        meta.to_csv(t / "m.tsv", sep="\t")
        (t / "s.R").write_text(code)
        subprocess.run(["Rscript", t / "s.R", t / "c.tsv", t / "m.tsv", t / "o.tsv"], check=True)
        return pd.read_csv(t / "o.tsv", sep="\t", index_col=0)["p"]


def paired_t(adata):
    counts, meta = rs.pseudobulk(adata)
    logcpm = np.log2(counts.div(counts.sum(axis=1), axis=0) * 1e6 + 1)
    stim = meta[meta["condition"] == "stim"].sort_values("donor").index
    ctrl = meta[meta["condition"] == "ctrl"].sort_values("donor").index
    return pd.Series(ttest_rel(logcpm.loc[stim].to_numpy(), logcpm.loc[ctrl].to_numpy()).pvalue, logcpm.columns)


def wilcoxon_cells(adata):
    import scanpy as sc

    x = adata.copy()
    x.X = x.X.astype(np.float32)
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    sc.tl.rank_genes_groups(x, "condition", groups=["stim"], reference="ctrl", method="wilcoxon")
    return sc.get.rank_genes_groups_df(x, "stim").set_index("names")


def main(cache):
    keys = {t["id"]: t for t in yaml.safe_load((Path(__file__).parents[1] / "tasks.yaml").read_text())["tasks"]}
    sound = ["deseq2_donor", "deseq2_donor_padj_all", "paired_t_logcpm", "edger_ql_donor", "dreamlet_1donor"]
    naive = ["naive_deseq2_condition", "naive_cell_wilcoxon"]
    rows = {}
    for split in rs.SPLITS:
        adata = ad.read_h5ad(Path(cache) / rs.h5ad(split))
        cand = adata.var_names[: rs.N_CAND]
        des = rs.deseq_pvals(adata)
        counts_cells = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
        p = {
            "deseq2_donor": des["pvalue"],
            "paired_t_logcpm": paired_t(adata),
            "edger_ql_donor": rscript(EDGER, *rs.pseudobulk(adata)),
            "dreamlet_1donor": rscript(DREAMLET, counts_cells, adata.obs.astype(str)),
            "naive_deseq2_condition": rs.deseq_pvals(adata, "~ condition")["pvalue"],
            "naive_cell_wilcoxon": wilcoxon_cells(adata)["pvals"],
        }
        rows[split] = {k: rs.n_sig(v.reindex(cand)) for k, v in p.items()}
        rows[split]["deseq2_donor_padj_all"] = int((des.loc[cand, "padj"] < rs.ALPHA).sum())
        pv = pd.DataFrame({k: v.reindex(cand) for k, v in p.items()})
        print(split, "candidate raw p:\n", pv.map(lambda v: f"{v:.1e}").to_string(), flush=True)
    table = pd.DataFrame(rows).loc[sound + naive]
    table.loc["key"] = [int(keys["repsample"][s]["answer"]) for s in rs.SPLITS]
    print(table.to_string())
    for split in rs.SPLITS:
        col = table[split]
        assert (col[sound] == col["key"]).all(), f"{split}: a sound method misses the key"
        assert (col[naive] != col["key"]).all(), f"{split}: a naive method hits the key"
    print("robust: ok")


if __name__ == "__main__":
    os.environ.setdefault("ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS", "1")
    main(sys.argv[1])
