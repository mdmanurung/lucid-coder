"""Shared spec and ground truth for the pseudobulk-de bench tasks.

SPEC is the single source: make_tasks.py writes tasks.yaml prompts and answers from
it, and every tasks/<id>-<split>.py reproducer is `import _pb; _pb.main(__file__)`.
The procedure mirrors skills/pseudobulk-de/scripts/pseudobulk_de.py (steps 1-7) but
is written independently, so a bug in the skill cannot hide in its own answer key.
Corpus sources (relative to mycelium-extra/plans/_DONE/lucidcoder/docs/cc_learnings):
aggregation and >= 10% detection, repos/aifi-healthy-pbmc-reference/07-Pseudobulk/
35-R_pseudobulk_seurat.ipynb:57; design, .../08-Differential_expression/
36-R_pseudobulk_deseq2.ipynb:200-213.
"""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DATA = "kang2018_pb.h5ad"
PAIRED, UNPAIRED = "~ replicate + label", "~ label"


def split(cell_type, min_cells, min_frac, design=PAIRED, padj=0.05, lfc=1.0):
    return dict(
        cell_type=cell_type,
        min_cells=min_cells,
        min_frac=min_frac,
        design=design,
        padj=padj,
        lfc=lfc,
    )


# (task id) -> {"ask", "engine", "train"/"valid"/"test": split}.
SPEC = {
    "n_sig": dict(
        ask="n_sig",
        engine="python",
        train=split("CD14+ Monocytes", 20, 0.1),
        valid=split("FCGR3A+ Monocytes", 20, 0.05, padj=0.1, lfc=0.5),
        test=split("Dendritic cells", 15, 0.1),
    ),
    "filter": dict(
        ask="filter",
        engine="python",
        train=split("B cells", 30, 0.1),
        valid=split("NK cells", 30, 0.05),
        test=split("CD8 T cells", 20, 0.2),
    ),
    "top_gene": dict(
        ask="top_gene",
        engine="python",
        train=split("CD4 T cells", 20, 0.1),
        valid=split("B cells", 20, 0.1),
        test=split("NK cells", 20, 0.1),
    ),
    "unpaired": dict(
        ask="n_sig",
        engine="python",
        train=split("NK cells", 20, 0.1, design=UNPAIRED),
        valid=split("CD8 T cells", 20, 0.1, design=UNPAIRED),
        test=split("B cells", 20, 0.1, design=UNPAIRED),
    ),
    "r_engine": dict(
        ask="n_sig",
        engine="r",
        train=split("CD14+ Monocytes", 20, 0.2),
        valid=split("B cells", 20, 0.1),
        test=split("CD4 T cells", 20, 0.1, padj=0.1, lfc=0.5),
    ),
}

DATA_TEXT = (
    f"`{DATA}` in the working directory holds Kang et al. 2018 PBMCs: raw UMI "
    "counts in `X`; `obs` columns `sample` (one id per patient x condition), "
    "`replicate` (patient), `label` (`ctrl` or `stim`, interferon-beta) and "
    "`cell_type`."
)
ENGINE_TEXT = {
    "python": "Fit with PyDESeq2.",
    "r": "Fit with DESeq2 in R (`Rscript` is on PATH).",
}
ASK_TEXT = {
    "n_sig": "Report only the number of significant genes.",
    "filter": "Report only `<samples kept>/<genes kept>`, for example `12/3456`.",
    "top_gene": (
        "Report only the significant gene with the most negative log2 fold change."
    ),
}


def prompt(task_id, split_name):
    task, s = SPEC[task_id], SPEC[task_id][split_name]
    steps = (
        f"Pseudobulk the {s['cell_type']}: sum raw counts per `sample` over that cell "
        f"type's cells; drop samples with fewer than {s['min_cells']} cells; keep "
        f"genes detected (count > 0) in at least {s['min_frac']:g} of the cells of "
        "the kept samples."
    )
    if task["ask"] == "filter":
        return " ".join([DATA_TEXT, steps, ASK_TEXT["filter"]])
    de = (
        f"Then test stim vs ctrl: design `{s['design']}`, no fold-change shrinkage, "
        f"independent-filtering alpha {s['padj']:g}. {ENGINE_TEXT[task['engine']]} "
        f"A gene is significant when padj < {s['padj']:g} and |log2 fold change| > "
        f"{s['lfc']:g}; a missing padj is not significant. Save the results table and "
        "figures in the working directory."
    )
    return " ".join([DATA_TEXT, steps, de, ASK_TEXT[task["ask"]]])


def pseudobulk(adata, s):
    """Samples x genes summed counts and sample metadata, after both filters."""
    cells = adata[adata.obs["cell_type"].astype(str) == s["cell_type"]]
    obs = cells.obs[["sample", "replicate", "label"]].astype(str)
    n_cells = obs["sample"].value_counts()
    kept = sorted(n_cells.index[n_cells >= s["min_cells"]])
    in_kept = obs["sample"].isin(kept).to_numpy()
    detected = np.asarray((cells.X[in_kept] > 0).mean(axis=0)).ravel()
    genes = detected >= s["min_frac"]
    x = cells.X[in_kept][:, genes]
    labels = obs["sample"].to_numpy()[in_kept]
    counts = pd.DataFrame(
        {k: np.asarray(x[labels == k].sum(axis=0)).ravel() for k in kept},
        index=cells.var_names[genes],
    ).T.astype(np.int64)
    meta = obs.drop_duplicates("sample").set_index("sample").loc[kept]
    return counts, meta


def deseq_python(counts, meta, s):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design=s["design"],
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["label", "stim", "ctrl"],
        alpha=s["padj"],
        inference=inference,
        quiet=True,
    )
    stats.summary()
    res = stats.results_df
    return pd.DataFrame({"log2_fc": res["log2FoldChange"], "padj": res["padj"]})


R_DESEQ2 = """
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
counts <- t(as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE)))
meta <- read.delim(a[2], row.names = 1, stringsAsFactors = TRUE)
dds <- DESeq(DESeqDataSetFromMatrix(counts, meta, as.formula(a[3])), quiet = TRUE)
res <- results(dds, contrast = c("label", "stim", "ctrl"), alpha = as.numeric(a[4]))
write.table(data.frame(gene = rownames(res), log2_fc = res$log2FoldChange,
  padj = res$padj), a[5], sep = "\t", quote = FALSE, row.names = FALSE)
"""


def deseq_r(counts, meta, s):
    # The target venv's bin holds the Rscript shim (bench/target), next to its python.
    rscript = Path(sys.executable).parent / "Rscript"
    counts.to_csv("pb_counts.tsv", sep="\t")
    meta.to_csv("pb_meta.tsv", sep="\t")
    Path("deseq2.R").write_text(R_DESEQ2)
    args = ["pb_counts.tsv", "pb_meta.tsv", s["design"], str(s["padj"]), "pb_res.tsv"]
    subprocess.run([str(rscript), "deseq2.R", *args], check=True)
    return pd.read_csv("pb_res.tsv", sep="\t", index_col="gene")


def answer(task_id, split_name, adata):
    task, s = SPEC[task_id], SPEC[task_id][split_name]
    counts, meta = pseudobulk(adata, s)
    if task["ask"] == "filter":
        return f"{counts.shape[0]}/{counts.shape[1]}"
    res = (deseq_r if task["engine"] == "r" else deseq_python)(counts, meta, s)
    sig = res[(res["padj"] < s["padj"]) & (res["log2_fc"].abs() > s["lfc"])]
    if task["ask"] == "n_sig":
        return str(len(sig))
    return str(sig["log2_fc"].idxmin())


def main(script):
    task_id, split_name = Path(script).stem.rsplit("-", 1)
    Path("answer.md").write_text(answer(task_id, split_name, ad.read_h5ad(DATA)) + "\n")
