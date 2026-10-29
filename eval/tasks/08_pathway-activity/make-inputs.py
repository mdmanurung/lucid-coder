"""Purpose: build the frozen inputs and truth of task 08 once; the task never fetches.
Inputs:  <kang_h5ad>  pertpy's kang_2018.h5ad, https://exampledata.scverse.org/pertpy/kang_2018.h5ad
                     (Kang et al. 2018, GEO GSE96583; NCBI GEO places no restriction on use)
         PROGENy human top 100 per pathway, fetched 2026-10-29 from OmniPath by decoupler
         (PROGENy, saezlab/progeny, Apache-2.0; license="academic" and "commercial" both
         return the same 1398 rows)
Outputs: inputs/data/de.tsv       CD14+ monocyte pseudobulk DESeq2 table, stim vs ctrl
         inputs/data/progeny.tsv  PROGENy weights (source, target, weight, padj)
         truth/pathways.tsv       ULM on the Wald stat: score and padj per pathway
Run:     uv run --python 3.12 --with-requirements eval/envs/python-tasks.txt \
         python -I eval/tasks/08_pathway-activity/make-inputs.py <kang_h5ad>
Notes:   Mirrors sc-best-practices "Differential gene expression" (CD14+ monocytes, design
         ~ label, unpaired) and "Gene set enrichment and pathway activity" (ulm on stat).
         Truth here: JAK-STAT 16.73, NFkB 2.78. The book reports 17.15 and 2.78 with its
         own 1400-row PROGENy copy; the gate uses the frozen numbers.
"""

import sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

TASK = Path(__file__).resolve().parent
DATA = TASK / "inputs" / "data"

# %% Pseudobulk CD14+ monocytes per replicate x label (book cells 23-51)
adata = ad.read_h5ad(sys.argv[1])
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.obs["sample"] = (
    adata.obs["replicate"].astype(str) + "_" + adata.obs["label"].astype(str)
)
pdata = dc.pp.pseudobulk(
    adata, sample_col="sample", groups_col="cell_type", layer="counts", mode="sum"
)
dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000)
mono = pdata[pdata.obs["cell_type"] == "CD14+ Monocytes"].copy()
dc.pp.filter_by_expr(
    mono, group="label", min_count=10, min_total_count=15, large_n=10, min_prop=0.7
)
dc.pp.filter_by_prop(mono, min_prop=0.1, min_smpls=2)

# %% DESeq2, stim vs ctrl
dds = DeseqDataSet(adata=mono, design="~ label", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["label", "stim", "ctrl"], quiet=True)
stats.summary()
de = (
    stats.results_df.rename_axis("gene")
    .reset_index()
    .sort_values("pvalue")
    .round({"baseMean": 3, "log2FoldChange": 6, "lfcSE": 6, "stat": 6})
)

# %% PROGENy and the truth
progeny = dc.op.progeny(organism="human", top=100)
wald = de.set_index("gene")[["stat"]].T
score, padj = dc.mt.ulm(data=wald, net=progeny)
truth = pd.DataFrame(
    {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)

de.to_csv(DATA / "de.tsv", sep="\t", index=False)
progeny.to_csv(DATA / "progeny.tsv", sep="\t", index=False)
truth.sort_values("score").to_csv(
    TASK / "truth" / "pathways.tsv", sep="\t", index=False
)
print(mono.shape, truth.sort_values("padj").head(5), sep="\n")
