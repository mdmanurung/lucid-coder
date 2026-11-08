"""Purpose: build the frozen inputs and truth of task 08b once; the task never fetches.
Inputs:  <kang_h5ad>  pertpy's kang_2018.h5ad, https://exampledata.scverse.org/pertpy/kang_2018.h5ad
                     (Kang et al. 2018, GEO GSE96583; NCBI GEO places no restriction on use)
         ../08_pathway-activity/inputs/data/progeny.tsv  task 08's frozen PROGENy copy
         ../08_pathway-activity/inputs/data/de.tsv       task 08's gene set (2,345 genes)
Outputs: inputs/data/counts.tsv    CD14+ monocyte pseudobulk counts, `gene` then 16 samples
         inputs/data/samples.tsv   sample, replicate (donor), label (ctrl or stim)
         inputs/data/progeny.tsv   identical to task 08's
         truth/pathways.tsv        ULM on the PyDESeq2 Wald stat, design ~ label
Run:     uv run --python 3.12 --with-requirements eval/envs/python-tasks.txt \
         python -I eval/tasks/08b_pathway-from-counts/make-inputs.py <kang_h5ad>
Notes:   Same pseudobulk and gene filters as task 08 (sc-best-practices "Differential gene
         expression"); the agent runs the DE itself. The truth is fitted on the written
         counts, so it is exactly what an agent can reproduce. Truth: JAK-STAT 16.72,
         NFkB 2.78 (task 08's frozen de.tsv gives 16.73).
"""

import shutil
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
TASK_08 = TASK.parent / "08_pathway-activity" / "inputs" / "data"

# %% Pseudobulk CD14+ monocytes per replicate x label (task 08, book cells 23-51)
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
genes_08 = pd.read_csv(TASK_08 / "de.tsv", sep="\t")["gene"]
if set(mono.var_names) != set(genes_08):
    raise ValueError("gene set differs from task 08's de.tsv")

# %% Counts and samples tables
mono.obs_names = mono.obs["sample"].astype(str)
samples = (
    mono.obs[["sample", "replicate", "label"]]
    .astype(str)
    .sort_values(["label", "replicate"])
)
counts = (
    mono.to_df()
    .loc[samples["sample"]]
    .T.round()
    .astype(int)
    .loc[sorted(mono.var_names)]
    .rename_axis("gene")
)

# %% Truth: PyDESeq2 ~ label on the written counts, ULM on the Wald stat
dds = DeseqDataSet(
    counts=counts.T, metadata=samples.set_index("sample"), design="~ label", quiet=True
)
dds.deseq2()
stats = DeseqStats(dds, contrast=["label", "stim", "ctrl"], quiet=True)
stats.summary()
progeny = pd.read_csv(TASK_08 / "progeny.tsv", sep="\t")
score, padj = dc.mt.ulm(data=stats.results_df[["stat"]].T, net=progeny)
truth = pd.DataFrame(
    {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)

counts.to_csv(DATA / "counts.tsv", sep="\t")
samples.to_csv(DATA / "samples.tsv", sep="\t", index=False)
shutil.copyfile(TASK_08 / "progeny.tsv", DATA / "progeny.tsv")
truth.sort_values("score").to_csv(
    TASK / "truth" / "pathways.tsv", sep="\t", index=False
)
print(counts.shape, samples, truth.sort_values("padj").head(5), sep="\n")
