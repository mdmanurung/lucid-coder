"""Purpose: build the frozen inputs and truth of task 10 once; the task never fetches.
Inputs:  <counts_tsv_gz>  GEO GSE151251 HSCs_Ctrl.vs.HSCs_TGFb.counts.tsv.gz (NCBI GEO: no
                          restriction on use; the decoupler bulk tutorial's dataset)
         CollecTRI human, fetched 2026-10-29 from OmniPath by decoupler with
                          license="commercial" (only resources open to any use; the
                          saezlab/CollecTRI build code is GPL-3.0)
Outputs: inputs/data/de.tsv.gz        DESeq2 TGF-beta vs control, gene IDs mostly Ensembl
         inputs/data/genes.tsv.gz     annotation: gene_id, gene_name, chr
         inputs/data/collectri.tsv.gz CollecTRI (source, target, weight)
         truth/tfs.tsv                ULM on the Wald stat after mapping IDs to symbols
Run:     uv run --python 3.12 --with-requirements eval/envs/python-tasks.txt \
         python -I eval/tasks/10_gene-ids/make-inputs.py <counts_tsv_gz>
Notes:   The planted trap: SYMBOL_FRACTION of rows keep the symbol (a second annotation
         merged in), the rest carry Ensembl IDs, so an unmapped ULM scores a few TFs.
"""

import sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

SYMBOL_FRACTION = 0.15
SEED = 0
TASK = Path(__file__).resolve().parent
DATA = TASK / "inputs" / "data"

# %% Counts, filtered as in the decoupler bulk tutorial
raw = pd.read_csv(sys.argv[1], sep="\t")
genes = raw[["GeneID", "GeneName", "Chr"]].set_axis(
    ["gene_id", "gene_name", "chr"], axis=1
)
counts = raw.set_index("GeneID").iloc[:, 5:].T
adata = ad.AnnData(counts.astype(float))
adata.obs["condition"] = np.where(
    adata.obs_names.str.contains("TGFb"), "treatment", "control"
)
dc.pp.filter_by_expr(
    adata, group="condition", min_count=10, min_total_count=15, large_n=10, min_prop=0.7
)

# %% DESeq2, treatment vs control
dds = DeseqDataSet(adata=adata, design="~ condition", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["condition", "treatment", "control"], quiet=True)
stats.summary()
de = stats.results_df.rename_axis("gene_id").reset_index()

# %% Plant mixed IDs: some rows carry the symbol instead of the Ensembl ID
symbol_of = genes.set_index("gene_id")["gene_name"]
rng = np.random.default_rng(SEED)
as_symbol = rng.random(len(de)) < SYMBOL_FRACTION
de["gene"] = np.where(as_symbol, de["gene_id"].map(symbol_of), de["gene_id"])
de = (
    de[["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
    .sort_values("pvalue")
    .round({"baseMean": 2, "log2FoldChange": 4, "lfcSE": 4, "stat": 4})
)

# %% CollecTRI and the truth: map to symbols, collapse duplicates by max |stat|
collectri = dc.op.collectri(organism="human", license="commercial")[
    ["source", "target", "weight"]
]
mapped = de.assign(symbol=de["gene"].map(symbol_of).fillna(de["gene"]))
mapped = mapped.loc[mapped["stat"].abs().groupby(mapped["symbol"]).idxmax()]
score, padj = dc.mt.ulm(data=mapped.set_index("symbol")[["stat"]].T, net=collectri)
truth = pd.DataFrame(
    {"source": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
)
truth = truth.sort_values("score", ascending=False)

# %% Naive check: the unmapped table must still score a few TFs without an error
naive, _ = dc.mt.ulm(data=de.set_index("gene")[["stat"]].T, net=collectri)

de.to_csv(DATA / "de.tsv.gz", sep="\t", index=False)
genes[genes["gene_id"].isin(stats.results_df.index)].to_csv(
    DATA / "genes.tsv.gz", sep="\t", index=False
)
collectri.to_csv(DATA / "collectri.tsv.gz", sep="\t", index=False)
truth.to_csv(TASK / "truth" / "tfs.tsv", sep="\t", index=False)
print(
    adata.shape,
    f"truth TFs {len(truth)}, naive TFs {naive.shape[1]}",
    truth.head(8),
    sep="\n",
)
print("naive top:", naive.iloc[0].sort_values(ascending=False).head(8).to_dict())
