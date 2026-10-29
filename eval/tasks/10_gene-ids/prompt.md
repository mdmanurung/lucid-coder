I have bulk RNA-seq of hepatic stellate cells, three control and three TGF-beta treated samples, and a DESeq2 table for treated vs control.

- `data/de.tsv.gz`: one row per gene, with `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj`.
- `data/genes.tsv.gz`: the gene annotation, with `gene_id`, `gene_name` and `chr`.
- `data/collectri.tsv.gz`: the CollecTRI TF-target network, with `source` (the TF), `target` (the gene) and `weight` (+1 activation, -1 repression).

Which transcription factors does TGF-beta activate or repress? Score TF activity on the DESeq2 Wald statistic (`stat`) with decoupler's univariate linear model (`ulm`) and write `tf_activity.tsv` (tab-separated) into the run directory, one row per TF scored, with the columns `source`, `score` (positive means more active after TGF-beta) and `padj`.

Write one Python script at `scripts/10_score-tfs.py` that takes the run directory as its only argument and writes its outputs there. Run it once with `uv run --python 3.12 --with-requirements requirements.txt python scripts/10_score-tfs.py results/run1`. Use only the standard library and: pandas, numpy, scipy, matplotlib, decoupler.
