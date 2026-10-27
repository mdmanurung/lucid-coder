I have a raw count table from a bulk RNA-seq run that still contains the control features.

- `data/counts.tsv`: raw counts, genes × samples. The first column is `gene`; every other column is one sample, headed by its sample ID. Spike-ins and mitochondrial genes are in there too, named `ERCC-` and `MT-`.

Please drop the spike-in and mitochondrial genes and write `counts_filtered.tsv` (tab-separated, same columns as the input) into the run directory. Also write `filter_log.tsv` there with one row per control set, columns `gene_set` (ERCC, MT) and `n_removed`.

Write one R script at `scripts/06_filter-controls.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/06_filter-controls.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
