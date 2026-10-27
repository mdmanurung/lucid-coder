I have single-cell RNA-seq data from six donors: three healthy and three with the disease, about 200 cells per donor, 50 genes.

- `data/counts.tsv`: raw UMI counts, genes × cells. The first column is `gene`; every other column is one cell, headed by its cell ID.
- `data/cells.tsv`: one row per cell, with `cell_id` (matches the count column names), `donor_id` and `condition` (`healthy` or `disease`).

Can you test which genes differ between disease and healthy? Write `de.tsv` (tab-separated) into the run directory with one row per gene and the columns `gene`, `log2_fc` (disease vs healthy), `p_value` and `q_value` (Benjamini-Hochberg over all genes).

Write one R script at `scripts/01_compare-conditions.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/01_compare-conditions.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
