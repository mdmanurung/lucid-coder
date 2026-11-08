I have annotated single-cell RNA-seq of mouse small-intestine epithelium: eight mice, four controls and four infected with Salmonella, in six epithelial cell types.

- `data/cells.tsv`: one row per cell, with `cell_id`, `sample_id` (one per mouse), `condition` (`Control` or `Salmonella`) and `cell_type`.

Goblet cells are known to be stable across these conditions, so treat them as the reference. Which of the other cell types change with infection? Write `composition.tsv` (tab-separated) into the run directory with one row per cell type and the columns `cell_type` (as written in `cells.tsv`), `log_fc` (Salmonella vs Control), `p_value` and `q_value` (Benjamini-Hochberg over the five non-reference cell types). The Goblet row may hold `NA`.

Write one R script at `scripts/09b_test-composition.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/09b_test-composition.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR, crumblr, variancePartition.
