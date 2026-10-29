I have annotated single-cell RNA-seq of mouse small-intestine epithelium (Haber et al. 2017): six mice, four controls and two infected with Salmonella, about 5,000 cells in eight epithelial cell types.

- `data/cells.tsv`: one row per cell, with `cell_id`, `sample_id` (one per mouse), `condition` (`Control` or `Salmonella`) and `cell_type`.

Does the infection change the cell-type composition of the epithelium? Write `composition.tsv` (tab-separated) into the run directory with one row per cell type and the columns `cell_type` (as written in `cells.tsv`), `log_fc` (Salmonella vs Control), `p_value` and `q_value` (Benjamini-Hochberg over the cell types).

Write one R script at `scripts/09_test-composition.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/09_test-composition.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR, crumblr, variancePartition.
