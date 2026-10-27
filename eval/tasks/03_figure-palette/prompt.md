I have cell-type proportions for eight donors, four healthy and four with the disease.

- `data/proportions.tsv`: one row per donor and cell type, with `donor_id`, `condition` (`healthy` or `disease`), `cell_type` and `proportion` (the share of that donor's cells).
- `data/palette.tsv`: our lab colours, with `label` and `hex`. It has one row per cell type and one per condition.

Please make a figure of cell-type proportion by condition, with one point per donor, in our colours. Save it into the run directory as `proportions_by-condition.pdf` and `proportions_by-condition.png`, and save the plotted data as `proportions_by-condition.csv`.

Write one R script at `scripts/03_plot-proportions.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/03_plot-proportions.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
