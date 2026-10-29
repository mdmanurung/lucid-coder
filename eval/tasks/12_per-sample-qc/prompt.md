I have single-cell RNA-seq from three PBMC libraries, each run through Cell Ranger, about 380 cells each.

- `data/SRR28679756/`, `data/SRR28679757/`, `data/SRR28679759/`: one 10x filtered matrix per library (`matrix.mtx`, `genes.tsv`, `barcodes.tsv`), raw UMI counts, genes × cells.

Before clustering I want to drop low-quality cells. Can you run cell QC and write `cells_qc.tsv` (tab-separated) into the run directory with one row per cell and the columns `sample` (the library folder name), `barcode` and `keep` (`TRUE` or `FALSE`)? Also write `qc_summary.tsv` there with one row per library: `sample`, `n_cells` and `n_kept`.

Write one R script at `scripts/12_filter-cells.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/12_filter-cells.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, ggplot2, here, Matrix, SingleCellExperiment, scuttle, DropletUtils.
