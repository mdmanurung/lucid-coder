I have PBMCs from eight donors, each split into a control sample and one stimulated with IFN-β. The single-cell counts are already summed per donor, condition and cell type (pseudobulk) for CD14+ monocytes and CD4 T cells: 32 samples, 2002 genes.

- `data/counts.tsv`: summed raw UMI counts, genes × samples. The first column is `gene`; every other column is one sample, headed by its sample ID.
- `data/samples.tsv`: one row per sample, with `sample_id` (matches the count column names), `donor_id`, `cell_type` (`CD14+ Monocytes` or `CD4 T cells`), `condition` (`ctrl` or `stim`) and `n_cells` (cells summed into the sample).

Which genes respond to IFN-β differently in CD14+ monocytes than in CD4 T cells? Write `response_difference.tsv` (tab-separated) into the run directory with one row per gene tested and the columns `gene`, `log2_fc_difference` (the stim vs ctrl log2 fold change in CD14+ monocytes minus the one in CD4 T cells) and `differs` (`TRUE` or `FALSE`: does the response differ between the two cell types, at a 5% FDR over all genes tested).

Write one R script at `scripts/14_compare-responses.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/14_compare-responses.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
