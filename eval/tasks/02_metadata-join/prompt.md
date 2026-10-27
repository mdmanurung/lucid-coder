I have bulk RNA-seq counts for 12 samples and a sample sheet from the sequencing core.

- `data/counts.tsv`: raw counts, genes × samples. The first column is `gene`; every other column is one sample, headed by its sample ID.
- `data/metadata.tsv`: one row per sample, with `sample_id`, `condition` and `batch`.

Please attach condition and batch to each sample and write `sample_summary.tsv` (tab-separated) into the run directory, with one row per sample and the columns `sample_id`, `condition`, `batch` and `library_size` (total counts in that sample).

Write one R script at `scripts/02_summarise-samples.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/02_summarise-samples.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here.
