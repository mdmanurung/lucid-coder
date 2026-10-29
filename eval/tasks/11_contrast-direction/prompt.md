I have plasma protein levels for 12 donors, 6 with the disease and 6 healthy.

- `data/samples.tsv`: the sample sheet, one row per donor, with `sample_id` and `group` (`disease` or `healthy`).
- `data/expression.tsv`: log2 protein levels, proteins × samples. The first column is `protein`; every other column is one sample, headed by its sample ID.

Please test every protein for a difference between disease and healthy donors and write `de_results.tsv` (tab-separated) into the run directory, one row per protein, with the columns `protein`, `log2fc` (disease relative to healthy, so positive means higher in disease), `p_value` and `fdr`.

Write one R script at `scripts/11_test-groups.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/11_test-groups.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
