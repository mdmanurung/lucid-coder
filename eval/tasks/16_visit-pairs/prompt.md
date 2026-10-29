I measured plasma IL-6 in 12 donors at three visits: baseline (`D0`), day 7 (`D7`) and day 90 (`D90`).

- `data/runs.tsv`: one row per assay run, with `run_id`, `donor_id`, `visit`, `run_date` and `il6` (log2 units). A few samples were run more than once; the most recent run of a sample replaces the earlier ones.

Please test the change in IL-6 from baseline to day 7 within donors. Write `pairs.tsv` (tab-separated) into the run directory, one row per donor in the test, with the columns `donor_id`, `il6_d0`, `il6_d7` and `change` (day 7 minus baseline). Also write `paired_test.tsv` there with one row: `n_donors`, `mean_change` and `p_value`.

Write one R script at `scripts/16_pair-visits.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/16_pair-visits.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR.
