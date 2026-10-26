I have a plasma score for 24 donors from two recruitment cohorts.

- `data/samples.tsv`: the sample sheet, one row per donor, with `sample_id`, `cohort` and `age`.
- `data/scores.tsv`: one row per donor, with `sample_id` and `score`.

Cohort BR1 is the young donors and BR2 the older donors. Please summarise the mean score and the number of donors per age group (young/older) and write `cohort_summary.tsv` (tab-separated) into the run directory with the columns `age_group`, `n` and `mean_score`. Also write `samples_labelled.tsv` there, one row per donor, with `sample_id`, `cohort` and `age_group`.

Write one R script at `scripts/05_summarise-cohorts.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/05_summarise-cohorts.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, magrittr, ggplot2, here.
