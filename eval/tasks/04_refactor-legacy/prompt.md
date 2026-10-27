`scripts/00_legacy.R` is an old script from a former lab member. It reads `data/counts.tsv` (raw counts: a `gene` column, then one column per sample) and `data/samples.tsv` (`sample_id`, `group`), and writes four tables into the directory given as its argument: `Rscript scripts/00_legacy.R <dir>`. We still use those tables.

Please refactor it into a new script. The outputs must stay the same: same files, same rows and columns, same values. Leave `scripts/00_legacy.R` as it is.

Write one R script at `scripts/04_group-means.R` that takes the run directory as its only argument and writes its outputs there. Run it once with `Rscript scripts/04_group-means.R results/run1`. Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here.
