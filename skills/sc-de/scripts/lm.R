# Purpose: per-feature linear-model engine for sc_de.py --mode lm; never run by hand.
# Inputs:  long.tsv (feature, value and every formula column; one row per sample)
# Outputs: <out_tsv> with feature, estimate, se, stat (t), pvalue (unadjusted);
#          the R version on stdout, appended to environment.txt by the caller.
# Run:     Rscript lm.R long.tsv "~ status + age + sex" status disease healthy out.tsv
#          A numeric term passes "" for level and reference.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c("data", "formula", "column", "level", "reference", "out")

data <- read.delim(args[["data"]], check.names = FALSE, stringsAsFactors = TRUE)
coef_name <- args[["column"]]
if (nzchar(args[["reference"]])) {
  data[[coef_name]] <- relevel(
    factor(data[[coef_name]]),
    ref = args[["reference"]]
  )
  coef_name <- paste0(coef_name, args[["level"]])
}
form <- as.formula(paste("value", args[["formula"]]))

# An aliased or unestimable term gives NA, as lm drops it from the coefficient table.
fit_one <- function(rows) {
  est <- coef(summary(lm(form, data = rows)))
  if (!coef_name %in% rownames(est)) {
    return(tibble(
      estimate = NA_real_,
      se = NA_real_,
      stat = NA_real_,
      pvalue = NA_real_
    ))
  }
  tibble(
    estimate = est[coef_name, 1],
    se = est[coef_name, 2],
    stat = est[coef_name, 3],
    pvalue = est[coef_name, 4]
  )
}

out <- split(data, ~feature, drop = TRUE) |>
  map(fit_one) |>
  list_rbind(names_to = "feature")
write.table(out, args[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(R.version.string, "\n", sep = "")
