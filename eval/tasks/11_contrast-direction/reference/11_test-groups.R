# Purpose: test each protein for a disease vs healthy difference, healthy as the reference.
# Run:     Rscript scripts/11_test-groups.R results/run1
# Inputs:  data/samples.tsv     [gap] one row per donor: sample_id, group
#          data/expression.tsv  [gap] log2 protein levels, proteins x samples
# Outputs: <run_dir>/de_results.tsv    per protein: log2fc (disease - healthy), p_value, fdr (BH
#                                      over all proteins)
#          <run_dir>/session_info.txt  R session

# Parameters -------------------------------------------------------------------
# Reference level first, so log2fc reads disease minus healthy.
GROUP_LEVELS <- c("healthy", "disease")

# Setup ------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
})
run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir is new" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Load data --------------------------------------------------------------------
samples <- readr::read_tsv(here::here("data", "samples.tsv"), show_col_types = FALSE)
expression <- read.delim(
  here::here("data", "expression.tsv"),
  row.names = 1,
  check.names = FALSE
) |>
  as.matrix()
stopifnot(
  "samples rows match expression columns" = identical(samples$sample_id, colnames(expression))
)

# Fit limma --------------------------------------------------------------------
# fct() stops on a group label outside GROUP_LEVELS.
group <- forcats::fct(samples$group, levels = GROUP_LEVELS)
fit <- lmFit(expression, model.matrix(~group)) |>
  eBayes()
de_results <- topTable(fit, coef = "groupdisease", number = Inf, sort.by = "none") |>
  as_tibble(rownames = "protein") |>
  select(protein, log2fc = logFC, p_value = P.Value, fdr = adj.P.Val)

# Write outputs ----------------------------------------------------------------
readr::write_tsv(de_results, file.path(run_dir, "de_results.tsv"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
