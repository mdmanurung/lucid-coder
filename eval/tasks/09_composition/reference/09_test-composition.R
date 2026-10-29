# Purpose: test whether Salmonella infection shifts epithelial cell-type composition, with the
#          mouse as the unit of replication.
# Run:     Rscript scripts/09_test-composition.R results/run1
# Inputs:  data/cells.tsv  [gap: one row per annotated cell; mouse, condition, cell type]
# Outputs: <run_dir>/composition.tsv     one row per cell type: log_fc, p, BH q
#          <run_dir>/session_info.txt    R session
# Notes:   [Stochastic] none. Cells are counted per mouse; crumblr turns the counts into
#          centred log-ratios with precision weights, so the types are modelled jointly and a
#          rise in one does not read as a fall in the rest. log_fc is on the CLR scale.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(limma)
  library(crumblr)
})

# Parameters --------------------------------------------------------------------

# Reference first, so log_fc reads Salmonella minus Control.
CONDITION_LEVELS <- c("Control", "Salmonella")
# Multiple-testing family: the eight cell types for this one contrast.
P_ADJUST_METHOD <- "BH"

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " exists; pick a new name.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Count cells per mouse ---------------------------------------------------------

cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)
counts <- unclass(table(cells$sample_id, cells$cell_type))
samples <- cells |>
  distinct(sample_id, condition) |>
  mutate(condition = forcats::fct(condition, levels = CONDITION_LEVELS)) |>
  arrange(match(sample_id, rownames(counts)))
stopifnot(
  "samples match counts" = identical(samples$sample_id, rownames(counts))
)

# Fit the compositional model ---------------------------------------------------

clr <- crumblr(counts)
design <- model.matrix(~condition, samples)
fit <- eBayes(lmFit(clr, design))
composition <- topTable(
  fit,
  coef = "conditionSalmonella",
  number = Inf,
  adjust.method = P_ADJUST_METHOD
) |>
  as_tibble(rownames = "cell_type") |>
  select(cell_type, log_fc = logFC, p_value = P.Value, q_value = adj.P.Val)

# Write outputs -----------------------------------------------------------------

write_tsv(composition, file.path(run_dir, "composition.tsv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
