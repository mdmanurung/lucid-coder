# Purpose: test which epithelial cell types change under Salmonella, with Goblet cells as the
#          known-stable reference and the mouse as the unit of replication.
# Run:     Rscript scripts/09b_test-composition.R results/run1
# Inputs:  data/cells.tsv  [gap: one row per annotated cell; mouse, condition, cell type]
# Outputs: <run_dir>/composition.tsv     one row per cell type: log_fc, p, BH q (Goblet NA)
#          <run_dir>/session_info.txt    R session
# Notes:   [Stochastic] none. Cells are counted per mouse. edgeR QL models each non-reference
#          type's count with offset log(Goblet count) in place of the library size, so log_fc
#          is the log2 change relative to Goblet and a rise in one type does not read as a fall
#          in the rest.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(edgeR)
})

# Parameters --------------------------------------------------------------------

# Reference first, so log_fc reads Salmonella minus Control.
CONDITION_LEVELS <- c("Control", "Salmonella")
# The cell type known to be stable across conditions.
REFERENCE_TYPE <- "Goblet"
# Multiple-testing family: the five non-reference cell types for this one contrast.
P_ADJUST_METHOD <- "BH"

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " exists; pick a new name.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Count cells per mouse ---------------------------------------------------------

cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)
counts <- unclass(table(cells$cell_type, cells$sample_id))
samples <- cells |>
  distinct(sample_id, condition) |>
  mutate(condition = forcats::fct(condition, levels = CONDITION_LEVELS)) |>
  arrange(match(sample_id, colnames(counts)))
stopifnot(
  "samples match counts" = identical(samples$sample_id, colnames(counts))
)

# Fit counts against the reference ----------------------------------------------

design <- model.matrix(~condition, samples)
dge <- DGEList(counts[rownames(counts) != REFERENCE_TYPE, ])
# Offset log(Goblet count) per mouse replaces the library size.
dge$offset <- matrix(
  log(counts[REFERENCE_TYPE, ]),
  nrow(dge),
  ncol(dge),
  byrow = TRUE
)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
composition <- topTags(
  glmQLFTest(fit, coef = "conditionSalmonella"),
  n = Inf,
  adjust.method = P_ADJUST_METHOD
)$table |>
  as_tibble(rownames = "cell_type") |>
  select(cell_type, log_fc = logFC, p_value = PValue, q_value = FDR) |>
  bind_rows(tibble(cell_type = REFERENCE_TYPE))

# Write outputs -----------------------------------------------------------------

write_tsv(composition, file.path(run_dir, "composition.tsv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
