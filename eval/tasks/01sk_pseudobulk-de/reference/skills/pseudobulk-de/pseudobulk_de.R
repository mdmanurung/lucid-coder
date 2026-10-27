# Purpose: test each gene for a disease-vs-healthy difference, with the donor as the unit.
# Inputs:  <input>/counts.tsv, <input>/cells.tsv   [gap: provided with the data]
# Outputs: <output>/de.tsv, donors.tsv, result.json, reproducibility/commands.sh,
#          session_info.txt; report.ipynb renders report.md from them.
# Run:     Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1
# Notes:   Cells from one donor are not independent, so counts are summed per donor
#          (pseudobulk) and edgeR quasi-likelihood tests 3 vs 3 donors.

library(dplyr)
library(edgeR)

# Setup ----
# Reference = healthy, so log2_fc reads as disease minus healthy.
condition_levels <- c("healthy", "disease")
q_cutoff <- 0.05

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[match("--input", args) + 1]
run_dir <- args[match("--output", args) + 1]
stopifnot(
  "usage: pseudobulk_de.R --input <data_dir> --output <out_dir>" = !anyNA(c(input_dir, run_dir))
)
stopifnot("output already exists; pick a new name" = !dir.exists(run_dir))
dir.create(file.path(run_dir, "reproducibility"), recursive = TRUE)

# Load inputs ----
counts <- readr::read_tsv(file.path(input_dir, "counts.tsv"), show_col_types = FALSE)
cells <- readr::read_tsv(file.path(input_dir, "cells.tsv"), show_col_types = FALSE)
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts$gene
# Cell order in cells.tsv is not guaranteed to match the count columns.
cells <- cells[match(colnames(count_matrix), cells$cell_id), ]
stopifnot("every count column has a row in cells.tsv" = !anyNA(cells$cell_id))

# Sum counts per donor ----
donors <- cells |>
  count(donor_id, condition, name = "n_cells") |>
  mutate(condition = forcats::fct(condition, levels = condition_levels))
stopifnot("each donor has one condition" = !anyDuplicated(donors$donor_id))
donor_counts <- sapply(donors$donor_id, \(donor) {
  rowSums(count_matrix[, cells$donor_id == donor, drop = FALSE])
})

# Filter and fit ----
design <- model.matrix(~condition, data = donors)
dge <- DGEList(donor_counts)
# filterByExpr keeps genes with enough counts in at least the smaller group (3 donors).
keep_gene <- filterByExpr(dge, design = design)
dge <- normLibSizes(dge[keep_gene, , keep.lib.sizes = FALSE])
fit <- glmQLFit(estimateDisp(dge, design), design)
ql_test <- glmQLFTest(fit, coef = "conditiondisease")
de_results <- tibble(
  gene = rownames(ql_test),
  log2_fc = ql_test$table$logFC,
  p_value = ql_test$table$PValue,
  # Multiple-testing family: every gene that passed filterByExpr.
  q_value = p.adjust(p_value, method = "BH")
) |>
  arrange(p_value)

# Write outputs ----
readr::write_tsv(de_results, file.path(run_dir, "de.tsv"))
readr::write_tsv(donors, file.path(run_dir, "donors.tsv"))
writeLines(
  c(
    paste("Rscript skills/pseudobulk-de/pseudobulk_de.R --input", input_dir, "--output", run_dir),
    paste("papermill skills/pseudobulk-de/report.ipynb", file.path(run_dir, "report.ipynb"), "-p out_dir", run_dir),
    paste("jupyter nbconvert --to markdown --no-input", file.path(run_dir, "report.ipynb"))
  ),
  file.path(run_dir, "reproducibility", "commands.sh")
)
jsonlite::write_json(
  list(
    parameters = list(
      input = input_dir,
      condition_levels = condition_levels,
      method = "edgeR glmQLFit on donor sums",
      p_adjust = "BH",
      q_cutoff = q_cutoff
    ),
    counts = list(
      donors = nrow(donors),
      cells = nrow(cells),
      genes = nrow(count_matrix),
      genes_tested = sum(keep_gene)
    ),
    n_significant = sum(de_results$q_value < q_cutoff)
  ),
  file.path(run_dir, "result.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
