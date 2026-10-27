# Purpose: test each gene for a disease-vs-healthy difference, with the donor as the unit.
# Inputs:  <input>/counts.tsv, <input>/cells.tsv   [gap: provided with the data]
# Outputs: <output>/de.tsv, report.md, result.json, reproducibility/commands.sh,
#          session_info.txt
# Run:     Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1
# Notes:   [Stochastic] none. Cells from one donor are not independent, so counts are summed
#          per donor (pseudobulk) and edgeR quasi-likelihood tests 3 vs 3 donors.

library(dplyr)
library(edgeR)

# ---- Setup ----
# Reference = healthy, so log2_fc reads as disease minus healthy.
CONDITION_LEVELS <- c("healthy", "disease")
P_ADJUST_METHOD <- "BH"
Q_CUTOFF <- 0.05
N_TOP_GENES <- 10

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[match("--input", args) + 1]
run_dir <- args[match("--output", args) + 1]
if (length(args) != 4 || anyNA(c(input_dir, run_dir))) {
  stop(
    "usage: Rscript skills/pseudobulk-de/pseudobulk_de.R --input <data_dir> --output <out_dir>",
    call. = FALSE
  )
}
if (dir.exists(run_dir)) {
  stop("output ", run_dir, " already exists; pick a new name.", call. = FALSE)
}
dir.create(file.path(run_dir, "reproducibility"), recursive = TRUE)

# ---- Load inputs ----
counts <- readr::read_tsv(
  file.path(input_dir, "counts.tsv"),
  show_col_types = FALSE
)
cells <- readr::read_tsv(
  file.path(input_dir, "cells.tsv"),
  show_col_types = FALSE
)
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts$gene
# Cell order in cells.tsv is not guaranteed to match the count columns.
cells <- cells[match(colnames(count_matrix), cells$cell_id), ]
stopifnot("every count column has a row in cells.tsv" = !anyNA(cells$cell_id))

# ---- Sum counts per donor ----
donors <- cells |>
  count(donor_id, condition, name = "n_cells") |>
  mutate(condition = forcats::fct(condition, levels = CONDITION_LEVELS))
stopifnot("each donor has one condition" = !anyDuplicated(donors$donor_id))
donor_counts <- sapply(donors$donor_id, \(donor) {
  rowSums(count_matrix[, cells$donor_id == donor, drop = FALSE])
})

# ---- Filter and fit ----
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
  q_value = p.adjust(p_value, method = P_ADJUST_METHOD)
) |>
  arrange(p_value)
n_significant <- sum(de_results$q_value < Q_CUTOFF)

# ---- Write tables ----
readr::write_tsv(de_results, file.path(run_dir, "de.tsv"))
command <- paste(
  "Rscript skills/pseudobulk-de/pseudobulk_de.R --input",
  input_dir,
  "--output",
  run_dir
)
writeLines(command, file.path(run_dir, "reproducibility", "commands.sh"))

# jsonlite is not on the allowed package list, so the JSON is written by hand.
result_json <- sprintf(
  paste0(
    '{"parameters": {"input": "%s", "output": "%s", "condition_levels": ["%s", "%s"], ',
    '"method": "edgeR glmQLFit, pseudobulk sum per donor", "p_adjust": "%s", "q_cutoff": %s}, ',
    '"counts": {"donors": %d, "cells": %d, "genes": %d, "genes_tested": %d}, ',
    '"n_significant": %d}'
  ),
  input_dir,
  run_dir,
  CONDITION_LEVELS[[1]],
  CONDITION_LEVELS[[2]],
  P_ADJUST_METHOD,
  Q_CUTOFF,
  nrow(donors),
  nrow(cells),
  nrow(count_matrix),
  nrow(de_results),
  n_significant
)
writeLines(result_json, file.path(run_dir, "result.json"))

# ---- Write report ----
donor_table <- c(
  "| donor | condition | cells |",
  "|---|---|---|",
  sprintf("| %s | %s | %d |", donors$donor_id, donors$condition, donors$n_cells)
)
top_table <- c(
  "| gene | log2_fc | p_value | q_value |",
  "|---|---|---|---|",
  with(
    head(de_results, N_TOP_GENES),
    sprintf("| %s | %.2f | %.2g | %.2g |", gene, log2_fc, p_value, q_value)
  )
)
report <- c(
  "# Pseudobulk differential expression: disease vs healthy",
  "",
  "## Data",
  sprintf(
    "%d donors (%d healthy, %d disease), %d cells, %d genes.",
    nrow(donors),
    sum(donors$condition == "healthy"),
    sum(donors$condition == "disease"),
    nrow(cells),
    nrow(count_matrix)
  ),
  "",
  donor_table,
  "",
  "## Procedure",
  "- **Replication unit:** the donor. Condition varies between donors, not within them,",
  "  so cells of one donor are pseudoreplicates and a cell-level test would overstate n.",
  "- **Aggregation:** raw UMI counts summed over all cells of each donor (pseudobulk).",
  sprintf(
    "- **Filter:** edgeR `filterByExpr` on the donor sums kept %d of %d genes (removed %d).",
    sum(keep_gene),
    length(keep_gene),
    sum(!keep_gene)
  ),
  "- **Normalisation:** TMM library sizes (`normLibSizes`).",
  "- **Model:** edgeR quasi-likelihood negative binomial GLM, `~ condition`. It models",
  "  counts directly and accounts for dispersion estimated across genes, which suits 3 vs 3.",
  "- **Contrast:** `conditiondisease`, so log2_fc is disease minus healthy (log2 CPM scale).",
  sprintf(
    "- **Multiple testing:** %s over the %d genes tested; significant at q < %s.",
    P_ADJUST_METHOD,
    nrow(de_results),
    Q_CUTOFF
  ),
  "",
  "## Results",
  sprintf(
    "%d genes significant at q < %s. Top %d by p-value:",
    n_significant,
    Q_CUTOFF,
    N_TOP_GENES
  ),
  "",
  top_table,
  "",
  "No figure is produced.",
  "",
  "## Command",
  "",
  paste0("`", command, "`")
)
writeLines(report, file.path(run_dir, "report.md"))

# ---- Session info ----
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
