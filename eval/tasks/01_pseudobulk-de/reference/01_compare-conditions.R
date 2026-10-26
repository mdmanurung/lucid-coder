# Purpose: test each gene for a disease-vs-healthy difference, with the donor as the unit.
# Inputs:  data/counts.tsv   [gap: provided with the task]
#          data/cells.tsv    [gap: provided with the task]
# Outputs: <run_dir>/de.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/01_compare-conditions.R results/run1
# Notes:   [Stochastic] none. Inference unit: donor. Cells from one donor are not independent,
#          so counts are summed per donor (pseudobulk) before testing 3 vs 3 donors.

# ---- Setup ----
# Reference = healthy, so log2_fc reads as disease minus healthy.
CONDITION_LEVELS <- c("healthy", "disease")
# Multiple-testing family: every gene in counts.tsv, adjusted together.
P_ADJUST_METHOD <- "BH"
# Pseudocount 1 on CPM keeps log2 finite for a gene with zero counts in a donor.
CPM_PSEUDOCOUNT <- 1
REQUIRED_CELL_COLS <- c("cell_id", "donor_id", "condition")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/01_compare-conditions.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
if (dir.exists(run_dir)) {
  stop(
    "run_dir ",
    run_dir,
    " already exists; pick a new name so nothing is overwritten.",
    call. = FALSE
  )
}
dir.create(run_dir, recursive = TRUE)

# ---- Load inputs ----
counts <- read.delim(here::here("data", "counts.tsv"), check.names = FALSE)
cells <- read.delim(here::here("data", "cells.tsv"))
missing_cols <- setdiff(REQUIRED_CELL_COLS, names(cells))
if (length(missing_cols) > 0) {
  stop("cells.tsv lacks columns ", toString(missing_cols), ".", call. = FALSE)
}
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts$gene
# Cell order in cells.tsv is not guaranteed to match the count columns.
cells <- cells[match(colnames(count_matrix), cells$cell_id), ]
if (anyNA(cells$cell_id)) {
  stop("some count columns have no row in cells.tsv; fix the cell table.", call. = FALSE)
}

# ---- Sum counts per donor ----
donors <- unique(cells[c("donor_id", "condition")])
if (anyDuplicated(donors$donor_id)) {
  stop("a donor has cells in more than one condition; fix cells.tsv.", call. = FALSE)
}
donors$condition <- factor(donors$condition, levels = CONDITION_LEVELS)
if (anyNA(donors$condition)) {
  stop("cells.tsv has conditions outside CONDITION_LEVELS; recode them at load.", call. = FALSE)
}
donor_counts <- sapply(donors$donor_id, function(donor) {
  rowSums(count_matrix[, cells$donor_id == donor, drop = FALSE])
})
donor_log_cpm <- log2(t(t(donor_counts) / colSums(donor_counts)) * 1e6 + CPM_PSEUDOCOUNT)

# ---- Test each gene ----
# Linear model on donor log-CPM: one observation per donor, so df = donors - 2.
fit_gene <- function(log_cpm) {
  fit <- lm(log_cpm ~ condition, data = donors)
  disease_term <- coef(summary(fit))["conditiondisease", ]
  c(log2_fc = disease_term[["Estimate"]], p_value = disease_term[["Pr(>|t|)"]])
}
de_results <- data.frame(gene = rownames(donor_log_cpm), t(apply(donor_log_cpm, 1, fit_gene)))
de_results$q_value <- p.adjust(de_results$p_value, method = P_ADJUST_METHOD)
de_results <- de_results[order(de_results$p_value), ]

# ---- Write outputs ----
write.table(de_results, file.path(run_dir, "de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
