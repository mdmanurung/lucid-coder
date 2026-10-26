# Purpose: drop spike-in and mitochondrial control genes from the raw count table.
# Inputs:  data/counts.tsv   [gap: provided with the task]
# Outputs: <run_dir>/counts_filtered.tsv, <run_dir>/filter_log.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/06_filter-controls.R results/run1
# Notes:   [Stochastic] none. Filtering uses !grepl(), so a control set with zero matches
#          removes nothing instead of emptying the table.

# ---- Setup ----
# Anchored prefixes of the control gene sets, keyed by the name written to the log.
CONTROL_PATTERNS <- c(ERCC = "^ERCC-", MT = "^MT-")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/06_filter-controls.R <run_dir>", call. = FALSE)
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
if (names(counts)[[1]] != "gene" || anyDuplicated(counts$gene)) {
  stop("counts.tsv: first column must be gene with unique values; fix the input.", call. = FALSE)
}

# ---- Filter control genes ----
is_control <- sapply(CONTROL_PATTERNS, grepl, x = counts$gene)
filter_log <- data.frame(gene_set = names(CONTROL_PATTERNS), n_removed = colSums(is_control))
counts_filtered <- counts[!rowSums(is_control), ]
if (nrow(counts_filtered) == 0) {
  stop("counts_filtered: no genes left after removing controls; check the patterns.", call. = FALSE)
}

# ---- Write outputs ----
readr::write_tsv(counts_filtered, file.path(run_dir, "counts_filtered.tsv"))
readr::write_tsv(filter_log, file.path(run_dir, "filter_log.tsv"))

# ---- Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
