# Purpose: attach condition and batch to each sample and report its library size.
# Inputs:  data/counts.tsv     [gap: provided with the task]
#          data/metadata.tsv   [gap: sample sheet from the sequencing core]
# Outputs: <run_dir>/sample_summary.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/02_summarise-samples.R results/run1
# Notes:   [Stochastic] none. Count columns and sheet rows are in different orders,
#          so samples are joined by sample_id, never by position.

# ---- Setup ----
library(dplyr)

REQUIRED_METADATA_COLS <- c("sample_id", "condition", "batch")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/02_summarise-samples.R <run_dir>", call. = FALSE)
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
metadata <- read.delim(here::here("data", "metadata.tsv"))
missing_cols <- setdiff(REQUIRED_METADATA_COLS, names(metadata))
if (length(missing_cols) > 0) {
  stop("metadata.tsv lacks columns ", toString(missing_cols), ".", call. = FALSE)
}

# ---- Join metadata by sample_id ----
library_sizes <- data.frame(
  sample_id = names(counts)[-1],
  library_size = colSums(counts[, -1])
)
# one-to-one: a duplicated sample_id on either side stops the join.
sample_summary <- library_sizes |>
  left_join(metadata, by = "sample_id", relationship = "one-to-one") |>
  select(sample_id, condition, batch, library_size)
unmatched <- sample_summary$sample_id[is.na(sample_summary$condition)]
if (length(unmatched) > 0) {
  stop(
    "no metadata row for samples ",
    toString(unmatched),
    "; fix the sample sheet.",
    call. = FALSE
  )
}

# ---- Write outputs ----
write.table(
  sample_summary,
  file.path(run_dir, "sample_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
