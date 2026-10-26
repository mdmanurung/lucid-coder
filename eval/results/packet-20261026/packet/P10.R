# Purpose: Attach condition and batch to each sample and compute its library size.
# Inputs:  data/counts.tsv   raw counts, genes x samples, first column `gene` [sequencing core]
#          data/metadata.tsv one row per sample: sample_id, condition, batch [sequencing core]
# Outputs: <run_dir>/sample_summary.tsv  sample_id, condition, batch, library_size
#          <run_dir>/session_info.txt
# Run:     Rscript scripts/02_summarise-samples.R results/run1
# Notes:   Count columns are not in the same order as the metadata rows, so samples are
#          matched by ID with a join, never by position.

# Setup ------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/02_summarise-samples.R <run_dir>", call. = FALSE)
}
run_dir <- here::here(args[[1]])
if (dir.exists(run_dir)) {
  stop("run_dir: ", run_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Load inputs ------------------------------------------------------------------

counts_wide <- read_tsv(here::here("data", "counts.tsv"), show_col_types = FALSE)
metadata <- read_tsv(here::here("data", "metadata.tsv"), show_col_types = FALSE)

missing_cols <- setdiff(c("sample_id", "condition", "batch"), names(metadata))
if (length(missing_cols) > 0) {
  stop(
    "metadata: missing column(s) ", toString(missing_cols), "; check the sample sheet.",
    call. = FALSE
  )
}

# Library size per sample ------------------------------------------------------

library_sizes <- counts_wide |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "count") |>
  summarise(library_size = sum(count), .by = sample_id)

# Attach metadata --------------------------------------------------------------

# Every sample must appear exactly once on both sides, otherwise rows are lost or mislabelled.
count_ids <- library_sizes$sample_id
meta_ids <- metadata$sample_id
if (anyDuplicated(meta_ids) > 0 || !setequal(count_ids, meta_ids)) {
  stop(
    "metadata: sample IDs differ from counts columns or repeat (counts only: ",
    toString(setdiff(count_ids, meta_ids)), "; metadata only: ",
    toString(setdiff(meta_ids, count_ids)), "); fix the sample sheet.",
    call. = FALSE
  )
}

sample_summary <- metadata |>
  select(sample_id, condition, batch) |>
  inner_join(library_sizes, by = "sample_id", relationship = "one-to-one") |>
  arrange(sample_id)

# Write outputs ----------------------------------------------------------------

write_tsv(sample_summary, file.path(run_dir, "sample_summary.tsv"))

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
