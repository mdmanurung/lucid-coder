# Purpose: Attach condition and batch to each sample and report its library size.
# Inputs:  data/counts.tsv   raw counts, genes x samples, first column `gene` [gap: core pipeline]
#          data/metadata.tsv one row per sample: sample_id, condition, batch [sequencing core]
# Outputs: <out_dir>/sample_summary.tsv  sample_id, condition, batch, library_size
#          <out_dir>/session_info.txt
# Run:     Rscript scripts/02_summarise-samples.R results/run1
# Notes:   Count columns are not in metadata order, so samples are matched by ID, never position.

options(warn = 1)

library(dplyr)
library(tidyr)
library(readr)

# Setup ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/02_summarise-samples.R <out_dir>.", call. = FALSE)
}
out_dir <- here::here(args[[1]])
if (dir.exists(out_dir)) {
  stop("out_dir: ", out_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE)

# Load inputs ----------------------------------------------------------------

counts_wide <- read_tsv(here::here("data", "counts.tsv"), show_col_types = FALSE)
metadata <- read_tsv(
  here::here("data", "metadata.tsv"),
  col_types = cols(.default = col_character())
)

missing_meta_cols <- setdiff(c("sample_id", "condition", "batch"), names(metadata))
if (length(missing_meta_cols) > 0) {
  stop(
    "metadata.tsv: missing columns ", toString(missing_meta_cols), "; check the sample sheet.",
    call. = FALSE
  )
}

# Library sizes --------------------------------------------------------------

# Long format first so each count is tied to its sample ID by name.
counts_long <- counts_wide |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "count")

library_sizes <- counts_long |>
  group_by(sample_id) |>
  summarise(library_size = sum(count), .groups = "drop")

# Attach metadata ------------------------------------------------------------

# Every count sample must have exactly one metadata row and vice versa;
# a mismatch would silently drop or mislabel a sample.
unmatched_ids <- union(
  setdiff(library_sizes$sample_id, metadata$sample_id),
  setdiff(metadata$sample_id, library_sizes$sample_id)
)
if (anyDuplicated(metadata$sample_id) > 0 || length(unmatched_ids) > 0) {
  stop(
    "metadata.tsv: sample IDs duplicated or unmatched to counts (", toString(unmatched_ids),
    "); fix the sample sheet.",
    call. = FALSE
  )
}

sample_summary <- metadata |>
  select(sample_id, condition, batch) |>
  inner_join(library_sizes, by = "sample_id", relationship = "one-to-one") |>
  arrange(sample_id)

# Write outputs --------------------------------------------------------------

write_tsv(sample_summary, file.path(out_dir, "sample_summary.tsv"))

# Session info ---------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
