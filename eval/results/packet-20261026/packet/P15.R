#!/usr/bin/env Rscript
# Attach condition and batch to each sample and compute library sizes.
# Usage: Rscript scripts/02_summarise-samples.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/02_summarise-samples.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_tsv(here("data", "counts.tsv"),
                   col_types = cols(gene = col_character(), .default = col_double()))
metadata <- read_tsv(here("data", "metadata.tsv"),
                     col_types = cols(.default = col_character())) %>%
  mutate(across(everything(), trimws))

# Sanity checks: sample IDs must be unique and match 1:1 between files.
count_ids <- setdiff(names(counts), "gene")
stopifnot(
  !anyDuplicated(count_ids),
  !anyDuplicated(metadata$sample_id),
  !anyDuplicated(counts$gene)
)
missing_meta <- setdiff(count_ids, metadata$sample_id)
missing_counts <- setdiff(metadata$sample_id, count_ids)
if (length(missing_meta) || length(missing_counts)) {
  stop("Sample ID mismatch. In counts only: ", paste(missing_meta, collapse = ", "),
       "; in metadata only: ", paste(missing_counts, collapse = ", "), call. = FALSE)
}
if (any(is.na(counts[count_ids]))) stop("counts.tsv contains missing values", call. = FALSE)
if (any(counts[count_ids] < 0)) stop("counts.tsv contains negative values", call. = FALSE)

# Join by sample ID (column order in counts differs from metadata row order).
library_sizes <- counts %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "count") %>%
  group_by(sample_id) %>%
  summarise(library_size = sum(count), .groups = "drop")

sample_summary <- metadata %>%
  select(sample_id, condition, batch) %>%
  inner_join(library_sizes, by = "sample_id", relationship = "one-to-one") %>%
  arrange(sample_id)

stopifnot(nrow(sample_summary) == length(count_ids))

out_file <- file.path(run_dir, "sample_summary.tsv")
write_tsv(sample_summary, out_file)
message("Wrote ", nrow(sample_summary), " samples to ", out_file)
