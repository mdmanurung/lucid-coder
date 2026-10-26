#!/usr/bin/env Rscript
# Attach condition/batch to each sample and compute library size.
# Usage: Rscript scripts/02_summarise-samples.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/02_summarise-samples.R <run_dir>")
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_tsv(here("data", "counts.tsv"), col_types = cols(gene = col_character(), .default = col_double()))
meta <- read_tsv(here("data", "metadata.tsv"), col_types = cols(.default = col_character())) %>%
  mutate(across(everything(), trimws))

# Count columns are in a different order than metadata rows: join by ID, never by position.
count_ids <- trimws(setdiff(names(counts), "gene"))
stopifnot(
  !anyDuplicated(count_ids), !anyDuplicated(meta$sample_id),
  setequal(count_ids, meta$sample_id),
  !anyNA(counts[-1])
)

summary_tbl <- counts %>%
  rename_with(trimws) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "count") %>%
  group_by(sample_id) %>%
  summarise(library_size = sum(count), .groups = "drop") %>%
  inner_join(meta, by = "sample_id", relationship = "one-to-one") %>%
  select(sample_id, condition, batch, library_size) %>%
  arrange(sample_id)

stopifnot(nrow(summary_tbl) == length(count_ids))
write_tsv(summary_tbl, file.path(run_dir, "sample_summary.tsv"))
