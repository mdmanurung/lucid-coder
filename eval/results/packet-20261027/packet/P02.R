#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 3
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 2
#$ q3 change checks (yes/no; would you add, remove or change a check?): yes
#$ q3 which (one sentence): see comment below about checks
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Attach condition and batch to each sample and compute library sizes.
#
# Usage: Rscript scripts/02_summarise-samples.R <run_dir>
# Writes <run_dir>/sample_summary.tsv

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

counts <- read_tsv(
  here("data", "counts.tsv"),
  col_types = cols(gene = col_character(), .default = col_double())
)
metadata <- read_tsv(here("data", "metadata.tsv"), col_types = cols(.default = col_character()))

# Samples are matched by ID: the count columns and metadata rows are in different orders.
count_ids <- setdiff(names(counts), "gene")
stopifnot( #$ not necessary
  !anyDuplicated(count_ids),
  !anyDuplicated(metadata$sample_id),
  !anyNA(counts)
)
missing_meta <- setdiff(count_ids, metadata$sample_id)
missing_counts <- setdiff(metadata$sample_id, count_ids)
if (length(missing_meta) > 0 || length(missing_counts) > 0) { #$ not necessary
  stop(
    "Sample IDs differ between counts and metadata. ",
    "No metadata: ", toString(missing_meta), "; no counts: ", toString(missing_counts),
    call. = FALSE
  )
}

library_sizes <- counts |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "count") |>
  group_by(sample_id) |>
  summarise(library_size = sum(count), .groups = "drop")

sample_summary <- metadata |>
  select(sample_id, condition, batch) |>
  inner_join(library_sizes, by = "sample_id", relationship = "one-to-one") |>
  arrange(sample_id)

stopifnot(nrow(sample_summary) == length(count_ids)) #$ ok for sanity check

write_tsv(sample_summary, file.path(run_dir, "sample_summary.tsv"))
