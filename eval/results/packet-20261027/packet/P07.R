#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 3
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 2
#$ q3 change checks (yes/no; would you add, remove or change a check?): no
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Purpose: attach condition and batch to each sample and compute its library size.
# Inputs: data/counts.tsv (raw counts, genes x samples) [sequencing core],
#   data/metadata.tsv (sample_id, condition, batch) [sequencing core]
# Outputs: <run_dir>/sample_summary.tsv, <run_dir>/session_info.txt
# Run: Rscript scripts/02_summarise-samples.R results/run1

library(dplyr)
library(readr)
library(tidyr)
library(forcats)

run_dir <- here::here(commandArgs(trailingOnly = TRUE)[[1]])
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# Load data ----
counts <- read_tsv(here::here("data/counts.tsv"), show_col_types = FALSE)
metadata <- read_tsv(here::here("data/metadata.tsv"), show_col_types = FALSE)

# Summarise samples ----
library_sizes <- counts |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "count") |>
  summarise(library_size = sum(count), .by = sample_id)

# Sample order differs between the two files, so samples are matched by ID. A full join keeps
# samples missing from either file, and the check below stops on them instead of dropping them.
sample_summary <- metadata |>
  full_join(library_sizes, join_by(sample_id), relationship = "one-to-one") |>
  mutate(sample_id = fct(sample_id, levels = sort(sample_id))) |>
  arrange(sample_id) |>
  select(sample_id, condition, batch, library_size)

stopifnot(
  "every sample has metadata and counts" = !anyNA(sample_summary[c("condition", "library_size")])
)

write_tsv(sample_summary, file.path(run_dir, "sample_summary.tsv"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
