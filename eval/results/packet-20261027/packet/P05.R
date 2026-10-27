#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Purpose: per-group mean and SD of log2(CPM + 1), log fold changes vs ctrl, and sample QC.
#   Refactor of scripts/00_legacy.R; outputs are identical to it.
# Inputs: data/counts.tsv (raw counts, gene x sample) [gap], data/samples.tsv (sample_id, group) [gap]
# Outputs: <run_dir>/group_means.tsv, top_genes.tsv, sample_qc.tsv, group_sizes.tsv, session_info.txt
# Run: Rscript scripts/04_group-means.R results/run1

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)
library(here)

min_library_size <- 20000
min_gene_count <- 10
min_samples_expressed <- 3
lfc_cutoff <- 1
group_levels <- c("ctrl", "treatA", "treatB")

run_dir <- here(commandArgs(trailingOnly = TRUE)[[1]])
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

# write.table prints 15 significant digits, as the legacy tables do; write_tsv prints more
write_legacy_tsv <- function(df, file) {
  write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Load data ----
counts_raw <- read_tsv(here("data/counts.tsv"), show_col_types = FALSE)
counts <- as.matrix(select(counts_raw, -gene))
rownames(counts) <- counts_raw[["gene"]]

samples <- tibble(sample_id = colnames(counts)) |>
  left_join(
    read_tsv(here("data/samples.tsv"), show_col_types = FALSE),
    by = join_by(sample_id),
    relationship = "one-to-one"
  ) |>
  mutate(group = fct(group, levels = group_levels))

# Filter samples and genes ----
keep_samples <- colSums(counts) > min_library_size
counts <- counts[, keep_samples]
samples <- samples[keep_samples, ]
message(sprintf(
  "Samples: removed %d with library size <= %d",
  sum(!keep_samples),
  min_library_size
))

keep_genes <- rowSums(counts >= min_gene_count) >= min_samples_expressed
counts <- counts[keep_genes, ]
message(sprintf(
  "Genes: removed %d with count >= %d in fewer than %d samples",
  sum(!keep_genes),
  min_gene_count,
  min_samples_expressed
))

stopifnot("samples align with count columns" = identical(samples[["sample_id"]], colnames(counts)))

log_cpm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)

# Group means and fold changes ----
group_stats <- map(group_levels, \(g) {
  group_log_cpm <- log_cpm[, samples[["group"]] == g, drop = FALSE]
  tibble(
    "{g}_mean" := rowMeans(group_log_cpm),
    "{g}_sd" := apply(group_log_cpm, 1, sd)
  )
}) |>
  list_cbind()

lfc_direction <- function(lfc) {
  case_when(
    lfc > lfc_cutoff ~ "up",
    lfc < -lfc_cutoff ~ "down",
    .default = "ns"
  )
}

group_means <- bind_cols(gene = rownames(log_cpm), group_stats) |>
  mutate(
    lfcA = treatA_mean - ctrl_mean,
    lfcB = treatB_mean - ctrl_mean,
    dirA = lfc_direction(lfcA),
    dirB = lfc_direction(lfcB)
  ) |>
  arrange(desc(abs(lfcA)))
write_legacy_tsv(group_means, file.path(run_dir, "group_means.tsv"))

group_means |>
  filter(dirA != "ns" | dirB != "ns") |>
  select(gene, lfcA, lfcB, dirA, dirB) |>
  write_legacy_tsv(file.path(run_dir, "top_genes.tsv"))

# Sample QC ----
# Library size and detected genes are counted after gene filtering, as in the legacy script
samples |>
  mutate(
    lib = colSums(counts),
    ngenes = colSums(counts > 0)
  ) |>
  write_legacy_tsv(file.path(run_dir, "sample_qc.tsv"))

samples |>
  count(group, .drop = FALSE) |>
  write_legacy_tsv(file.path(run_dir, "group_sizes.tsv"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
