#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 4
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 5
#$ q3 change checks (yes/no; would you add, remove or change a check?): no
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Per-group mean/SD of log2-CPM expression and log fold changes vs ctrl.
# Refactor of scripts/00_legacy.R; outputs are identical to the legacy script.
#
# Usage: Rscript scripts/04_group-means.R <run_dir>
#
# Writes to <run_dir>:
#   group_means.tsv  per-gene group mean/SD, lfcA/lfcB, direction calls
#   top_genes.tsv    genes called up/down in treatA or treatB
#   sample_qc.tsv    retained samples with library size and detected genes
#   group_sizes.tsv  number of retained samples per group

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(here)
})

# ---- Parameters ---------------------------------------------------------------

min_lib_size <- 20000 # samples need a total count above this
min_count <- 10 # a gene counts as expressed in a sample at >= this
min_samples <- 3 # ... and must be expressed in at least this many samples
lfc_cutoff <- 1 # |log2 fold change| above this is called up/down
groups <- c("ctrl", "treatA", "treatB")
reference <- "ctrl"

# ---- Arguments ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/04_group-means.R <run_dir>", call. = FALSE)
}
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Base write.table keeps the legacy number formatting (15 significant digits).
write_out <- function(df, name) {
  write.table(
    df,
    file.path(out_dir, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ---- Read ---------------------------------------------------------------------

counts_df <- read.delim(here("data", "counts.tsv"), check.names = FALSE)
samples <- read.delim(here("data", "samples.tsv"))

counts <- as.matrix(counts_df[, -1])
rownames(counts) <- counts_df$gene
samples <- samples[match(colnames(counts), samples$sample_id), ]

# ---- Filter -------------------------------------------------------------------

keep_samples <- colSums(counts) > min_lib_size
counts <- counts[, keep_samples, drop = FALSE]
samples <- samples[keep_samples, ]

counts <- counts[rowSums(counts >= min_count) >= min_samples, , drop = FALSE]

# ---- Normalise ----------------------------------------------------------------

cpm <- t(t(counts) / colSums(counts)) * 1e6
log_cpm <- log2(cpm + 1)

# ---- Group summaries ----------------------------------------------------------

group_stats <- function(g) {
  m <- log_cpm[, samples$group == g, drop = FALSE]
  tibble(
    "{g}_mean" := rowMeans(m),
    "{g}_sd" := apply(m, 1, sd)
  )
}

call_direction <- function(lfc) {
  case_when(
    lfc > lfc_cutoff ~ "up",
    lfc < -lfc_cutoff ~ "down",
    .default = "ns"
  )
}

group_means <- bind_cols(
  tibble(gene = rownames(log_cpm)),
  map(groups, group_stats)
) |>
  mutate(
    lfcA = treatA_mean - .data[[paste0(reference, "_mean")]],
    lfcB = treatB_mean - .data[[paste0(reference, "_mean")]],
    dirA = call_direction(lfcA),
    dirB = call_direction(lfcB)
  ) |>
  arrange(desc(abs(lfcA)))

top_genes <- group_means |>
  filter(dirA != "ns" | dirB != "ns") |>
  select(gene, lfcA, lfcB, dirA, dirB)

# ---- Sample QC ----------------------------------------------------------------

sample_qc <- tibble(
  sample_id = colnames(counts),
  group = samples$group,
  lib = colSums(counts),
  ngenes = colSums(counts > 0)
)

group_sizes <- tibble(
  group = groups,
  n = map_int(groups, \(g) sum(samples$group == g))
)

# ---- Write --------------------------------------------------------------------

write_out(group_means, "group_means.tsv")
write_out(top_genes, "top_genes.tsv")
write_out(sample_qc, "sample_qc.tsv")
write_out(group_sizes, "group_sizes.tsv")
