# Purpose: Per-group mean and SD of log2(CPM + 1) per gene, log2 fold changes of each treatment
#   versus ctrl, and per-sample QC. Refactor of scripts/00_legacy.R with identical outputs.
# Inputs:
#   data/counts.tsv  - raw counts: `gene`, then one column per sample [gap: producer unknown]
#   data/samples.tsv - `sample_id`, `group` [gap: producer unknown]
# Outputs (in <out_dir>):
#   group_means.tsv - one row per kept gene, sorted by |lfcA| descending
#   top_genes.tsv   - genes with |lfc| > LFC_CUTOFF in either treatment
#   sample_qc.tsv   - library size and detected genes per kept sample (after gene filter)
#   group_sizes.tsv - number of kept samples per group
# Run: Rscript scripts/04_group-means.R results/<run_dir>
# Notes:
#   - Library sizes in sample_qc.tsv and CPM denominators use counts after the gene filter,
#     as in the legacy script.
#   - Tables are written with write.table() so numbers keep the legacy 15-significant-digit text.

options(warn = 1)

# Setup ----------------------------------------------------------------------------------------

library(dplyr)
library(purrr)

MIN_LIB_SIZE <- 20000 # samples at or below this total count are dropped (legacy cutoff)
MIN_GENE_COUNT <- 10 # a gene counts as expressed in a sample at >= 10 reads
MIN_EXPRESSED_SAMPLES <- 3 # keep genes expressed in at least 3 kept samples
LFC_CUTOFF <- 1 # |log2 fold change| above this is called up/down
REFERENCE_GROUP <- "ctrl"
GROUP_LEVELS <- c("ctrl", "treatA", "treatB")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("args: expected one run directory; run `Rscript scripts/04_group-means.R <out_dir>`.",
    call. = FALSE
  )
}
out_dir <- here::here(args[[1]])
if (dir.exists(out_dir)) {
  stop("out_dir: ", out_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE)

# Load inputs ----------------------------------------------------------------------------------

counts_table <- read.delim(
  here::here("data", "counts.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
samples <- read.delim(here::here("data", "samples.tsv"), stringsAsFactors = FALSE)
if (!all(c("sample_id", "group") %in% names(samples))) {
  stop("samples.tsv: needs columns sample_id and group; fix the sample sheet.", call. = FALSE)
}

counts <- as.matrix(counts_table[, setdiff(names(counts_table), "gene")])
rownames(counts) <- counts_table$gene

# Align sample sheet to count columns by ID, not position
samples <- samples[match(colnames(counts), samples$sample_id), ]
if (!identical(samples$sample_id, colnames(counts))) {
  stop("samples.tsv: sample IDs do not match count columns one-to-one; fix the sheet.",
    call. = FALSE
  )
}
if (!all(samples$group %in% GROUP_LEVELS)) {
  stop("samples.tsv: group outside ", toString(GROUP_LEVELS), "; update GROUP_LEVELS.",
    call. = FALSE
  )
}

# Filter samples by library size ---------------------------------------------------------------

keep_sample <- colSums(counts) > MIN_LIB_SIZE
message(
  "Removed ", sum(!keep_sample), " sample(s) with library size <= ", MIN_LIB_SIZE, ": ",
  toString(colnames(counts)[!keep_sample])
)
counts <- counts[, keep_sample, drop = FALSE]
samples <- samples[keep_sample, ]

# Filter lowly expressed genes -----------------------------------------------------------------

keep_gene <- rowSums(counts >= MIN_GENE_COUNT) >= MIN_EXPRESSED_SAMPLES
message(
  "Removed ", sum(!keep_gene), " of ", length(keep_gene), " gene(s) with >= ", MIN_GENE_COUNT,
  " counts in fewer than ", MIN_EXPRESSED_SAMPLES, " samples"
)
counts <- counts[keep_gene, , drop = FALSE]

# Normalise ------------------------------------------------------------------------------------

cpm <- t(t(counts) / colSums(counts)) * 1e6
log_cpm <- log2(cpm + 1)

# Group summaries ------------------------------------------------------------------------------

# rowMeans() and apply(sd) on the matrix, as in the legacy script, so values match to the last bit
group_stats <- map(GROUP_LEVELS, function(group) {
  group_log_cpm <- log_cpm[, samples$group == group, drop = FALSE]
  stats <- data.frame(rowMeans(group_log_cpm), apply(group_log_cpm, 1, sd))
  names(stats) <- paste0(group, c("_mean", "_sd"))
  stats
})

call_direction <- function(lfc, cutoff) {
  case_when(lfc > cutoff ~ "up", lfc < -cutoff ~ "down", .default = "ns")
}

group_means <- bind_cols(data.frame(gene = rownames(log_cpm)), group_stats) |>
  mutate(
    lfcA = treatA_mean - .data[[paste0(REFERENCE_GROUP, "_mean")]],
    lfcB = treatB_mean - .data[[paste0(REFERENCE_GROUP, "_mean")]],
    dirA = call_direction(lfcA, LFC_CUTOFF),
    dirB = call_direction(lfcB, LFC_CUTOFF)
  ) |>
  arrange(desc(abs(lfcA)))

top_genes <- group_means |>
  filter(dirA != "ns" | dirB != "ns") |>
  select(gene, lfcA, lfcB, dirA, dirB)

sample_qc <- data.frame(
  sample_id = colnames(counts),
  group = samples$group,
  lib = colSums(counts),
  ngenes = colSums(counts > 0)
)

group_sizes <- data.frame(
  group = GROUP_LEVELS,
  n = map_int(GROUP_LEVELS, \(group) sum(samples$group == group))
)

# Write outputs --------------------------------------------------------------------------------

output_tables <- list(
  group_means = group_means,
  top_genes = top_genes,
  sample_qc = sample_qc,
  group_sizes = group_sizes
)
iwalk(output_tables, function(table, name) {
  if (nrow(table) == 0) {
    warning(name, ": table is empty.", call. = FALSE)
  }
  write.table(
    table, file.path(out_dir, paste0(name, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
})

# Session info ---------------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
