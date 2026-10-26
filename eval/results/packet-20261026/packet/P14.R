# Purpose: Per-group mean and SD of log2(CPM + 1) per gene, log2 fold changes of each
#          treatment versus ctrl, and per-sample QC. Refactor of scripts/00_legacy.R with
#          identical outputs.
# Inputs:  data/counts.tsv  [gap: producer unknown] raw counts; `gene`, then one column per sample
#          data/samples.tsv [gap: producer unknown] `sample_id`, `group`
# Outputs: <run_dir>/group_means.tsv  one row per kept gene, sorted by |lfcA| descending
#          <run_dir>/top_genes.tsv    genes with |lfc| > 1 in either treatment
#          <run_dir>/sample_qc.tsv    kept samples: library size and detected genes
#          <run_dir>/group_sizes.tsv  kept samples per group
#          <run_dir>/session_info.txt
# Run:     Rscript scripts/04_group-means.R results/<run_dir>
# Notes:   To match the legacy tables byte for byte:
#          - CPM library sizes and the sample_qc `lib` and `ngenes` columns are computed
#            AFTER gene filtering, as in the legacy script, not on the raw library.
#          - Tables are written with utils::write.table() (15 significant digits), not
#            readr::write_tsv(), which prints doubles with more digits.
#          - Means use rowMeans() and SDs apply(, 1, sd) so floating-point results match.

# Setup ------------------------------------------------------------------------

library(dplyr)
library(purrr)

MIN_LIBRARY_SIZE <- 20000 # samples with total counts <= this are dropped (legacy cutoff)
MIN_GENE_COUNT <- 10 # a gene counts as expressed in a sample at >= this many reads
MIN_EXPRESSED_SAMPLES <- 3 # keep genes expressed in at least this many samples
LFC_CUTOFF <- 1 # |log2 fold change| above this is called "up" / "down"
CPM_SCALE <- 1e6
PSEUDOCOUNT <- 1 # added to CPM before log2
REFERENCE_GROUP <- "ctrl" # fold changes are treatment minus this group
GROUPS <- c("ctrl", "treatA", "treatB") # column and row order in the outputs

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Arguments: expected one run directory; run `Rscript scripts/04_group-means.R <dir>`.",
    call. = FALSE
  )
}
run_dir <- here::here(args[[1]])
if (dir.exists(run_dir)) {
  stop("run_dir: ", run_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

write_table <- function(table, file_name) {
  utils::write.table(table, file.path(run_dir, file_name),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# Load inputs ------------------------------------------------------------------

counts_table <- utils::read.delim(here::here("data", "counts.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
samples <- utils::read.delim(here::here("data", "samples.tsv"), stringsAsFactors = FALSE)

missing_samples <- setdiff(colnames(counts_table)[-1], samples$sample_id)
if (length(missing_samples) > 0) {
  stop("samples.tsv: no group for ", toString(missing_samples), "; add them.", call. = FALSE)
}

counts <- as.matrix(counts_table[, colnames(counts_table) != "gene"])
rownames(counts) <- counts_table$gene
sample_groups <- samples$group[match(colnames(counts), samples$sample_id)]

# Filter samples by library size -----------------------------------------------

keep_samples <- colSums(counts) > MIN_LIBRARY_SIZE
message(
  "Removed ", sum(!keep_samples), " sample(s) with library size <= ", MIN_LIBRARY_SIZE,
  ": ", toString(colnames(counts)[!keep_samples])
)
counts <- counts[, keep_samples, drop = FALSE]
sample_groups <- sample_groups[keep_samples]

# Filter lowly expressed genes -------------------------------------------------

keep_genes <- rowSums(counts >= MIN_GENE_COUNT) >= MIN_EXPRESSED_SAMPLES
message(
  "Removed ", sum(!keep_genes), " of ", length(keep_genes), " gene(s) with >= ",
  MIN_GENE_COUNT, " counts in fewer than ", MIN_EXPRESSED_SAMPLES, " samples"
)
counts <- counts[keep_genes, , drop = FALSE]

# Normalise --------------------------------------------------------------------

library_sizes <- colSums(counts)
log_cpm <- log2(t(t(counts) / library_sizes) * CPM_SCALE + PSEUDOCOUNT)

# Group means and fold changes -------------------------------------------------

group_stats <- map(GROUPS, function(group) {
  group_log_cpm <- log_cpm[, sample_groups == group, drop = FALSE]
  tibble(
    "{group}_mean" := rowMeans(group_log_cpm),
    "{group}_sd" := apply(group_log_cpm, 1, sd)
  )
})

call_direction <- function(lfc) {
  case_when(
    lfc > LFC_CUTOFF ~ "up",
    lfc < -LFC_CUTOFF ~ "down",
    .default = "ns"
  )
}

ref_mean <- paste0(REFERENCE_GROUP, "_mean")
group_means <- bind_cols(tibble(gene = rownames(log_cpm)), group_stats) |>
  mutate(
    lfcA = treatA_mean - .data[[ref_mean]],
    lfcB = treatB_mean - .data[[ref_mean]],
    dirA = call_direction(lfcA),
    dirB = call_direction(lfcB)
  ) |>
  arrange(desc(abs(lfcA)))

top_genes <- group_means |>
  filter(dirA != "ns" | dirB != "ns") |>
  select(gene, lfcA, lfcB, dirA, dirB)

# Sample QC --------------------------------------------------------------------

sample_qc <- data.frame(
  sample_id = colnames(counts),
  group = sample_groups,
  lib = library_sizes,
  ngenes = colSums(counts > 0)
)

group_sizes <- data.frame(
  group = GROUPS,
  n = map_int(GROUPS, \(group) sum(sample_groups == group))
)

# Write outputs ----------------------------------------------------------------

if (nrow(group_means) == 0) {
  stop("group_means: no genes passed filtering; check cutoffs.", call. = FALSE)
}
write_table(group_means, "group_means.tsv")
write_table(top_genes, "top_genes.tsv")
write_table(sample_qc, "sample_qc.tsv")
write_table(group_sizes, "group_sizes.tsv")

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
