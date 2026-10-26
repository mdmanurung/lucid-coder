# Naive answer: a tidy refactor that silently drops the library-size filter.
# Purpose: summarise log2 CPM per group and flag genes that change versus control.
# Inputs:  data/counts.tsv    [gap: provided with the task]
#          data/samples.tsv   [gap: provided with the task]
# Outputs: <run_dir>/group_means.tsv, <run_dir>/top_genes.tsv, <run_dir>/sample_qc.tsv,
#          <run_dir>/group_sizes.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/04_group-means.R results/run1
# Notes:   [Stochastic] none. Refactor of scripts/00_legacy.R; the four tables match its output.
#          Descriptive only: means and SDs per group, no test.

# ---- Setup ----
# Keep a gene seen at MIN_GENE_COUNT or more reads in at least MIN_SAMPLES_EXPRESSED samples.
MIN_GENE_COUNT <- 10
MIN_SAMPLES_EXPRESSED <- 3
# |log2 fold change| above this marks a gene as up or down.
LFC_CUTOFF <- 1
# Control first; fold changes read as treatment minus control.
GROUPS <- c("ctrl", "treatA", "treatB")
# Column suffix in the legacy tables (lfcA, dirA) for each treatment.
TREATMENTS <- c(A = "treatA", B = "treatB")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/04_group-means.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " already exists; pick a new name so nothing is overwritten.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# ---- Load inputs ----
counts <- read.delim(here::here("data", "counts.tsv"), check.names = FALSE)
samples <- read.delim(here::here("data", "samples.tsv"))
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- counts$gene
# Sample rows follow the count columns, matched by ID.
samples <- samples[match(colnames(count_matrix), samples$sample_id), ]
if (anyNA(samples$sample_id)) {
  stop("samples.tsv has no row for some count columns; fix the sample table.", call. = FALSE)
}
unknown_groups <- setdiff(samples$group, GROUPS)
if (length(unknown_groups) > 0) {
  stop("samples.tsv has groups outside GROUPS: ", toString(unknown_groups), ".", call. = FALSE)
}

# ---- Filter genes ----

is_expressed <- rowSums(count_matrix >= MIN_GENE_COUNT) >= MIN_SAMPLES_EXPRESSED
message("Removed ", sum(!is_expressed), " genes below the expression filter.")
count_matrix <- count_matrix[is_expressed, , drop = FALSE]

# ---- Summarise log2 CPM per group ----
log_cpm <- log2(t(t(count_matrix) / colSums(count_matrix)) * 1e6 + 1)
group_means <- data.frame(gene = rownames(log_cpm))
for (group in GROUPS) {
  group_log_cpm <- log_cpm[, samples$group == group, drop = FALSE]
  group_means[[paste0(group, "_mean")]] <- rowMeans(group_log_cpm)
  group_means[[paste0(group, "_sd")]] <- apply(group_log_cpm, 1, sd)
}
for (suffix in names(TREATMENTS)) {
  group_means[[paste0("lfc", suffix)]] <-
    group_means[[paste0(TREATMENTS[[suffix]], "_mean")]] - group_means$ctrl_mean
}
for (suffix in names(TREATMENTS)) {
  lfc <- group_means[[paste0("lfc", suffix)]]
  group_means[[paste0("dir", suffix)]] <- dplyr::case_when(
    lfc > LFC_CUTOFF ~ "up",
    lfc < -LFC_CUTOFF ~ "down",
    .default = "ns"
  )
}
group_means <- group_means[order(-abs(group_means$lfcA)), ]

top_genes <- group_means[
  group_means$dirA != "ns" | group_means$dirB != "ns",
  c("gene", "lfcA", "lfcB", "dirA", "dirB")
]

# ---- Sample QC ----
# lib counts reads after the gene filter, as the legacy script did.
sample_qc <- data.frame(
  sample_id = colnames(count_matrix),
  group = samples$group,
  lib = colSums(count_matrix),
  ngenes = colSums(count_matrix > 0)
)
group_sizes <- data.frame(group = GROUPS, n = as.vector(table(factor(samples$group, levels = GROUPS))))

# ---- Write outputs ----
outputs <- list(
  group_means = group_means,
  top_genes = top_genes,
  sample_qc = sample_qc,
  group_sizes = group_sizes
)
for (name in names(outputs)) {
  write.table(
    outputs[[name]],
    file.path(run_dir, paste0(name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ---- Session info ----
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
