# Purpose: freeze the task 12 inputs and truth from three nf-core scdownstream PBMC samples:
#          plant low-quality cells in every sample and thin one sample to plant a depth gap.
# Inputs:  <download_dir>/SRR28679756_filtered_matrix.h5ad      [nf-core/test-datasets, branch
#          <download_dir>/SRR28679757_filtered_matrix.sce.rds    scdownstream, pbmc/; MIT]
#          <download_dir>/SRR28679759_filtered_matrix.h5ad
# Outputs: eval/tasks/12_per-sample-qc/inputs/data/<sample>/{matrix.mtx,genes.tsv,barcodes.tsv}
#          eval/tasks/12_per-sample-qc/truth/cells.tsv   (sample, barcode, planted)
# Run:     Rscript eval/tasks/12_per-sample-qc/make-inputs.R <download_dir>
# Notes:   [Stochastic] one seed. Not wired into eval/fixtures/make-fixtures.R; rerun by hand.
#          Each sample is 384 cells x 3000 genes. LOW_QUALITY_PER_SAMPLE cells per sample are
#          binomially thinned to about LOW_QUALITY_FRACTION of the sample's median UMIs; then
#          every cell of THIN_SAMPLE is thinned to THIN_FRACTION. A global cutoff then either
#          keeps the deep samples' planted cells or drops most of THIN_SAMPLE; per-sample MAD
#          does neither.

# ---- Setup ----
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

SEED <- 12
LOW_QUALITY_PER_SAMPLE <- 24
# Planted cells are thinned to 15% of their sample's median UMIs: an outlier within the sample,
# yet deeper than most cells of THIN_SAMPLE, so no single cutoff separates both.
LOW_QUALITY_FRACTION <- 0.15
# The thinned sample keeps 10% of its counts: its median falls to about 14 UMIs.
THIN_SAMPLE <- "SRR28679757"
THIN_FRACTION <- 0.1

download_dir <- commandArgs(trailingOnly = TRUE)[[1]]
task_dir <- here::here("eval", "tasks", "12_per-sample-qc")

# ---- Load samples ----
# h5ad X is CSC (cells x genes): indptr runs over genes, indices are cells.
read_h5ad_counts <- function(path) {
  x <- rhdf5::h5read(path, "X")
  cells <- sub(".*::", "", as.character(rhdf5::h5read(path, "obs/_index")))
  genes <- as.character(rhdf5::h5read(path, "var/_index"))
  counts <- sparseMatrix(
    i = as.integer(x$indices) + 1L,
    p = as.integer(x$indptr),
    x = as.numeric(x$data),
    dims = c(length(cells), length(genes)),
    dimnames = list(cells, genes)
  )
  t(counts)
}
samples <- list(
  SRR28679756 = read_h5ad_counts(file.path(
    download_dir,
    "SRR28679756_filtered_matrix.h5ad"
  )),
  SRR28679757 = as(
    assay(
      readRDS(file.path(download_dir, "SRR28679757_filtered_matrix.sce.rds")),
      "X"
    ),
    "CsparseMatrix"
  ),
  SRR28679759 = read_h5ad_counts(file.path(
    download_dir,
    "SRR28679759_filtered_matrix.h5ad"
  ))
)
stopifnot(
  "samples share genes in one order" = identical(
    rownames(samples$SRR28679756),
    rownames(samples$SRR28679757)
  ) &&
    identical(rownames(samples$SRR28679756), rownames(samples$SRR28679759))
)

# ---- Plant low-quality cells and the depth gap ----
thin <- function(counts, fraction) {
  counts@x <- as.numeric(rbinom(length(counts@x), counts@x, fraction))
  drop0(counts)
}
set.seed(SEED)
truth <- list()
for (sample_id in names(samples)) {
  counts <- samples[[sample_id]]
  low_quality <- sample(ncol(counts), LOW_QUALITY_PER_SAMPLE)
  target <- LOW_QUALITY_FRACTION * median(colSums(counts))
  for (j in low_quality) {
    counts[, j] <- thin(
      counts[, j, drop = FALSE],
      min(1, target / sum(counts[, j]))
    )
  }
  if (sample_id == THIN_SAMPLE) {
    counts <- thin(counts, THIN_FRACTION)
  }
  planted <- ifelse(
    seq_len(ncol(counts)) %in% low_quality,
    "low_quality",
    "cell"
  )
  # A droplet thinned to zero UMIs never reaches a filtered matrix, so it leaves the fixture.
  has_counts <- colSums(counts) > 0
  message(sprintf(
    "%s: dropped %d zero-UMI planted cells",
    sample_id,
    sum(!has_counts)
  ))
  samples[[sample_id]] <- counts[, has_counts]
  truth[[sample_id]] <- data.frame(
    sample = sample_id,
    barcode = colnames(counts)[has_counts],
    planted = planted[has_counts]
  )
}
truth <- do.call(rbind, truth)

# ---- Write inputs ----
# 10x v2 layout, uncompressed: read10xCounts() and sc.read_10x_mtx() both read it.
for (sample_id in names(samples)) {
  out_dir <- file.path(task_dir, "inputs", "data", sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  counts <- samples[[sample_id]]
  writeMM(counts, file.path(out_dir, "matrix.mtx"))
  write.table(
    data.frame(rownames(counts), rownames(counts)),
    file.path(out_dir, "genes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  writeLines(colnames(counts), file.path(out_dir, "barcodes.tsv"))
}

# ---- Write truth ----
dir.create(file.path(task_dir, "truth"), showWarnings = FALSE)
write.table(
  truth,
  file.path(task_dir, "truth", "cells.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
