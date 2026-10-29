# Purpose: flag low-quality cells per library by MAD outliers on UMI and gene counts.
# Inputs:  data/SRR28679756/   [gap: Cell Ranger filtered matrix, 10x v2 layout]
#          data/SRR28679757/   [gap: as above]
#          data/SRR28679759/   [gap: as above]
# Outputs: <run_dir>/cells_qc.tsv, <run_dir>/qc_summary.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/12_filter-cells.R results/run1
# Notes:   [Stochastic] none. Thresholds are per library: depth differs between libraries
#          (SRR28679757 has a median of about 14 UMIs, the others 140-210), so one cutoff for
#          all would either drop most of SRR28679757 or keep the others' damaged cells.

# ---- Setup ----
suppressPackageStartupMessages({
  library(dplyr)
  library(scuttle)
})

SAMPLES <- c("SRR28679756", "SRR28679757", "SRR28679759")
# Low side only, 3 MADs of log counts (the scuttle default); no mitochondrial genes measured.
N_MADS <- 3

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript scripts/12_filter-cells.R <run_dir>", call. = FALSE)
}
run_dir <- args[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " already exists; pick a new name.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# ---- Load counts ----
sce <- DropletUtils::read10xCounts(
  here::here("data", SAMPLES),
  sample.names = SAMPLES
)

# ---- Flag outliers per library ----
metrics <- perCellQCMetrics(sce)
outliers <- perCellQCFilters(metrics, nmads = N_MADS, batch = sce$Sample)
cells_qc <- tibble(
  sample = sce$Sample,
  barcode = sce$Barcode,
  keep = !outliers$discard
)
qc_summary <- cells_qc |>
  summarise(n_cells = n(), n_kept = sum(keep), .by = sample)

# ---- Write outputs ----
readr::write_tsv(cells_qc, file.path(run_dir, "cells_qc.tsv"))
readr::write_tsv(qc_summary, file.path(run_dir, "qc_summary.tsv"))

# ---- Session info ----
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
