# Purpose: gate for task 12. In every sample the QC must drop the planted low-quality cells
#          and keep the rest; one global cutoff cannot do both once one sample is shallow.
# Inputs:  <run_repo>/results/run1/cells_qc.tsv, <task_dir>/truth/cells.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule, per sample: keep >= MIN_KEPT of the real cells and drop >= MIN_DROPPED of the
#          planted low-quality cells. Inputs are built by make-inputs.R: SRR28679757 is thinned
#          to 10% depth (median about 14 UMIs); 24 cells per sample are thinned to about 15% of
#          their sample's median (23 left in SRR28679757), so the deep samples' planted cells
#          sit above most real cells of SRR28679757.
#          Checked on the frozen inputs: global MAD (3 or 5, scaled or not) drops 0-17 of 24
#          planted cells; fixed UMI floors 10-20 drop 0-11 of 24 or keep under half of
#          SRR28679757. Per-sample MAD (3 scaled lower, 5 scaled lower, 3 or 5 unscaled
#          two-sided) passes.
#          qc_summary.tsv (sample, n_cells, n_kept) must list every sample once, with n_cells
#          equal to the cells of that sample in truth and n_kept equal to the TRUE rows of the
#          script's own cells_qc.tsv.

# ---- Setup ----
# Per-sample 3-MAD two-sided keeps 82% of real cells; the floor leaves room for that.
MIN_KEPT <- 0.75
# At ~14 UMIs the shallow sample's planted cells overlap its real ones: per-sample 5 scaled
# MADs drops only 15 of its 23 (65%). Global rules drop at most 71% in a deep sample and
# fail MIN_KEPT there or in SRR28679757.
MIN_DROPPED <- 0.6

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  message("usage: Rscript check.R <run_repo> <task_dir>")
  quit(status = 1)
}
run_repo <- args[[1]]
task_dir <- args[[2]]

gate <- function(pass, detail) {
  verdict <- list(gate = if (pass) "pass" else "fail", detail = detail)
  cat(jsonlite::toJSON(verdict, auto_unbox = TRUE), "\n", sep = "")
  quit(status = 0)
}

# ---- Load truth ----
truth <- read.delim(file.path(task_dir, "truth", "cells.tsv"))

# ---- Check cells_qc.tsv ----
qc_path <- file.path(run_repo, "results", "run1", "cells_qc.tsv")
if (!file.exists(qc_path)) {
  gate(FALSE, "results/run1/cells_qc.tsv is missing")
}
qc <- tryCatch(
  read.delim(qc_path, colClasses = "character"),
  error = function(e) NULL
)
if (is.null(qc) || !all(c("sample", "barcode", "keep") %in% names(qc))) {
  gate(
    FALSE,
    "cells_qc.tsv is not a tab-separated table with sample, barcode, keep"
  )
}
keep_text <- toupper(trimws(qc$keep))
if (!all(keep_text %in% c("TRUE", "FALSE", "T", "F", "1", "0"))) {
  gate(FALSE, "cells_qc.tsv: keep must be TRUE/FALSE")
}
qc_keep <- setNames(
  keep_text %in% c("TRUE", "T", "1"),
  paste(qc$sample, qc$barcode)
)
keep <- qc_keep[paste(truth$sample, truth$barcode)]
if (anyNA(keep)) {
  gate(
    FALSE,
    sprintf("cells_qc.tsv lacks %d of %d cells", sum(is.na(keep)), nrow(truth))
  )
}

# ---- Score each sample ----
per_sample <- lapply(split(data.frame(truth, keep), truth$sample), function(x) {
  real <- x$planted == "cell"
  data.frame(
    sample = x$sample[[1]],
    kept = sum(x$keep[real]),
    n_real = sum(real),
    dropped = sum(!x$keep[!real]),
    n_planted = sum(!real)
  )
})
per_sample <- do.call(rbind, per_sample)
ok <- per_sample$kept >= MIN_KEPT * per_sample$n_real &
  per_sample$dropped >= MIN_DROPPED * per_sample$n_planted
summary <- sprintf(
  "%s kept %d/%d cells, dropped %d/%d planted",
  per_sample$sample,
  per_sample$kept,
  per_sample$n_real,
  per_sample$dropped,
  per_sample$n_planted
)

# ---- Check qc_summary.tsv ----
# The summary must agree with truth on cells per sample and with the script's own keep column.
summary_path <- file.path(run_repo, "results", "run1", "qc_summary.tsv")
qc_summary <- if (file.exists(summary_path)) {
  tryCatch(read.delim(summary_path), error = function(e) NULL)
}
summary_problem <- if (is.null(qc_summary)) {
  "qc_summary.tsv is missing or unreadable"
} else if (!all(c("sample", "n_cells", "n_kept") %in% names(qc_summary))) {
  "qc_summary.tsv lacks sample, n_cells or n_kept"
} else {
  expected <- data.frame(
    sample = sort(unique(truth$sample)),
    n_cells = as.vector(table(truth$sample)[sort(unique(truth$sample))]),
    n_kept = as.vector(tapply(qc_keep, qc$sample, sum)[sort(unique(truth$sample))])
  )
  got <- qc_summary[match(expected$sample, qc_summary$sample), ]
  if (
    nrow(qc_summary) != nrow(expected) ||
      anyNA(got$sample) ||
      any(got$n_cells != expected$n_cells) ||
      any(got$n_kept != expected$n_kept)
  ) {
    "qc_summary.tsv does not match the cells per sample and cells_qc.tsv keep counts"
  }
}
gate(
  all(ok) && is.null(summary_problem),
  paste(
    c(
      if (!all(ok)) paste("failed:", toString(per_sample$sample[!ok])),
      summary_problem,
      summary
    ),
    collapse = "; "
  )
)
