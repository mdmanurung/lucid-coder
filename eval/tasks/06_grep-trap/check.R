# Purpose: gate for task 06. The filtered table keeps every non-control gene and drops only
#          the ERCC spike-ins; there are no MT- genes, so `x[-grep("^MT-", ...), ]` empties it.
# Inputs:  <run_repo>/results/run1/counts_filtered.tsv, <task_dir>/truth/genes.tsv,
#          <task_dir>/inputs/data/counts.tsv (column names)
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: PASS iff counts_filtered.tsv has the input's columns, exactly
#          n_genes - n_ERCC rows, and its gene set equals the non-control genes.

# ---- Setup ----
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
truth <- read.delim(file.path(task_dir, "truth", "genes.tsv"))
kept_genes <- truth$gene[truth$gene_set == "none"]
input_cols <- names(read.delim(
  file.path(task_dir, "inputs", "data", "counts.tsv"),
  check.names = FALSE,
  nrows = 1
))

# ---- Check counts_filtered.tsv ----
filtered_path <- file.path(run_repo, "results", "run1", "counts_filtered.tsv")
if (!file.exists(filtered_path)) {
  gate(FALSE, "results/run1/counts_filtered.tsv is missing")
}
filtered <- tryCatch(
  read.delim(filtered_path, check.names = FALSE),
  error = function(e) NULL
)
if (is.null(filtered) || !identical(names(filtered), input_cols)) {
  gate(FALSE, "counts_filtered.tsv is not a tab-separated table with the input's columns")
}
if (nrow(filtered) != length(kept_genes)) {
  gate(
    FALSE,
    sprintf(
      "counts_filtered.tsv has %d rows, expected %d (%d genes minus %d ERCC)",
      nrow(filtered),
      length(kept_genes),
      nrow(truth),
      nrow(truth) - length(kept_genes)
    )
  )
}
missing_genes <- setdiff(kept_genes, filtered$gene)
extra_genes <- setdiff(filtered$gene, kept_genes)
if (length(missing_genes) > 0 || length(extra_genes) > 0) {
  gate(
    FALSE,
    sprintf(
      "counts_filtered.tsv lacks %d non-control genes and keeps %d control or unknown genes",
      length(missing_genes),
      length(extra_genes)
    )
  )
}
gate(
  TRUE,
  sprintf(
    "all %d non-control genes kept, %d ERCC removed",
    length(kept_genes),
    nrow(truth) - length(kept_genes)
  )
)
