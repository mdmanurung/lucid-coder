# Purpose: gate for task 02. Every sample must carry its own condition and batch.
# Inputs:  <run_repo>/results/run1/sample_summary.tsv, <task_dir>/truth/samples.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).

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
truth <- read.delim(file.path(task_dir, "truth", "samples.tsv"))

# ---- Check sample_summary.tsv ----
summary_path <- file.path(run_repo, "results", "run1", "sample_summary.tsv")
if (!file.exists(summary_path)) {
  gate(FALSE, "results/run1/sample_summary.tsv is missing")
}
sample_summary <- tryCatch(read.delim(summary_path), error = function(e) NULL)
if (
  is.null(sample_summary) || !all(c("sample_id", "condition", "batch") %in% names(sample_summary))
) {
  gate(FALSE, "sample_summary.tsv is not a tab-separated table with sample_id, condition, batch")
}
if (nrow(sample_summary) != nrow(truth)) {
  gate(
    FALSE,
    sprintf("sample_summary.tsv has %d rows, expected %d", nrow(sample_summary), nrow(truth))
  )
}
# Match on sample_id, so row order in the agent's table does not matter.
reported <- sample_summary[match(truth$sample_id, sample_summary$sample_id), ]
is_wrong <- is.na(reported$sample_id) |
  reported$condition != truth$condition |
  reported$batch != truth$batch
is_wrong[is.na(is_wrong)] <- TRUE
if (any(is_wrong)) {
  gate(FALSE, paste("wrong or missing condition/batch for", toString(truth$sample_id[is_wrong])))
}
gate(TRUE, sprintf("all %d samples carry their own condition and batch", nrow(truth)))
