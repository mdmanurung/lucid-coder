# Purpose: gate for task 03. The figure and its source data exist; the data has one row per
#          donor x cell type.
# Inputs:  <run_repo>/results/run1/proportions_by-condition.{pdf,png,csv},
#          <task_dir>/truth/donor_cell_types.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).

# ---- Setup ----
PANEL <- "proportions_by-condition"

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
truth <- read.delim(file.path(task_dir, "truth", "donor_cell_types.tsv"))

# ---- Check the panel files ----
panel_files <- paste0(PANEL, c(".pdf", ".png", ".csv"))
panel_paths <- file.path(run_repo, "results", "run1", panel_files)
# file.size is NA for a missing file.
is_absent <- is.na(file.size(panel_paths)) | file.size(panel_paths) == 0
if (any(is_absent)) {
  gate(FALSE, paste("missing or empty in results/run1:", toString(panel_files[is_absent])))
}
source_data <- tryCatch(read.csv(panel_paths[[3]]), error = function(e) NULL)
if (is.null(source_data)) {
  gate(FALSE, paste(panel_files[[3]], "is not a readable CSV"))
}
if (nrow(source_data) != nrow(truth)) {
  gate(
    FALSE,
    sprintf(
      "%s has %d rows, expected %d (one per donor x cell type)",
      panel_files[[3]],
      nrow(source_data),
      nrow(truth)
    )
  )
}
gate(TRUE, sprintf("PDF, PNG and CSV present; CSV has %d rows", nrow(truth)))
