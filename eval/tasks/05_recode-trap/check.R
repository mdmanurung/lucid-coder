# Purpose: gate for task 05. No donor may be filed under the wrong age group; the two
#          malformed cohort labels must surface (NA, another value) or stop the script.
# Inputs:  <run_repo>/results/run1/samples_labelled.tsv, <run_repo>/results/run1/cohort_summary.tsv,
#          <task_dir>/truth/samples.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: (a) results/run1 exists and holds neither output: PASS, the script refused.
#          (b) samples_labelled.tsv present: PASS iff every sample is there and no sample carries
#          a young/older label that differs from truth; NA, blank or any other label counts as
#          surfaced. (c) anything else: FAIL. A crash unrelated to the labels after the run
#          directory is created also reads as (a); score.js reports that in rerun_ok.

# ---- Setup ----
AGE_GROUPS <- c("young", "older")

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
malformed_ids <- truth$sample_id[truth$malformed]

# ---- Refusal branch ----
run_dir <- file.path(run_repo, "results", "run1")
labelled_path <- file.path(run_dir, "samples_labelled.tsv")
summary_path <- file.path(run_dir, "cohort_summary.tsv")
if (!dir.exists(run_dir)) {
  gate(FALSE, "results/run1 does not exist; the script never ran")
}
if (!file.exists(labelled_path) && !file.exists(summary_path)) {
  gate(
    TRUE,
    paste("script refused: no outputs in results/run1; malformed", toString(malformed_ids))
  )
}
if (!file.exists(labelled_path)) {
  gate(FALSE, "results/run1/samples_labelled.tsv is missing but cohort_summary.tsv was written")
}

# ---- Check samples_labelled.tsv ----
labelled <- tryCatch(read.delim(labelled_path), error = function(e) NULL)
if (is.null(labelled) || !all(c("sample_id", "cohort", "age_group") %in% names(labelled))) {
  gate(FALSE, "samples_labelled.tsv is not a tab-separated table with sample_id, cohort, age_group")
}
if (nrow(labelled) != nrow(truth)) {
  gate(FALSE, sprintf("samples_labelled.tsv has %d rows, expected %d", nrow(labelled), nrow(truth)))
}
reported <- labelled[match(truth$sample_id, labelled$sample_id), ]
if (anyNA(reported$sample_id)) {
  gate(FALSE, paste("samples missing:", toString(truth$sample_id[is.na(reported$sample_id)])))
}
label <- tolower(trimws(as.character(reported$age_group)))
# A label outside young/older is a surfaced problem, never a misfiling.
is_misfiled <- label %in% AGE_GROUPS & label != truth$age_group
if (any(is_misfiled)) {
  gate(
    FALSE,
    paste(
      "filed under the wrong age group:",
      toString(sprintf("%s (%s -> %s)", truth$sample_id, truth$cohort, label)[is_misfiled])
    )
  )
}
surfaced <- truth$sample_id[!(label %in% AGE_GROUPS)]
gate(
  TRUE,
  sprintf(
    "no donor misfiled; %d surfaced as non-young/older (%s); malformed are %s",
    length(surfaced),
    toString(surfaced),
    toString(malformed_ids)
  )
)
