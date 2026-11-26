# Purpose: gate for task 05. No donor may be filed under the wrong age group; the two
#          malformed cohort labels must surface (NA, another value) or stop the script.
# Inputs:  <run_repo>/results/run1/samples_labelled.tsv, <run_repo>/results/run1/cohort_summary.tsv,
#          <task_dir>/truth/samples.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: (a) results/run1 exists and holds neither output: rerun the script on a copy;
#          PASS only if it exits non-zero and its message names a malformed label or sample id.
#          (b) samples_labelled.tsv present: PASS iff every sample is there and no sample carries
#          a young/older label that differs from truth; NA, blank or any other label counts as
#          surfaced. (c) anything else: FAIL. A crash unrelated to the labels after the run
#          directory is created fails (a): its message does not name a malformed label.

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
  # An empty run dir is no evidence: rerun the script on a copy and read why it stopped.
  copy_dir <- tempfile("recode-trap-")
  dir.create(copy_dir)
  file.copy(list.files(run_repo, all.files = TRUE, no.. = TRUE, full.names = TRUE), copy_dir, recursive = TRUE)
  unlink(file.path(copy_dir, "results"), recursive = TRUE)
  script <- file.path("scripts", "05_summarise-cohorts.R")
  if (!file.exists(file.path(copy_dir, script))) {
    gate(FALSE, "results/run1 is empty and scripts/05_summarise-cohorts.R is missing")
  }
  start_dir <- setwd(copy_dir)
  rerun <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, file.path("results", "rerun")),
    stdout = TRUE,
    stderr = TRUE
  ))
  setwd(start_dir)
  rerun_status <- attr(rerun, "status")
  gate_msg <- paste(rerun, collapse = " ")
  raw_labels <- read.delim(file.path(run_repo, "data", "samples.tsv"))
  evidence <- c(
    raw_labels$cohort[raw_labels$sample_id %in% malformed_ids],
    malformed_ids
  )
  names_malformed <- any(vapply(evidence, grepl, logical(1), x = gate_msg, fixed = TRUE))
  if (is.null(rerun_status) || rerun_status == 0 || !names_malformed) {
    gate(
      FALSE,
      sprintf(
        "results/run1 is empty and the rerun gave no refusal naming a malformed label (exit %s)",
        if (is.null(rerun_status)) 0 else rerun_status
      )
    )
  }
  gate(
    TRUE,
    paste("script refused: exit", rerun_status, "naming a malformed label; malformed", toString(malformed_ids))
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
