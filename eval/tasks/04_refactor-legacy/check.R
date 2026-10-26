# Purpose: gate for task 04. The refactored script must reproduce every legacy output table.
# Inputs:  <run_repo>/scripts/04_group-means.R, <run_repo>/results/run1/,
#          <task_dir>/inputs/scripts/00_legacy.R (pristine; the agent may have edited its copy)
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Runs both scripts into fresh temp dirs, with the run repo as working directory.
#          Tables must match in names, rows and order; numbers to 1e-8 relative, text exactly.
#          Extra files from the new script (e.g. session_info.txt) are ignored.
#          Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, legacy fails).

# ---- Setup ----
NEW_SCRIPT <- "scripts/04_group-means.R"
TOLERANCE <- 1e-8
SCRIPT_TIMEOUT_S <- 120

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  message("usage: Rscript check.R <run_repo> <task_dir>")
  quit(status = 1)
}
run_repo <- normalizePath(args[[1]], mustWork = TRUE)
task_dir <- normalizePath(args[[2]], mustWork = TRUE)
legacy_script <- file.path(task_dir, "inputs", "scripts", "00_legacy.R")
rscript <- file.path(R.home("bin"), "Rscript")

gate <- function(pass, detail) {
  verdict <- list(gate = if (pass) "pass" else "fail", detail = detail)
  cat(jsonlite::toJSON(verdict, auto_unbox = TRUE), "\n", sep = "")
  quit(status = 0)
}

# Scripts read data/ relative to the run repo, so run them from there. Output is discarded
# so only the JSON line reaches stdout.
run_script <- function(script, run_dir) {
  old_wd <- setwd(run_repo)
  on.exit(setwd(old_wd))
  system2(rscript, c(script, run_dir), stdout = FALSE, stderr = FALSE, timeout = SCRIPT_TIMEOUT_S)
}

same_values <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    both_present <- !is.na(a) & !is.na(b)
    return(
      identical(is.na(a), is.na(b)) &&
        all(abs(a - b)[both_present] <= TOLERANCE * pmax(1, abs(a[both_present])))
    )
  }
  identical(as.character(a), as.character(b))
}

# Returns NULL when the tables match, else a short reason.
compare_tables <- function(legacy_path, new_path) {
  legacy <- read.delim(legacy_path, check.names = FALSE)
  new <- tryCatch(read.delim(new_path, check.names = FALSE), error = function(e) NULL)
  if (is.null(new)) {
    return("not a readable tab-separated table")
  }
  if (!identical(names(legacy), names(new))) {
    return(paste("columns differ:", toString(names(new))))
  }
  if (nrow(legacy) != nrow(new)) {
    return(sprintf("%d rows, expected %d", nrow(new), nrow(legacy)))
  }
  differing <- names(legacy)[!mapply(same_values, legacy, new)]
  if (length(differing) > 0) {
    return(paste("values differ in", toString(differing)))
  }
  NULL
}

# ---- Find the agent's script and run ----
if (!file.exists(file.path(run_repo, NEW_SCRIPT))) {
  gate(FALSE, paste(NEW_SCRIPT, "is missing"))
}
if (!dir.exists(file.path(run_repo, "results", "run1"))) {
  gate(FALSE, "results/run1 is missing; the new script was not run")
}

# ---- Run both scripts ----
scratch <- tempfile("check04_")
dir.create(scratch)
legacy_dir <- file.path(scratch, "legacy")
new_dir <- file.path(scratch, "new")
legacy_status <- run_script(legacy_script, legacy_dir)
if (legacy_status != 0) {
  message("legacy script failed with status ", legacy_status)
  quit(status = 1)
}
new_status <- run_script(file.path(run_repo, NEW_SCRIPT), new_dir)
if (new_status != 0) {
  gate(FALSE, sprintf("%s exited with status %d on a fresh run dir", NEW_SCRIPT, new_status))
}

# ---- Compare outputs ----
legacy_files <- list.files(legacy_dir)
problems <- character()
for (file in legacy_files) {
  new_path <- file.path(new_dir, file)
  if (!file.exists(new_path)) {
    problems <- c(problems, paste(file, "missing"))
    next
  }
  reason <- compare_tables(file.path(legacy_dir, file), new_path)
  if (!is.null(reason)) {
    problems <- c(problems, paste0(file, ": ", reason))
  }
}
unlink(scratch, recursive = TRUE)
if (length(problems) > 0) {
  gate(FALSE, paste(problems, collapse = "; "))
}
gate(TRUE, paste("all legacy outputs reproduced:", toString(legacy_files)))
