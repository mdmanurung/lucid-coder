# Purpose: gate for task 15. The edge signs must reach the scores: repressor T2 has a negative
#          mean score in group B, and activator T1 a positive one in group A.
# Inputs:  <run_repo>/results/run1/tf_scores.tsv, <task_dir>/inputs/data/cells.tsv,
#          <task_dir>/truth/group_means.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Truth (ulm, signed): T2 in B -3.42, T1 in A 5.17. decoupler 2.3.1 reads only a
#          `weight` column; with `mor` left as is every edge scores +1 and T2 in B reads +3.20.

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
cells <- read.delim(file.path(task_dir, "inputs", "data", "cells.tsv"))
truth <- read.delim(file.path(task_dir, "truth", "group_means.tsv"))
truth_t2_b <- truth$mean_score[truth$group == "B" & truth$source == "T2"]

# ---- Check tf_scores.tsv ----
res_path <- file.path(run_repo, "results", "run1", "tf_scores.tsv")
if (!file.exists(res_path)) {
  gate(FALSE, "results/run1/tf_scores.tsv is missing")
}
res <- tryCatch(read.delim(res_path), error = function(e) NULL)
if (is.null(res) || !all(c("cell_id", "source", "score") %in% names(res))) {
  gate(
    FALSE,
    "tf_scores.tsv is not a tab-separated table with cell_id, source, score"
  )
}
group <- cells$group[match(res$cell_id, cells$cell_id)]
mean_of <- function(tf, grp) {
  mean(suppressWarnings(as.numeric(res$score[
    res$source == tf & group %in% grp
  ])))
}
t2_b <- mean_of("T2", "B")
t1_a <- mean_of("T1", "A")
summary <- sprintf(
  "T2 mean score in B = %.2f (signed truth %.2f); T1 mean score in A = %.2f",
  t2_b,
  truth_t2_b,
  t1_a
)
problems <- c(
  if (is.na(t2_b) || t2_b >= 0) "repressor T2 not negative in B",
  if (is.na(t1_a) || t1_a <= 0) "activator T1 not positive in A"
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
