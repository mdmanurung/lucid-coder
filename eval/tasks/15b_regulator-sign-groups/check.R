# Purpose: gate for task 15b. The edge signs must reach the scores: repressor T2 has a negative
#          mean score in group B, and activator T1 a positive one in group A. The A vs B test
#          must use the sample as replicate: each gated TF's p-value within a tolerance of a
#          paired t-test on the 3 per-sample group means, and its difference of the same sign.
# Inputs:  <run_repo>/results/run1/tf_scores.tsv, <run_repo>/results/run1/tf_group_test.tsv,
#          <task_dir>/inputs/data/cells.tsv, <task_dir>/truth/group_means.tsv,
#          <task_dir>/truth/tf_group_test.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Truth (ulm, signed): T2 in B -3.42, T1 in A 5.17. decoupler 2.3.1 reads only a
#          `weight` column; with `mor` left as is every edge scores +1 and T2 in B reads +3.20.
#          Group test truth (paired t, n = 3 samples): T1 9.75e-05, T2 2.45e-05, T3 1.76e-04,
#          T4 3.07e-03, T5 0.189. Tolerance, measured on these TFs as |log10 p - log10 truth|:
#          sample-level tests stay within 4.01 on T1-T4 (unpaired t on the 6 means 3.65,
#          Welch on the 6 means 2.81, Wilcoxon signed-rank 4.01, floored at p = 0.25 with 3
#          pairs); cell-level tests on the 30 cells move T1-T4 by 10.90 at least (Welch 10.90
#          to 27.06, Student 11.08 to 27.16, OLS with sample as covariate 10.72 to 24.73).
#          Gated on T1-T4; T5 is null and moves under 1.07 whatever the test. Known soft
#          spot: cell-level Mann-Whitney passes (gaps 0.86 to 2.96), because scipy's
#          asymptotic p at full separation, 3.39e-06, lands near the paired truth; these data
#          cannot gate it. A limma-style test was skipped: 6 means, nothing to moderate.

# ---- Setup ----
LOG10_P_TOL <- 6
GATED_TFS <- c("T1", "T2", "T3", "T4") # T5 is null; every test agrees on it

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
truth_test <- read.delim(file.path(task_dir, "truth", "tf_group_test.tsv"))
truth_test <- truth_test[match(GATED_TFS, truth_test$source), ]

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

# ---- Check tf_group_test.tsv ----
test_path <- file.path(run_repo, "results", "run1", "tf_group_test.tsv")
if (!file.exists(test_path)) {
  gate(
    FALSE,
    paste(
      c(problems, "results/run1/tf_group_test.tsv is missing"),
      collapse = "; "
    )
  )
}
test <- tryCatch(read.delim(test_path), error = function(e) NULL)
if (
  is.null(test) || !all(c("source", "difference", "p_value") %in% names(test))
) {
  gate(
    FALSE,
    paste(
      c(
        problems,
        "tf_group_test.tsv is not a tab-separated table with source, difference, p_value"
      ),
      collapse = "; "
    )
  )
}
reported <- test[match(GATED_TFS, test$source), ]
p_value <- suppressWarnings(as.numeric(reported$p_value))
difference <- suppressWarnings(as.numeric(reported$difference))
shift <- abs(log10(p_value) - log10(truth_test$p_value))
is_off <- is.na(shift) | shift > LOG10_P_TOL
is_flipped <- is.na(difference) |
  sign(difference) != sign(truth_test$difference)
test_summary <- sprintf(
  "group test |log10 p shift| on %s max %.2f (tolerance %.2f)",
  toString(GATED_TFS),
  max(shift, na.rm = TRUE),
  LOG10_P_TOL
)
problems <- c(
  problems,
  if (any(is_off)) {
    paste(
      "group p-value off the paired sample-level test:",
      toString(sprintf(
        "%s (%.2g vs %.2g)",
        GATED_TFS,
        p_value,
        truth_test$p_value
      )[is_off])
    )
  },
  if (any(is_flipped)) {
    paste("difference sign wrong or missing:", toString(GATED_TFS[is_flipped]))
  }
)
gate(
  length(problems) == 0,
  paste(c(problems, summary, test_summary), collapse = "; ")
)
