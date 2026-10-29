# Purpose: gate for task 09. Only Enterocyte rises under Salmonella (scCODA, Haber 2017); the
#          mouse, not the cell, is the unit of replication.
# Inputs:  <run_repo>/results/run1/composition.tsv, <task_dir>/truth/cell_types.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: (a) no unchanged type with q < 0.01 (closure: other types only look depleted);
#          (b) no q below 1e-4, out of reach for 4 vs 2 mice, so cells were the unit;
#          (c) the up type has log_fc > 0 and q < 0.25. No smallest-p rule: a per-mouse
#          quasi-binomial GLM ranks TA above Enterocyte. Stem/TA direction is not gated (no
#          ground truth). Rule (a) sits at 0.01 and (c) at 0.25 so every replicate-aware
#          method passes: crumblr (with and without eBayes(robust = TRUE)), edgeR QL on
#          per-mouse counts, per-mouse t-tests on proportions, per-mouse quasi-binomial GLM
#          (Enterocyte q 0.127), limma on per-mouse logit proportions with robust eBayes.
#          Cells as the unit (pooled prop.test) fails rule (a) and (b).
#          Known soft spots: per-mouse tests that ignore closure (t-tests, quasi-binomial) pass,
#          since 2 infected mice give too little power to call the depleted types. log_fc is
#          gated for sign only, not scale. The q < 0.25 limit also lets a weak rise pass.

# ---- Setup ----
Q_CUTOFF <- 0.01
Q_UP_CUTOFF <- 0.25
Q_FLOOR <- 1e-4

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
truth <- read.delim(file.path(task_dir, "truth", "cell_types.tsv"))
up_type <- truth$cell_type[truth$role == "up"]

# ---- Check composition.tsv ----
out_path <- file.path(run_repo, "results", "run1", "composition.tsv")
if (!file.exists(out_path)) {
  gate(FALSE, "results/run1/composition.tsv is missing")
}
out <- tryCatch(read.delim(out_path), error = function(e) NULL)
if (is.null(out) || !all(c("cell_type", "log_fc", "p_value", "q_value") %in% names(out))) {
  gate(
    FALSE,
    "composition.tsv is not a tab-separated table with cell_type, log_fc, p_value, q_value"
  )
}
reported <- out[match(truth$cell_type, out$cell_type), ]
if (anyNA(reported$cell_type)) {
  gate(FALSE, paste("cell types missing:", toString(truth$cell_type[is.na(reported$cell_type)])))
}
q_value <- setNames(suppressWarnings(as.numeric(reported$q_value)), truth$cell_type)
log_fc <- setNames(suppressWarnings(as.numeric(reported$log_fc)), truth$cell_type)
false_calls <- names(q_value)[truth$role != "up" & !is.na(q_value) & q_value < Q_CUTOFF]
summary <- sprintf(
  "%s log_fc = %.2f, q = %.3g; smallest other q = %.3g",
  up_type,
  log_fc[[up_type]],
  q_value[[up_type]],
  min(q_value[names(q_value) != up_type], na.rm = TRUE)
)
problems <- c(
  if (length(false_calls) > 0) paste("unchanged types called:", toString(false_calls)),
  if (any(q_value < Q_FLOOR, na.rm = TRUE)) {
    sprintf("q below %g with 6 mice: cells were the unit", Q_FLOOR)
  },
  if (
    is.na(q_value[[up_type]]) ||
      q_value[[up_type]] >= Q_UP_CUTOFF ||
      !isTRUE(log_fc[[up_type]] > 0)
  ) {
    paste(up_type, "not increased with q <", Q_UP_CUTOFF)
  }
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
