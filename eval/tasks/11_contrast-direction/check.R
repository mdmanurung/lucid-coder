# Purpose: gate for task 11. Every planted protein must carry its planted sign of disease minus
#          healthy; "disease" sorts first, so a default factor makes it the reference and flips
#          every sign.
# Inputs:  <run_repo>/results/run1/de_results.tsv, <task_dir>/truth/proteins.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: PASS iff de_results.tsv has protein, log2fc, p_value and fdr, one row per
#          protein, sign(log2fc) equals the planted sign for all 20 planted proteins, p_value
#          and fdr lie in [0, 1], and fdr equals BH of p_value (1e-6 relative). Data: make-data.R.

# ---- Setup ----
FDR_REL_TOL <- 1e-6

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
truth <- read.delim(file.path(task_dir, "truth", "proteins.tsv"))
planted <- truth[truth$direction != "none", ]

# ---- Check de_results.tsv ----
results_path <- file.path(run_repo, "results", "run1", "de_results.tsv")
if (!file.exists(results_path)) {
  gate(FALSE, "results/run1/de_results.tsv is missing")
}
results <- tryCatch(read.delim(results_path), error = function(e) NULL)
if (is.null(results) || !all(c("protein", "log2fc", "p_value", "fdr") %in% names(results))) {
  gate(FALSE, "de_results.tsv is not a tab-separated table with protein, log2fc, p_value, fdr")
}
if (anyDuplicated(results$protein) || !setequal(results$protein, truth$protein)) {
  gate(
    FALSE,
    sprintf("de_results.tsv has %d rows, expected one per protein (%d)", nrow(results), nrow(truth))
  )
}
reported <- results$log2fc[match(planted$protein, results$protein)]
is_flipped <- is.na(reported) | sign(reported) != sign(planted$planted_log2fc)
if (any(is_flipped)) {
  gate(
    FALSE,
    sprintf(
      "sign flipped or missing for %d of %d planted proteins: %s",
      sum(is_flipped),
      nrow(planted),
      toString(sprintf(
        "%s (%s, log2fc %.2f)",
        planted$protein,
        planted$direction,
        reported
      )[is_flipped])
    )
  )
}
p_value <- suppressWarnings(as.numeric(results$p_value))
fdr <- suppressWarnings(as.numeric(results$fdr))
is_prob <- function(x) !anyNA(x) && all(x >= 0 & x <= 1)
if (!is_prob(p_value) || !is_prob(fdr)) {
  gate(FALSE, "p_value and fdr must be numbers in [0, 1]")
}
if (!isTRUE(all(abs(fdr - p.adjust(p_value, "BH")) <= FDR_REL_TOL * p.adjust(p_value, "BH")))) {
  gate(FALSE, "fdr is not BH of p_value over all proteins")
}
gate(TRUE, sprintf("all %d planted proteins carry the planted sign", nrow(planted)))
