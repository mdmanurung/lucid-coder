# Purpose: gate for task 07b. Task 07 with planted missing TF scores: each TF's time p-value
#          must match the donor mixed model fitted on that TF's own scored samples, and n_samples
#          must say how many samples that was; dropping every incomplete sample misses both.
# Inputs:  <run_repo>/results/run1/tf_time.tsv, <task_dir>/truth/tf_time.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Missing: HOXB2 in 17 samples, RUNX1 in 2, TET1 in 1, so 20 of 69 samples are
#          incomplete (make-inputs.R). Truth is task 07's ML lmer on each TF's available rows.
#          Tolerance as in task 07; on these inputs REML stays within 0.26 log10 units and lm
#          without the donor moves every trap TF by 1.19 to 1.99. A wide drop of incomplete
#          samples (drop_na, na.omit or complete.cases over all TFs, 49 samples left) moves 10
#          TFs past the tolerance, PLAGL1 by 2.55, ELF4 by 1.72 and TGIF2 by 1.60, so it fails
#          on p-values as well as on n_samples. Code that pivots long and lets lmer drop each
#          TF's own NA rows is already right and passes; that is a fair outcome. Such code still
#          fails if it counts n_samples before lmer drops the NAs (69 for HOXB2, not 52). q_value
#          must equal BH over the script's own p_value (1e-6 relative), so a missing, unadjusted,
#          Bonferroni or BY q fails.

# ---- Setup ----
LOG10_P_TOL <- 0.75
Q_REL_TOL <- 1e-6 # q is arithmetic on the script's p, not a model choice

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
truth <- read.delim(file.path(task_dir, "truth", "tf_time.tsv"))

# ---- Check tf_time.tsv ----
out_path <- file.path(run_repo, "results", "run1", "tf_time.tsv")
if (!file.exists(out_path)) {
  gate(FALSE, "results/run1/tf_time.tsv is missing")
}
out <- tryCatch(read.delim(out_path), error = function(e) NULL)
if (is.null(out) || !all(c("tf", "estimate", "p_value", "q_value", "n_samples") %in% names(out))) {
  gate(
    FALSE,
    paste(
      "tf_time.tsv is not a tab-separated table with columns",
      "tf, estimate, p_value, q_value, n_samples"
    )
  )
}
reported <- out[match(truth$tf, out$tf), ]
if (anyNA(reported$tf)) {
  gate(FALSE, paste("TFs missing:", toString(truth$tf[is.na(reported$tf)])))
}
p_value <- suppressWarnings(as.numeric(reported$p_value))
estimate <- suppressWarnings(as.numeric(reported$estimate))
shift <- abs(log10(p_value) - log10(truth$p_value))
is_off <- is.na(shift) | shift > LOG10_P_TOL
# q: BH over the script's own 30 p-values, so its choice of model does not move the check.
q_value <- suppressWarnings(as.numeric(reported$q_value))
q_expected <- p.adjust(p_value, method = "BH")
is_q_ok <- abs(q_value - q_expected) <= Q_REL_TOL * q_expected
is_q_off <- !(is_q_ok %in% TRUE)
# n_samples: rows used per TF; a wide drop of incomplete samples gives 49 for every TF.
n_samples <- suppressWarnings(as.numeric(reported$n_samples))
is_n_off <- !((n_samples == truth$n_samples) %in% TRUE)
# A slope with the wrong sign is wrong whatever its p-value.
is_flipped <- !is.na(estimate) & sign(estimate) != sign(truth$estimate) & truth$p_value < 0.05
trap <- truth$role == "trap"
summary <- sprintf(
  "trap TFs |log10 p shift| median %.2f, max %.2f (tolerance %.2f)",
  median(shift[trap], na.rm = TRUE),
  max(shift[trap], na.rm = TRUE),
  LOG10_P_TOL
)
problems <- c(
  if (any(is_off)) {
    paste(
      "time p-value off the per-TF mixed model:",
      toString(sprintf("%s (%.2g vs %.2g)", truth$tf, p_value, truth$p_value)[is_off])
    )
  },
  if (any(is_n_off)) {
    paste(
      "n_samples off:",
      toString(sprintf("%s (%s vs %d)", truth$tf, n_samples, truth$n_samples)[is_n_off])
    )
  },
  if (any(is_flipped)) paste("time slope sign flipped:", toString(truth$tf[is_flipped])),
  if (any(is_q_off)) {
    paste("q_value not BH of its own p_value:", toString(truth$tf[is_q_off]))
  }
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
