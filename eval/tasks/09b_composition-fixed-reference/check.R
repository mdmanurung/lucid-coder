# Purpose: gate for task 09b. Only Enterocyte rises (10-fold) under Salmonella and Goblet is a
#          known-stable reference; the other four types do not change, so a test that ignores
#          the reference reads the Enterocyte rise as their depletion (closure).
# Inputs:  <run_repo>/results/run1/composition.tsv, <task_dir>/truth/cell_types.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: the Goblet row is ignored; (a) no fixed type (Stem, TA, Endocrine, Tuft) with
#          q < 0.05 and log_fc < 0; (b) Enterocyte log_fc > 0 and q < 0.05; (c) all five
#          non-reference types are reported. The fixture is one simulated draw (make-inputs.R):
#          6 types, 4 vs 4 mice, negative binomial phi 0.1, Enterocyte x10. Over 200 simulated
#          draws of this design the share of draws with a fixed type called down at q < 0.05
#          was 0.76 for the per-mouse proportion t-test, 0.92 for the quasi-binomial GLM and
#          0.73 for crumblr + limma, but 0.08 for edgeR QL with offset log(Goblet count) and
#          0.03 for the t-test on log(count / Goblet count). So q < 0.05 separates methods that
#          use the reference from those that ignore it; the seed was chosen so this draw shows
#          that split for all five methods. Only depletion of a fixed type is gated: closure
#          pushes them down, never up. log_fc is gated for sign only, not scale, so any log base
#          and a CLR or reference-ratio scale pass.
#          Known soft spots: a closure-unaware method can pass on a lucky draw, and a rare
#          reference-aware one can fail; this draw was picked so neither happens for the five
#          methods above. edgeR QL with default TMM factors and no reference passes this draw
#          (smallest fixed q 0.68); TMM absorbs much of the closure (0.17 of simulated draws
#          call a fixed type down), so the gate does not catch it.

# ---- Setup ----
QVALUE_CUTOFF <- 0.05

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
truth <- truth[truth$role != "reference", ]
up_type <- truth$cell_type[truth$role == "up"]

# ---- Check composition.tsv ----
out_path <- file.path(run_repo, "results", "run1", "composition.tsv")
if (!file.exists(out_path)) {
  gate(FALSE, "results/run1/composition.tsv is missing")
}
out <- tryCatch(read.delim(out_path), error = function(e) NULL)
if (
  is.null(out) ||
    !all(c("cell_type", "log_fc", "p_value", "q_value") %in% names(out))
) {
  gate(
    FALSE,
    "composition.tsv is not a tab-separated table with cell_type, log_fc, p_value, q_value"
  )
}
reported <- out[match(truth$cell_type, out$cell_type), ]
if (anyNA(reported$cell_type)) {
  gate(
    FALSE,
    paste(
      "cell types missing:",
      toString(truth$cell_type[is.na(reported$cell_type)])
    )
  )
}
q_value <- setNames(
  suppressWarnings(as.numeric(reported$q_value)),
  truth$cell_type
)
log_fc <- setNames(
  suppressWarnings(as.numeric(reported$log_fc)),
  truth$cell_type
)
fixed <- truth$cell_type[truth$role == "unchanged"]
false_calls <- fixed[
  !is.na(q_value[fixed]) &
    q_value[fixed] < QVALUE_CUTOFF &
    !is.na(log_fc[fixed]) &
    log_fc[fixed] < 0
]
summary <- sprintf(
  "%s log_fc = %.2f, q = %.3g; smallest fixed-type q = %.3g",
  up_type,
  log_fc[[up_type]],
  q_value[[up_type]],
  min(q_value[fixed], na.rm = TRUE)
)
problems <- c(
  if (length(false_calls) > 0) {
    paste("fixed types called down:", toString(false_calls))
  },
  if (
    is.na(q_value[[up_type]]) ||
      q_value[[up_type]] >= QVALUE_CUTOFF ||
      !isTRUE(log_fc[[up_type]] > 0)
  ) {
    paste(up_type, "not increased with q <", QVALUE_CUTOFF)
  }
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
