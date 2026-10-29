# Purpose: gate for task 16. Each donor with a D0 and a D7 sample appears once, paired from the
#          latest run of each visit; the donor without a D0 run is not paired.
# Inputs:  <run_repo>/results/run1/pairs.tsv, <task_dir>/truth/pairs.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: PASS iff pairs.tsv has one row per paired donor and no other donor, il6_d0 and
#          il6_d7 equal the latest runs' values and change equals il6_d7 - il6_d0 (within 0.005),
#          and paired_test.tsv has n_donors, mean_change and a p_value within 0.1 log10 units of
#          t.test(truth change), so a flipped sign or an unpaired test fails. Every problem is
#          listed. Data: make-data.R.

# ---- Setup ----
TOLERANCE <- 0.005
LOG10_P_TOL <- 0.1 # a Welch test on the same donors moves p by more than this

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
truth <- read.delim(file.path(task_dir, "truth", "pairs.tsv"))
paired <- truth[truth$paired, ]

# ---- Check pairs.tsv ----
pairs_path <- file.path(run_repo, "results", "run1", "pairs.tsv")
if (!file.exists(pairs_path)) {
  gate(FALSE, "results/run1/pairs.tsv is missing")
}
pairs <- tryCatch(read.delim(pairs_path), error = function(e) NULL)
if (is.null(pairs) || !all(c("donor_id", "il6_d0", "il6_d7", "change") %in% names(pairs))) {
  gate(FALSE, "pairs.tsv is not a tab-separated table with donor_id, il6_d0, il6_d7, change")
}

problems <- character()
row_counts <- table(pairs$donor_id)
repeated <- row_counts[row_counts > 1]
if (length(repeated) > 0) {
  problems <- c(
    problems,
    paste("donors with more than one row:", toString(sprintf("%s (%d)", names(repeated), repeated)))
  )
}
extra <- setdiff(pairs$donor_id, paired$donor_id)
if (length(extra) > 0) {
  problems <- c(problems, paste("extra donors (no D0 or D7 run):", toString(extra)))
}
missing <- setdiff(paired$donor_id, pairs$donor_id)
if (length(missing) > 0) {
  problems <- c(problems, paste("paired donors missing:", toString(missing)))
}
single <- pairs[pairs$donor_id %in% names(row_counts)[row_counts == 1], ]
expected <- paired[match(single$donor_id, paired$donor_id), ]
# NA (a blank value) counts as wrong; extra donors are already reported above.
is_match <- abs(single$il6_d0 - expected$il6_d0) < TOLERANCE &
  abs(single$il6_d7 - expected$il6_d7) < TOLERANCE
is_wrong <- !is.na(expected$donor_id) & !(is_match %in% TRUE)
if (any(is_wrong)) {
  problems <- c(
    problems,
    paste(
      "values not from the latest D0 and D7 runs:",
      toString(sprintf(
        "%s (%.2f, %.2f; expected %.2f, %.2f)",
        single$donor_id,
        single$il6_d0,
        single$il6_d7,
        expected$il6_d0,
        expected$il6_d7
      )[is_wrong])
    )
  )
}
is_change <- abs(pairs$change - (pairs$il6_d7 - pairs$il6_d0)) < TOLERANCE
wrong_change <- !(is_change %in% TRUE)
if (any(wrong_change)) {
  problems <- c(
    problems,
    paste("change is not il6_d7 - il6_d0:", toString(pairs$donor_id[wrong_change]))
  )
}

# ---- Check paired_test.tsv ----
test_path <- file.path(run_repo, "results", "run1", "paired_test.tsv")
test_out <- if (file.exists(test_path)) tryCatch(read.delim(test_path), error = function(e) NULL)
if (is.null(test_out) || !all(c("n_donors", "mean_change", "p_value") %in% names(test_out))) {
  problems <- c(problems, "paired_test.tsv is missing or lacks n_donors, mean_change, p_value")
} else {
  true_change <- paired$il6_d7 - paired$il6_d0
  n_donors <- suppressWarnings(as.numeric(test_out$n_donors[[1]]))
  mean_change <- suppressWarnings(as.numeric(test_out$mean_change[[1]]))
  p_value <- suppressWarnings(as.numeric(test_out$p_value[[1]]))
  p_shift <- abs(log10(p_value) - log10(t.test(true_change)$p.value))
  if (!isTRUE(n_donors == nrow(paired))) {
    problems <- c(problems, sprintf("paired_test n_donors %s, expected %d", n_donors, nrow(paired)))
  }
  if (!isTRUE(abs(mean_change - mean(true_change)) < TOLERANCE)) {
    problems <- c(
      problems,
      sprintf("paired_test mean_change %.3f, expected %.3f", mean_change, mean(true_change))
    )
  }
  if (!isTRUE(p_shift < LOG10_P_TOL)) {
    problems <- c(
      problems,
      sprintf("paired_test p_value %.3g, expected the paired t-test %.3g", p_value, t.test(true_change)$p.value)
    )
  }
}
if (length(problems) > 0) {
  gate(FALSE, paste(problems, collapse = "; "))
}
gate(
  TRUE,
  sprintf(
    "%d donors paired from the latest runs; %s left out (no D0)",
    nrow(paired),
    toString(truth$donor_id[!truth$paired])
  )
)
