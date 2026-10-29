# Purpose: gate for task 10. Gene IDs must be mapped to symbols before scoring: most CollecTRI
#          TFs are scored and SRF ranks in the top 5.
# Inputs:  <run_repo>/results/run1/tf_activity.tsv, <task_dir>/truth/tfs.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Truth (decoupler 2.3.1, ulm, tmin 5, IDs mapped): 655 TFs, top SRF, MYOCD, SMAD3.
#          Unmapped, 15% of rows carry symbols, so ulm scores 229 TFs with no error and SRF
#          drops out of the top. MIN_TFS sits between the two and leaves room for a stricter tmin.

# ---- Setup ----
MIN_TFS <- 400
TOP_N <- 5
MARKER_TF <- "SRF"

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
truth <- read.delim(file.path(task_dir, "truth", "tfs.tsv"))

# ---- Check tf_activity.tsv ----
res_path <- file.path(run_repo, "results", "run1", "tf_activity.tsv")
if (!file.exists(res_path)) {
  gate(FALSE, "results/run1/tf_activity.tsv is missing")
}
res <- tryCatch(read.delim(res_path), error = function(e) NULL)
if (is.null(res) || !all(c("source", "score") %in% names(res))) {
  gate(
    FALSE,
    "tf_activity.tsv is not a tab-separated table with source and score"
  )
}
res <- res[!is.na(suppressWarnings(as.numeric(res$score))), ]
top <- head(res$source[order(-as.numeric(res$score))], TOP_N)
summary <- sprintf(
  "%d TFs scored (truth %d); top %d: %s",
  nrow(res),
  nrow(truth),
  TOP_N,
  toString(top)
)
problems <- c(
  if (nrow(res) < MIN_TFS) sprintf("fewer than %d TFs scored", MIN_TFS),
  if (!MARKER_TF %in% top) sprintf("%s not in the top %d", MARKER_TF, TOP_N)
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
