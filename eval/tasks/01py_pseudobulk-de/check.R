# Purpose: gate for task 01. The donor-level-null gene must not be called; the true gene must be.
# Inputs:  <run_repo>/results/run1/de.tsv, <task_dir>/truth/genes.tsv,
#          <task_dir>/inputs/data/counts.tsv (gene count)
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).

# ---- Setup ----
Q_CUTOFF <- 0.05

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
truth <- read.delim(file.path(task_dir, "truth", "genes.tsv"))
pseudorep_gene <- truth$gene[truth$role == "pseudoreplicated"]
true_gene <- truth$gene[truth$role == "true_effect"]
n_genes <- nrow(read.delim(file.path(task_dir, "inputs", "data", "counts.tsv")))

# ---- Check de.tsv ----
de_path <- file.path(run_repo, "results", "run1", "de.tsv")
if (!file.exists(de_path)) {
  gate(FALSE, "results/run1/de.tsv is missing")
}
de <- tryCatch(read.delim(de_path), error = function(e) NULL)
if (is.null(de) || !all(c("gene", "q_value") %in% names(de))) {
  gate(FALSE, "de.tsv is not a tab-separated table with columns gene and q_value")
}
if (nrow(de) != n_genes) {
  gate(FALSE, sprintf("de.tsv has %d rows, expected %d", nrow(de), n_genes))
}
q_value <- setNames(suppressWarnings(as.numeric(de$q_value)), de$gene)
q_pseudorep <- q_value[pseudorep_gene]
q_true <- q_value[true_gene]
summary <- sprintf(
  "%s (no donor-level effect) q = %.3g; %s (true effect) q = %.3g",
  pseudorep_gene,
  q_pseudorep,
  true_gene,
  q_true
)
problems <- c(
  if (is.na(q_pseudorep) || q_pseudorep < Q_CUTOFF) paste(pseudorep_gene, "called or missing"),
  if (is.na(q_true) || q_true >= Q_CUTOFF) paste(true_gene, "not called or missing")
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
