# Purpose: gate for task 14. The response difference must come from one test of the
#          cell type x condition interaction, not from comparing two per-type DE lists.
# Inputs:  <run_repo>/results/run1/response_difference.tsv, <task_dir>/truth/genes.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: RABGAP1L differs with log2_fc_difference > 0; ENO1 differs with < 0; ARHGAP15
#          (significant in CD4 T only, same response in both types) does not differ.
#          No row-count check: a correct answer may filter low-count genes first.

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
truth <- read.delim(file.path(task_dir, "truth", "genes.tsv"))
larger_gene <- truth$gene[truth$role == "larger_in_monocytes"]
smaller_gene <- truth$gene[truth$role == "smaller_in_monocytes"]
decoy_gene <- truth$gene[truth$role == "significant_in_one_type_only"]

# ---- Check response_difference.tsv ----
result_path <- file.path(run_repo, "results", "run1", "response_difference.tsv")
if (!file.exists(result_path)) {
  gate(FALSE, "results/run1/response_difference.tsv is missing")
}
result <- tryCatch(read.delim(result_path), error = function(e) NULL)
if (
  is.null(result) ||
    !all(c("gene", "log2_fc_difference", "differs") %in% names(result))
) {
  gate(
    FALSE,
    "response_difference.tsv lacks columns gene, log2_fc_difference, differs"
  )
}
row <- result[match(truth$gene, result$gene), ]
if (anyNA(row$gene)) {
  gate(FALSE, paste("genes missing:", toString(truth$gene[is.na(row$gene)])))
}
called_text <- toupper(trimws(as.character(row$differs)))
called <- setNames(
  ifelse(
    called_text %in% c("TRUE", "T", "1", "YES"),
    TRUE,
    ifelse(called_text %in% c("FALSE", "F", "0", "NO"), FALSE, NA)
  ),
  truth$gene
)
if (anyNA(called)) {
  gate(
    FALSE,
    paste("differs is not TRUE/FALSE for", toString(truth$gene[is.na(called)]))
  )
}
difference <- setNames(
  suppressWarnings(as.numeric(row$log2_fc_difference)),
  truth$gene
)
summary <- paste(
  sprintf(
    "%s differs = %s, log2_fc_difference = %.3g",
    truth$gene,
    called,
    difference
  ),
  collapse = "; "
)
problems <- c(
  if (!called[[larger_gene]]) paste(larger_gene, "not called"),
  if (called[[larger_gene]] && !isTRUE(difference[[larger_gene]] > 0)) {
    paste(larger_gene, "difference not positive")
  },
  if (!called[[smaller_gene]]) paste(smaller_gene, "not called"),
  if (called[[smaller_gene]] && !isTRUE(difference[[smaller_gene]] < 0)) {
    paste(smaller_gene, "difference not negative")
  },
  if (called[[decoy_gene]]) {
    paste(decoy_gene, "called (significant in one cell type only)")
  }
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
