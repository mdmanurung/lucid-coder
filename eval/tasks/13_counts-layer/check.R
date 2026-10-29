# Purpose: gate for task 13. Library size and genes detected per cell must come from raw
#          counts: X in lab_a.h5ad, layers["counts"] in lab_b.h5ad (whose X is transformed).
# Inputs:  <run_repo>/results/run1/cells.tsv, <task_dir>/truth/cells.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Rule: every (lab, cell_id) present; library_size and n_genes equal truth exactly.
#          Reading lab_b's X gives float library sizes and n_genes = 128 (X has no zeros).
#          Inputs are verbatim nf-core scdownstream fixtures; make-inputs.py builds the truth.
#          Arms that load the plugin see the skills/sc-normalize/SKILL.md gotcha on checking an
#          unknown h5ad's format, so this task partly measures whether that note is used.
#          Known soft spot: the prompt also asks for combined.h5ad (log1p CP10k in X, counts
#          in layers["counts"]), which is not checked.

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
truth <- read.delim(file.path(task_dir, "truth", "cells.tsv"))

# ---- Check cells.tsv ----
cells_path <- file.path(run_repo, "results", "run1", "cells.tsv")
if (!file.exists(cells_path)) {
  gate(FALSE, "results/run1/cells.tsv is missing")
}
cells <- tryCatch(read.delim(cells_path), error = function(e) NULL)
if (
  is.null(cells) ||
    !all(c("lab", "cell_id", "library_size", "n_genes") %in% names(cells))
) {
  gate(
    FALSE,
    "cells.tsv is not a tab-separated table with lab, cell_id, library_size, n_genes"
  )
}
reported <- cells[
  match(paste(truth$lab, truth$cell_id), paste(cells$lab, cells$cell_id)),
]
if (anyNA(reported$cell_id)) {
  gate(
    FALSE,
    sprintf(
      "cells.tsv lacks %d of %d cells",
      sum(is.na(reported$cell_id)),
      nrow(truth)
    )
  )
}

# ---- Compare per lab ----
wrong_size <- abs(
  suppressWarnings(as.numeric(reported$library_size)) - truth$library_size
) >
  1e-6
wrong_genes <- suppressWarnings(as.numeric(reported$n_genes)) != truth$n_genes
wrong_size[is.na(wrong_size)] <- TRUE
wrong_genes[is.na(wrong_genes)] <- TRUE
summary <- vapply(
  split(seq_len(nrow(truth)), truth$lab),
  function(i) {
    sprintf(
      "%s: %d/%d library_size wrong, %d/%d n_genes wrong (%s: %s vs %d)",
      truth$lab[i[[1]]],
      sum(wrong_size[i]),
      length(i),
      sum(wrong_genes[i]),
      length(i),
      truth$cell_id[i[[1]]],
      reported$library_size[i[[1]]],
      truth$library_size[i[[1]]]
    )
  },
  ""
)
gate(!any(wrong_size | wrong_genes), paste(summary, collapse = "; "))
