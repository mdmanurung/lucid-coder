# Purpose: gate for task 08. Pathway scores must come from the signed Wald statistic: JAK-STAT
#          matches the ULM score on `stat`, and NFkB is positive and significant.
# Inputs:  <run_repo>/results/run1/pathways.tsv, <task_dir>/truth/pathways.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Truth: decoupler 2.3.1 ulm on `stat` gives JAK-STAT 16.73 and NFkB 2.78 (padj 0.036);
#          sc-best-practices reports 17.15 and 2.78 with an older PROGENy copy. log2FoldChange
#          gives JAK-STAT 23.47 (NFkB not significant), |stat| 17.85, -log10 p 18.65 and signed
#          -log10 p 18.97, so JAK_TOLERANCE separates the Wald stat from every other vector.

# ---- Setup ----
JAK_TOLERANCE <- 0.5
PADJ_CUTOFF <- 0.05

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
truth <- read.delim(file.path(task_dir, "truth", "pathways.tsv"))
truth_jak <- truth$score[truth$pathway == "JAK-STAT"]

# ---- Check pathways.tsv ----
res_path <- file.path(run_repo, "results", "run1", "pathways.tsv")
if (!file.exists(res_path)) {
  gate(FALSE, "results/run1/pathways.tsv is missing")
}
res <- tryCatch(read.delim(res_path, check.names = FALSE), error = function(e) {
  NULL
})
if (is.null(res) || !all(c("pathway", "score", "padj") %in% names(res))) {
  gate(
    FALSE,
    "pathways.tsv is not a tab-separated table with pathway, score, padj"
  )
}
score <- setNames(suppressWarnings(as.numeric(res$score)), res$pathway)
padj <- setNames(suppressWarnings(as.numeric(res$padj)), res$pathway)
jak <- score["JAK-STAT"]
nfkb <- score["NFkB"]
summary <- sprintf(
  "JAK-STAT score %.2f (Wald stat gives %.2f); NFkB score %.2f, padj %.3g",
  jak,
  truth_jak,
  nfkb,
  padj["NFkB"]
)
problems <- c(
  if (is.na(jak) || abs(jak - truth_jak) > JAK_TOLERANCE) {
    "JAK-STAT not scored on the Wald stat"
  },
  if (
    is.na(nfkb) ||
      nfkb <= 0 ||
      is.na(padj["NFkB"]) ||
      padj["NFkB"] >= PADJ_CUTOFF
  ) {
    "NFkB not positive and significant"
  }
)
gate(length(problems) == 0, paste(c(problems, summary), collapse = "; "))
