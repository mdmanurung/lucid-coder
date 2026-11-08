# Purpose: gate for task 08b. The agent runs the DE from counts; pathway scores must come
#          from the signed Wald statistic: JAK-STAT inside the band, NFkB positive and
#          significant.
# Inputs:  <run_repo>/results/run1/pathways.tsv, <task_dir>/truth/pathways.tsv
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Truth: PyDESeq2 0.5.4 ~ label, decoupler 2.3.1 ulm on `stat` gives JAK-STAT 16.72
#          and NFkB 2.78 (padj 0.036). The prompt does not name the design, so the band
#          [JAK_LOW, JAK_HIGH] covers both defensible fits: unpaired ~ label (16.72) and
#          donor-paired ~ replicate + label (14.89).
#          Audit of wrong vectors (JAK-STAT; NFkB padj): ~ label gives log2FoldChange 23.47
#          (0.152), |stat| 17.82, -log10 p 18.61, signed -log10 p 18.93 (0.046); log-CPM
#          mean difference 24.82, raw-count mean difference 9.71. Under ~ replicate + label,
#          |stat| 15.07 (0.155), -log10 p 15.92 (0.195) and signed -log10 p 16.58 (0.095)
#          fall inside the band, so the NFkB clause is what fails them: no wrong vector
#          passes both clauses.

# ---- Setup ----
JAK_LOW <- 14.3
JAK_HIGH <- 17.3
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
  "JAK-STAT score %.2f (band %.1f-%.1f, Wald stat gives %.2f); NFkB score %.2f, padj %.3g",
  jak,
  JAK_LOW,
  JAK_HIGH,
  truth_jak,
  nfkb,
  padj["NFkB"]
)
problems <- c(
  if (is.na(jak) || jak < JAK_LOW || jak > JAK_HIGH) {
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
