# Purpose: test which genes respond to IFN-beta differently in CD14+ monocytes than in CD4 T
#          cells, with one test of the cell type x condition interaction per gene.
# Run:     Rscript scripts/14_compare-responses.R results/run1
# Inputs:  data/counts.tsv    [gap: pseudobulk counts provided with the task]
#          data/samples.tsv   [gap: sample sheet provided with the task]
# Outputs: <run_dir>/response_difference.tsv   one row per gene tested: interaction log2 FC,
#                                              p, BH q over all genes tested, differs
#          <run_dir>/session_info.txt          R session info

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(forcats)
})

# Parameters -------------------------------------------------------------------
# Reference levels first, so the interaction reads (monocyte response) - (CD4 T response).
CELL_TYPE_LEVELS <- c("CD4 T cells", "CD14+ Monocytes")
CONDITION_LEVELS <- c("ctrl", "stim")
# Donor as a block: each donor gives all four samples, so donor effects are paired out.
DESIGN_FORMULA <- ~ donor_id + cell_type * condition
FDR_CUTOFF <- 0.05

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
if (dir.exists(run_dir)) {
  stop("run_dir ", run_dir, " already exists; pick a new name.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Load data --------------------------------------------------------------------
# Samples aligned to the count columns by name.
counts <- read.delim(here::here("data", "counts.tsv"), check.names = FALSE)
count_matrix <- as.matrix(counts[-1])
rownames(count_matrix) <- counts$gene
samples <- read.delim(here::here("data", "samples.tsv")) |>
  mutate(
    cell_type = fct(cell_type, levels = CELL_TYPE_LEVELS),
    condition = fct(condition, levels = CONDITION_LEVELS)
  )
samples <- samples[match(colnames(count_matrix), samples$sample_id), ]
stopifnot(
  "samples match count columns" = identical(
    samples$sample_id,
    colnames(count_matrix)
  )
)

# Filter genes -----------------------------------------------------------------
# edgeR's own filter for this design; the count removed is reported when nonzero.
design <- model.matrix(DESIGN_FORMULA, samples)
dge <- DGEList(count_matrix)
keep <- filterByExpr(dge, design)
if (any(!keep)) {
  message("filterByExpr removed ", sum(!keep), " of ", length(keep), " genes")
}
dge <- normLibSizes(dge[keep, , keep.lib.sizes = FALSE])

# Test the interaction ---------------------------------------------------------
# The last coefficient is cell_typeCD14+ Monocytes:conditionstim. BH over all genes tested.
interaction <- dge |>
  estimateDisp(design) |>
  glmQLFit(design) |>
  glmQLFTest(coef = ncol(design)) |>
  topTags(n = Inf)
response_difference <- interaction$table |>
  as_tibble(rownames = "gene") |>
  transmute(
    gene,
    log2_fc_difference = logFC,
    p_value = PValue,
    q_value = FDR,
    differs = q_value < FDR_CUTOFF
  )

# Write outputs ----------------------------------------------------------------
readr::write_tsv(
  response_difference,
  file.path(run_dir, "response_difference.tsv")
)

# Session info -----------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(run_dir, "session_info.txt")
)
