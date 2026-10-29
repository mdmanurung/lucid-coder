# Purpose: freeze task 14 inputs (a seeded gene subset of the Kang 2018 pseudobulk) and the
#          gate truth from an edgeR interaction fit on exactly those inputs.
# Run:     Rscript eval/tasks/14_interaction-contrast/make-inputs/02_freeze-inputs.R <work_dir>
# Inputs:  <work_dir>/pseudobulk_counts.tsv.gz   [01_pseudobulk.py]
#          <work_dir>/pseudobulk_samples.tsv     [01_pseudobulk.py]
# Outputs: eval/tasks/14_interaction-contrast/inputs/data/counts.tsv    genes x 32 samples
#          eval/tasks/14_interaction-contrast/inputs/data/samples.tsv   sample sheet
#          eval/tasks/14_interaction-contrast/truth/genes.tsv           gate genes and their role
# Notes:   Source: Kang et al. 2018 Nat Biotechnol 36:89, GEO GSE96583, via the pertpy mirror.
#          Neither the mirror nor pertpy states a data licence; NCBI puts no restriction on GEO
#          data. What ships is an aggregate (gene sums per donor x condition x cell type).
#          Reruns are byte-identical (seed 14).

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(forcats)
})

# Parameters -------------------------------------------------------------------
# Keeps counts.tsv well under the 1 MB fixture limit; the known genes are added on top.
N_GENES <- 2000
# Responds more (RABGAP1L) and less (ENO1) in CD14 monocytes than CD4 T (sc-best-practices).
KNOWN_GENES <- c(
  larger_in_monocytes = "RABGAP1L",
  smaller_in_monocytes = "ENO1"
)
FDR_CUTOFF <- 0.05

work_dir <- commandArgs(trailingOnly = TRUE)[[1]]
task_dir <- here::here("eval", "tasks", "14_interaction-contrast")

# Load pseudobulk --------------------------------------------------------------
# Sample IDs short and readable; the sheet keeps patient, cell type and condition.
counts_all <- read.delim(
  file.path(work_dir, "pseudobulk_counts.tsv.gz"),
  check.names = FALSE
)
samples <- read.delim(file.path(work_dir, "pseudobulk_samples.tsv")) |>
  mutate(
    sample_id = paste(
      sub("patient_", "p", replicate),
      label,
      if_else(cell_type == "CD4 T cells", "cd4t", "cd14mono"),
      sep = "_"
    ),
    cell_type = fct(cell_type, levels = c("CD4 T cells", "CD14+ Monocytes")),
    condition = fct(label, levels = c("ctrl", "stim"))
  ) |>
  select(
    sample_id,
    donor_id = replicate,
    cell_type,
    condition,
    n_cells,
    sample
  )
count_matrix <- as.matrix(counts_all[samples$sample])
dimnames(count_matrix) <- list(counts_all$gene, samples$sample_id)

# Pick genes -------------------------------------------------------------------
# A seeded sample of the genes edgeR keeps, so the shipped genes are all testable.
design <- model.matrix(~ donor_id + cell_type * condition, samples)
expressed <- rownames(count_matrix)[filterByExpr(DGEList(count_matrix), design)]
set.seed(14)
genes <- sort(union(sample(expressed, N_GENES), KNOWN_GENES))
count_matrix <- count_matrix[genes, ]

# Fit interaction and per-type responses ---------------------------------------
# Each fit filters with filterByExpr first, as an analyst would; a filtered gene gets NA.
fit_edger <- function(count_matrix, design, coef) {
  dge <- DGEList(count_matrix)
  dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE] |>
    normLibSizes() |>
    estimateDisp(design) |>
    glmQLFit(design) |>
    glmQLFTest(coef = coef) |>
    topTags(n = Inf) |>
    as.data.frame() |>
    as_tibble(rownames = "gene") |>
    select(gene, log2_fc = logFC, log_cpm = logCPM, q = FDR)
}
# Interaction: the last coefficient, (stim - ctrl in monocytes) - (stim - ctrl in CD4 T).
interaction <- fit_edger(count_matrix, design, ncol(design))
# The naive route: one stim-vs-ctrl fit per cell type, then compare which genes are significant.
per_type <- split(samples, samples$cell_type) |>
  lapply(\(type_samples) {
    fit_edger(
      count_matrix[, type_samples$sample_id],
      model.matrix(~ donor_id + condition, type_samples),
      "conditionstim"
    )
  })
gene_stats <- tibble(gene = genes) |>
  left_join(
    select(interaction, gene, interaction_log2_fc = log2_fc, interaction_q = q),
    by = join_by(gene),
    relationship = "one-to-one"
  ) |>
  left_join(
    rename_with(
      per_type[["CD14+ Monocytes"]],
      \(x) paste0("monocyte_", x),
      -gene
    ),
    by = join_by(gene),
    relationship = "one-to-one"
  ) |>
  left_join(
    rename_with(per_type[["CD4 T cells"]], \(x) paste0("cd4t_", x), -gene),
    by = join_by(gene),
    relationship = "one-to-one"
  )

# Pick the decoy ---------------------------------------------------------------
# Significant in one cell type only, kept by every filter, but the responses do not differ: the
# gene a comparison of two DE lists calls and an interaction test does not.
decoy <- gene_stats |>
  filter(
    (monocyte_q < 1e-4 & cd4t_q > 0.25) | (cd4t_q < 1e-4 & monocyte_q > 0.25),
    !is.na(interaction_q)
  ) |>
  slice_max(interaction_q, n = 1, with_ties = FALSE)
truth <- bind_rows(
  tibble(gene = KNOWN_GENES, role = names(KNOWN_GENES)),
  tibble(gene = decoy$gene, role = "significant_in_one_type_only")
) |>
  left_join(gene_stats, by = join_by(gene), relationship = "one-to-one")
n_naive_calls <- sum(
  xor(
    gene_stats$monocyte_q < FDR_CUTOFF,
    gene_stats$cd4t_q < FDR_CUTOFF
  ),
  na.rm = TRUE
)
message(
  "task 14: ",
  sum(gene_stats$interaction_q < FDR_CUTOFF, na.rm = TRUE),
  " genes with an interaction at FDR 0.05; ",
  n_naive_calls,
  " significant in exactly one cell type"
)
print(as.data.frame(truth), digits = 3)

# Write fixtures ---------------------------------------------------------------
write.table(
  data.frame(gene = genes, count_matrix, check.names = FALSE),
  file.path(task_dir, "inputs", "data", "counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  select(samples, -sample),
  file.path(task_dir, "inputs", "data", "samples.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mutate(truth, across(where(is.double), \(x) signif(x, 4))),
  file.path(task_dir, "truth", "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
