options(warn = 1)

# Purpose: Test which genes differ between disease and healthy donors, using
#          donor-level pseudobulk counts.
# Inputs:  data/counts.tsv  raw UMI counts, genes x cells (first column `gene`) [gap]
#          data/cells.tsv   cell_id, donor_id, condition (healthy / disease) [gap]
# Outputs: <out_dir>/de.tsv            gene, log2_fc (disease vs healthy), p_value, q_value
#          <out_dir>/pseudobulk_counts.tsv  summed counts per gene and donor
#          <out_dir>/session_info.txt
# Run:     Rscript scripts/01_compare-conditions.R results/run1
# Notes:   The donor is the replication unit (3 vs 3), so cells are summed per donor
#          before testing; a cell-level test would treat ~200 cells per donor as
#          independent and give far too small p-values. Each gene is fitted with a
#          quasi-Poisson GLM, count ~ condition + offset(log(library size)), because
#          edgeR / DESeq2 are not available; quasi-Poisson lets the variance exceed
#          the Poisson variance across donors. p-values come from the Wald t-test on
#          the condition coefficient (4 residual df). BH family: all genes in the input.

# Setup ----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# Reference level first so the coefficient is disease vs healthy.
CONDITION_LEVELS <- c("healthy", "disease")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Arguments: expected one run directory; run as `Rscript <script> results/<run>`.", call. = FALSE)
}
out_dir <- here::here(args[[1]])
if (dir.exists(out_dir)) {
  stop("out_dir: ", out_dir, " already exists; choose a new run directory.", call. = FALSE)
}

# Load inputs ----------------------------------------------------------------

counts_wide <- read_tsv(here::here("data", "counts.tsv"), show_col_types = FALSE)
cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)

missing_cols <- setdiff(c("cell_id", "donor_id", "condition"), names(cells))
if (length(missing_cols) > 0) {
  stop("cells.tsv: missing columns ", toString(missing_cols), "; check the input.", call. = FALSE)
}

# Pseudobulk per donor -------------------------------------------------------

counts_long <- counts_wide |>
  pivot_longer(-gene, names_to = "cell_id", values_to = "count")

# inner_join would silently drop cells without metadata, so join left and require a match.
counts_annotated <- counts_long |>
  left_join(cells, by = "cell_id", relationship = "many-to-one")
if (nrow(counts_annotated) != nrow(counts_long) || anyNA(counts_annotated$condition)) {
  stop("counts.tsv: some cell IDs are missing from cells.tsv; fix the metadata.", call. = FALSE)
}

pseudobulk <- counts_annotated |>
  summarise(count = sum(count), .by = c(gene, donor_id, condition)) |>
  mutate(
    library_size = sum(count),
    condition = factor(condition, levels = CONDITION_LEVELS),
    .by = donor_id
  )

# A gene with zero counts in a whole condition has an infinite fold change.
genes_zero_in_group <- pseudobulk |>
  summarise(group_count = sum(count), .by = c(gene, condition)) |>
  filter(group_count == 0) |>
  pull(gene)
if (length(genes_zero_in_group) > 0) {
  stop("pseudobulk: zero counts in a condition for ", toString(genes_zero_in_group),
       "; filter these genes first.", call. = FALSE)
}

# Test each gene -------------------------------------------------------------

fit_gene <- function(gene_counts) {
  fit <- glm(
    count ~ condition + offset(log(library_size)),
    family = quasipoisson(),
    data = gene_counts
  )
  coef_row <- summary(fit)$coefficients["conditiondisease", ]
  tibble(log2_fc = coef_row[["Estimate"]] / log(2), p_value = coef_row[["Pr(>|t|)"]])
}

de_results <- pseudobulk |>
  nest(gene_counts = -gene) |>
  mutate(test = map(gene_counts, fit_gene)) |>
  select(gene, test) |>
  unnest(test) |>
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

# Write outputs --------------------------------------------------------------

if (nrow(de_results) == 0) {
  stop("de_results: no genes were tested; check counts.tsv.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE)
write_tsv(de_results, file.path(out_dir, "de.tsv"))
write_tsv(
  select(pseudobulk, gene, donor_id, condition, count, library_size),
  file.path(out_dir, "pseudobulk_counts.tsv")
)

# Session info ---------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
