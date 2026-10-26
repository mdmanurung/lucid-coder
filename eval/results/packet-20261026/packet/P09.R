# Purpose: Test each gene for differential expression between disease and healthy donors.
#   Cells from one donor are not independent, so counts are summed per donor (pseudobulk)
#   and the donor is the replication unit (3 vs 3).
# Inputs:  data/counts.tsv  raw UMI counts, genes x cells, first column `gene` [gap: producer]
#          data/cells.tsv   cell_id, donor_id, condition [gap: producer]
# Outputs: <run_dir>/de.tsv               gene, log2_fc (disease vs healthy), p_value, q_value
#          <run_dir>/pseudobulk_counts.tsv summed counts per gene and donor
#          <run_dir>/session_info.txt
# Run:     Rscript scripts/01_compare-conditions.R results/run1
# Notes:   Model per gene: quasi-Poisson GLM on donor sums with log(donor library size) as
#          offset; the quasi dispersion absorbs between-donor variation. p-values from the
#          Wald t-test on the condition coefficient (residual df = 6 donors - 2 = 4).
#          Multiple-testing family: all genes in counts.tsv, Benjamini-Hochberg.

# Setup ------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)

CONDITION_LEVELS <- c("healthy", "disease") # healthy is the reference, so log2_fc is disease vs healthy
COUNTS_PATH <- here::here("data", "counts.tsv")
CELLS_PATH <- here::here("data", "cells.tsv")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Arguments: expected exactly one run directory; pass e.g. results/run1.", call. = FALSE)
}
run_dir <- here::here(args[[1]])
if (dir.exists(run_dir)) {
  stop("run_dir: ", run_dir, " already exists; choose a new run directory.", call. = FALSE)
}
dir.create(run_dir, recursive = TRUE)

# Load inputs ------------------------------------------------------------------

counts_wide <- read_tsv(COUNTS_PATH, col_types = cols(gene = col_character(), .default = col_double()))
cells <- read_tsv(CELLS_PATH, col_types = cols(.default = col_character()))

missing_cell_cols <- setdiff(c("cell_id", "donor_id", "condition"), names(cells))
if (length(missing_cell_cols) > 0) {
  stop("cells.tsv: missing columns ", toString(missing_cell_cols), "; check the input.", call. = FALSE)
}

# Every count column must have metadata, and vice versa, or cells would silently drop
count_cell_ids <- setdiff(names(counts_wide), "gene")
if (!setequal(count_cell_ids, cells$cell_id) || anyDuplicated(cells$cell_id) > 0) {
  stop("cells.tsv: cell IDs do not match count columns one-to-one; fix the metadata.", call. = FALSE)
}

# Pseudobulk -------------------------------------------------------------------

counts_long <- counts_wide |>
  pivot_longer(-gene, names_to = "cell_id", values_to = "umi_count")

counts_annotated <- counts_long |>
  inner_join(cells, by = "cell_id", relationship = "many-to-one")
if (nrow(counts_annotated) != nrow(counts_long)) {
  stop("Join: rows lost joining counts to cells.tsv; check cell IDs.", call. = FALSE)
}

pseudobulk <- counts_annotated |>
  group_by(gene, donor_id, condition) |>
  summarise(umi_count = sum(umi_count), n_cells = n(), .groups = "drop")

donor_library_sizes <- pseudobulk |>
  group_by(donor_id) |>
  summarise(library_size = sum(umi_count), .groups = "drop")

pseudobulk <- pseudobulk |>
  left_join(donor_library_sizes, by = "donor_id") |>
  mutate(condition = factor(condition, levels = CONDITION_LEVELS))

# Differential expression -------------------------------------------------------

fit_gene <- function(gene_counts) {
  fit <- glm(
    umi_count ~ condition + offset(log(library_size)),
    family = quasipoisson(),
    data = gene_counts
  )
  coefficient <- summary(fit)$coefficients["conditiondisease", ]
  tibble(
    log2_fc = coefficient[["Estimate"]] / log(2),
    p_value = coefficient[["Pr(>|t|)"]]
  )
}

de_results <- pseudobulk |>
  group_by(gene) |>
  group_modify(\(gene_counts, key) fit_gene(gene_counts)) |>
  ungroup() |>
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

# Write outputs ----------------------------------------------------------------

if (nrow(de_results) != nrow(counts_wide)) {
  stop("de_results: not one row per gene; check the model fits.", call. = FALSE)
}
write_tsv(de_results, file.path(run_dir, "de.tsv"))
write_tsv(
  select(pseudobulk, gene, donor_id, condition, n_cells, umi_count, library_size),
  file.path(run_dir, "pseudobulk_counts.tsv")
)

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
