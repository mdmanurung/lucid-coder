#!/usr/bin/env Rscript
# Differential expression between disease and healthy donors.
#
# Donors, not cells, are the independent replicates (3 vs 3), so cells are
# summed into one pseudobulk profile per donor before testing. Testing at the
# cell level would treat ~200 cells per donor as independent and give badly
# anti-conservative p-values.
#
# Per gene: quasi-Poisson GLM on the donor pseudobulk counts,
#   counts ~ condition + offset(log(size factor)),
# with an F-test on the condition term (dispersion estimated from
# donor-to-donor variation, 4 residual df). Size factors use the
# median-of-ratios method, which is robust to a minority of truly changing
# genes when there are only 50 genes.
#
# Usage: Rscript scripts/01_compare-conditions.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(magrittr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/01_compare-conditions.R <run_dir>")
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_tsv(here("data", "counts.tsv"), show_col_types = FALSE)
cells  <- read_tsv(here("data", "cells.tsv"),
                   col_types = cols(.default = col_character()))

# ---- Checks -----------------------------------------------------------------
cell_ids <- setdiff(names(counts), "gene")
stopifnot(
  !anyDuplicated(counts$gene),
  !anyDuplicated(cells$cell_id),
  setequal(cell_ids, cells$cell_id),
  all(cells$condition %in% c("healthy", "disease")),
  # each donor belongs to exactly one condition
  all(table(distinct(cells, donor_id, condition)$donor_id) == 1)
)
count_mat <- as.matrix(counts[, cell_ids])
rownames(count_mat) <- counts$gene
stopifnot(!anyNA(count_mat), all(count_mat >= 0), all(count_mat == round(count_mat)))

# ---- Pseudobulk: sum counts per donor ---------------------------------------
donor_of_cell <- cells$donor_id[match(cell_ids, cells$cell_id)]
pb <- t(rowsum(t(count_mat), group = donor_of_cell))   # genes x donors

donors <- cells %>%
  distinct(donor_id, condition) %>%
  mutate(condition = factor(condition, levels = c("healthy", "disease"))) %>%
  arrange(match(donor_id, colnames(pb)))
stopifnot(identical(donors$donor_id, colnames(pb)))

# ---- Size factors (median of ratios) ----------------------------------------
expressed <- rowSums(pb == 0) == 0
log_geo_mean <- rowMeans(log(pb[expressed, , drop = FALSE]))
size_factors <- apply(pb[expressed, , drop = FALSE], 2,
                      function(x) exp(median(log(x) - log_geo_mean)))
size_factors <- size_factors / exp(mean(log(size_factors)))

# ---- Per-gene quasi-Poisson GLM ---------------------------------------------
fit_gene <- function(y) {
  if (all(y == 0)) return(tibble(log2_fc = NA_real_, p_value = NA_real_))
  d <- data.frame(y = y, condition = donors$condition, sf = size_factors)
  fit  <- glm(y ~ condition + offset(log(sf)), family = quasipoisson(), data = d)
  fit0 <- glm(y ~ 1 + offset(log(sf)), family = quasipoisson(), data = d)
  tibble(
    log2_fc = unname(coef(fit)["conditiondisease"]) / log(2),
    p_value = anova(fit0, fit, test = "F")[2, "Pr(>F)"]
  )
}

de <- map_dfr(seq_len(nrow(pb)), function(i) fit_gene(pb[i, ])) %>%
  mutate(gene = rownames(pb), .before = 1) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_tsv(de, file.path(run_dir, "de.tsv"), na = "NA")

# Supporting output: donor pseudobulk counts with size factors
pseudobulk_long <- as_tibble(pb, rownames = "gene") %>%
  pivot_longer(-gene, names_to = "donor_id", values_to = "count") %>%
  left_join(mutate(donors, size_factor = size_factors), by = "donor_id")
write_tsv(pseudobulk_long, file.path(run_dir, "pseudobulk_counts.tsv"))

message("Genes tested: ", nrow(de),
        "; q < 0.05: ", sum(de$q_value < 0.05, na.rm = TRUE))
