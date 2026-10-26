# Pseudobulk DE: disease vs healthy, donor is the unit of replication.
# Per-cell counts are summed per donor; each gene is fit with a quasi-Poisson
# GLM with a log library-size offset (3 vs 3 donors, df = 4).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

out_dir <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_tsv(here::here("data", "counts.tsv"), show_col_types = FALSE)
cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)

mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$gene
stopifnot(setequal(colnames(mat), cells$cell_id))
cells <- cells[match(colnames(mat), cells$cell_id), ]

# Sum cells -> donors (genes x donors)
pb <- t(rowsum(t(mat), cells$donor_id))
donors <- distinct(cells, donor_id, condition)
donors <- donors[match(colnames(pb), donors$donor_id), ]
condition <- factor(donors$condition, levels = c("healthy", "disease"))
lib <- colSums(pb)

fit_gene <- function(y) {
  if (all(y == 0)) return(c(NA_real_, NA_real_))
  co <- summary(glm(y ~ condition + offset(log(lib)), family = quasipoisson()))$coefficients
  c(co["conditiondisease", "Estimate"] / log(2), co["conditiondisease", "Pr(>|t|)"])
}

res <- t(apply(pb, 1, fit_gene))
de <- tibble(gene = rownames(pb), log2_fc = res[, 1], p_value = res[, 2]) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

write_tsv(de, file.path(out_dir, "de.tsv"))
