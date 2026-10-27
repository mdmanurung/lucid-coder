#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 4
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 2
#$ q3 change checks (yes/no; would you add, remove or change a check?): no
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Differential expression between disease and healthy donors.
#
# Cells from the same donor are not independent replicates, so testing at the
# cell level would inflate significance. Instead, counts are summed per donor
# (pseudobulk), giving 3 healthy vs 3 disease samples, and each gene is fitted
# with a quasi-Poisson GLM with a log library-size offset:
#   counts ~ condition + offset(log(library size))
# The quasi-Poisson dispersion absorbs between-donor variability; the condition
# effect is tested with an F-test (residual df = 4).
#
# Usage: Rscript scripts/01_compare-conditions.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/01_compare-conditions.R <run_dir>")
}
run_dir <- args[[1]]
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read_tsv(here::here("data", "counts.tsv"), show_col_types = FALSE)
cells <- read_tsv(here::here("data", "cells.tsv"), show_col_types = FALSE)

# Every count column must have exactly one metadata row, and vice versa.
cell_ids <- setdiff(names(counts), "gene")
stopifnot(
  !anyDuplicated(cells$cell_id),
  setequal(cell_ids, cells$cell_id),
  all(cells$condition %in% c("healthy", "disease")),
  # each donor belongs to a single condition
  all(count(distinct(cells, donor_id, condition), donor_id)$n == 1)
)

# Pseudobulk: sum raw UMI counts per gene per donor.
pseudobulk <- counts |>
  pivot_longer(-gene, names_to = "cell_id", values_to = "count") |>
  inner_join(cells, by = "cell_id") |>
  group_by(gene, donor_id, condition) |>
  summarise(count = sum(count), .groups = "drop")

donors <- pseudobulk |>
  group_by(donor_id, condition) |>
  summarise(lib_size = sum(count), .groups = "drop") |>
  mutate(condition = factor(condition, levels = c("healthy", "disease")))

message("Donors per condition:")
print(table(donors$condition))

fit_gene <- function(df) {
  df <- inner_join(select(df, donor_id, count), donors, by = "donor_id")
  fit <- glm(
    count ~ condition + offset(log(lib_size)),
    family = quasipoisson(),
    data = df
  )
  test <- anova(fit, test = "F")
  tibble(
    log2_fc = unname(coef(fit)["conditiondisease"]) / log(2),
    p_value = test[["Pr(>F)"]][2]
  )
}

de <- pseudobulk |>
  group_by(gene) |>
  group_modify(~ fit_gene(.x)) |>
  ungroup() |>
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  # keep genes in input order
  slice(match(counts$gene, gene)) |>
  select(gene, log2_fc, p_value, q_value)

write_tsv(de, file.path(run_dir, "de.tsv"))

message(sprintf(
  "Wrote %d genes to %s; %d with q < 0.05",
  nrow(de), file.path(run_dir, "de.tsv"), sum(de$q_value < 0.05, na.rm = TRUE)
))
