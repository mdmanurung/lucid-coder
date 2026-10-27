#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ ---- your #$ comments: on their own line above the code, or after the code on its line ----
# Purpose: Test each gene for differential expression, disease vs healthy, with donor as the unit.
# Inputs:  data/counts.tsv (raw UMI, genes x cells), data/cells.tsv (cell_id, donor_id, condition) [gap]
# Outputs: <run_dir>/de.tsv, <run_dir>/pseudobulk_counts.tsv, <run_dir>/session_info.txt
# Run:     Rscript scripts/01_compare-conditions.R results/run1
# Notes:   Cells from one donor are not independent, so counts are summed per donor (pseudobulk)
#          and each gene gets a quasi-Poisson GLM on the 6 donor totals with a log library-size
#          offset; the quasi dispersion absorbs between-donor variability (F-test, 4 residual df).

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(forcats)

condition_levels <- c("healthy", "disease") # healthy is the reference: log2_fc = disease / healthy

run_dir <- commandArgs(trailingOnly = TRUE)[[1]]
stopifnot("run dir already exists; refusing to overwrite" = !dir.exists(run_dir))
dir.create(run_dir, recursive = TRUE)

# Load data ----
counts_wide <- read_tsv(here::here("data/counts.tsv"), show_col_types = FALSE)
cells <- read_tsv(here::here("data/cells.tsv"), show_col_types = FALSE)

count_mat <- as.matrix(select(counts_wide, -gene))
rownames(count_mat) <- counts_wide[["gene"]]
cells <- cells[match(colnames(count_mat), cells[["cell_id"]]), ]
stopifnot(
  "cell metadata matches count columns" = identical(cells[["cell_id"]], colnames(count_mat))
)

# Pseudobulk per donor ----
donor_mat <- t(rowsum(t(count_mat), group = cells[["donor_id"]]))

donors <- cells |>
  distinct(donor_id, condition) |>
  mutate(condition = fct(condition, levels = condition_levels))

pseudobulk <- donor_mat |>
  as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "donor_id", values_to = "count") |>
  mutate(lib_size = sum(count), .by = donor_id) |>
  left_join(donors, join_by(donor_id), relationship = "many-to-one")

write_tsv(pseudobulk, file.path(run_dir, "pseudobulk_counts.tsv"))

# Test per gene ----
de <- pseudobulk |>
  split(~gene) |>
  map(\(gene_data) {
    fit <- glm(
      count ~ condition + offset(log(lib_size)),
      family = quasipoisson(),
      data = gene_data
    )
    tibble(
      log2_fc = coef(fit)[["conditiondisease"]] / log(2),
      p_value = anova(fit, test = "F")["condition", "Pr(>F)"]
    )
  }) |>
  list_rbind(names_to = "gene") |>
  # Multiple-testing family: all 50 genes in this one contrast
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

write_tsv(de, file.path(run_dir, "de.tsv"))

writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info.txt"))
