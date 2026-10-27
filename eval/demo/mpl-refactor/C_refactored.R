# Purpose: Test MPL-stimulated subcluster frequencies (per lineage) for a linear and quadratic
#          trend across the rural-urban gradient with limma, adjusting for age, sex and batch.
# Inputs:  data/MPL_data.qs        per-cell markers and subclusters [Zenodo 14852065]
#          data/MPL_metadata.rds   per-patient group2, age, sex [gap: producer not in repo]
# Outputs: data/MPL_cluster_median_markers.qs   per-subcluster median marker expression
#          outputs/MPL_limma_trend_test.rds     trend statistics per subcluster; read by Fig4-6
# Run:     Rscript src/_MPL_differential_analysis_refactored.R   (from the repository root)
# Notes:   Output paths are fixed because downstream figure scripts read them; both are
#          overwritten on rerun. The median table is written into data/ as in the original.

# Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
})

source(here::here("src/weightedcontrasts.R"))

PATHS <- list(
  data = here::here("data/MPL_data.qs"),
  meta = here::here("data/MPL_metadata.rds"),
  cluster_medians = here::here("data/MPL_cluster_median_markers.qs"),
  limma_results = here::here("outputs/MPL_limma_trend_test.rds")
)

# Threads for decompressing the per-cell table; affects speed only
N_THREADS <- 6
# Reference (control) sample run alongside patients; see _PMA_differential_analysis.R
REFERENCE_SAMPLE <- "LD413"
# Iterations for robust (M-estimation) fitting in lmFit(); raised from the rlm() default of 20
ROBUST_MAXIT <- 100


# Load inputs -------------------------------------------------------------

# Subcluster names are "<lineage>_<annotation>"; the lineage is taken from that prefix
cells <- qs::qread(PATHS$data, nthreads = N_THREADS) |>
  separate(subclusterPG_merged, c("subset", "anno"), "_", remove = FALSE, convert = TRUE) |>
  mutate(lineage = subset) |>
  as.data.table()

# group2 is the gradient used in the model. [gap] Its level order comes from the metadata file;
# the linear (.L) term assumes that order follows the rural-urban gradient.
patient_meta <- readRDS(PATHS$meta) |>
  select(patient_id, group2, age, sex)


# Median marker expression per subcluster ---------------------------------

marker_cols <- colnames(cells[, CD45:IL10])
cluster_medians <- cells[, lapply(.SD, median), .SDcols = marker_cols, by = subclusterPG_merged]


# Subcluster frequencies per lineage --------------------------------------

lineage_counts <- dplyr::count(cells, patient_id, lineage, name = "n_lineage")
cluster_counts <- dplyr::count(cells, patient_id, subclusterPG_merged, name = "n_cluster")
cluster_lineage <- cells |>
  select(lineage, subclusterPG_merged) |>
  distinct()
cluster_counts_lineage <- left_join(cluster_lineage, cluster_counts, by = "subclusterPG_merged")
sample_meta <- cells |>
  select(patient_id, batch, group3) |>
  distinct()

counts_joined <- patient_meta |>
  left_join(sample_meta, by = "patient_id") |>
  left_join(lineage_counts, by = "patient_id") |>
  left_join(cluster_counts_lineage, by = c("patient_id", "lineage")) |>
  mutate(freq_cluster_perLineage = n_cluster / n_lineage * 100) |>
  # Arcsine square-root transform stabilises the variance of proportions for the linear model
  mutate(freq_cluster_perLineage_trans = asin(sqrt(freq_cluster_perLineage / 100)))

# Rows with any missing value (patient absent from cells or missing covariates) are dropped
counts_complete <- drop_na(counts_joined)
message(sprintf(
  "Dropped %d rows with missing values; %d -> %d patients.",
  nrow(counts_joined) - nrow(counts_complete),
  n_distinct(counts_joined$patient_id),
  n_distinct(counts_complete$patient_id)
))

counts <- filter(counts_complete, patient_id != REFERENCE_SAMPLE)
message(sprintf(
  "Removed reference sample %s (%d rows).",
  REFERENCE_SAMPLE,
  nrow(counts_complete) - nrow(counts)
))

# Subclusters a patient lacks get 0, which is also their transformed frequency
freq_matrix <- counts |>
  select(patient_id, subclusterPG_merged, freq_cluster_perLineage_trans) |>
  pivot_wider(
    names_from = patient_id,
    values_from = freq_cluster_perLineage_trans,
    values_fill = 0
  ) |>
  column_to_rownames(var = "subclusterPG_merged")

design_meta <- counts |>
  select(patient_id, group = group2, age, sex, batch) |>
  distinct()

# Design rows are matched to matrix columns by position
if (!identical(as.character(design_meta$patient_id), colnames(freq_matrix))) {
  stop("design_meta: rows do not match freq_matrix columns; check metadata.", call. = FALSE)
}


# limma trend test --------------------------------------------------------

# Polynomial contrasts weighted by group sizes, so the linear term is orthogonal to the intercept
contrasts(design_meta$group) <- contr.poly.weighted(design_meta$group)
design_trend <- model.matrix(~ group + age + sex + batch, data = design_meta)

fit <- lmFit(freq_matrix, design_trend, method = "robust", maxit = ROBUST_MAXIT)
efit <- eBayes(fit, robust = TRUE)

# BH-adjusted across all subclusters, separately for each trend coefficient
trend_table <- function(efit, coef, suffix) {
  topTable(efit, coef = coef, number = Inf) |>
    rownames_to_column(var = "cluster") |>
    as_tibble() |>
    dplyr::rename(!!paste0("pval_", suffix) := P.Value, !!paste0("fdr_", suffix) := adj.P.Val)
}

res_linear <- trend_table(efit, coef = "group.L", suffix = "L")
res_quadratic <- trend_table(efit, coef = "group.Q", suffix = "Q")

trend_results <- cluster_lineage |>
  dplyr::rename(cluster = subclusterPG_merged) |>
  left_join(res_linear, by = "cluster") |>
  left_join(select(res_quadratic, cluster, pval_Q, fdr_Q), by = "cluster")


# Write outputs -----------------------------------------------------------

qs::qsave(cluster_medians, PATHS$cluster_medians)
saveRDS(trend_results, PATHS$limma_results)


# Session info ------------------------------------------------------------

sessionInfo()
