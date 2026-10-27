# Purpose: Test rural-to-urban linear and quadratic trends in MPL subcluster frequencies (limma).
# Inputs:  data/MPL_data.qs [gap], data/MPL_metadata.rds [gap], src/weightedcontrasts.R
# Outputs: data/MPL_cluster_median_markers.qs, outputs/MPL_limma_trend_test.rds
# Run:     Rscript src/MPL_trend_test.R

options(warn = 1)

library(tidyverse)
library(qs)
library(limma)
library(here)

source(here("src", "weightedcontrasts.R"))

excluded_patients <- "LD413"
trend_formula <- ~ group2 + age + sex + batch
out_dir <- here("outputs")

# Load data ----
cells <- qread(here("data", "MPL_data.qs"))
metadata <- readRDS(here("data", "MPL_metadata.rds"))

# Median marker intensity per subcluster ----
cluster_medians <- cells |>
  summarise(across(CD45:IL10, median), .by = subclusterPG_merged)
qsave(cluster_medians, here("data", "MPL_cluster_median_markers.qs"))

# Subcluster frequencies within lineage ----
# batch is a patient-level covariate: a patient in two batches duplicates the key and fails the join
patient_info <- cells |>
  distinct(patient_id, batch) |>
  left_join(
    metadata |>
      select(patient_id, group2, age, sex),
    join_by(patient_id),
    relationship = "one-to-one"
  )

cluster_freq <- cells |>
  count(patient_id, subclusterPG_merged, name = "n_cells") |>
  separate_wider_delim(
    subclusterPG_merged,
    delim = "_",
    names = c("lineage", "cluster"),
    too_many = "merge",
    cols_remove = FALSE
  ) |>
  mutate(freq = n_cells / sum(n_cells) * 100, .by = c(patient_id, lineage)) |>
  mutate(asin_freq = asin(sqrt(freq / 100))) |>
  left_join(patient_info, join_by(patient_id), relationship = "many-to-one")

n_patients_all <- n_distinct(cluster_freq$patient_id)
cluster_freq <- cluster_freq |>
  filter(!patient_id %in% excluded_patients)
n_patients_kept <- n_distinct(cluster_freq$patient_id)
message(sprintf("Excluded patients: %d -> %d patients", n_patients_all, n_patients_kept))
cluster_freq <- cluster_freq |>
  drop_na()
message(sprintf(
  "Missing covariates: %d -> %d patients",
  n_patients_kept,
  n_distinct(cluster_freq$patient_id)
))

# Subcluster x patient matrix; a subcluster absent in a patient is a frequency of 0
freq_mat <- cluster_freq |>
  select(subclusterPG_merged, patient_id, asin_freq) |>
  pivot_wider(names_from = patient_id, values_from = asin_freq, values_fill = 0) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# limma trend test ----
design_data <- patient_info |>
  slice(match(colnames(freq_mat), patient_id)) |>
  mutate(batch = fct(as.character(batch)))
stopifnot(
  "design rows match matrix columns" = identical(design_data$patient_id, colnames(freq_mat))
)

# Weights come from the analysed patients, so set contrasts after all filtering
contrasts(design_data$group2) <- contr.poly.weighted(design_data$group2)
design <- model.matrix(trend_formula, data = design_data)

fit <- lmFit(freq_mat, design, method = "robust", maxit = 100) |>
  eBayes(robust = TRUE)

# BH FDR across all subclusters, separately for each trend term
trend_linear <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "none") |>
  rownames_to_column("subclusterPG_merged") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)
trend_quadratic <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") |>
  rownames_to_column("subclusterPG_merged") |>
  select(subclusterPG_merged, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_results <- trend_linear |>
  left_join(trend_quadratic, join_by(subclusterPG_merged), relationship = "one-to-one") |>
  separate_wider_delim(
    subclusterPG_merged,
    delim = "_",
    names = c("lineage", "cluster"),
    too_many = "merge"
  ) |>
  arrange(pval_L)
saveRDS(trend_results, file.path(out_dir, "MPL_limma_trend_test.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
