# Purpose: Test for linear and quadratic rural-to-urban trends in MPL subcluster frequencies (limma).
# Inputs:  data/MPL_data.qs [gap], data/MPL_metadata.rds [gap], src/weightedcontrasts.R
# Outputs: data/MPL_cluster_median_markers.qs, outputs/MPL_limma_trend_test.rds
# Run:     Rscript src/MPL_trend_test.R
options(warn = 1)

library(tidyverse)
library(limma)
library(qs)
library(here)

source(here("src", "weightedcontrasts.R"))

excluded_patients <- "LD413"
model_formula <- ~ group2 + age + sex + batch

out_dir <- here("outputs")
dir.create(out_dir, showWarnings = FALSE)

# Load data ----
cells <- qread(here("data", "MPL_data.qs")) |>
  rename(cluster = subclusterPG_merged) |>
  mutate(lineage = str_extract(cluster, "^[^_]+"))
metadata <- readRDS(here("data", "MPL_metadata.rds"))

# Median marker intensity per subcluster ----
cells |>
  summarise(across(CD45:IL10, median), .by = cluster) |>
  qsave(here("data", "MPL_cluster_median_markers.qs"))

# Subcluster frequencies within lineage ----
# batch is recorded per cell; the one-to-one join fails if a patient spans two batches
patient_covariates <- cells |>
  distinct(patient_id, batch) |>
  left_join(
    metadata |> select(patient_id, group2, age, sex),
    join_by(patient_id),
    relationship = "one-to-one"
  )

n_patients_all <- nrow(patient_covariates)
patient_covariates <- drop_na(patient_covariates)
n_patients_complete <- nrow(patient_covariates)
patient_covariates <- filter(patient_covariates, !patient_id %in% excluded_patients)
message(sprintf(
  "Patients: %d with cells, %d after dropping missing covariates, %d after excluding %s",
  n_patients_all,
  n_patients_complete,
  nrow(patient_covariates),
  str_flatten_comma(excluded_patients)
))

# Every subcluster x patient pair gets a count, 0 if absent; a patient with no cells
# in a whole lineage has 0/0 there, which is set to frequency 0 as specified
cluster_freqs <- cells |>
  filter(patient_id %in% patient_covariates$patient_id) |>
  count(patient_id, lineage, cluster) |>
  complete(patient_id, nesting(lineage, cluster), fill = list(n = 0)) |>
  mutate(n_lineage = sum(n), .by = c(patient_id, lineage)) |>
  mutate(
    freq = if_else(n_lineage == 0, 0, n / n_lineage * 100),
    freq_asin = asin(sqrt(freq / 100))
  )

freq_mat <- cluster_freqs |>
  select(cluster, patient_id, freq_asin) |>
  pivot_wider(names_from = patient_id, values_from = freq_asin) |>
  column_to_rownames("cluster") |>
  as.matrix()

# limma trend test ----
design_data <- patient_covariates |>
  mutate(batch = fct(as.character(batch))) |>
  column_to_rownames("patient_id")
design_data <- design_data[colnames(freq_mat), ]
stopifnot("design rows match matrix columns" = identical(rownames(design_data), colnames(freq_mat)))

# Weights come from the patients in the model, so .L/.Q are orthogonal for this sample
contrasts(design_data$group2) <- contr.poly.weighted(design_data$group2)
design <- model.matrix(model_formula, data = design_data)

fit <- lmFit(freq_mat, design, method = "robust", maxit = 100) |>
  eBayes(robust = TRUE)

linear_res <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)
quadratic_res <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

# FDR is Benjamini-Hochberg across all subclusters, separately per trend term
trend_res <- linear_res |>
  left_join(quadratic_res, join_by(cluster), relationship = "one-to-one") |>
  mutate(lineage = str_extract(cluster, "^[^_]+"), .before = cluster) |>
  arrange(pval_L)

saveRDS(trend_res, file.path(out_dir, "MPL_limma_trend_test.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
