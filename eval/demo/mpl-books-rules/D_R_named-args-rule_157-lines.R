# Purpose: Test rural-to-urban linear and quadratic trends in subcluster frequencies after MPL stimulation (limma).
# Inputs: data/MPL_data.qs [gap], data/MPL_metadata.rds [gap], src/weightedcontrasts.R
# Outputs: data/MPL_cluster_median_markers.qs, outputs/MPL_limma_trend_test.rds, outputs/session_info.txt
# Run: Rscript src/MPL_trend_test.R
# Notes: frequencies are per lineage (subcluster cells / lineage cells), not per total cells.
#   Group contrasts are polynomial, weighted by group size among the modelled patients.
#   Multiple-testing family: all subclusters, BH, separately for the linear and quadratic term.

options(warn = 1)

library(tidyverse)
library(limma)
library(here)

source(here("src", "weightedcontrasts.R"))

excluded_patients <- "LD413"
covariates <- c("group2", "age", "sex", "batch")

out_dir <- here("outputs")
dir.create(out_dir, showWarnings = FALSE)

# Load data ----
cells <- qs::qread(here("data", "MPL_data.qs"))
metadata <- readRDS(here("data", "MPL_metadata.rds"))

# Median marker intensity per subcluster ----
cluster_median_markers <- cells |>
  summarise(
    across(CD45:IL10, median),
    .by = subclusterPG_merged
  )
qs::qsave(cluster_median_markers, here("data", "MPL_cluster_median_markers.qs"))

# Per-lineage subcluster frequencies ----
patient_batch <- distinct(cells, patient_id, batch)
stopifnot("each patient is in one batch" = !anyDuplicated(patient_batch$patient_id))

# A subcluster absent in a patient has proportion 0, and asin(sqrt(0)) = 0,
# so filling 0 in the matrix below matches the transformed scale.
cluster_freq <- cells |>
  count(patient_id, subclusterPG_merged, name = "n_cells") |>
  separate_wider_delim(
    subclusterPG_merged,
    delim = "_",
    names = c("lineage", "annotation"),
    too_many = "merge",
    cols_remove = FALSE
  ) |>
  mutate(
    n_lineage = sum(n_cells),
    .by = c(patient_id, lineage)
  ) |>
  mutate(
    freq = n_cells / n_lineage * 100,
    asin_freq = asin(sqrt(freq / 100))
  ) |>
  left_join(
    patient_batch,
    join_by(patient_id),
    relationship = "many-to-one"
  ) |>
  left_join(
    select(metadata, patient_id, group2, age, sex),
    join_by(patient_id),
    relationship = "many-to-one"
  )

n_before <- n_distinct(cluster_freq$patient_id)
cluster_freq <- drop_na(cluster_freq, all_of(covariates))
message(sprintf(
  "Missing covariates: dropped %d of %d patients",
  n_before - n_distinct(cluster_freq$patient_id),
  n_before
))

n_before <- n_distinct(cluster_freq$patient_id)
cluster_freq <- filter(cluster_freq, !patient_id %in% excluded_patients)
message(sprintf(
  "Excluded patients: dropped %d of %d patients",
  n_before - n_distinct(cluster_freq$patient_id),
  n_before
))

# Subcluster x patient matrix ----
patients <- cluster_freq |>
  distinct(patient_id, group2, age, sex, batch) |>
  mutate(batch = fct_drop(as_factor(batch)))

freq_mat <- cluster_freq |>
  select(subclusterPG_merged, patient_id, asin_freq) |>
  pivot_wider(
    names_from = patient_id,
    values_from = asin_freq,
    values_fill = 0
  ) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()
freq_mat <- freq_mat[, patients$patient_id]

# limma trend test ----
# Weights come from the modelled patients, so they are set after filtering.
contrasts(patients$group2) <- contr.poly.weighted(patients$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = patients)

fit <- freq_mat |>
  lmFit(
    design = design,
    method = "robust",
    maxit = 100
  ) |>
  eBayes(robust = TRUE)

trend_linear <- topTable(
  fit,
  coef = "group2.L",
  number = Inf,
  sort.by = "none"
) |>
  rownames_to_column("cluster") |>
  rename(
    pval_L = P.Value,
    fdr_L = adj.P.Val
  )

trend_quadratic <- topTable(
  fit,
  coef = "group2.Q",
  number = Inf,
  sort.by = "none"
) |>
  rownames_to_column("cluster") |>
  select(
    cluster,
    pval_Q = P.Value,
    fdr_Q = adj.P.Val
  )

cluster_lineage <- distinct(cluster_freq, cluster = subclusterPG_merged, lineage)

trend_results <- trend_linear |>
  left_join(
    trend_quadratic,
    join_by(cluster),
    relationship = "one-to-one"
  ) |>
  left_join(
    cluster_lineage,
    join_by(cluster),
    relationship = "one-to-one"
  ) |>
  relocate(lineage, cluster) |>
  arrange(pval_L)

saveRDS(trend_results, file.path(out_dir, "MPL_limma_trend_test.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
