# Purpose: Test rural-to-urban linear and quadratic trends in MPL-stimulated subcluster frequencies.
# Inputs:  data/MPL_data.qs [gap], data/MPL_metadata.rds [gap], src/weightedcontrasts.R
# Outputs: data/MPL_cluster_median_markers.qs, outputs/MPL_limma_trend_test.rds
# Run:     Rscript src/MPL_trend_test.R
# Notes:   Frequencies are within lineage (subcluster cells / lineage cells), so a lineage shift
#          does not appear as a change in every subcluster of that lineage.

options(warn = 1)

library(dplyr)
library(tidyr)
library(tibble)
library(limma)

source(here::here("src", "weightedcontrasts.R"))

EXCLUDED_PATIENTS <- "LD413"
DESIGN_FORMULA <- ~ group2 + age + sex + batch

out_dir <- here::here("outputs")
dir.create(out_dir, showWarnings = FALSE)

cells <- qs::qread(here::here("data", "MPL_data.qs"))
metadata <- readRDS(here::here("data", "MPL_metadata.rds"))

# Median marker intensity per subcluster -----------------------------------
cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median)) |>
  qs::qsave(here::here("data", "MPL_cluster_median_markers.qs"))

# Within-lineage frequencies per patient ------------------------------------
# One batch per patient; the many-to-one join below fails if a patient spans batches.
# factor() so a numeric batch ID gets one term per batch, not a slope.
patient_batch <- cells |>
  distinct(patient_id, batch) |>
  mutate(batch = factor(batch))

patient_freq <- cells |>
  count(patient_id, subclusterPG_merged, name = "n_cells") |>
  mutate(lineage = sub("_.*", "", subclusterPG_merged)) |>
  mutate(n_lineage = sum(n_cells), .by = c(patient_id, lineage)) |>
  mutate(
    freq = n_cells / n_lineage * 100,
    asin_freq = asin(sqrt(freq / 100))
  ) |>
  left_join(patient_batch, by = "patient_id", relationship = "many-to-one") |>
  left_join(
    select(metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )

n_before <- n_distinct(patient_freq$patient_id)
patient_freq <- patient_freq |>
  drop_na() |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)
message(sprintf(
  "Patients: %d before, %d after dropping missing covariates and %s",
  n_before,
  n_distinct(patient_freq$patient_id),
  paste(EXCLUDED_PATIENTS, collapse = ", ")
))

# Subcluster x patient matrix; a patient with no cells in a subcluster has frequency 0.
freq_mat <- patient_freq |>
  select(subclusterPG_merged, patient_id, asin_freq) |>
  pivot_wider(names_from = patient_id, values_from = asin_freq, values_fill = 0) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# limma trend test ---------------------------------------------------------
patients <- patient_freq |>
  distinct(patient_id, batch, group2, age, sex) |>
  column_to_rownames("patient_id")
patients <- patients[colnames(freq_mat), ]

# Weights come from group sizes among the analysed patients, so set contrasts after filtering.
patients$group2 <- droplevels(patients$group2)
contrasts(patients$group2) <- contr.poly.weighted(patients$group2)

design <- model.matrix(DESIGN_FORMULA, data = patients)
fit <- lmFit(freq_mat, design, method = "robust", maxit = 100) |>
  eBayes(robust = TRUE)

trend_linear <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

trend_quadratic <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

# FDR is Benjamini-Hochberg across all subclusters, separately for each trend term.
trend_results <- trend_linear |>
  inner_join(trend_quadratic, by = "cluster", relationship = "one-to-one") |>
  mutate(lineage = sub("_.*", "", cluster), .before = cluster) |>
  arrange(pval_L)

saveRDS(trend_results, file.path(out_dir, "MPL_limma_trend_test.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
