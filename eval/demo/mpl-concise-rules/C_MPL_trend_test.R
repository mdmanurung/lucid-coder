# Purpose: Test linear and quadratic rural-to-urban trends in MPL subcluster frequencies (limma).
# Inputs:  data/MPL_data.qs [gap], data/MPL_metadata.rds [gap], src/weightedcontrasts.R
# Outputs: data/MPL_cluster_median_markers.qs, outputs/MPL_limma_trend_test.rds,
#          outputs/session_info.txt
# Run:     Rscript src/MPL_trend_test.R

library(dplyr)
library(tidyr)
library(purrr)
library(limma)

source(here::here("src/weightedcontrasts.R"))

EXCLUDED_PATIENTS <- "LD413"
TREND_TERMS <- c(L = "group2.L", Q = "group2.Q")

cells <- qs::qread(here::here("data/MPL_data.qs"))
metadata <- readRDS(here::here("data/MPL_metadata.rds"))

# Marker medians per subcluster, over all cells ----------------------------
cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median)) |>
  qs::qsave(here::here("data/MPL_cluster_median_markers.qs"))

# Subcluster frequencies within lineage, per patient -----------------------
subcluster_freq_all <- cells |>
  count(patient_id, batch, subclusterPG_merged) |>
  mutate(lineage = sub("_.*", "", subclusterPG_merged)) |>
  group_by(patient_id, lineage) |>
  # Percent within lineage is n / sum(n) * 100; the transform takes the proportion
  mutate(asin_freq = asin(sqrt(n / sum(n)))) |>
  ungroup() |>
  left_join(
    select(metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )

subcluster_freq <- subcluster_freq_all |>
  drop_na() |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)
message(sprintf(
  "Removed %d of %d patient-subcluster rows (missing values or excluded patients)",
  nrow(subcluster_freq_all) - nrow(subcluster_freq),
  nrow(subcluster_freq_all)
))

# A subcluster absent in a patient is a proportion of 0, and asin(sqrt(0)) = 0
freq_mat <- subcluster_freq |>
  select(subclusterPG_merged, patient_id, asin_freq) |>
  pivot_wider(names_from = patient_id, values_from = asin_freq, values_fill = 0) |>
  tibble::column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# Design -------------------------------------------------------------------
patient_meta <- subcluster_freq |>
  distinct(patient_id, batch, group2, age, sex) |>
  mutate(across(c(batch, sex), factor))
stopifnot("one batch per patient" = !anyDuplicated(patient_meta$patient_id))
patient_meta <- patient_meta[match(colnames(freq_mat), patient_meta$patient_id), ]
stopifnot("design rows match matrix columns" = identical(patient_meta$patient_id, colnames(freq_mat)))

# Weights come from the patients actually modelled, so set after exclusions
contrasts(patient_meta$group2) <- contr.poly.weighted(patient_meta$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = patient_meta)

fit <- lmFit(freq_mat, design, method = "robust", maxit = 100) |>
  eBayes(robust = TRUE)

# Trend tables -------------------------------------------------------------
# BH FDR across all subclusters, separately for each trend term
trend_tables <- map(TREND_TERMS, \(term) {
  topTable(fit, coef = term, number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("subclusterPG_merged")
})

trend_test <- trend_tables[["L"]] |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val) |>
  left_join(
    trend_tables[["Q"]] |>
      select(subclusterPG_merged, pval_Q = P.Value, fdr_Q = adj.P.Val),
    by = "subclusterPG_merged",
    relationship = "one-to-one"
  ) |>
  mutate(
    lineage = sub("_.*", "", subclusterPG_merged),
    cluster = sub("^[^_]+_", "", subclusterPG_merged),
    .before = 1
  ) |>
  select(-subclusterPG_merged) |>
  arrange(pval_L)

dir.create(here::here("outputs"), showWarnings = FALSE)
saveRDS(trend_test, here::here("outputs/MPL_limma_trend_test.rds"))

writeLines(capture.output(sessionInfo()), here::here("outputs/session_info.txt"))
