# Purpose: Test for linear and quadratic rural-to-urban trends in subcluster frequencies
#          (as % of lineage) in MPL-stimulated samples, using robust limma.
# Inputs:  data/MPL_data.qs        one row per cell [gap: producer script]
#          data/MPL_metadata.rds   one row per patient [gap: producer script]
#          src/weightedcontrasts.R contr.poly.weighted()
# Outputs: data/MPL_cluster_median_markers.qs  median marker intensity per subcluster
#          outputs/MPL_limma_trend_test.rds    limma linear + quadratic trend statistics
# Run:     Rscript src/MPL_trend_test.R   (from the project root)
# Notes:   One limma "sample" is one patient. Contrasts are polynomial, weighted by the number
#          of patients per group2 level after filtering. BH FDR is per term, across all
#          subclusters.

# Setup ----------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(limma)

source(here::here("src", "weightedcontrasts.R"))

# Reason for exclusion not recorded here [gap]; listed as requested.
EXCLUDED_PATIENTS <- c("LD413")
# Subcluster names are "<lineage>_<annotation>"; lineage is the part before the first "_".
LINEAGE_PATTERN <- "_.*$"
ANNOTATION_PATTERN <- "^[^_]*_"
# Coefficient names produced by model.matrix() for the weighted polynomial contrasts.
COEF_LINEAR <- "group2.L"
COEF_QUADRATIC <- "group2.Q"

# Load inputs ----------------------------------------------------------------------------------

mpl_cells <- qs::qread(here::here("data", "MPL_data.qs"))
patient_metadata <- readRDS(here::here("data", "MPL_metadata.rds"))

required_cell_cols <- c("patient_id", "batch", "group3", "subclusterPG_merged", "CD45", "IL10")
missing_cell_cols <- setdiff(required_cell_cols, colnames(mpl_cells))
if (length(missing_cell_cols) > 0) {
  stop(
    sprintf("MPL_data.qs: missing columns %s; check the input.", toString(missing_cell_cols)),
    call. = FALSE
  )
}
required_meta_cols <- c("patient_id", "group2", "age", "sex")
missing_meta_cols <- setdiff(required_meta_cols, colnames(patient_metadata))
if (length(missing_meta_cols) > 0) {
  stop(
    sprintf("MPL_metadata.rds: missing columns %s; check the input.", toString(missing_meta_cols)),
    call. = FALSE
  )
}

# Median marker intensity per subcluster -------------------------------------------------------

cluster_median_markers <- mpl_cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median), .groups = "drop")

# Subcluster frequency per patient -------------------------------------------------------------

# Batch is a patient-level covariate; many-to-one below fails if a patient spans two batches.
patient_batch <- mpl_cells |>
  distinct(patient_id, batch)

subcluster_freq <- mpl_cells |>
  count(patient_id, subclusterPG_merged, name = "n_cells") |>
  mutate(lineage = sub(LINEAGE_PATTERN, "", subclusterPG_merged)) |>
  group_by(patient_id, lineage) |>
  mutate(n_lineage_cells = sum(n_cells)) |>
  ungroup() |>
  mutate(
    freq_pct = n_cells / n_lineage_cells * 100,
    freq_asin = asin(sqrt(freq_pct / 100))
  )

n_rows_before_join <- nrow(subcluster_freq)
subcluster_freq <- subcluster_freq |>
  left_join(patient_batch, by = "patient_id", relationship = "many-to-one") |>
  left_join(
    select(patient_metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )
if (nrow(subcluster_freq) != n_rows_before_join) {
  stop("subcluster_freq: row count changed in metadata join; check patient_id keys.", call. = FALSE)
}

incomplete_patients <- subcluster_freq |>
  filter(!complete.cases(pick(everything()))) |>
  distinct(patient_id) |>
  pull(patient_id)
message(sprintf(
  "Removed %d patient(s) with missing values (incl. no metadata): %s",
  length(incomplete_patients),
  toString(incomplete_patients)
))
excluded_present <- intersect(EXCLUDED_PATIENTS, subcluster_freq$patient_id)
message(sprintf(
  "Removed %d excluded patient(s): %s",
  length(excluded_present),
  toString(excluded_present)
))

subcluster_freq <- subcluster_freq |>
  drop_na() |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)

# Subcluster x patient matrix ------------------------------------------------------------------

# 0 where a patient has no cells in a subcluster (also when the whole lineage is absent).
freq_matrix <- subcluster_freq |>
  pivot_wider(
    id_cols = subclusterPG_merged,
    names_from = patient_id,
    values_from = freq_asin,
    values_fill = 0
  ) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# Design matrix --------------------------------------------------------------------------------

patient_design_data <- subcluster_freq |>
  distinct(patient_id, group2, age, sex, batch)
patient_design_data <- patient_design_data[
  match(colnames(freq_matrix), patient_design_data$patient_id),
]
if (!identical(patient_design_data$patient_id, colnames(freq_matrix))) {
  stop(
    "patient_design_data: patient order differs from freq_matrix columns; check IDs.",
    call. = FALSE
  )
}

# An empty group2 level gives zero weight, so the weighted contrasts would be NaN.
empty_levels <- levels(patient_design_data$group2)[table(patient_design_data$group2) == 0]
if (length(empty_levels) > 0) {
  stop(
    sprintf("group2: no patients left in level(s) %s; check filtering.", toString(empty_levels)),
    call. = FALSE
  )
}

patient_design_data <- patient_design_data |>
  mutate(
    sex = factor(sex),
    batch = factor(batch)
  )
contrasts(patient_design_data$group2) <- contr.poly.weighted(patient_design_data$group2)

design <- model.matrix(~ group2 + age + sex + batch, data = patient_design_data)

# Robust limma trend test ----------------------------------------------------------------------

trend_fit <- lmFit(freq_matrix, design, method = "robust", maxit = 100) |>
  eBayes(robust = TRUE)

linear_trend <- topTable(trend_fit, coef = COEF_LINEAR, number = Inf, sort.by = "none") |>
  rownames_to_column("subcluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

quadratic_trend <- topTable(trend_fit, coef = COEF_QUADRATIC, number = Inf, sort.by = "none") |>
  rownames_to_column("subcluster") |>
  select(subcluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_results <- linear_trend |>
  left_join(quadratic_trend, by = "subcluster", relationship = "one-to-one") |>
  mutate(
    lineage = sub(LINEAGE_PATTERN, "", subcluster),
    cluster = sub(ANNOTATION_PATTERN, "", subcluster),
    .before = 1
  ) |>
  select(-subcluster) |>
  arrange(pval_L)

# Write outputs --------------------------------------------------------------------------------

qs::qsave(cluster_median_markers, here::here("data", "MPL_cluster_median_markers.qs"))
saveRDS(trend_results, here::here("outputs", "MPL_limma_trend_test.rds"))

# Session info ---------------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), here::here("outputs", "MPL_trend_test_session_info.txt"))
