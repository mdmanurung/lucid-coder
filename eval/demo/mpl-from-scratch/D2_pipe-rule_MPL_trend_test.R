options(warn = 1)

# Purpose: Test for linear and quadratic rural-to-urban trends (group2) in subcluster
#          frequencies of MPL-stimulated cells, using robust limma on arcsine-square-root
#          transformed within-lineage frequencies, adjusted for age, sex and batch.
# Inputs:  data/MPL_data.qs        one row per cell (patient_id, batch, group3,
#                                  subclusterPG_merged, markers CD45..IL10) [gap: producer]
#          data/MPL_metadata.rds   one row per patient (patient_id, group2, age, sex) [gap: producer]
#          src/weightedcontrasts.R contr.poly.weighted()
# Outputs: data/MPL_cluster_median_markers.qs  median marker intensity per subcluster
#          outputs/MPL_limma_trend_test.rds    linear-trend statistics + quadratic p/FDR
# Run:     Rscript src/MPL_trend_test.R   (from the project root)
# Notes:   Frequency denominator is the patient's cells in the same lineage, not all cells.
#          Group contrasts are weighted by the number of patients per group2 level, so they
#          are computed on the patient-level factor after all exclusions.

# Setup ----------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(limma)

source(here::here("src", "weightedcontrasts.R"))

# Patient excluded from the trend test at the analyst's request.
EXCLUDED_PATIENTS <- c("LD413")
# limma::lmFit robust regression needs more iterations than the default 20 to converge.
ROBUST_MAXIT <- 100
# Coefficient names produced by contr.poly.weighted() on group2.
COEF_LINEAR <- "group2.L"
COEF_QUADRATIC <- "group2.Q"

median_markers_path <- here::here("data", "MPL_cluster_median_markers.qs")
trend_test_path <- here::here("outputs", "MPL_limma_trend_test.rds")

existing_outputs <- c(median_markers_path, trend_test_path)[
  file.exists(c(median_markers_path, trend_test_path))
]
if (length(existing_outputs) > 0) {
  stop(
    "Output already exists: ",
    paste(existing_outputs, collapse = ", "),
    "; move it or delete it first.",
    call. = FALSE
  )
}

# Load inputs ----------------------------------------------------------------------------

cells <- qs::qread(here::here("data", "MPL_data.qs"))
patient_metadata <- readRDS(here::here("data", "MPL_metadata.rds"))

required_cell_cols <- c("patient_id", "batch", "subclusterPG_merged", "CD45", "IL10")
missing_cell_cols <- setdiff(required_cell_cols, colnames(cells))
if (length(missing_cell_cols) > 0) {
  stop("MPL_data: missing columns ", toString(missing_cell_cols), "; check input.", call. = FALSE)
}
required_meta_cols <- c("patient_id", "group2", "age", "sex")
missing_meta_cols <- setdiff(required_meta_cols, colnames(patient_metadata))
if (length(missing_meta_cols) > 0) {
  stop(
    "MPL_metadata: missing columns ",
    toString(missing_meta_cols),
    "; check input.",
    call. = FALSE
  )
}

# Median marker intensity per subcluster -------------------------------------------------

# Markers are the contiguous block CD45..IL10 in the cell table.
cluster_median_markers <- cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median), .groups = "drop")

# Subcluster frequencies per patient -----------------------------------------------------

# Lineage is the part of the subcluster name before the first underscore (CD4T_Th1 -> CD4T).
cluster_counts <- cells |>
  count(patient_id, batch, subclusterPG_merged, name = "n_cells") |>
  mutate(lineage = sub("_.*$", "", subclusterPG_merged)) |>
  group_by(patient_id, lineage) |>
  mutate(n_lineage_cells = sum(n_cells)) |>
  ungroup() |>
  mutate(
    freq_pct = n_cells / n_lineage_cells * 100,
    freq_asin = asin(sqrt(n_cells / n_lineage_cells))
  )

# Batch is a patient-level covariate; a patient in two batches would duplicate design rows.
patient_batch <- distinct(cells, patient_id, batch)
if (anyDuplicated(patient_batch$patient_id) > 0) {
  stop("MPL_data: some patients span more than one batch; resolve before modelling.", call. = FALSE)
}

n_rows_before_join <- nrow(cluster_counts)
cluster_freq <- cluster_counts |>
  left_join(
    select(patient_metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )
if (nrow(cluster_freq) != n_rows_before_join) {
  stop("Metadata join changed the row count; check duplicate patient_id.", call. = FALSE)
}

patients_with_na <- cluster_freq |>
  filter(if_any(c(group2, age, sex, batch), is.na)) |>
  distinct(patient_id) |>
  pull(patient_id)
message(
  "Removed ", length(patients_with_na), " patient(s) with missing group2/age/sex/batch: ",
  toString(patients_with_na)
)
message("Removed patient(s) by exclusion list: ", toString(EXCLUDED_PATIENTS))

cluster_freq <- cluster_freq |>
  drop_na(group2, age, sex, batch) |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)

# Subcluster x patient matrix ------------------------------------------------------------

# A patient with no cells in a subcluster has proportion 0, and asin(sqrt(0)) = 0.
freq_matrix <- cluster_freq |>
  select(subclusterPG_merged, patient_id, freq_asin) |>
  pivot_wider(names_from = patient_id, values_from = freq_asin, values_fill = 0) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# limma robust trend test ----------------------------------------------------------------

patient_design_data <- cluster_freq |>
  distinct(patient_id, group2, age, sex, batch) |>
  mutate(
    group2 = droplevels(group2),
    sex = factor(sex),
    batch = factor(batch)
  )
design_row_order <- match(colnames(freq_matrix), patient_design_data$patient_id)
patient_design_data <- patient_design_data[design_row_order, ]
if (!identical(patient_design_data$patient_id, colnames(freq_matrix))) {
  stop("Design rows do not match matrix columns; check patient_id.", call. = FALSE)
}

# Weighted by patients per group2 level, so .L and .Q stay orthogonal for unbalanced groups.
contrasts(patient_design_data$group2) <- contr.poly.weighted(patient_design_data$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = patient_design_data)

fit <- lmFit(freq_matrix, design, method = "robust", maxit = ROBUST_MAXIT)
fit <- eBayes(fit, robust = TRUE)

# Benjamini-Hochberg FDR within each coefficient, across all subclusters.
linear_trend <- topTable(fit, coef = COEF_LINEAR, number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)
quadratic_trend <- topTable(fit, coef = COEF_QUADRATIC, number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_test <- linear_trend |>
  inner_join(quadratic_trend, by = "cluster", relationship = "one-to-one") |>
  mutate(lineage = sub("_.*$", "", cluster)) |>
  relocate(lineage, cluster) |>
  arrange(pval_L)

# Write outputs --------------------------------------------------------------------------

dir.create(here::here("outputs"), showWarnings = FALSE)
qs::qsave(cluster_median_markers, median_markers_path)
saveRDS(trend_test, trend_test_path)

# Session info ---------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  here::here("outputs", "MPL_trend_test_session_info.txt")
)
