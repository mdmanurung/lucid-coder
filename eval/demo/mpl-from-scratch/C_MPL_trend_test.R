# Purpose: Test for rural-to-urban trends in subcluster frequencies of MPL-stimulated cells,
#          using robust limma on arcsine-square-root frequencies with weighted polynomial
#          contrasts for group2 (linear and quadratic terms).
# Inputs:  data/MPL_data.qs         one row per cell [gap: producer script not recorded]
#          data/MPL_metadata.rds    one row per patient [gap: producer script not recorded]
#          src/weightedcontrasts.R  contr.poly.weighted()
# Outputs: data/MPL_cluster_median_markers.qs  median marker intensity per subcluster
#          outputs/MPL_limma_trend_test.rds     linear-trend statistics + quadratic p/FDR
#          outputs/MPL_trend_test_session_info.txt
# Run:     Rscript src/MPL_trend_test.R   (from the project root)
# Notes:   - Frequency denominator: all cells of the same lineage in the same patient.
#          - batch is taken from the cell table; each patient must belong to one batch.
#          - Column names are kept as-is (no janitor::clean_names()) because the marker range
#            CD45:IL10 and the names asked for downstream are the original ones.
#          - The median table is written to data/ as requested, although data/ is otherwise
#            treated as read-only.

# Setup ------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(here)

source(here("src", "weightedcontrasts.R"))

# Excluded at the analyst's request; reason not recorded here.
EXCLUDED_PATIENTS <- c("LD413")
# limma::lmFit(method = "robust") passes maxit to MASS::rlm(); the default 20 can fail to converge.
RLM_MAXIT <- 100
# Coefficient names produced by contr.poly.weighted() on group2.
COEF_LINEAR <- "group2.L"
COEF_QUADRATIC <- "group2.Q"

cell_path <- here("data", "MPL_data.qs")
metadata_path <- here("data", "MPL_metadata.rds")
medians_path <- here("data", "MPL_cluster_median_markers.qs")
out_dir <- here("outputs")
results_path <- file.path(out_dir, "MPL_limma_trend_test.rds")

dir.create(out_dir, showWarnings = FALSE)

# Load inputs ------------------------------------------------------------------------------------

cells <- qs::qread(cell_path)
missing_cell_cols <- setdiff(
  c("patient_id", "batch", "subclusterPG_merged", "CD45", "IL10"),
  colnames(cells)
)
if (length(missing_cell_cols) > 0) {
  stop(
    "cells: missing columns ",
    paste(missing_cell_cols, collapse = ", "),
    "; check MPL_data.qs.",
    call. = FALSE
  )
}

metadata <- readRDS(metadata_path)
missing_meta_cols <- setdiff(c("patient_id", "group2", "age", "sex"), colnames(metadata))
if (length(missing_meta_cols) > 0) {
  stop(
    "metadata: missing columns ",
    paste(missing_meta_cols, collapse = ", "),
    "; check MPL_metadata.rds.",
    call. = FALSE
  )
}

# Character IDs so joins and the column-order check below compare like with like.
cells <- mutate(cells, patient_id = as.character(patient_id))
metadata <- mutate(metadata, patient_id = as.character(patient_id))

# Median marker intensity per subcluster ---------------------------------------------------------

cluster_medians <- cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median), .groups = "drop")

# Subcluster frequencies per patient -------------------------------------------------------------

# Lineage is the prefix before the first "_" in "<lineage>_<annotation>".
cluster_counts <- cells |>
  count(patient_id, subclusterPG_merged, name = "n_cells") |>
  mutate(lineage = sub("_.*$", "", subclusterPG_merged)) |>
  group_by(patient_id, lineage) |>
  mutate(n_lineage = sum(n_cells)) |>
  ungroup() |>
  mutate(
    freq = n_cells / n_lineage * 100,
    freq_asin = asin(sqrt(freq / 100))
  )

patient_batch <- cells |>
  distinct(patient_id, batch)
if (anyDuplicated(patient_batch$patient_id) > 0) {
  stop("patient_batch: a patient spans several batches; fix batch labels.", call. = FALSE)
}

# many-to-one: errors if metadata has duplicate patient_id rows.
cluster_freq <- cluster_counts |>
  left_join(patient_batch, by = "patient_id", relationship = "many-to-one") |>
  left_join(
    select(metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )

# Patients absent from the metadata get NA here and are removed with other incomplete rows.
incomplete_patients <- cluster_freq |>
  filter(if_any(c(batch, group2, age, sex), is.na)) |>
  distinct(patient_id) |>
  pull(patient_id)
message(
  "Removed ",
  length(incomplete_patients),
  " patient(s) with missing batch/group2/age/sex: ",
  paste(incomplete_patients, collapse = ", ")
)
message(
  "Removed ",
  sum(cluster_freq$patient_id %in% EXCLUDED_PATIENTS),
  " row(s) for excluded patient(s): ",
  paste(EXCLUDED_PATIENTS, collapse = ", ")
)

cluster_freq <- cluster_freq |>
  drop_na() |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)

# Subcluster x patient matrix --------------------------------------------------------------------

# A patient with no cells in a subcluster has frequency 0, and asin(sqrt(0)) = 0.
freq_matrix <- cluster_freq |>
  select(subclusterPG_merged, patient_id, freq_asin) |>
  pivot_wider(names_from = patient_id, values_from = freq_asin, values_fill = 0) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

# limma trend test -------------------------------------------------------------------------------

# One row per patient, in the same order as the matrix columns.
patient_info <- cluster_freq |>
  distinct(patient_id, batch, group2, age, sex) |>
  mutate(
    group2 = droplevels(group2),
    sex = factor(sex),
    batch = factor(batch)
  )
patient_info <- patient_info[match(colnames(freq_matrix), patient_info$patient_id), ]
if (!identical(patient_info$patient_id, colnames(freq_matrix))) {
  stop("patient_info: order differs from freq_matrix columns; check the pivot.", call. = FALSE)
}

# Weights come from the patients actually in the model, so set contrasts after all filtering.
contrasts(patient_info$group2) <- contr.poly.weighted(patient_info$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = patient_info)

fit <- lmFit(freq_matrix, design, method = "robust", maxit = RLM_MAXIT)
fit <- eBayes(fit, robust = TRUE)

# BH FDR is computed per term, across all subclusters.
linear_table <- topTable(fit, coef = COEF_LINEAR, number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

quadratic_table <- topTable(fit, coef = COEF_QUADRATIC, number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_results <- linear_table |>
  left_join(quadratic_table, by = "cluster", relationship = "one-to-one") |>
  left_join(
    distinct(cluster_counts, subclusterPG_merged, lineage),
    by = c("cluster" = "subclusterPG_merged"),
    relationship = "one-to-one"
  ) |>
  relocate(lineage, .before = cluster) |>
  arrange(pval_L)

# Write outputs ----------------------------------------------------------------------------------

if (nrow(cluster_medians) == 0 || nrow(trend_results) == 0) {
  stop("outputs: median or limma table is empty; check the filtering steps.", call. = FALSE)
}

qs::qsave(cluster_medians, medians_path)
saveRDS(trend_results, results_path)

# Session info -----------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "MPL_trend_test_session_info.txt")
)
