options(warn = 1)

# Purpose: Test for rural-to-urban linear and quadratic trends in subcluster frequencies of
#   MPL-stimulated cells, adjusted for age, sex and batch, using robust limma.
# Inputs:
#   data/MPL_data.qs         one row per cell [gap: producer not recorded]
#   data/MPL_metadata.rds    one row per patient; group2 ordered rural -> urban [gap]
#   src/weightedcontrasts.R  contr.poly.weighted()
# Outputs:
#   data/MPL_cluster_median_markers.qs  median intensity of every marker per subcluster
#   outputs/MPL_limma_trend_test.rds    linear-trend statistics plus pval_Q, fdr_Q per subcluster
#   outputs/session_info.txt
# Run: Rscript src/MPL_trend_test.R  (from the project root)
# Notes:
#   - Frequency = subcluster cells / cells of the same lineage in that patient x 100, then
#     asin(sqrt(proportion)).
#   - Marker medians use every cell, including patients later excluded from the trend test.
#   - FDR is Benjamini-Hochberg across all subclusters, separately for the linear and the
#     quadratic term.
#   - Column names are kept as delivered (no janitor::clean_names()) so marker and subcluster
#     names match the panel.

# Setup -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(here)

source(here("src", "weightedcontrasts.R"))

# LD413 is excluded from the trend test by study decision
EXCLUDED_PATIENTS <- c("LD413")
# Lineage is the subcluster name up to the first underscore, e.g. CD4T_Th1 -> CD4T
LINEAGE_PATTERN <- "_.*$"
# Iterations for MASS::rlm() via lmFit(method = "robust"); the default 20 may not converge
ROBUST_MAXIT <- 100

out_dir <- here("outputs")
median_markers_path <- here("data", "MPL_cluster_median_markers.qs")
trend_test_path <- file.path(out_dir, "MPL_limma_trend_test.rds")

existing_outputs <- c(median_markers_path, trend_test_path)
existing_outputs <- existing_outputs[file.exists(existing_outputs)]
if (length(existing_outputs) > 0) {
  stop(
    "Outputs already exist: ",
    paste(existing_outputs, collapse = ", "),
    "; delete them to rerun.",
    call. = FALSE
  )
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load inputs -----------------------------------------------------------------

# patient_id as character in both inputs so it matches the matrix column names exactly
mpl_cells <- qs::qread(here("data", "MPL_data.qs")) |>
  mutate(patient_id = as.character(patient_id))
cell_cols <- c("patient_id", "batch", "subclusterPG_merged", "CD45", "IL10")
missing_cell_cols <- setdiff(cell_cols, colnames(mpl_cells))
if (length(missing_cell_cols) > 0) {
  stop("MPL_data.qs: missing columns ", toString(missing_cell_cols), ".", call. = FALSE)
}

patient_metadata <- readRDS(here("data", "MPL_metadata.rds")) |>
  mutate(patient_id = as.character(patient_id))
metadata_cols <- c("patient_id", "group2", "age", "sex")
missing_metadata_cols <- setdiff(metadata_cols, colnames(patient_metadata))
if (length(missing_metadata_cols) > 0) {
  stop("MPL_metadata.rds: missing columns ", toString(missing_metadata_cols), ".", call. = FALSE)
}

# Median marker intensity per subcluster ---------------------------------------

cluster_median_markers <- mpl_cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median), .groups = "drop")

# Subcluster frequencies within lineage ----------------------------------------

# Batch enters the model per patient, so each patient must come from a single batch
batches_per_patient <- mpl_cells |> distinct(patient_id, batch) |> count(patient_id)
if (any(batches_per_patient$n > 1)) {
  stop(
    "MPL_data.qs: patients in more than one batch: ",
    toString(batches_per_patient$patient_id[batches_per_patient$n > 1]),
    "; resolve batch per patient.",
    call. = FALSE
  )
}

subcluster_freq <- mpl_cells |>
  count(patient_id, batch, subclusterPG_merged, name = "n_cells") |>
  mutate(lineage = sub(LINEAGE_PATTERN, "", subclusterPG_merged)) |>
  group_by(patient_id, lineage) |>
  mutate(n_lineage_cells = sum(n_cells)) |>
  ungroup() |>
  mutate(
    freq = n_cells / n_lineage_cells * 100,
    asin_freq = asin(sqrt(freq / 100))
  ) |>
  left_join(
    patient_metadata |> select(all_of(metadata_cols)),
    by = "patient_id",
    relationship = "many-to-one"
  )

complete_freq <- drop_na(subcluster_freq)
dropped_na <- setdiff(subcluster_freq$patient_id, complete_freq$patient_id)
message(
  "Removed ", nrow(subcluster_freq) - nrow(complete_freq), " rows with missing values; ",
  "patients lost entirely: ", length(dropped_na), " (", toString(dropped_na), ")."
)

analysis_freq <- complete_freq |> filter(!patient_id %in% EXCLUDED_PATIENTS)
message(
  "Removed ", nrow(complete_freq) - nrow(analysis_freq), " rows from excluded patients (",
  toString(intersect(EXCLUDED_PATIENTS, complete_freq$patient_id)), ")."
)

# Subcluster x patient matrix --------------------------------------------------

# A patient without cells in a subcluster has frequency 0, and asin(sqrt(0)) = 0
freq_mat <- analysis_freq |>
  select(subclusterPG_merged, patient_id, asin_freq) |>
  pivot_wider(names_from = patient_id, values_from = asin_freq, values_fill = 0) |>
  column_to_rownames("subclusterPG_merged") |>
  as.matrix()

patient_info <- analysis_freq |>
  distinct(patient_id, batch, group2, age, sex) |>
  mutate(batch = factor(batch))
patient_info <- patient_info[match(colnames(freq_mat), patient_info$patient_id), ]
if (!identical(patient_info$patient_id, colnames(freq_mat))) {
  stop("patient_info: rows do not match matrix columns; check patient_id.", call. = FALSE)
}

# Robust limma trend test -------------------------------------------------------

# Weighted polynomial contrasts make group2.L and group2.Q the linear and quadratic trends
# across rural -> urban, weighted by the number of patients per group in this analysis
contrasts(patient_info$group2) <- contr.poly.weighted(patient_info$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = patient_info)

trend_fit <- lmFit(freq_mat, design, method = "robust", maxit = ROBUST_MAXIT)
trend_fit <- eBayes(trend_fit, robust = TRUE)

linear_trend <- topTable(trend_fit, coef = "group2.L", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

quadratic_trend <- topTable(trend_fit, coef = "group2.Q", number = Inf, sort.by = "none") |>
  rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_test <- linear_trend |>
  left_join(quadratic_trend, by = "cluster", relationship = "one-to-one") |>
  mutate(lineage = sub(LINEAGE_PATTERN, "", cluster)) |>
  relocate(lineage, cluster) |>
  arrange(pval_L)

# Write outputs ----------------------------------------------------------------

qs::qsave(cluster_median_markers, median_markers_path)
saveRDS(trend_test, trend_test_path)

# Session info -----------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
