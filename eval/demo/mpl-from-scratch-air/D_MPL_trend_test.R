# Purpose: Test for linear and quadratic rural-to-urban trends in subcluster frequencies
#          (within lineage) of MPL-stimulated cells, using robust limma.
# Inputs:  data/MPL_data.qs          one row per cell [gap: producer not recorded]
#          data/MPL_metadata.rds     one row per patient [gap: producer not recorded]
#          src/weightedcontrasts.R   contr.poly.weighted()
# Outputs: data/MPL_cluster_median_markers.qs      median marker intensity per subcluster
#          outputs/MPL_limma_trend_test.rds         linear-trend statistics plus quadratic p/FDR
#          outputs/MPL_trend_test_session_info.txt
# Run:     Rscript src/MPL_trend_test.R   (from the project root)
# Notes:   Output paths are fixed by the request; the script refuses to overwrite them.
#          Column names are kept as delivered (no clean_names()) because marker and
#          subcluster names are carried into the outputs.

options(warn = 1)

# Setup -------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(limma)

source(here::here("src", "weightedcontrasts.R"))

# Patient flagged for exclusion from the trend test.
EXCLUDED_PATIENTS <- c("LD413")
# Lineage is the subcluster prefix before the first underscore, e.g. "CD4T" in "CD4T_Th1".
LINEAGE_PATTERN <- "_.*$"
# rlm() defaults to 20 iterations; raised as specified so robust fits converge.
LMFIT_MAXIT <- 100

out_median_file <- here::here("data", "MPL_cluster_median_markers.qs")
out_limma_file <- here::here("outputs", "MPL_limma_trend_test.rds")
out_session_file <- here::here("outputs", "MPL_trend_test_session_info.txt")

output_files <- c(out_median_file, out_limma_file)
existing_outputs <- output_files[file.exists(output_files)]
if (length(existing_outputs) > 0) {
  stop(
    sprintf("%s: already exists; delete it to rerun.", paste(existing_outputs, collapse = ", ")),
    call. = FALSE
  )
}

# Load inputs -------------------------------------------------------------------------------------

cells <- qs::qread(here::here("data", "MPL_data.qs"))
metadata <- readRDS(here::here("data", "MPL_metadata.rds"))

missing_cell_cols <- setdiff(
  c("patient_id", "batch", "group3", "subclusterPG_merged", "CD45", "IL10"),
  colnames(cells)
)
if (length(missing_cell_cols) > 0) {
  stop(
    sprintf("MPL_data.qs: missing columns %s.", paste(missing_cell_cols, collapse = ", ")),
    call. = FALSE
  )
}

missing_meta_cols <- setdiff(c("patient_id", "group2", "age", "sex"), colnames(metadata))
if (length(missing_meta_cols) > 0) {
  stop(
    sprintf("MPL_metadata.rds: missing columns %s.", paste(missing_meta_cols, collapse = ", ")),
    call. = FALSE
  )
}

# Median marker intensity per subcluster ----------------------------------------------------------

# Markers are the contiguous block CD45 ... IL10.
cluster_median_markers <- cells |>
  group_by(subclusterPG_merged) |>
  summarise(across(CD45:IL10, median), .groups = "drop")

# Subcluster frequencies per patient --------------------------------------------------------------

# One batch per patient; a patient spread over batches would duplicate rows in the design.
patient_batch <- cells |>
  distinct(patient_id, batch)
if (anyDuplicated(patient_batch$patient_id) > 0) {
  stop(
    "MPL_data.qs: some patients occur in more than one batch; resolve batch per patient.",
    call. = FALSE
  )
}

subcluster_freq <- cells |>
  mutate(lineage = sub(LINEAGE_PATTERN, "", subclusterPG_merged)) |>
  count(patient_id, lineage, subclusterPG_merged, name = "n_cells") |>
  group_by(patient_id, lineage) |>
  mutate(n_lineage_cells = sum(n_cells)) |>
  ungroup() |>
  # Denominator: the patient's cells in the same lineage.
  mutate(
    prop = n_cells / n_lineage_cells,
    freq = prop * 100,
    asin_sqrt_freq = asin(sqrt(prop))
  )

# relationship = "many-to-one" errors on duplicated patient_id in the metadata.
subcluster_freq <- subcluster_freq |>
  left_join(
    select(metadata, patient_id, group2, age, sex),
    by = "patient_id",
    relationship = "many-to-one"
  )

dropped_na_rows <- subcluster_freq |>
  filter(if_any(everything(), is.na))
message(sprintf(
  "Removed %d row(s) with missing values across %d patient(s): %s",
  nrow(dropped_na_rows),
  n_distinct(dropped_na_rows$patient_id),
  paste(unique(dropped_na_rows$patient_id), collapse = ", ")
))

excluded_present <- intersect(EXCLUDED_PATIENTS, subcluster_freq$patient_id)
message(sprintf(
  "Removed %d excluded patient(s): %s",
  length(excluded_present),
  paste(excluded_present, collapse = ", ")
))

subcluster_freq <- subcluster_freq |>
  drop_na() |>
  filter(!patient_id %in% EXCLUDED_PATIENTS)

# Subcluster x patient matrix ---------------------------------------------------------------------

# A patient with no cells in a subcluster has proportion 0, and asin(sqrt(0)) = 0.
freq_wide <- subcluster_freq |>
  select(patient_id, subclusterPG_merged, asin_sqrt_freq) |>
  pivot_wider(names_from = patient_id, values_from = asin_sqrt_freq, values_fill = 0)

freq_matrix <- freq_wide |>
  select(-subclusterPG_merged) |>
  as.matrix()
rownames(freq_matrix) <- freq_wide$subclusterPG_merged

# limma trend test --------------------------------------------------------------------------------

design_data <- subcluster_freq |>
  distinct(patient_id, group2, age, sex) |>
  left_join(patient_batch, by = "patient_id", relationship = "one-to-one") |>
  mutate(batch = factor(batch)) |>
  arrange(match(patient_id, colnames(freq_matrix)))

if (!identical(design_data$patient_id, colnames(freq_matrix))) {
  stop(
    "design_data: patient order differs from freq_matrix columns; check patient IDs.",
    call. = FALSE
  )
}

# Weights come from the analysed patients, so trends are orthogonal for this sample's group sizes.
contrasts(design_data$group2) <- contr.poly.weighted(design_data$group2)
design <- model.matrix(~ group2 + age + sex + batch, data = design_data)

fit <- lmFit(freq_matrix, design, method = "robust", maxit = LMFIT_MAXIT)
fit <- eBayes(fit, robust = TRUE)

# Results -----------------------------------------------------------------------------------------

# FDR family: all subclusters, separately for each trend term.
linear_trend <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("cluster") |>
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

quadratic_trend <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("cluster") |>
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

trend_results <- linear_trend |>
  left_join(quadratic_trend, by = "cluster", relationship = "one-to-one") |>
  mutate(lineage = sub(LINEAGE_PATTERN, "", cluster)) |>
  relocate(lineage, cluster)

# Write outputs -----------------------------------------------------------------------------------

if (nrow(cluster_median_markers) == 0 || nrow(trend_results) == 0) {
  stop("outputs: a result table is empty; check filtering upstream.", call. = FALSE)
}

qs::qsave(cluster_median_markers, out_median_file)
saveRDS(trend_results, out_limma_file)

# Session info ------------------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), out_session_file)
