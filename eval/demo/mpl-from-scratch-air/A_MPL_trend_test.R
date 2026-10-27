# Trend test of MPL-stimulated subcluster frequencies across group2
# (rural -> urban) using limma with weighted polynomial contrasts.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(qs)
  library(limma)
})

source("src/weightedcontrasts.R")

dir.create("outputs", showWarnings = FALSE)

# Data -------------------------------------------------------------------------

cells <- qread("data/MPL_data.qs")
meta <- readRDS("data/MPL_metadata.rds")

# Lineage is everything before the first underscore, e.g. "CD4T_Th1" -> "CD4T".
cells <- cells %>%
  mutate(
    subcluster = as.character(subclusterPG_merged),
    lineage = sub("_.*$", "", subcluster)
  )

# 1. Median marker intensity per subcluster ------------------------------------

cluster_medians <- cells %>%
  group_by(subclusterPG_merged) %>%
  summarise(across(CD45:IL10, \(x) median(x, na.rm = TRUE)), .groups = "drop")

qsave(cluster_medians, "data/MPL_cluster_median_markers.qs")

# 2. Per-patient subcluster frequencies within lineage -------------------------

# Batch is recorded per cell; collapse to one value per patient.
patient_batch <- cells %>%
  distinct(patient_id, batch)

stopifnot(!anyDuplicated(patient_batch$patient_id))

freq <- cells %>%
  count(patient_id, lineage, subcluster, name = "n_cells") %>%
  group_by(patient_id, lineage) %>%
  mutate(
    n_lineage = sum(n_cells),
    freq = n_cells / n_lineage * 100,
    asin_freq = asin(sqrt(freq / 100))
  ) %>%
  ungroup() %>%
  left_join(patient_batch, by = "patient_id") %>%
  left_join(
    meta %>% select(patient_id, group2, age, sex),
    by = "patient_id"
  ) %>%
  drop_na() %>%
  filter(patient_id != "LD413")

# 3. Subcluster x patient matrix -----------------------------------------------

freq_mat <- freq %>%
  select(subcluster, patient_id, asin_freq) %>%
  pivot_wider(
    names_from = patient_id,
    values_from = asin_freq,
    values_fill = 0
  ) %>%
  column_to_rownames("subcluster") %>%
  as.matrix()

# 4. limma robust fit with weighted polynomial contrasts -----------------------

design_df <- freq %>%
  distinct(patient_id, group2, age, sex, batch) %>%
  mutate(
    group2 = droplevels(group2),
    sex = factor(sex),
    batch = factor(batch)
  )

stopifnot(!anyDuplicated(design_df$patient_id))

design_df <- design_df[match(colnames(freq_mat), design_df$patient_id), ]
stopifnot(identical(design_df$patient_id, colnames(freq_mat)))

contrasts(design_df$group2) <- contr.poly.weighted(design_df$group2)

design <- model.matrix(~ group2 + age + sex + batch, data = design_df)
stopifnot(all(c("group2.L", "group2.Q") %in% colnames(design)))

fit <- lmFit(freq_mat, design, method = "robust", maxit = 100)
fit <- eBayes(fit, robust = TRUE)

# 5. Linear and quadratic trend tables -----------------------------------------

res_L <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "none") %>%
  rownames_to_column("subcluster") %>%
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

res_Q <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") %>%
  rownames_to_column("subcluster") %>%
  select(subcluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

res <- res_L %>%
  left_join(res_Q, by = "subcluster") %>%
  mutate(lineage = sub("_.*$", "", subcluster)) %>%
  select(lineage, cluster = subcluster, everything()) %>%
  arrange(pval_L)

saveRDS(res, "outputs/MPL_limma_trend_test.rds")
