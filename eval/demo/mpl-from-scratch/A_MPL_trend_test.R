# Trend test (rural -> urban) of subcluster frequencies in MPL-stimulated samples
#
# Inputs:
#   data/MPL_data.qs          one row per cell
#   data/MPL_metadata.rds     one row per patient
# Outputs:
#   data/MPL_cluster_median_markers.qs
#   outputs/MPL_limma_trend_test.rds

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(limma)
})

source("src/weightedcontrasts.R")

dir.create("outputs", showWarnings = FALSE)

cells <- qs::qread("data/MPL_data.qs")
meta <- readRDS("data/MPL_metadata.rds")

# 1. Median marker intensity per subcluster ----------------------------------
cluster_medians <- cells %>%
  group_by(subclusterPG_merged) %>%
  summarise(across(CD45:IL10, ~ median(.x, na.rm = TRUE)), .groups = "drop")

qs::qsave(cluster_medians, "data/MPL_cluster_median_markers.qs")

# 2. Per-patient subcluster frequencies within lineage -----------------------
# Subcluster names are "<lineage>_<annotation>"; the annotation itself may
# contain underscores, so split on the first one only.
cells <- cells %>%
  mutate(lineage = sub("_.*$", "", subclusterPG_merged))

patient_batch <- cells %>%
  distinct(patient_id, batch)
stopifnot(!anyDuplicated(patient_batch$patient_id))

freq <- cells %>%
  count(patient_id, lineage, subclusterPG_merged, name = "n_cluster") %>%
  group_by(patient_id, lineage) %>%
  mutate(n_lineage = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(
    freq = n_cluster / n_lineage * 100,
    asin_freq = asin(sqrt(freq / 100))
  ) %>%
  left_join(patient_batch, by = "patient_id") %>%
  left_join(meta %>% select(patient_id, group2, age, sex), by = "patient_id") %>%
  drop_na() %>%
  filter(patient_id != "LD413")

# 3. Subcluster x patient matrix (0 where a patient has no cells) -------------
mat <- freq %>%
  select(subclusterPG_merged, patient_id, asin_freq) %>%
  pivot_wider(names_from = patient_id, values_from = asin_freq,
              values_fill = 0) %>%
  column_to_rownames("subclusterPG_merged") %>%
  as.matrix()

# Patient-level design data, aligned to the matrix columns
design_df <- freq %>%
  distinct(patient_id, batch, group2, age, sex) %>%
  mutate(
    group2 = droplevels(group2),
    sex = factor(sex),
    batch = factor(batch)
  )
design_df <- design_df[match(colnames(mat), design_df$patient_id), ]
stopifnot(identical(design_df$patient_id, colnames(mat)))
stopifnot(nlevels(design_df$group2) == 3)

# 4. limma with weighted polynomial contrasts on group2 ----------------------
# Weights are computed from the patients actually in the model.
contrasts(design_df$group2) <- contr.poly.weighted(design_df$group2)

design <- model.matrix(~ group2 + age + sex + batch, data = design_df)

fit <- lmFit(mat, design, method = "robust", maxit = 100)
fit <- eBayes(fit, robust = TRUE)

# 5. Linear and quadratic trend results --------------------------------------
res_L <- topTable(fit, coef = "group2.L", number = Inf, sort.by = "none") %>%
  rownames_to_column("cluster") %>%
  rename(pval_L = P.Value, fdr_L = adj.P.Val)

res_Q <- topTable(fit, coef = "group2.Q", number = Inf, sort.by = "none") %>%
  rownames_to_column("cluster") %>%
  select(cluster, pval_Q = P.Value, fdr_Q = adj.P.Val)

res <- res_L %>%
  left_join(res_Q, by = "cluster") %>%
  mutate(lineage = sub("_.*$", "", cluster)) %>%
  select(lineage, cluster, everything()) %>%
  arrange(pval_L)

saveRDS(res, "outputs/MPL_limma_trend_test.rds")
