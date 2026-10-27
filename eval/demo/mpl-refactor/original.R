# functions --------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(qs)
    library(scattermore)
    library(ggbeeswarm)
    library(cowplot)
    library(ggpubr)
    library(data.table)
    library(ComplexHeatmap)
    library(patchwork)
    library(limma)
    library(ggforce)

    source("src/utils.R")
    source("src/utils_plot.R")
    source("src/weightedcontrasts.R")

    conflicted::conflict_prefer("count", "dplyr")
    conflicted::conflict_prefer("filter", "dplyr")
    conflicted::conflict_prefer("select", "dplyr")

})

PATHS <- list(
    data = "data/MPL_data.qs",
    meta = "data/MPL_metadata.rds",
    results = "outputs/"
)


# prepare data ------------------------------------------------------------

dat <- qs::qread(PATHS$data, nthreads=6)
dat <- dat %>%
    separate(subclusterPG_merged, c("subset", "anno"), "_", remove = F, convert = TRUE) %>%
    mutate(lineage = subset) %>%
    as.data.table()

markersAll <- colnames(dat[,CD45:IL10])
cluster.medians <- dat[, lapply(.SD, median), .SDcols = markersAll, by = subclusterPG_merged]
qs::qsave(cluster.medians, "data/MPL_cluster_median_markers.qs")


countsLineage <- dplyr::count(dat, patient_id, lineage, name = "n_lineage")
countsCluster <- dplyr::count(dat, patient_id, subclusterPG_merged, name = "n_cluster")
metaLineage <- dat %>%
    select(lineage, subclusterPG_merged) %>%
    distinct()
metaSample <- dat %>%
    select(patient_id, batch, group3) %>%
    distinct()

meta <- readRDS(PATHS$meta) %>%
    select(patient_id, group2, age, sex)
counts <- left_join(metaLineage, countsCluster, by = "subclusterPG_merged") %>%
    left_join(countsLineage, ., by = c("patient_id", "lineage")) %>%
    left_join(metaSample, ., by = "patient_id") %>%
    left_join(meta, ., by = "patient_id") %>%
    as.data.table() %>%
    mutate(freq_cluster_perLineage = n_cluster / n_lineage * 100) %>%
    mutate(freq_cluster_perLineage_trans = asin(sqrt(freq_cluster_perLineage / 100))) %>%
    drop_na() %>%
    filter(patient_id != "LD413")


matCounts <- counts %>%
    select(patient_id, subclusterPG_merged, freq_cluster_perLineage_trans) %>%
    pivot_wider(names_from = patient_id, values_from = freq_cluster_perLineage_trans, values_fill = 0) %>%
    column_to_rownames(var = "subclusterPG_merged")

matMeta <- counts %>%
    select(patient_id, group3 = group2, groupCat = group3, age, sex, batch) %>%
    distinct()


# limma differential analysis ---------------------------------------------

contrasts(matMeta$group3) <- contr.poly.weighted(matMeta$group3)
designTrend <- model.matrix(~ group3 + age + sex + batch, data = matMeta)

fit <- lmFit(matCounts, designTrend, method = "robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.L <- topTable(efit, coef = 2, number = Inf) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cluster") %>%
    as_tibble() %>%
    dplyr::rename(pval_L = P.Value, fdr_L = adj.P.Val)
res.Q <- topTable(efit, coef = 3, number = Inf) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cluster") %>%
    as_tibble() %>%
    dplyr::rename(pval_Q = P.Value, fdr_Q = adj.P.Val)
res.LQ <- left_join(res.L,
                    res.Q %>% select(cluster, pval_Q, fdr_Q),
                    by="cluster") %>%
    left_join(metaLineage %>% dplyr::rename(cluster = subclusterPG_merged), .,
              by = "cluster")

saveRDS(res.LQ, "outputs/MPL_limma_trend_test.rds")
