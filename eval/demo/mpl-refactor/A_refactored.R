# MPL differential abundance analysis
#
# Tests for a linear and quadratic trend across the ordered levels of group2
# in the per-lineage frequency of each subcluster, using limma with weighted
# orthogonal polynomial contrasts, adjusted for age, sex and batch.
#
# Outputs:
#   data/MPL_cluster_median_markers.qs   median marker expression per subcluster
#   outputs/MPL_limma_trend_test.rds     limma linear (L) and quadratic (Q) trend results


# setup -------------------------------------------------------------------

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
    source("src/weightedcontrasts.R")  # contr.poly.weighted()

    conflicted::conflict_prefer("count", "dplyr")
    conflicted::conflict_prefer("filter", "dplyr")
    conflicted::conflict_prefer("select", "dplyr")
})

PATHS <- list(
    data = "data/MPL_data.qs",
    meta = "data/MPL_metadata.rds",
    medians = "data/MPL_cluster_median_markers.qs",
    results = "outputs/"
)

REFERENCE_SAMPLE <- "LD413"  # excluded from the analysis


# helpers -----------------------------------------------------------------

# limma topTable for one coefficient as a tibble, with p-value and FDR
# columns suffixed (e.g. pval_L, fdr_L).
tidy_toptable <- function(efit, coef, suffix) {
    topTable(efit, coef = coef, number = Inf) %>%
        as.data.frame() %>%
        rownames_to_column(var = "cluster") %>%
        as_tibble() %>%
        dplyr::rename(!!paste0("pval_", suffix) := P.Value,
                      !!paste0("fdr_", suffix) := adj.P.Val)
}


# read data ---------------------------------------------------------------

# subclusterPG_merged is "<lineage>_<annotation>"
dat <- qs::qread(PATHS$data, nthreads = 6) %>%
    separate(subclusterPG_merged, c("subset", "anno"), "_", remove = F, convert = TRUE) %>%
    mutate(lineage = subset) %>%
    as.data.table()


# cluster median marker expression ----------------------------------------

markersAll <- colnames(dat[, CD45:IL10])
cluster.medians <- dat[, lapply(.SD, median), .SDcols = markersAll, by = subclusterPG_merged]
qs::qsave(cluster.medians, PATHS$medians)


# count tables ------------------------------------------------------------

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

# Join cluster counts -> lineage counts -> sample metadata -> clinical metadata.
# The join order determines the row order of `counts`, which in turn fixes the
# column order of matCounts and row order of matMeta (see below).
clusterCounts <- left_join(metaLineage, countsCluster, by = "subclusterPG_merged")
clusterLineageCounts <- left_join(countsLineage, clusterCounts, by = c("patient_id", "lineage"))
sampleCounts <- left_join(metaSample, clusterLineageCounts, by = "patient_id")

counts <- left_join(meta, sampleCounts, by = "patient_id") %>%
    as.data.table() %>%
    mutate(
        freq_cluster_perLineage = n_cluster / n_lineage * 100,
        # arcsine square-root transform of the proportion
        freq_cluster_perLineage_trans = asin(sqrt(freq_cluster_perLineage / 100))
    ) %>%
    drop_na() %>%
    filter(patient_id != REFERENCE_SAMPLE)


# limma inputs ------------------------------------------------------------

# NOTE: matCounts columns (patients) and matMeta rows are aligned by position
# only: both follow the order in which patients first appear in `counts`.
# lmFit() does not match them by name, so do not reorder either one.

# subclusters x patients
matCounts <- counts %>%
    select(patient_id, subclusterPG_merged, freq_cluster_perLineage_trans) %>%
    pivot_wider(names_from = patient_id, values_from = freq_cluster_perLineage_trans, values_fill = 0) %>%
    column_to_rownames(var = "subclusterPG_merged")

# The model uses group2 from the clinical metadata, renamed to `group3`;
# the original group3 from the data is kept as `groupCat`.
matMeta <- counts %>%
    select(patient_id, group3 = group2, groupCat = group3, age, sex, batch) %>%
    distinct()


# limma differential analysis ---------------------------------------------

contrasts(matMeta$group3) <- contr.poly.weighted(matMeta$group3)
designTrend <- model.matrix(~ group3 + age + sex + batch, data = matMeta)

fit <- lmFit(matCounts, designTrend, method = "robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)

# coef 2 = group3.L (linear trend), coef 3 = group3.Q (quadratic trend)
res.L <- tidy_toptable(efit, coef = 2, suffix = "L")
res.Q <- tidy_toptable(efit, coef = 3, suffix = "Q")

res.LQ <- left_join(res.L, res.Q %>% select(cluster, pval_Q, fdr_Q), by = "cluster")
res.LQ <- left_join(metaLineage %>% dplyr::rename(cluster = subclusterPG_merged),
                    res.LQ,
                    by = "cluster")

saveRDS(res.LQ, paste0(PATHS$results, "MPL_limma_trend_test.rds"))
