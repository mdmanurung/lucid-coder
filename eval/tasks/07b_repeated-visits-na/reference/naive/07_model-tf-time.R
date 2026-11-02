# Naive answer: the right mixed model, but every sample with any missing TF score is dropped
# first (a wide drop_na), so the TFs with no gaps lose those samples too.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lmerTest)
})

out <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/samples.tsv")
samples$donor_id <- as.character(samples$donor_id)
scores <- read.delim("data/tf_scores.tsv", check.names = FALSE) |>
  drop_na()
df <- inner_join(samples, scores, by = "sample_id")

res <- do.call(
  rbind,
  lapply(setdiff(names(scores), "sample_id"), function(tf) {
    df$score <- df[[tf]]
    fit <- lmer(score ~ age + sex + time + time_sqrd + (1 | donor_id), df, REML = FALSE)
    cf <- coef(summary(fit))
    data.frame(tf = tf, estimate = cf["time", 1], p_value = cf["time", 5], n_samples = nobs(fit))
  })
)
res$q_value <- p.adjust(res$p_value, method = "BH")
write.table(res, file.path(out, "tf_time.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
