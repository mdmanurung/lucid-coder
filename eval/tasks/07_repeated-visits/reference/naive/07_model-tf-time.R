# Naive answer: every visit is treated as an independent sample (lm, no donor term).
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/samples.tsv")
scores <- read.delim("data/tf_scores.tsv", check.names = FALSE)
df <- merge(samples, scores, by = "sample_id")

res <- do.call(
  rbind,
  lapply(setdiff(names(scores), "sample_id"), function(tf) {
    fit <- lm(df[[tf]] ~ age + sex + time + time_sqrd, data = df)
    cf <- coef(summary(fit))
    data.frame(tf = tf, estimate = cf["time", 1], p_value = cf["time", 4])
  })
)
res$q_value <- p.adjust(res$p_value, method = "BH")
write.table(res, file.path(out, "tf_time.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
