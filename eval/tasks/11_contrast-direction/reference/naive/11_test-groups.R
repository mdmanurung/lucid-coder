# Naive answer: a default factor sorts "disease" first, so coefficient 2 is healthy - disease.
library(limma)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/samples.tsv")
expr <- as.matrix(read.delim("data/expression.tsv", row.names = 1, check.names = FALSE))

design <- model.matrix(~ factor(samples$group))
fit <- eBayes(lmFit(expr, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
res <- data.frame(
  protein = rownames(tt),
  log2fc = tt$logFC,
  p_value = tt$P.Value,
  fdr = tt$adj.P.Val
)

write.table(res, file.path(out, "de_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
