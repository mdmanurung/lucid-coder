# Naive answer: a sound stim-vs-ctrl fit per cell type, then "differs" = significant in one
# cell type but not the other. A difference in significance is not a significant difference.
library(edgeR)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim("data/counts.tsv", check.names = FALSE)
samples <- read.delim("data/samples.tsv")
mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$gene
samples <- samples[match(colnames(mat), samples$sample_id), ]

de <- lapply(c(mono = "CD14+ Monocytes", cd4t = "CD4 T cells"), function(type) {
  s <- samples[samples$cell_type == type, ]
  s$condition <- factor(s$condition, levels = c("ctrl", "stim"))
  design <- model.matrix(~ donor_id + condition, s)
  y <- DGEList(mat[, s$sample_id])
  y <- normLibSizes(y[filterByExpr(y, design), , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  topTags(glmQLFTest(glmQLFit(y, design), coef = "conditionstim"), n = Inf)$table
})

genes <- intersect(rownames(de$mono), rownames(de$cd4t))
sig_mono <- de$mono[genes, "FDR"] < 0.05
sig_cd4t <- de$cd4t[genes, "FDR"] < 0.05
res <- data.frame(
  gene = genes,
  log2_fc_difference = de$mono[genes, "logFC"] - de$cd4t[genes, "logFC"],
  differs = xor(sig_mono, sig_cd4t)
)
write.table(
  res,
  file.path(out, "response_difference.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
