# Naive answer: treats every cell as a replicate (t-test on log-normalised cells).
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim("data/counts.tsv", check.names = FALSE)
cells <- read.delim("data/cells.tsv")
mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$gene
cells <- cells[match(colnames(mat), cells$cell_id), ]

norm <- log1p(t(t(mat) / colSums(mat)) * 1e4)
disease <- cells$condition == "disease"

res <- do.call(rbind, lapply(rownames(norm), function(g) {
  tt <- t.test(norm[g, disease], norm[g, !disease])
  data.frame(
    gene = g,
    log2_fc = log2(mean(expm1(norm[g, disease])) + 1) - log2(mean(expm1(norm[g, !disease])) + 1),
    p_value = tt$p.value
  )
}))
res$q_value <- p.adjust(res$p_value, method = "BH")
write.table(res, file.path(out, "de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
