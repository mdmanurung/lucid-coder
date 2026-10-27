# Naive answer: treats every cell as a replicate (t-test on log-normalised cells).
args <- commandArgs(trailingOnly = TRUE)
inp <- args[which(args == "--input") + 1]
out <- args[which(args == "--output") + 1]
dir.create(file.path(out, "reproducibility"), recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(file.path(inp, "counts.tsv"), check.names = FALSE)
cells <- read.delim(file.path(inp, "cells.tsv"))
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
writeLines(c("# DE", "Welch t-test per gene over cells."), file.path(out, "report.md"))
writeLines(sprintf('{"n_significant": %d}', sum(res$q_value < 0.05)), file.path(out, "result.json"))
writeLines(paste("Rscript skills/pseudobulk-de/pseudobulk_de.R", paste(args, collapse = " ")),
  file.path(out, "reproducibility", "commands.sh"))
