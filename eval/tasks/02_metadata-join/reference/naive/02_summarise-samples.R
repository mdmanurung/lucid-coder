# Naive answer: binds the metadata rows onto the count columns by position.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim("data/counts.tsv", check.names = FALSE)
meta <- read.delim("data/metadata.tsv")

res <- data.frame(
  sample_id = colnames(counts)[-1],
  condition = meta$condition,
  batch = meta$batch,
  library_size = colSums(counts[, -1])
)
write.table(res, file.path(out, "sample_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
