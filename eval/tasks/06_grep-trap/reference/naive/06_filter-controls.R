# Naive answer: x[-grep(), ] per control set; a set with zero matches empties the table.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim("data/counts.tsv", check.names = FALSE)
n_ercc <- length(grep("^ERCC-", counts$gene))
n_mt <- length(grep("^MT-", counts$gene))

counts <- counts[-grep("^ERCC-", counts$gene), ]
counts <- counts[-grep("^MT-", counts$gene), ]

write.table(counts, file.path(out, "counts_filtered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(gene_set = c("ERCC", "MT"), n_removed = c(n_ercc, n_mt)),
  file.path(out, "filter_log.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
