# Per-group log2-CPM means/SDs and log fold changes vs ctrl.
# Refactor of scripts/00_legacy.R; outputs are identical.
# Usage: Rscript scripts/04_group-means.R <run_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/04_group-means.R <run_dir>")
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

groups <- c("ctrl", "treatA", "treatB")
min_lib_size <- 20000 # samples with total counts <= this are dropped
min_count <- 10       # a gene is kept if >= min_samples samples have >= min_count
min_samples <- 3
lfc_cutoff <- 1

write_out <- function(df, name) {
  write.table(df, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

counts <- read.delim(here::here("data/counts.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.delim(here::here("data/samples.tsv"), stringsAsFactors = FALSE)

x <- as.matrix(counts[, -1])
rownames(x) <- counts$gene
samples <- samples[match(colnames(x), samples$sample_id), ]

# Drop low-depth samples, then lowly expressed genes
keep <- colSums(x) > min_lib_size
x <- x[, keep]
samples <- samples[keep, ]
x <- x[rowSums(x >= min_count) >= min_samples, ]

log_cpm <- log2(t(t(x) / colSums(x)) * 1e6 + 1)

res <- data.frame(gene = rownames(log_cpm), stringsAsFactors = FALSE)
for (g in groups) {
  m <- log_cpm[, samples$group == g, drop = FALSE]
  res[[paste0(g, "_mean")]] <- rowMeans(m)
  res[[paste0(g, "_sd")]] <- apply(m, 1, sd)
}

direction <- function(lfc) ifelse(lfc > lfc_cutoff, "up", ifelse(lfc < -lfc_cutoff, "down", "ns"))
res$lfcA <- res$treatA_mean - res$ctrl_mean
res$lfcB <- res$treatB_mean - res$ctrl_mean
res$dirA <- direction(res$lfcA)
res$dirB <- direction(res$lfcB)
res <- res[order(-abs(res$lfcA)), ]
write_out(res, "group_means.tsv")

top <- res[res$dirA != "ns" | res$dirB != "ns", c("gene", "lfcA", "lfcB", "dirA", "dirB")]
write_out(top, "top_genes.tsv")

# QC on the filtered gene set (as in the legacy script)
write_out(data.frame(sample_id = colnames(x), group = samples$group,
                     lib = colSums(x), ngenes = colSums(x > 0)), "sample_qc.tsv")
write_out(data.frame(group = groups, n = sapply(groups, function(g) sum(samples$group == g))),
          "group_sizes.tsv")
