# Purpose: DESeq2 engine for pseudobulk_de.py; never run by hand.
# Inputs:  pseudobulk_counts.tsv (genes x samples), sample_meta.tsv (one row per sample)
# Outputs: <out_tsv> with gene, base_mean, log2_fc, lfc_se, stat, pvalue, padj;
#          R and DESeq2 versions on stdout, appended to environment.txt by the caller.
# Run:     Rscript deseq2.R counts.tsv meta.tsv "~ patient + condition" \
#            condition stim ctrl 0.05 0 out.tsv

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "counts",
  "meta",
  "design",
  "factor",
  "numerator",
  "denominator",
  "alpha",
  "seed",
  "out"
)
# DESeq2's default path draws no random numbers; the seed is set for parity with Python.
set.seed(as.integer(args[["seed"]]))

counts <- as.matrix(read.delim(
  args[["counts"]],
  row.names = 1,
  check.names = FALSE
))
meta <- read.delim(
  args[["meta"]],
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = TRUE
)
stopifnot(identical(colnames(counts), rownames(meta)))

dds <- DESeqDataSetFromMatrix(
  counts,
  meta,
  design = as.formula(args[["design"]])
)
dds <- DESeq(dds, quiet = TRUE)
# alpha = --padj-max, so independent filtering targets the user's cutoff, as in PyDESeq2.
res <- results(
  dds,
  contrast = c(args[["factor"]], args[["numerator"]], args[["denominator"]]),
  alpha = as.numeric(args[["alpha"]])
)

out <- data.frame(
  gene = rownames(res),
  base_mean = res$baseMean,
  log2_fc = res$log2FoldChange,
  lfc_se = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj
)
write.table(out, args[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
# Reference level R chose for each factor term (alphabetical; DESeq2 does not relevel).
factor_terms <- intersect(
  all.vars(as.formula(args[["design"]])),
  colnames(meta)
)
factor_terms <- factor_terms[sapply(meta[factor_terms], is.factor)]
write.table(
  data.frame(
    term = factor_terms,
    reference = sapply(meta[factor_terms], \(x) levels(x)[1])
  ),
  paste0(args[["out"]], ".levels"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat(R.version.string, "\n", sep = "")
cat("DESeq2 ", as.character(packageVersion("DESeq2")), "\n", sep = "")
