# Purpose: dreamlet engine for pseudobulk_de.py; never run by hand.
# Inputs:  pseudobulk_counts.tsv (genes x samples), sample_meta.tsv (one row per
#          sample, design terms plus n_cells)
# Outputs: <out_tsv> with gene, base_mean, log2_fc, lfc_se, stat, pvalue, padj;
#          R and package versions on stdout, appended to environment.txt by the caller.
# Run:     Rscript dreamlet.R counts.tsv meta.tsv "~ condition + (1|patient)" \
#            condition stim ctrl 0.05 0 out.tsv

suppressPackageStartupMessages({
  library(dreamlet)
  library(SingleCellExperiment)
})

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
# dreamlet's default path draws no random numbers; the seed is set for the record.
# alpha is unused: limma/dream has no independent filtering.
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
factor <- args[["factor"]]
meta[[factor]] <- relevel(factor(meta[[factor]]), ref = args[["denominator"]])

# processAssays reads three things that aggregateToPseudoBulk writes:
# int_colData(x)$n_cells, metadata(x)$agg_pars$by and metadata(x)$aggr_means.
# Our counts are already summed, so we rebuild a cell-level SCE whose cells sum
# back to them (each sample: its pseudobulk vector plus n_cells - 1 zero cells)
# and let the public aggregateToPseudoBulk set all three, instead of copying
# its private slot layout by hand.
n_cells <- meta$n_cells
cell_sample <- rep(rownames(meta), n_cells)
first_cell <- cumsum(c(1, head(n_cells, -1)))
cells <- Matrix::sparseMatrix(
  i = integer(0),
  j = integer(0),
  x = numeric(0),
  dims = c(nrow(counts), sum(n_cells)),
  dimnames = list(rownames(counts), paste0("cell", seq_len(sum(n_cells))))
)
cells[, first_cell] <- counts
cd <- meta[cell_sample, setdiff(colnames(meta), "n_cells"), drop = FALSE]
cd$sample_id <- cell_sample
cd$cluster_id <- "cells"
rownames(cd) <- colnames(cells)
sce <- SingleCellExperiment(list(counts = cells), colData = cd)
pb <- aggregateToPseudoBulk(
  sce,
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)
stopifnot(all(cellCounts(pb)[rownames(meta), "cells"] == n_cells))

# Our --min-cells and --min-frac already filtered samples and genes, so every
# dreamlet filter is switched off: min.cells = 0 keeps every sample,
# min.count = 0, min.prop = 0 and min.total.count = 0 make filterByExpr keep
# every gene, and min.samples = 2 never binds because the caller requires
# 2 samples per contrast level. edgeR's DGEList still drops all-zero genes.
formula <- as.formula(args[["design"]])
proc <- processAssays(
  pb,
  formula,
  min.cells = 0,
  min.count = 0,
  min.samples = 2,
  min.prop = 0,
  min.total.count = 0,
  quiet = TRUE
)
fit <- dreamlet(proc, formula, quiet = TRUE)
coef <- paste0(factor, args[["numerator"]])
if (!coef %in% coefNames(fit)) {
  stop("coefficient ", coef, " not in dreamlet fit: ", toString(coefNames(fit)))
}
res <- as.data.frame(topTable(fit, coef = coef, number = Inf))

# base_mean is AveExpr (mean log2 CPM), not DESeq2's mean normalised count.
out <- data.frame(
  gene = res$ID,
  base_mean = res$AveExpr,
  log2_fc = res$logFC,
  lfc_se = res$logFC / res$t,
  stat = res$t,
  pvalue = res$P.Value,
  padj = res$adj.P.Val
)
write.table(out, args[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(R.version.string, "\n", sep = "")
for (pkg in c("dreamlet", "variancePartition", "limma", "edgeR", "lme4")) {
  cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")
}
