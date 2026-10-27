# Purpose: variancePartition engine for variance_partition.py; never run by hand.
# Inputs:  counts.tsv (genes x samples, summed raw counts), meta.tsv (one row per
#          sample, formula columns)
# Outputs: <out_dir>/cancor.tsv always; <out_dir>/variance_fractions.tsv only when
#          every canonical correlation is at most max_cancor (the caller fails
#          otherwise). R and package versions on stdout.
# Run:     Rscript varpart.R counts.tsv meta.tsv "~ (1|donor) + (1|condition) + rin" \
#            0.9 0 out_dir

suppressPackageStartupMessages({
  library(variancePartition)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c("counts", "meta", "formula", "max_cancor", "seed", "out")
# lme4 fits draw no random numbers; the seed is set for the record.
set.seed(as.integer(args[["seed"]]))

counts <- as.matrix(read.delim(
  args[["counts"]],
  row.names = 1,
  check.names = FALSE
))
meta <- read.delim(args[["meta"]], row.names = 1, check.names = FALSE)
stopifnot(
  "meta rows match counts columns" = identical(rownames(meta), colnames(counts))
)
form <- as.formula(args[["formula"]])
# A (1|x) term is categorical even when its values look numeric (donor 1, 2, 3).
random_terms <- vapply(reformulas::findbars(form), \(bar) deparse(bar[[3]]), "")
meta[random_terms] <- lapply(meta[random_terms], factor)

# canCorPairs reads the variables, not the random-effect syntax.
cancor <- canCorPairs(form, meta)
write.table(
  data.frame(term = rownames(cancor), cancor, check.names = FALSE),
  file.path(args[["out"]], "cancor.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat(R.version.string, "\n", sep = "")
for (pkg in c("variancePartition", "limma", "edgeR", "lme4")) {
  cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")
}
if (max(cancor[upper.tri(cancor)]) > as.numeric(args[["max_cancor"]])) {
  quit(status = 0)
}

# Genes and samples are already filtered by the caller; TMM, then voom weights
# from the same formula (dreamlet processAssays does the same per cell type).
dge <- calcNormFactors(DGEList(counts))
vobj <- voomWithDreamWeights(dge, form, meta, BPPARAM = SerialParam())
vp <- fitExtractVarPartModel(vobj, form, meta, BPPARAM = SerialParam())
write.table(
  data.frame(gene = rownames(vp), as.data.frame(vp), check.names = FALSE),
  file.path(args[["out"]], "variance_fractions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
