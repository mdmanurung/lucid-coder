# Purpose: fgsea engine for enrichment.py; never run by hand.
# Inputs:  ranks.tsv (contrast, gene, rank), sets.tsv (set, gene; long)
# Outputs: <out_tsv> with contrast, set, size, nes, pvalue, leading_edge;
#          R and fgsea versions on stdout, appended to environment.txt by the caller.
# Run:     Rscript fgsea.R ranks.tsv sets.tsv 10 500 10000 0 out.tsv
# Notes:   Exit 3 with the package name on stderr when fgsea is missing.

if (!requireNamespace("fgsea", quietly = TRUE)) {
  cat("fgsea\n", file = stderr())
  quit(status = 3)
}

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "ranks",
  "sets",
  "min_size",
  "max_size",
  "n_perm_simple",
  "seed",
  "out"
)

ranks <- read.delim(
  args[["ranks"]],
  colClasses = c("character", "character", "numeric")
)
sets <- read.delim(args[["sets"]], colClasses = "character")
pathways <- split(sets$gene, sets$set)

fits <- lapply(split(ranks, ranks$contrast), \(r) {
  # Reseed per contrast so each result is independent of the other contrasts.
  set.seed(as.integer(args[["seed"]]))
  fit <- fgsea::fgsea(
    pathways = pathways,
    stats = setNames(r$rank, r$gene),
    minSize = as.integer(args[["min_size"]]),
    maxSize = as.integer(args[["max_size"]]),
    nPermSimple = as.integer(args[["n_perm_simple"]]),
    BPPARAM = BiocParallel::SerialParam()
  )
  data.frame(
    contrast = r$contrast[[1]],
    set = fit$pathway,
    size = fit$size,
    nes = fit$NES,
    pvalue = fit$pval,
    leading_edge = vapply(fit$leadingEdge, paste, character(1), collapse = ",")
  )
})
write.table(
  do.call(rbind, fits),
  args[["out"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat(R.version.string, "\n", sep = "")
cat(
  "fgsea ",
  as.character(packageVersion("fgsea")),
  " (BPPARAM = SerialParam())\n",
  sep = ""
)
