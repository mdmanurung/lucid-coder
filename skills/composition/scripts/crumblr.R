# Purpose: crumblr engine for composition.py; never run by hand.
# Inputs:  counts.tsv (samples x labels, first column sample), meta.tsv (one row per
#          sample, formula terms)
# Outputs: <out_da> with label, logFC, se, t, ci_low, ci_high, pvalue, padj;
#          <out_clr> with crumblr's CLR values (samples x labels);
#          R and package versions on stdout, appended to environment.txt by the caller.
# Run:     Rscript crumblr.R counts.tsv meta.tsv "~ condition + (1|donor)" \
#            conditionstim 1 0 da.tsv clr.tsv
# Notes:   Calls follow ~/.claude/skills/crumblr-compositional/references/
#          getting-started.md:15-28,66-76 (crumblr, dream, eBayes, topTable).

suppressPackageStartupMessages({
  library(crumblr)
  library(variancePartition)
})

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "counts",
  "meta",
  "formula",
  "coef",
  "tau",
  "seed",
  "out_da",
  "out_clr"
)
# crumblr, dream and eBayes draw no random numbers; the seed is set for the record.
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
stopifnot(
  "meta rows match counts" = identical(rownames(meta), rownames(counts))
)

# crumblr takes samples as rows, labels as columns; $E comes back labels x samples.
cobj <- crumblr(counts, tau = as.numeric(args[["tau"]]))
fit <- dream(
  cobj,
  as.formula(args[["formula"]]),
  meta,
  BPPARAM = BiocParallel::SerialParam(),
  quiet = TRUE
)
fit <- eBayes(fit)
coef <- args[["coef"]]
if (!coef %in% colnames(coef(fit))) {
  stop(
    "coefficient ",
    coef,
    " not in the fit; choose one of: ",
    toString(colnames(coef(fit)))
  )
}
# confint = TRUE adds CI.L / CI.R, the 95% CI on the moderated t (dream's df).
res <- topTable(
  fit,
  coef = coef,
  number = Inf,
  sort.by = "none",
  confint = TRUE
)
out <- data.frame(
  label = rownames(res),
  logFC = res$logFC,
  se = res$logFC / res$t,
  t = res$t,
  ci_low = res$CI.L,
  ci_high = res$CI.R,
  pvalue = res$P.Value,
  padj = res$adj.P.Val
)
write.table(out, args[["out_da"]], sep = "\t", quote = FALSE, row.names = FALSE)
clr <- data.frame(sample = colnames(cobj$E), t(cobj$E), check.names = FALSE)
write.table(
  clr,
  args[["out_clr"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat(R.version.string, "\n", sep = "")
for (pkg in c("crumblr", "variancePartition", "limma", "lme4")) {
  cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")
}
