# Purpose: MAST hurdle-model engine for sc_de.py --mode mast; never run by hand.
# Inputs:  expr.tsv (genes x cells, log2 expression), meta.tsv (one row per cell)
# Outputs: <out_tsv> with feature, log2_fc, se, stat, pvalue (hurdle LRT, unadjusted);
#          R and MAST versions on stdout, appended to environment.txt by the caller.
# Run:     Rscript mast.R expr.tsv meta.tsv "~ condition + cngeneson" \
#            condition stim ctrl 0 out.tsv
#          A numeric term passes "" for level and reference.

suppressPackageStartupMessages(library(MAST))

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "expr",
  "meta",
  "formula",
  "column",
  "level",
  "reference",
  "seed",
  "out"
)
set.seed(as.integer(args[["seed"]]))
# zlm fans out over mc.cores; one core keeps runs identical.
options(mc.cores = 1)

expr <- as.matrix(read.delim(
  args[["expr"]],
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
  "meta rows match expr columns" = identical(rownames(meta), colnames(expr))
)

coef_name <- args[["column"]]
if (nzchar(args[["reference"]])) {
  meta[[coef_name]] <- relevel(
    factor(meta[[coef_name]]),
    ref = args[["reference"]]
  )
  coef_name <- paste0(coef_name, args[["level"]])
}

sca <- FromMatrix(
  expr,
  cData = meta,
  fData = data.frame(primerid = rownames(expr), row.names = rownames(expr)),
  check_sanity = FALSE
)
form <- as.formula(args[["formula"]])
# A (1 | donor) term needs lme4's glmer; empirical Bayes is not defined for it.
fit <- if (grepl("|", args[["formula"]], fixed = TRUE)) {
  zlm(form, sca, method = "glmer", ebayes = FALSE, silent = TRUE)
} else {
  zlm(form, sca, silent = TRUE)
}
stopifnot(
  "term is an estimable coefficient" = coef_name %in% colnames(coef(fit, "D"))
)

lrt <- lrTest(fit, CoefficientHypothesis(coef_name))
lfc <- as.data.frame(getLogFC(fit))
lfc <- lfc[lfc$contrast == coef_name, ]
lfc <- lfc[match(dimnames(lrt)[[1]], lfc$primerid), ]

out <- data.frame(
  feature = dimnames(lrt)[[1]],
  log2_fc = lfc$logFC,
  se = sqrt(lfc$varLogFC),
  stat = lrt[, "hurdle", "lambda"],
  pvalue = lrt[, "hurdle", "Pr(>Chisq)"]
)
write.table(out, args[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(R.version.string, "\n", sep = "")
cat("MAST ", as.character(packageVersion("MAST")), "\n", sep = "")
cat("lme4 ", as.character(packageVersion("lme4")), "\n", sep = "")
