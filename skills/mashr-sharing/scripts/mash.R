# Purpose: mashr engine for mashr_sharing.py; never run by hand.
# Inputs:  bhat.tsv, shat.tsv (genes x cell types; NA = gene not measured there)
# Outputs: <outdir>/{pm,psd,lfsr}.tsv (NA where the gene was not measured),
#          <outdir>/sharing.tsv (cell type x cell type), <outdir>/null_corr.tsv (V);
#          stdout: "n_strong <n>", then R, mashr and ashr versions.
# Run:     Rscript mash.R bhat.tsv shat.tsv canonical+ed identity 0.05 0.05 0.5 0 outdir
# Exit 3: mashr missing. Exit 4: fewer strong signals than cell types for ED.

for (pkg in c("mashr", "ashr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("R package not installed: ", pkg)
    quit(status = 3)
  }
}
suppressPackageStartupMessages(library(mashr))

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "bhat",
  "shat",
  "covariances",
  "null_corr",
  "strong_lfsr",
  "lfsr_max",
  "factor",
  "seed",
  "outdir"
)
seed <- as.integer(args[["seed"]])
set.seed(seed)

read_matrix <- function(path) {
  as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
}
bhat <- read_matrix(args[["bhat"]])
shat <- read_matrix(args[["shat"]])
measured <- !is.na(bhat)

# mash_set_data turns NA into Bhat = 0, Shat = 1e6 (no information).
data <- mash_set_data(bhat, shat)
U <- cov_canonical(data)
# em: as dreamlet's run_mash, V by EM over the canonical covariances, before ED.
V <- diag(ncol(bhat))
if (args[["null_corr"]] == "em") {
  V <- mash_estimate_corr_em(data, U, details = TRUE)$V
  data <- mash_update_data(data, V = V)
}
dimnames(V) <- list(colnames(bhat), colnames(bhat))
n_strong <- NA
if (args[["covariances"]] == "canonical+ed") {
  strong <- get_significant_results(
    mash_1by1(data),
    thresh = as.numeric(args[["strong_lfsr"]])
  )
  n_strong <- length(strong)
  if (n_strong < ncol(bhat)) {
    message("n_strong ", n_strong)
    quit(status = 4)
  }
  # npc = all cell types: no PCs are dropped, so npc is no choice to make.
  U_pca <- cov_pca(data, ncol(bhat), subset = strong)
  U <- c(U, cov_ed(data, U_pca, subset = strong))
}
m <- mash(data, U, verbose = FALSE, seed = seed)

# As dreamlet's run_mash: no posterior where the gene was not measured.
for (field in c("PosteriorMean", "PosteriorSD", "lfsr")) {
  m$result[[field]][!measured] <- NA
}

# mashr's own sharing, per pair, over genes measured in both cell types.
k <- ncol(bhat)
sharing <- diag(k)
dimnames(sharing) <- list(colnames(bhat), colnames(bhat))
for (i in seq_len(k - 1)) {
  for (j in (i + 1):k) {
    both <- measured[, i] & measured[, j]
    m_ij <- m
    m_ij$result <- lapply(
      m$result[c("PosteriorMean", "PosteriorSD", "lfsr")],
      function(x) x[both, c(i, j), drop = FALSE]
    )
    sharing[i, j] <- sharing[j, i] <- get_pairwise_sharing(
      m_ij,
      factor = as.numeric(args[["factor"]]),
      lfsr_thresh = as.numeric(args[["lfsr_max"]])
    )[1, 2]
  }
}

write_matrix <- function(x, name) {
  out <- data.frame(gene = rownames(x), x, check.names = FALSE)
  path <- file.path(args[["outdir"]], name)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_matrix(get_pm(m), "pm.tsv")
write_matrix(get_psd(m), "psd.tsv")
write_matrix(get_lfsr(m), "lfsr.tsv")
write_matrix(sharing, "sharing.tsv")
write_matrix(V, "null_corr.tsv")

cat("n_strong ", n_strong, "\n", sep = "")
cat(R.version.string, "\n", sep = "")
for (pkg in c("mashr", "ashr")) {
  cat(pkg, " ", as.character(packageVersion(pkg)), "\n", sep = "")
}
