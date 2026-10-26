# Per-group log2-CPM means/SDs, log fold changes vs ctrl, and sample QC.
# Refactor of scripts/00_legacy.R; outputs are identical to the legacy script.
#
# Usage: Rscript scripts/04_group-means.R <run_dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/04_group-means.R <run_dir>")
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters (unchanged from the legacy script)
min_lib_size  <- 20000  # samples need library size > this
min_count     <- 10     # a gene needs >= min_count reads ...
min_samples   <- 3      # ... in >= min_samples samples
lfc_threshold <- 1      # |log2 FC| > this is called up/down
groups        <- c("ctrl", "treatA", "treatB")
ref_group     <- "ctrl"
contrasts     <- c(A = "treatA", B = "treatB")

# Base write.table keeps the legacy number formatting (15 significant digits, NA as "NA")
write_out <- function(df, name) {
  write.table(df, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Read ----
counts_df <- read_tsv(here("data", "counts.tsv"), col_types = cols(gene = col_character(), .default = col_double()))
samples   <- read_tsv(here("data", "samples.tsv"), col_types = cols(.default = col_character()))

counts <- counts_df %>% select(-gene) %>% as.matrix()
rownames(counts) <- counts_df$gene
samples <- samples[match(colnames(counts), samples$sample_id), ]

# ---- Filter samples, then genes ----
keep_samples <- colSums(counts) > min_lib_size
counts  <- counts[, keep_samples, drop = FALSE]
samples <- samples[keep_samples, ]

counts <- counts[rowSums(counts >= min_count) >= min_samples, , drop = FALSE]

# ---- Normalise ----
log_cpm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)

# ---- Per-group mean and SD ----
group_stats <- map(groups, function(g) {
  m <- log_cpm[, samples$group == g, drop = FALSE]
  tibble(mean = rowMeans(m), sd = apply(m, 1, sd)) %>%
    set_names(paste0(g, c("_mean", "_sd")))
})

call_direction <- function(lfc) {
  # nested if_else (not case_when) so NA lfc stays NA, as in the legacy ifelse()
  if_else(lfc > lfc_threshold, "up", if_else(lfc < -lfc_threshold, "down", "ns"))
}

stat_cols <- paste0(rep(groups, each = 2), c("_mean", "_sd"))  # ctrl_mean, ctrl_sd, treatA_mean, ...

group_means <- bind_cols(tibble(gene = rownames(log_cpm)), group_stats)
for (k in names(contrasts)) {
  lfc <- group_means[[paste0(contrasts[[k]], "_mean")]] - group_means[[paste0(ref_group, "_mean")]]
  group_means[[paste0("lfc", k)]] <- lfc
  group_means[[paste0("dir", k)]] <- call_direction(lfc)
}
group_means <- group_means %>%
  select(gene, all_of(stat_cols), starts_with("lfc"), starts_with("dir")) %>%
  arrange(desc(abs(lfcA)))

top_genes <- group_means %>%
  filter(dirA != "ns" | dirB != "ns") %>%
  select(gene, lfcA, lfcB, dirA, dirB)

# ---- Sample QC (computed on the gene-filtered counts, as in the legacy script) ----
sample_qc <- tibble(
  sample_id = colnames(counts),
  group     = samples$group,
  lib       = colSums(counts),
  ngenes    = colSums(counts > 0)
)

group_sizes <- tibble(group = groups, n = map_int(groups, ~ sum(samples$group == .x)))

# ---- Write ----
write_out(group_means, "group_means.tsv")
write_out(top_genes,   "top_genes.tsv")
write_out(sample_qc,   "sample_qc.tsv")
write_out(group_sizes, "group_sizes.tsv")
