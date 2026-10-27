#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ q4 procedure sound (yes/no; design, replication unit, method, filters, contrast): 
#$ q4 which (one sentence): 
#$ ---- your comments: #$ for code, #$P for procedure; own line above the code, or after it ----
options(warn = 1)
# Purpose: donor-level pseudobulk differential expression, disease vs healthy (edgeR quasi-likelihood).
# Inputs:  <input>/counts.tsv (raw UMI, genes x cells, first column `gene`) [user-supplied]
#          <input>/cells.tsv (cell_id, donor_id, condition) [user-supplied]
# Outputs: <output>/de.tsv, report.md, result.json, top_genes_donor_logcpm.{pdf,png,tsv},
#          reproducibility/{commands.sh,session_info.txt}
# Run:     Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(edgeR)
})

# Parameters ----
# First level is the reference, so log2_fc > 0 means higher in disease.
condition_levels <- c("healthy", "disease")
min_cells_per_donor <- 10
fdr_cutoff <- 0.05
n_top_genes <- 10
n_plot_genes <- 6
condition_colours <- c(healthy = "#0072B2", disease = "#D55E00")

# Paths are taken as given, relative to the caller's working directory.
cli_args <- commandArgs(trailingOnly = TRUE)
input_dir <- cli_args[match("--input", cli_args) + 1]
out_dir <- cli_args[match("--output", cli_args) + 1]
usage <- "Usage: Rscript pseudobulk_de.R --input <data_dir> --output <out_dir>"
if (is.na(input_dir) || is.na(out_dir)) {
  stop(usage, call. = FALSE)
}
if (file.exists(out_dir)) {
  stop("Output path already exists, refusing to overwrite: ", out_dir, call. = FALSE)
}

# Load data ----
counts <- read_tsv(file.path(input_dir, "counts.tsv"), show_col_types = FALSE)
cells <- read_tsv(file.path(input_dir, "cells.tsv"), show_col_types = FALSE)

counts_mat <- as.matrix(select(counts, -gene))
rownames(counts_mat) <- counts[["gene"]]
cells <- cells[match(colnames(counts_mat), cells[["cell_id"]]), ]
stopifnot(
  "cells.tsv rows match count columns" = identical(cells[["cell_id"]], colnames(counts_mat))
)

# Aggregate cells to donors ----
# The donor is the replication unit: cells from one donor are not independent samples.
donors <- cells |>
  count(donor_id, condition, name = "n_cells") |>
  mutate(condition = fct(condition, levels = condition_levels))
stopifnot("each donor has exactly one condition" = !anyDuplicated(donors[["donor_id"]]))

kept_donors <- filter(donors, n_cells >= min_cells_per_donor)
kept_cells <- cells[["donor_id"]] %in% kept_donors[["donor_id"]]

# Sum raw UMI counts over all of a donor's cells: one library per donor.
pseudobulk <- t(rowsum(t(counts_mat[, kept_cells]), cells[["donor_id"]][kept_cells]))
kept_donors <- kept_donors[match(colnames(pseudobulk), kept_donors[["donor_id"]]), ]

# Fit edgeR quasi-likelihood model ----
dge <- DGEList(pseudobulk, group = kept_donors[["condition"]])
expressed <- filterByExpr(dge, group = kept_donors[["condition"]])
dge <- normLibSizes(dge[expressed, , keep.lib.sizes = FALSE])

design <- model.matrix(~condition, data = kept_donors)
fit <- glmQLFit(dge, design)
contrast_coef <- paste0("condition", condition_levels[2])
qlf <- glmQLFTest(fit, coef = contrast_coef)

# Multiple-testing family: every gene that passes filterByExpr, this one contrast.
tested <- tibble(
  gene = rownames(qlf[["table"]]),
  log2_fc = qlf[["table"]][["logFC"]],
  p_value = qlf[["table"]][["PValue"]],
  q_value = p.adjust(p_value, method = "BH")
)
# Genes removed by the expression filter keep a row with NA statistics.
de <- tibble(gene = rownames(counts_mat)) |>
  left_join(tested, by = join_by(gene), relationship = "one-to-one") |>
  arrange(p_value)

n_significant <- sum(de[["q_value"]] < fdr_cutoff, na.rm = TRUE)

# Write outputs ----
dir.create(file.path(out_dir, "reproducibility"), recursive = TRUE)
write_tsv(de, file.path(out_dir, "de.tsv"), na = "NA")

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
writeLines(
  c(
    "#!/usr/bin/env bash",
    paste("# Run from:", getwd()),
    paste(c("Rscript", script_path, cli_args), collapse = " ")
  ),
  file.path(out_dir, "reproducibility", "commands.sh")
)

# Figure: donor-level expression of the top genes ----
plot_genes <- head(de[["gene"]], n_plot_genes)
log_cpm <- cpm(dge, log = TRUE)
plot_data <- log_cpm[plot_genes, , drop = FALSE] |>
  as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "donor_id", values_to = "log2_cpm") |>
  left_join(kept_donors, by = join_by(donor_id), relationship = "many-to-one") |>
  mutate(gene = fct(gene, levels = plot_genes))

theme_set(theme_bw(base_size = 8))
p <- ggplot(plot_data, aes(x = condition, y = log2_cpm, colour = condition)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "grey40", linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(~gene, nrow = 2) +
  scale_colour_manual(values = condition_colours, guide = "none") +
  labs(x = NULL, y = "TMM-normalised expression (log2 CPM)")

figure_stem <- file.path(out_dir, "top_genes_donor_logcpm")
ggsave(
  paste0(figure_stem, ".pdf"),
  plot = p,
  width = 120,
  height = 80,
  units = "mm",
  device = cairo_pdf
)
ggsave(paste0(figure_stem, ".png"), plot = p, width = 120, height = 80, units = "mm", dpi = 200)
write_tsv(p[["data"]], paste0(figure_stem, ".tsv"))

# Report and result.json ----
n_genes <- nrow(counts_mat)
n_expressed <- sum(expressed)
donors_per_condition <- count(kept_donors, condition, .drop = FALSE)[["n"]] |>
  set_names(levels(kept_donors[["condition"]]))
removed_donors <- anti_join(donors, kept_donors, by = join_by(donor_id))

top_table <- head(de, n_top_genes) |>
  mutate(
    row = sprintf(
      "| %s | %.3f | %s | %s |",
      gene,
      log2_fc,
      formatC(p_value, format = "e", digits = 2),
      formatC(q_value, format = "e", digits = 2)
    )
  )

report <- c(
  "# Pseudobulk differential expression: disease vs healthy",
  "",
  "## Data",
  "",
  sprintf("- Input: `%s` (counts.tsv, cells.tsv)", input_dir),
  sprintf(
    "- %d donors (%s), %d cells, %d genes.",
    nrow(donors),
    paste(
      sprintf(
        "%d %s",
        count(donors, condition, .drop = FALSE)[["n"]],
        levels(donors[["condition"]])
      ),
      collapse = ", "
    ),
    ncol(counts_mat),
    n_genes
  ),
  sprintf("- Cells per donor: %d to %d.", min(donors[["n_cells"]]), max(donors[["n_cells"]])),
  "",
  "## Replication unit and aggregation",
  "",
  "- The replication unit is the **donor**. Cells from one donor share genotype, environment and",
  "  sample handling, so testing cells as replicates would understate variance and inflate false",
  "  positives (pseudoreplication; Squair et al. 2021).",
  "- Raw UMI counts were **summed** over all cells of each donor, giving one pseudobulk library per",
  "  donor. No cell-type annotation was supplied, so all cells were pooled; a shift in cell-type",
  "  composition between conditions would therefore appear as differential expression.",
  "",
  "## Filters",
  "",
  sprintf(
    "1. Donors with fewer than %d cells: removed %d donor(s) (%d cells); %d donors kept.",
    min_cells_per_donor,
    nrow(removed_donors),
    sum(removed_donors[["n_cells"]]),
    nrow(kept_donors)
  ),
  sprintf(
    "2. Genes failing `edgeR::filterByExpr` (default thresholds, grouped by condition): removed %d of %d; %d genes tested.",
    n_genes - n_expressed,
    n_genes,
    n_expressed
  ),
  "",
  "## Model",
  "",
  "- Library sizes normalised with TMM (`edgeR::normLibSizes`). TMM assumes most genes are not",
  sprintf(
    "  differentially expressed; with only %d genes this assumption rests on few genes. Normalisation factors: %s.",
    n_expressed,
    paste(sprintf("%s %.3f", colnames(dge), dge[["samples"]][["norm.factors"]]), collapse = ", ")
  ),
  "- Negative binomial GLM with quasi-likelihood F-tests (`edgeR::glmQLFit`, `glmQLFTest`), design",
  sprintf(
    "  `~ condition`, %d disease vs %d healthy donors.",
    donors_per_condition[[condition_levels[2]]],
    donors_per_condition[[condition_levels[1]]]
  ),
  "  Why: pseudobulk counts are overdispersed counts with few replicates; the QL framework",
  "  shares information across genes for the dispersion and propagates its uncertainty into the",
  "  test, which keeps error rates controlled at small sample sizes (Lund et al. 2012; Chen et al. 2016).",
  sprintf(
    "- Empirical Bayes prior degrees of freedom for the QL dispersions: %.1f. A large value means the",
    median(fit[["df.prior"]])
  ),
  "  gene-wise dispersions were shrunk almost entirely onto the fitted trend; with few genes the",
  "  trend itself is estimated from little information.",
  "",
  "## Contrast",
  "",
  sprintf(
    "- Coefficient `%s`: %s minus %s on the log2 scale. `log2_fc > 0` means higher in %s.",
    contrast_coef,
    condition_levels[2],
    condition_levels[1],
    condition_levels[2]
  ),
  "- edgeR adds a small prior count (0.125) when computing log2 fold changes, so values for",
  "  near-zero genes are shrunk towards 0.",
  "",
  "## Multiple testing",
  "",
  sprintf(
    "- Family: the %d genes tested (all genes passing filter 2) for this single contrast. Benjamini-Hochberg",
    n_expressed
  ),
  "  adjustment gives `q_value`; genes removed by the filter have `NA` statistics in `de.tsv`.",
  sprintf("- Genes with q < %.2f: **%d** of %d tested.", fdr_cutoff, n_significant, n_expressed),
  "",
  sprintf("## Top %d genes by p-value", n_top_genes),
  "",
  "| gene | log2_fc | p_value | q_value |",
  "|---|---:|---:|---:|",
  top_table[["row"]],
  "",
  "## Figure",
  "",
  "![Donor-level expression of the top genes](top_genes_donor_logcpm.png)",
  "",
  sprintf(
    "**Figure 1. Donor-level expression of the %d genes with the smallest p-values.** Each point is one",
    n_plot_genes
  ),
  sprintf(
    "donor's pseudobulk library (n = %d healthy, %d disease donors); the grey bar marks the mean of",
    donors_per_condition[[condition_levels[1]]],
    donors_per_condition[[condition_levels[2]]]
  ),
  "donor values in each condition. y-axis: TMM-normalised log2 counts per million (edgeR `cpm`,",
  "prior count 2). Panels share the y-axis. Source data: `top_genes_donor_logcpm.tsv`; vector",
  "version: `top_genes_donor_logcpm.pdf`.",
  "",
  "## Reproducibility",
  "",
  "- Command: `reproducibility/commands.sh`; R session: `reproducibility/session_info.txt`."
)
writeLines(report, file.path(out_dir, "report.md"))

result <- list(
  parameters = list(
    input = input_dir,
    output = out_dir,
    replication_unit = "donor",
    aggregation = "sum of raw UMI counts per donor",
    reference_level = condition_levels[1],
    test_level = condition_levels[2],
    contrast = contrast_coef,
    min_cells_per_donor = min_cells_per_donor,
    gene_filter = "edgeR::filterByExpr (defaults, group = condition)",
    normalisation = "TMM",
    method = "edgeR glmQLFit + glmQLFTest",
    multiple_testing = "Benjamini-Hochberg over all tested genes",
    fdr_cutoff = fdr_cutoff
  ),
  counts = list(
    donors = nrow(donors),
    donors_kept = nrow(kept_donors),
    donors_per_condition = as.list(donors_per_condition),
    cells = ncol(counts_mat),
    cells_kept = sum(kept_cells),
    genes = n_genes,
    genes_tested = n_expressed
  ),
  n_significant = n_significant
)
jsonlite::write_json(
  result,
  file.path(out_dir, "result.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = NA
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "reproducibility", "session_info.txt")
)
