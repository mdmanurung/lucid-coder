#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ q4 procedure sound (yes/no; design, replication unit, method, filters, contrast): 
#$ q4 which (one sentence): 
#$ ---- your comments: #$ for code, #$P for procedure; own line above the code, or after it ----
#!/usr/bin/env Rscript
# Purpose: Pseudobulk differential expression, disease vs healthy, with the donor as replication unit.
# Inputs:  <input>/counts.tsv (raw UMI, genes x cells), <input>/cells.tsv (cell_id, donor_id,
#          condition) [gap: upstream single-cell QC]
# Outputs: <output>/de.tsv, report.md, result.json, top_genes.png, top_genes_source.tsv,
#          reproducibility/commands.sh, reproducibility/session_info.txt
# Run:     Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(edgeR)
  library(jsonlite)
})

# Parameters ----
condition_levels <- c("healthy", "disease") # first level is the reference
condition_colours <- c(healthy = "#4477AA", disease = "#CC6677")
min_cells_per_donor <- 10
fdr_cutoff <- 0.05
n_top_table <- 10
n_top_plot <- 6

theme_set(theme_bw(base_size = 8))

# Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[match("--input", args) + 1]
out_dir <- args[match("--output", args) + 1]
if (anyNA(c(input_dir, out_dir)) || any(str_starts(c(input_dir, out_dir), "--"))) {
  stop("Usage: Rscript pseudobulk_de.R --input <dir> --output <new_dir>", call. = FALSE)
}
if (file.exists(out_dir)) {
  stop(sprintf("Output '%s' already exists; choose a new directory.", out_dir), call. = FALSE)
}

# Load data ----
cells <- read_tsv(file.path(input_dir, "cells.tsv"), col_types = cols(.default = col_character()))
counts_tbl <- read_tsv(
  file.path(input_dir, "counts.tsv"),
  col_types = cols(gene = col_character(), .default = col_double())
)
counts <- as.matrix(select(counts_tbl, -gene))
rownames(counts) <- counts_tbl$gene

# fct() errors on any condition outside condition_levels, so no label is silently dropped.
cells <- cells[match(colnames(counts), cells$cell_id), ] |>
  mutate(condition = fct(condition, levels = condition_levels))
stopifnot("cell metadata matches count columns" = identical(cells$cell_id, colnames(counts)))

# Aggregate cells to donors ----
donors <- cells |>
  count(donor_id, condition, name = "n_cells")
stopifnot("each donor has exactly one condition" = !anyDuplicated(donors$donor_id))

pseudobulk <- t(rowsum(t(counts), cells$donor_id))[, donors$donor_id]
donors <- mutate(donors, total_umi = colSums(pseudobulk))

# Filter donors and genes ----
donors_kept <- filter(donors, n_cells >= min_cells_per_donor)
dge <- DGEList(counts = pseudobulk[, donors_kept$donor_id], group = donors_kept$condition)
keep_genes <- filterByExpr(dge, group = donors_kept$condition)
dge <- calcNormFactors(dge[keep_genes, , keep.lib.sizes = FALSE], method = "TMM")

filter_summary <- tibble(
  step = c(
    sprintf("Donors with < %d cells", min_cells_per_donor),
    "Genes failing edgeR::filterByExpr"
  ),
  before = c(nrow(donors), nrow(pseudobulk)),
  removed = c(nrow(donors) - nrow(donors_kept), sum(!keep_genes))
) |>
  mutate(after = before - removed)
message(sprintf(
  "%s: removed %d of %d\n",
  filter_summary$step,
  filter_summary$removed,
  filter_summary$before
))

# Fit model ----
design <- model.matrix(~condition, data = donors_kept)
contrast_coef <- paste0("condition", condition_levels[[2]])
fit <- glmQLFit(dge, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = contrast_coef)

# BH family: every gene that passed the filter; filtered genes stay in de.tsv with NA.
de_tested <- tibble(
  gene = rownames(qlf$table),
  log2_fc = qlf$table$logFC,
  p_value = qlf$table$PValue
) |>
  mutate(q_value = p.adjust(p_value, method = "BH"))
de <- tibble(gene = rownames(counts)) |>
  left_join(de_tested, join_by(gene), relationship = "one-to-one") |>
  arrange(p_value)

n_significant <- sum(de$q_value < fdr_cutoff, na.rm = TRUE)
n_significant_up <- sum(de$q_value < fdr_cutoff & de$log2_fc > 0, na.rm = TRUE)

# Figure: per-donor expression of the top genes ----
top_genes <- head(de_tested$gene[order(de_tested$p_value)], n_top_plot)
top_gene_cpm <- cpm(dge, log = TRUE, prior.count = 2)[top_genes, , drop = FALSE] |>
  as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "donor_id", values_to = "log2_cpm") |>
  left_join(donors_kept, join_by(donor_id), relationship = "many-to-one") |>
  mutate(gene = fct(gene, levels = top_genes))

p_top_genes <- ggplot(top_gene_cpm, aes(x = condition, y = log2_cpm, colour = condition)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.5) +
  facet_wrap(~gene, nrow = 2) +
  scale_colour_manual(values = condition_colours) +
  labs(x = NULL, y = "Expression (log2 CPM, TMM)") +
  theme(legend.position = "none")

# Write outputs ----
dir.create(file.path(out_dir, "reproducibility"), recursive = TRUE)
write_tsv(de, file.path(out_dir, "de.tsv"), na = "NA")
ggsave(
  file.path(out_dir, "top_genes.png"),
  plot = p_top_genes,
  width = 120,
  height = 75,
  units = "mm",
  dpi = 300
)
write_tsv(p_top_genes$data, file.path(out_dir, "top_genes_source.tsv"))

n_donors_by_condition <- donors_kept |>
  count(condition)
result <- list(
  skill = "pseudobulk-de",
  parameters = list(
    input = input_dir,
    output = out_dir,
    condition_levels = condition_levels,
    reference_level = condition_levels[[1]],
    contrast = sprintf("%s vs %s", condition_levels[[2]], condition_levels[[1]]),
    aggregation = "sum of raw UMI counts per gene per donor",
    min_cells_per_donor = min_cells_per_donor,
    gene_filter = "edgeR::filterByExpr, defaults, group = condition",
    normalisation = "TMM",
    model = "edgeR quasi-likelihood negative binomial GLM, design ~ condition, robust = TRUE",
    test = sprintf("glmQLFTest, coefficient %s", contrast_coef),
    multiple_testing = "Benjamini-Hochberg over all genes passing the gene filter",
    fdr_cutoff = fdr_cutoff
  ),
  counts = list(
    n_donors_input = nrow(donors),
    n_donors_tested = nrow(donors_kept),
    n_donors_by_condition = as.list(set_names(
      n_donors_by_condition$n,
      n_donors_by_condition$condition
    )),
    n_cells_input = ncol(counts),
    n_cells_tested = sum(donors_kept$n_cells),
    n_genes_input = nrow(counts),
    n_genes_tested = nrow(de_tested)
  ),
  n_significant = n_significant,
  n_significant_up = n_significant_up,
  n_significant_down = n_significant - n_significant_up
)
write_json(result, file.path(out_dir, "result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# Report ----
md_table <- \(df) {
  c(
    paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    paste("|", do.call(paste, c(df, sep = " | ")), "|")
  )
}

top_table <- de |>
  filter(!is.na(p_value)) |>
  head(n_top_table) |>
  mutate(
    log2_fc = sprintf("%.2f", log2_fc),
    across(c(p_value, q_value), \(x) formatC(x, format = "e", digits = 2))
  )
donor_table <- donors |>
  mutate(
    tested = if_else(donor_id %in% donors_kept$donor_id, "yes", "no"),
    tmm_factor = if_else(
      tested == "yes",
      sprintf("%.3f", dge$samples[donor_id, "norm.factors"]),
      "NA"
    )
  )
n_group <- sprintf("%d %s", n_donors_by_condition$n, n_donors_by_condition$condition)

report <- c(
  "# Pseudobulk differential expression: disease vs healthy",
  "",
  sprintf(
    "Generated %s by `pseudobulk_de.R`; exact command in `reproducibility/commands.sh`.",
    Sys.Date()
  ),
  "",
  "## Data summary",
  "",
  sprintf(
    "Input `%s`: %d donors (%s), %d cells, %d genes. Counts are raw UMIs.",
    input_dir,
    nrow(donors),
    paste(
      sprintf("%d %s", count(donors, condition)$n, count(donors, condition)$condition),
      collapse = ", "
    ),
    ncol(counts),
    nrow(counts)
  ),
  "",
  md_table(donor_table),
  "",
  "## Replication unit",
  "",
  sprintf(
    "The donor. Cells from one donor share genotype, environment and sample handling, so they are not independent replicates of the condition; treating cells as replicates inflates false positives (pseudoreplication). Inference therefore rests on %s donors.",
    paste(n_group, collapse = " vs ")
  ),
  "",
  "## Aggregation",
  "",
  "Raw UMI counts were summed per gene across all cells of each donor, giving one pseudobulk count profile per donor (genes x donors). Summing raw counts keeps the data as counts, which the negative binomial model requires; no cell-level normalisation was applied before summing. All cells of a donor were pooled: no cell-type stratification was done.",
  "",
  "## Filters",
  "",
  md_table(filter_summary),
  "",
  sprintf(
    "- Donors: kept if they contribute at least %d cells, so that each pseudobulk profile averages over enough cells.",
    min_cells_per_donor
  ),
  "- Genes: `edgeR::filterByExpr` with `group = condition` and default thresholds (min.count = 10, min.total.count = 15); a gene is kept if it reaches the CPM equivalent of 10 counts in at least as many donors as the smallest condition group. Low-count genes carry little information and destabilise dispersion estimation.",
  if (any(!keep_genes)) {
    sprintf("- Genes removed: %s.", paste(rownames(pseudobulk)[!keep_genes], collapse = ", "))
  },
  "",
  "## Normalisation and model",
  "",
  "Library sizes were scaled by TMM normalisation factors (`edgeR::calcNormFactors`; factors in the donor table). Each gene was fitted with a negative binomial GLM, design `~ condition`, using edgeR's quasi-likelihood pipeline (`glmQLFit(robust = TRUE)`, `glmQLFTest`).",
  "",
  "Why: pseudobulk counts are overdispersed counts with biological variation between donors; the negative binomial GLM models both. With few replicates per group, the quasi-likelihood F-test, which shares information across genes through empirical Bayes moderation of the dispersion, controls the type I error rate better than likelihood-ratio or cell-level tests. Pseudobulk methods with edgeR or limma are among the best-calibrated approaches for multi-donor single-cell DE (Squair et al. 2021; Crowell et al. 2020).",
  "",
  "## Contrast",
  "",
  sprintf(
    "Coefficient `%s`: log2 fold change of %s relative to %s (reference level). Positive `log2_fc` means higher expression in %s.",
    contrast_coef,
    condition_levels[[2]],
    condition_levels[[1]],
    condition_levels[[2]]
  ),
  "",
  "## Multiple testing",
  "",
  sprintf(
    "Family: the %d genes that passed the gene filter, one test per gene for one contrast. P-values were adjusted with Benjamini-Hochberg (`q_value`). Genes removed by the filter appear in `de.tsv` with `NA`.",
    nrow(de_tested)
  ),
  "",
  sprintf(
    "**%d of %d tested genes have q < %.2f** (%d higher and %d lower in %s).",
    n_significant,
    nrow(de_tested),
    fdr_cutoff,
    n_significant_up,
    n_significant - n_significant_up,
    condition_levels[[2]]
  ),
  "",
  "## Top results",
  "",
  sprintf(
    "The %d genes with the smallest p-values; the full table is in `de.tsv`.",
    nrow(top_table)
  ),
  "",
  md_table(top_table),
  "",
  "![Per-donor expression of the top genes](top_genes.png)",
  "",
  sprintf(
    "**Figure 1.** Expression of the %d genes with the smallest p-values. Each point is one donor (%s); grey bars mark the condition mean. Expression is log2 counts per million after TMM normalisation (prior count 2). Panels share the y axis. Source data: `top_genes_source.tsv`.",
    length(top_genes),
    paste(n_group, collapse = ", ")
  ),
  "",
  "## Limitations",
  "",
  sprintf(
    "- %s donors per condition: residual degrees of freedom are few, so power is limited and per-gene dispersion leans heavily on the shared prior.",
    paste(n_donors_by_condition$n, collapse = " and ")
  ),
  sprintf(
    "- TMM assumes most genes are not differentially expressed. With a %d-gene panel, a shift shared by many genes can be partly absorbed by the normalisation factors; check the factors in the donor table.",
    nrow(counts)
  ),
  "- Cells were pooled across cell types within each donor; a difference in cell-type composition between conditions would appear as differential expression."
)
writeLines(report, file.path(out_dir, "report.md"))

# Reproducibility ----
script_path <- str_remove(str_subset(commandArgs(trailingOnly = FALSE), "^--file="), "^--file=")
command_words <- c("Rscript", script_path, args)
command_words <- if_else(
  str_detect(command_words, "^[A-Za-z0-9_./=-]+$"),
  command_words,
  shQuote(command_words)
)
writeLines(
  c(
    "#!/usr/bin/env bash",
    "# Command that produced this run; rerun with a new --output (existing directories are refused).",
    paste("cd", shQuote(getwd())),
    paste(command_words, collapse = " ")
  ),
  file.path(out_dir, "reproducibility", "commands.sh")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "reproducibility", "session_info.txt"))
