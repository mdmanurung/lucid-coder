#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====
#$ q1 one read (1-5; can you say what this script does after one read?): 
#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): 
#$ q3 change checks (yes/no; would you add, remove or change a check?): 
#$ q3 which (one sentence): 
#$ q4 procedure sound (yes/no; design, replication unit, method, filters, contrast): 
#$ q4 which (one sentence): 
#$ ---- your comments: #$ for code, #$P for procedure; own line above the code, or after it ----
#!/usr/bin/env Rscript
# pseudobulk-de: donor-level pseudobulk differential expression, disease vs healthy.
# Usage: Rscript skills/pseudobulk-de/pseudobulk_de.R --input <dir> --output <new_dir> [--min-cells 10]

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(jsonlite)
})

usage <- "Usage: Rscript pseudobulk_de.R --input <dir with counts.tsv, cells.tsv> --output <new dir> [--min-cells 10]"
die <- function(...) {
  message("Error: ", ...)
  quit(save = "no", status = 1)
}

# ---- arguments ----
args <- commandArgs(trailingOnly = TRUE)
opt <- list(input = NULL, output = NULL, min_cells = "10")
i <- 1
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!startsWith(args[i], "--") || !key %in% names(opt) || i == length(args)) {
    die("unrecognised or incomplete argument '", args[i], "'\n", usage)
  }
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
if (is.null(opt$input) || is.null(opt$output)) die("--input and --output are both required\n", usage)
min_cells <- suppressWarnings(as.integer(opt$min_cells))
if (is.na(min_cells) || min_cells < 1) die("--min-cells must be a positive integer")
if (file.exists(opt$output)) die("output path already exists, refusing to overwrite: ", opt$output)

ref_level <- "healthy"
test_level <- "disease"
fdr <- 0.05

# ---- read and validate inputs ----
counts_path <- file.path(opt$input, "counts.tsv")
cells_path <- file.path(opt$input, "cells.tsv")
for (p in c(counts_path, cells_path)) if (!file.exists(p)) die("missing input file: ", p)

raw <- read.delim(counts_path, check.names = FALSE, colClasses = c(gene = "character"))
cells <- read.delim(cells_path, colClasses = "character", na.strings = c("", "NA"))

if (colnames(raw)[1] != "gene") die("first column of counts.tsv must be 'gene'")
miss <- setdiff(c("cell_id", "donor_id", "condition"), colnames(cells))
if (length(miss)) die("cells.tsv lacks column(s): ", paste(miss, collapse = ", "))
if (anyDuplicated(raw$gene)) die("duplicated gene names in counts.tsv")
if (anyDuplicated(colnames(raw))) die("duplicated cell IDs in counts.tsv header")
if (anyDuplicated(cells$cell_id)) die("duplicated cell_id in cells.tsv")

m <- as.matrix(raw[, -1, drop = FALSE])
rownames(m) <- raw$gene
if (!is.numeric(m) || anyNA(m) || any(m < 0) || any(m != round(m))) {
  die("counts.tsv must hold non-negative integer counts with no missing values")
}
bad <- setdiff(na.omit(cells$condition), c(ref_level, test_level))
if (length(bad)) die("unexpected condition value(s): ", paste(bad, collapse = ", "))
n_genes_in <- nrow(m)
n_cells_in <- ncol(m)

# ---- cell filter 1: cells with complete metadata ----
meta <- cells[match(colnames(m), cells$cell_id), ]
has_meta <- !is.na(meta$cell_id) & !is.na(meta$donor_id) & !is.na(meta$condition)
m <- m[, has_meta, drop = FALSE]
meta <- meta[has_meta, ]
n_cells_no_meta <- sum(!has_meta)
n_meta_rows_unused <- sum(!cells$cell_id %in% colnames(m))

donor_cond <- tapply(meta$condition, meta$donor_id, unique, simplify = FALSE)
mixed <- names(donor_cond)[lengths(donor_cond) > 1]
if (length(mixed)) die("donor(s) appear under both conditions: ", paste(mixed, collapse = ", "))
donor_cond <- unlist(donor_cond)
n_donors_in <- length(donor_cond)

# ---- donor filter 2: minimum cells per donor ----
cells_per_donor <- table(meta$donor_id)
low <- names(cells_per_donor)[cells_per_donor < min_cells]
keep_cell <- !meta$donor_id %in% low
n_cells_low_donor <- sum(!keep_cell)
m <- m[, keep_cell, drop = FALSE]
meta <- meta[keep_cell, ]
donor_cond <- donor_cond[!names(donor_cond) %in% low]
n_per_cond <- table(factor(donor_cond, levels = c(ref_level, test_level)))
if (any(n_per_cond < 2)) {
  die("need >= 2 donors per condition after filtering; have healthy=", n_per_cond[[ref_level]],
      ", disease=", n_per_cond[[test_level]])
}

# ---- aggregate: sum raw counts per donor ----
pb <- t(rowsum(t(m), meta$donor_id))
group <- factor(donor_cond[colnames(pb)], levels = c(ref_level, test_level))

# ---- gene filter 3 + edgeR quasi-likelihood ----
y <- DGEList(pb, group = group)
keep_gene <- filterByExpr(y, group = group)
genes_removed <- rownames(y)[!keep_gene]
y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = paste0("group", test_level))
tt <- topTags(qlf, n = Inf, sort.by = "none")$table

de <- data.frame(gene = rownames(pb), log2_fc = NA_real_, p_value = NA_real_)
idx <- match(rownames(tt), de$gene)
de$log2_fc[idx] <- tt$logFC
de$p_value[idx] <- tt$PValue
de$q_value <- p.adjust(de$p_value, method = "BH") # NAs (filtered genes) are outside the family
de <- de[order(de$p_value, na.last = TRUE), ]
sig <- !is.na(de$q_value) & de$q_value < fdr
n_sig <- sum(sig)
n_up <- sum(sig & de$log2_fc > 0)
n_down <- sum(sig & de$log2_fc < 0)

# ---- write outputs ----
out <- opt$output
dir.create(file.path(out, "reproducibility"), recursive = TRUE)
write.table(de, file.path(out, "de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

tested <- de[!is.na(de$p_value), ]
tested$label <- ifelse(tested$q_value < fdr, tested$gene, "")
tested$class <- ifelse(tested$q_value < fdr, ifelse(tested$log2_fc > 0, "Up in disease", "Down in disease"),
                       "Not significant")
volcano <- ggplot(tested, aes(log2_fc, -log10(p_value), colour = class)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 2) +
  geom_text(aes(label = label), size = 3, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_colour_manual(values = c("Up in disease" = "#C0392B", "Down in disease" = "#2471A3",
                                 "Not significant" = "grey60"), name = NULL) +
  labs(x = "log2 fold change (disease / healthy)", y = "-log10 p-value") +
  theme_bw()
ggsave(file.path(out, "volcano.png"), volcano, width = 6.5, height = 4.5, dpi = 150)

q <- function(x) ifelse(grepl("^[A-Za-z0-9_./=:@%+-]+$", x), x, shQuote(x))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
cmd <- paste(q(c("Rscript", script_path, args)), collapse = " ")
writeLines(c("#!/usr/bin/env bash",
             paste("# Generated by pseudobulk-de on", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("# Working directory:", getwd()),
             cmd),
           file.path(out, "reproducibility", "commands.sh"))

pkgs <- c("edgeR", "limma", "ggplot2", "jsonlite")
params <- list(
  input = opt$input, output = out, min_cells_per_donor = min_cells,
  replication_unit = "donor", aggregation = "sum of raw UMI counts per donor",
  gene_filter = "edgeR::filterByExpr(group = condition), default thresholds",
  normalisation = "TMM", design = "~ condition", reference_level = ref_level,
  test_level = test_level, test = "edgeR quasi-likelihood F-test (glmQLFit, robust = TRUE)",
  multiple_testing = "Benjamini-Hochberg over all genes passing the gene filter",
  significance_threshold_q = fdr
)
counts <- list(
  genes_input = n_genes_in, genes_removed_by_filter = length(genes_removed), genes_tested = nrow(tt),
  cells_input = n_cells_in, cells_removed_missing_metadata = n_cells_no_meta,
  cells_removed_low_cell_donors = n_cells_low_donor, cells_used = ncol(m),
  donors_input = n_donors_in, donors_removed_low_cells = length(low), donors_used = ncol(pb),
  donors_healthy = n_per_cond[[ref_level]], donors_disease = n_per_cond[[test_level]]
)
write_json(list(skill = "pseudobulk-de", parameters = params, counts = counts,
                n_significant_q05 = n_sig, n_significant_q05_up = n_up, n_significant_q05_down = n_down,
                software = c(list(R = R.version.string),
                             setNames(lapply(pkgs, \(p) as.character(packageVersion(p))), pkgs))),
           file.path(out, "result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ---- report ----
md_table <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, \(r) paste("|", paste(r, collapse = " | "), "|")))
}
fmt_p <- function(x) formatC(x, format = "g", digits = 3)
donor_tab <- data.frame(
  donor = colnames(pb), condition = as.character(group),
  cells = as.integer(cells_per_donor[colnames(pb)]),
  pseudobulk_UMIs = format(colSums(pb), big.mark = ","),
  TMM_factor = sprintf("%.3f", y$samples$norm.factors)
)
top <- head(tested, 10)
top_tab <- data.frame(gene = top$gene, log2_fc = sprintf("%.2f", top$log2_fc),
                      p_value = fmt_p(top$p_value), q_value = fmt_p(top$q_value))
res_df <- fit$df.residual[1]

report <- c(
  "# Pseudobulk differential expression: disease vs healthy", "",
  paste0("Generated by the `pseudobulk-de` skill on ", format(Sys.Date()), ". Command: `", cmd, "`"), "",
  "## Data summary", "",
  md_table(data.frame(
    item = c("Donors", "Cells", "Genes"),
    input = c(n_donors_in, n_cells_in, n_genes_in),
    used_in_model = c(ncol(pb), ncol(m), nrow(tt))
  )), "",
  sprintf("Donors per condition in the model: %d healthy, %d disease.", n_per_cond[[ref_level]],
          n_per_cond[[test_level]]), "",
  md_table(donor_tab), "",
  "## Replication unit", "",
  paste("The **donor** is the unit of replication. Cells from the same donor share genotype, environment",
        "and sample handling, so they are not independent observations of the disease effect; treating",
        "cells as replicates (pseudoreplication) badly inflates false positives. All inference is",
        sprintf("therefore made across %d donors, leaving %d residual degrees of freedom.", ncol(pb), res_df)),
  "",
  "## Aggregation", "",
  paste("Raw UMI counts were **summed** over all cells of each donor, giving a genes x donors count",
        "matrix (one pseudobulk sample per donor). Sums, not means, keep the data on the count scale",
        "the negative-binomial model expects; differing cell numbers per donor are absorbed by the",
        "library-size normalisation. The input is treated as a single cell population: if it mixes",
        "cell types, run the skill once per cell type."), "",
  "## Filters", "",
  sprintf("1. **Cells without metadata** (no row in `cells.tsv`, or missing `donor_id`/`condition`): removed %d of %d cells. %d `cells.tsv` rows had no matching count column and were ignored.",
          n_cells_no_meta, n_cells_in, n_meta_rows_unused),
  sprintf("2. **Donors with fewer than %d cells**: removed %d donor(s)%s, %d cells.", min_cells, length(low),
          if (length(low)) paste0(" (", paste(low, collapse = ", "), ")") else "", n_cells_low_donor),
  sprintf("3. **Low-expression genes** (`edgeR::filterByExpr` on pseudobulk counts, grouped by condition, default thresholds `min.count = 10`, `min.total.count = 15`): removed %d of %d genes%s. Removed genes appear in `de.tsv` with `NA` statistics.",
          length(genes_removed), n_genes_in,
          if (length(genes_removed)) paste0(" (", paste(genes_removed, collapse = ", "), ")") else ""),
  "",
  "No cell-level QC (mitochondrial fraction, doublets, ambient RNA) is done here; it is assumed upstream.",
  "",
  "## Normalisation", "",
  "TMM normalisation factors (`edgeR::normLibSizes`) correct for pseudobulk library size and composition; see the TMM column above.",
  "",
  "## Model and why", "",
  paste("Negative-binomial generalised linear model with design `~ condition` (healthy as reference),",
        "fitted with edgeR's quasi-likelihood pipeline: `estimateDisp`, then `glmQLFit(robust = TRUE)`",
        "and a QL F-test. The QL F-test accounts for uncertainty in the dispersion estimates and gives",
        "the best type-I error control of the edgeR tests at small sample sizes; robust empirical Bayes",
        "shrinks gene-wise dispersions toward a trend, stabilising them with only a handful of donors.",
        "Pseudobulk + edgeR/limma is among the best-calibrated approaches for multi-donor single-cell",
        "DE in independent benchmarks (Crowell et al. 2020; Squair et al. 2021).",
        sprintf("Common biological coefficient of variation between donors: %.3f.", sqrt(y$common.dispersion))),
  "",
  "## Contrast", "",
  paste0("Coefficient `group", test_level, "`: `log2_fc` = log2(disease / healthy) of normalised expression. ",
         "Positive values mean higher expression in disease."), "",
  "## Multiple testing", "",
  sprintf(paste("Benjamini-Hochberg FDR over the family of all %d genes that passed the gene filter",
                "(%d of %d input genes). At q < %.2f, **%d gene(s) significant** (%d up, %d down in disease)."),
          nrow(tt), nrow(tt), n_genes_in, fdr, n_sig, n_up, n_down), "",
  "## Top results", "",
  sprintf("Top %d genes by p-value (full table: `de.tsv`).", nrow(top_tab)), "",
  md_table(top_tab), "",
  "![Volcano plot](volcano.png)", "",
  sprintf(paste("**Figure 1.** Volcano plot of the %d tested genes. x: log2 fold change, disease vs healthy",
                "(edgeR QL fit on donor pseudobulks, %d healthy vs %d disease donors); y: -log10 QL F-test",
                "p-value. Red/blue: q < %.2f (BH), higher/lower in disease, labelled by gene; grey: not",
                "significant."),
          nrow(tested), n_per_cond[[ref_level]], n_per_cond[[test_level]], fdr), "",
  "## Caveats", "",
  sprintf("- With %d vs %d donors, power is limited; non-significant genes are not evidence of no effect.",
          n_per_cond[[ref_level]], n_per_cond[[test_level]]),
  "- Each donor belongs to one condition, so any donor-level batch aligned with condition (processing date, site) cannot be separated from disease.",
  "- If cell-type composition differs between conditions, pseudobulk of a mixed population conflates expression change with composition shift.",
  "",
  "## Files", "",
  "- `de.tsv`: gene, log2_fc, p_value, q_value for every input gene",
  "- `result.json`: parameters, counts, number of significant genes, software versions",
  "- `volcano.png`: Figure 1",
  "- `reproducibility/commands.sh`: command that produced these outputs"
)
writeLines(report, file.path(out, "report.md"))
message(sprintf("Done: %d genes tested, %d significant at q < %.2f. Outputs in %s", nrow(tt), n_sig, fdr, out))
