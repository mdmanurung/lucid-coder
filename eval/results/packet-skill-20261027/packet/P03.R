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
#
# Usage: Rscript skills/pseudobulk-de/pseudobulk_de.R --input <dir> --output <out_dir>
#
# <dir> must contain counts.tsv (genes x cells, first column `gene`) and cells.tsv
# (`cell_id`, `donor_id`, `condition` in {healthy, disease}). <out_dir> must not exist.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(jsonlite)
})

SKILL_VERSION <- "1.0.0"
MIN_CELLS_PER_DONOR <- 10
MIN_DONORS_PER_CONDITION <- 2
FDR_THRESHOLD <- 0.05
N_TOP <- 10
LEVELS <- c("healthy", "disease") # first level is the reference
COLORS <- c(healthy = "#2a78d6", disease = "#eb6834")
INK <- "#0b0b0b"
INK_2 <- "#52514e"
MUTED <- "#b5b4ae"
SURFACE <- "#fcfcfb"

usage <- "Usage: Rscript pseudobulk_de.R --input <dir> --output <out_dir>"

die <- function(...) {
  message("ERROR: ", ...)
  quit(save = "no", status = 1)
}

# ---- Arguments ---------------------------------------------------------------

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(save = "no", status = 0)
    }
    if (!a %in% c("--input", "--output")) die("unknown argument '", a, "'\n", usage)
    if (i == length(args) || startsWith(args[i + 1], "--")) die(a, " needs a value\n", usage)
    out[[sub("^--", "", a)]] <- args[i + 1]
    i <- i + 2
  }
  for (k in c("input", "output")) if (is.null(out[[k]])) die("--", k, " is required\n", usage)
  out
}

sh_quote <- function(x) ifelse(grepl("^[A-Za-z0-9_./=:+-]+$", x), x, shQuote(x))

raw_args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(raw_args)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
command_line <- paste(c("Rscript", sh_quote(script_path), sh_quote(raw_args)), collapse = " ")
work_dir <- getwd()

if (file.exists(opt$output)) {
  die("output path '", opt$output, "' already exists; refusing to overwrite. ",
      "Choose a new --output directory.")
}
counts_path <- file.path(opt$input, "counts.tsv")
cells_path <- file.path(opt$input, "cells.tsv")
for (p in c(counts_path, cells_path)) if (!file.exists(p)) die("missing input file '", p, "'")

# ---- Read and validate -------------------------------------------------------

counts_df <- read_tsv(counts_path, col_types = cols(.default = col_character()),
                      name_repair = "minimal", progress = FALSE)
if (names(counts_df)[1] != "gene") die("first column of counts.tsv must be 'gene'")
genes <- counts_df$gene
if (anyNA(genes) || any(genes == "")) die("counts.tsv has empty gene IDs")
if (anyDuplicated(genes)) die("duplicated gene IDs in counts.tsv: ",
                              paste(unique(genes[duplicated(genes)]), collapse = ", "))
cell_cols <- names(counts_df)[-1]
if (length(cell_cols) == 0) die("counts.tsv has no cell columns")
if (anyDuplicated(cell_cols)) die("duplicated cell IDs in counts.tsv header")

counts <- suppressWarnings(apply(as.matrix(counts_df[, -1]), 2, as.numeric))
counts <- matrix(counts, nrow = length(genes), dimnames = list(genes, cell_cols))
if (anyNA(counts)) die("counts.tsv contains missing or non-numeric values")
if (any(counts < 0) || any(counts != round(counts))) {
  die("counts.tsv must hold raw UMI counts (non-negative integers)")
}

cells <- read_tsv(cells_path, col_types = cols(.default = col_character()), progress = FALSE)
missing_cols <- setdiff(c("cell_id", "donor_id", "condition"), names(cells))
if (length(missing_cols)) die("cells.tsv is missing column(s): ", paste(missing_cols, collapse = ", "))
cells <- cells %>% select(cell_id, donor_id, condition)
if (any(is.na(as.matrix(cells)) | as.matrix(cells) == "")) die("cells.tsv has empty values")
if (anyDuplicated(cells$cell_id)) die("duplicated cell_id in cells.tsv")
bad_cond <- setdiff(unique(cells$condition), LEVELS)
if (length(bad_cond)) die("condition must be 'healthy' or 'disease'; found: ", paste(bad_cond, collapse = ", "))

only_counts <- setdiff(cell_cols, cells$cell_id)
only_meta <- setdiff(cells$cell_id, cell_cols)
if (length(only_counts) || length(only_meta)) {
  die(length(only_counts), " cell(s) in counts.tsv lack metadata and ", length(only_meta),
      " cell(s) in cells.tsv lack counts (e.g. ", paste(head(c(only_counts, only_meta), 3), collapse = ", "),
      "); cell IDs must match exactly")
}
cells <- cells[match(cell_cols, cells$cell_id), ]

donor_conditions <- cells %>% distinct(donor_id, condition) %>% count(donor_id) %>% filter(n > 1)
if (nrow(donor_conditions)) {
  die("donor(s) assigned to more than one condition: ", paste(donor_conditions$donor_id, collapse = ", "))
}

n_genes_input <- length(genes)
n_cells_input <- ncol(counts)
n_donors_input <- n_distinct(cells$donor_id)

# ---- Filters -------------------------------------------------------------------

filter_log <- list()
log_filter <- function(step, unit, rule, before, after, removed_ids = character()) {
  filter_log[[length(filter_log) + 1]] <<- list(
    step = step, unit = unit, rule = rule, before = before, removed = before - after,
    kept = after, removed_ids = as.list(removed_ids)
  )
}

# 1. Cells with no UMIs in the panel carry no information and would inflate cell counts.
cell_total <- colSums(counts)
keep_cell <- cell_total > 0
log_filter("empty cells", "cells", "total UMI count = 0", length(keep_cell), sum(keep_cell))
counts <- counts[, keep_cell, drop = FALSE]
cells <- cells[keep_cell, ]

# 2. Donors with too few cells give noisy pseudobulk profiles.
cells_per_donor <- cells %>% count(donor_id, condition, name = "n_cells")
small_donors <- cells_per_donor$donor_id[cells_per_donor$n_cells < MIN_CELLS_PER_DONOR]
log_filter("low-cell donors", "donors", paste0("fewer than ", MIN_CELLS_PER_DONOR, " cells"),
           nrow(cells_per_donor), nrow(cells_per_donor) - length(small_donors), small_donors)
keep_cell <- !cells$donor_id %in% small_donors
counts <- counts[, keep_cell, drop = FALSE]
cells <- cells[keep_cell, ]

donors_per_condition <- cells %>% distinct(donor_id, condition) %>% count(condition)
for (lv in LEVELS) {
  n_lv <- sum(donors_per_condition$n[donors_per_condition$condition == lv])
  if (n_lv < MIN_DONORS_PER_CONDITION) {
    die("only ", n_lv, " ", lv, " donor(s) after filtering; need at least ",
        MIN_DONORS_PER_CONDITION, " per condition (donors are the replicates)")
  }
}

# ---- Aggregate to donor pseudobulk -------------------------------------------

donors <- cells %>%
  count(donor_id, condition, name = "n_cells") %>%
  mutate(condition = factor(condition, levels = LEVELS)) %>%
  arrange(condition, donor_id)
pseudobulk <- t(rowsum(t(counts), group = cells$donor_id))[, donors$donor_id, drop = FALSE]

y <- DGEList(counts = pseudobulk, group = donors$condition)

# 3. Genes too lowly expressed in the pseudobulk samples to test.
keep_gene <- filterByExpr(y, group = donors$condition)
log_filter("low-expression genes", "genes",
           "edgeR::filterByExpr defaults (min.count = 10, min.total.count = 15, by condition)",
           length(keep_gene), sum(keep_gene), names(keep_gene)[!keep_gene])
if (sum(keep_gene) < 2) die("fewer than 2 genes pass the expression filter")
y <- y[keep_gene, , keep.lib.sizes = FALSE]

# ---- Model ---------------------------------------------------------------------

y <- normLibSizes(y, method = "TMM")
design <- model.matrix(~condition, data = donors)
rownames(design) <- donors$donor_id
coef_name <- "conditiondisease"

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = coef_name)
tested <- topTags(qlf, n = Inf, sort.by = "none")$table
disp <- estimateDisp(y, design) # diagnostic only: between-donor variability (BCV)

# BH over the full gene family: genes removed by the expression filter count as untested
# (p = 1), so q-values are never smaller than BH over the tested subset alone.
de <- data.frame(gene = genes, log2_fc = NA_real_, p_value = NA_real_)
idx <- match(rownames(tested), de$gene)
de$log2_fc[idx] <- tested$logFC
de$p_value[idx] <- tested$PValue
de$q_value <- NA_real_
de$q_value[idx] <- p.adjust(tested$PValue, method = "BH", n = n_genes_input)
de <- de[order(de$p_value, na.last = TRUE), ]
rownames(de) <- NULL
sig <- de %>% filter(!is.na(q_value), q_value < FDR_THRESHOLD)

# Sensitivity: library-size normalization only (no TMM composition correction).
y_ls <- y
y_ls$samples$norm.factors <- 1
q_ls <- p.adjust(glmQLFTest(glmQLFit(y_ls, design), coef = coef_name)$table$PValue,
                 method = "BH", n = n_genes_input)
sig_ls <- rownames(y_ls)[q_ls < FDR_THRESHOLD]

# ---- Write outputs -------------------------------------------------------------

if (file.exists(opt$output)) die("output path '", opt$output, "' appeared during the run; refusing to overwrite")
dir.create(opt$output, recursive = TRUE)
dir.create(file.path(opt$output, "reproducibility"))
dir.create(file.path(opt$output, "figures"))
out <- function(...) file.path(opt$output, ...)

write_tsv(de, out("de.tsv"), na = "NA")

writeLines(c(
  "#!/usr/bin/env bash",
  "# Exact command that produced the outputs in this directory.",
  paste0("# Working directory: ", work_dir),
  paste0("# Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "# The script refuses an existing --output directory; change --output to re-run.",
  command_line
), out("reproducibility", "commands.sh"))
Sys.chmod(out("reproducibility", "commands.sh"), mode = "0755")
writeLines(capture.output(sessionInfo()), out("reproducibility", "sessionInfo.txt"))

# Figures
theme_report <- theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = SURFACE, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
    text = element_text(colour = INK),
    axis.text = element_text(colour = INK_2),
    plot.title = element_text(face = "bold"),
    legend.position = "top", legend.justification = "left"
  )

mds <- plotMDS(y, plot = FALSE)
var_pct <- round(100 * mds$var.explained[1:2])
mds_df <- data.frame(donor_id = colnames(y), condition = donors$condition, x = mds$x, y = mds$y)
p_mds <- ggplot(mds_df, aes(x, y, colour = condition)) +
  geom_point(size = 3.2) +
  geom_text(aes(label = donor_id), colour = INK_2, size = 3.2, vjust = -1.1, show.legend = FALSE) +
  scale_colour_manual(values = COLORS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(title = "Pseudobulk samples (one point per donor)",
       x = paste0("Leading log2 fold change, dim 1 (", var_pct[1], "%)"),
       y = paste0("Leading log2 fold change, dim 2 (", var_pct[2], "%)")) +
  theme_report
ggsave(out("figures", "mds.png"), p_mds, width = 6, height = 4.5, dpi = 150, bg = SURFACE)

vol_df <- de %>%
  filter(!is.na(p_value)) %>%
  mutate(significant = !is.na(q_value) & q_value < FDR_THRESHOLD, nlp = -log10(p_value))
p_vol <- ggplot(vol_df, aes(log2_fc, nlp)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_point(data = filter(vol_df, !significant), colour = MUTED, size = 2.4) +
  geom_point(data = filter(vol_df, significant), colour = COLORS[["disease"]], size = 2.8) +
  geom_text(data = filter(vol_df, significant), aes(label = gene), colour = INK,
            size = 3.2, hjust = 1.15) +
  labs(title = paste0("Disease vs healthy: ", nrow(sig), " of ", n_genes_input,
                      " genes at q < ", FDR_THRESHOLD),
       x = "log2 fold change (disease / healthy)", y = "-log10 p-value") +
  theme_report
ggsave(out("figures", "volcano.png"), p_vol, width = 6, height = 4.5, dpi = 150, bg = SURFACE)

# result.json
pkg_versions <- sapply(c("edgeR", "limma", "dplyr", "readr", "ggplot2", "jsonlite"),
                       function(p) as.character(packageVersion(p)), simplify = FALSE)
n_tested <- nrow(tested)
result <- list(
  skill = "pseudobulk-de",
  version = SKILL_VERSION,
  command = command_line,
  parameters = list(
    input = opt$input,
    output = opt$output,
    replication_unit = "donor_id",
    aggregation = "sum of raw UMI counts over all cells of a donor",
    min_cells_per_donor = MIN_CELLS_PER_DONOR,
    min_donors_per_condition = MIN_DONORS_PER_CONDITION,
    gene_filter = "edgeR::filterByExpr(group = condition), default thresholds",
    normalization = "TMM (edgeR::normLibSizes)",
    method = "edgeR quasi-likelihood negative binomial GLM (glmQLFit + glmQLFTest)",
    design = "~ condition",
    reference_level = LEVELS[1],
    contrast = "disease - healthy (log2 scale)",
    multiple_testing = "Benjamini-Hochberg",
    bh_family_size = n_genes_input,
    fdr_threshold = FDR_THRESHOLD
  ),
  counts = list(
    genes_input = n_genes_input,
    genes_tested = n_tested,
    cells_input = n_cells_input,
    cells_used = ncol(counts),
    donors_input = n_donors_input,
    donors_used = nrow(donors),
    donors_per_condition = as.list(table(donors$condition)),
    cells_per_donor = as.list(setNames(donors$n_cells, donors$donor_id)),
    filters = filter_log
  ),
  n_significant_q05 = sum(de$q_value < 0.05, na.rm = TRUE),
  significant_genes = sig$gene,
  diagnostics = list(
    ql_prior_df = unname(median(fit$df.prior)),
    common_bcv = sqrt(disp$common.dispersion),
    tmm_norm_factors = as.list(setNames(y$samples$norm.factors, colnames(y))),
    sensitivity_no_tmm_n_significant = length(sig_ls)
  ),
  software = c(list(R = R.version.string), pkg_versions)
)
write_json(result, out("result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ---- report.md -----------------------------------------------------------------

md_table <- function(df) {
  df <- as.data.frame(lapply(df, as.character), check.names = FALSE)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
fmt_p <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

n_by <- table(donors$condition)
filter_tbl <- bind_rows(lapply(filter_log, function(f) data.frame(
  Step = f$step, Unit = f$unit, Rule = f$rule, Before = f$before, Removed = f$removed,
  Kept = f$kept,
  `Removed IDs` = if (length(f$removed_ids)) paste(unlist(f$removed_ids), collapse = ", ") else "none",
  check.names = FALSE
)))
donor_tbl <- data.frame(
  Donor = donors$donor_id, Condition = as.character(donors$condition),
  Cells = donors$n_cells, `Pseudobulk UMIs` = fmt_n(colSums(pseudobulk)),
  `TMM factor` = sprintf("%.3f", y$samples$norm.factors), check.names = FALSE
)
logcpm <- tested$logCPM[match(de$gene, rownames(tested))]
top_tbl <- head(data.frame(
  Rank = seq_len(nrow(de)), Gene = de$gene, `log2 FC` = sprintf("%.3f", de$log2_fc),
  `Avg log2 CPM` = sprintf("%.2f", logcpm), `p-value` = fmt_p(de$p_value),
  `q-value` = fmt_p(de$q_value), check.names = FALSE
), N_TOP)

sig_text <- if (nrow(sig)) {
  paste0(nrow(sig), " of ", n_genes_input, " genes differ at q < ", FDR_THRESHOLD, ": ",
         paste0(sig$gene, " (log2 FC ", sprintf("%+.2f", sig$log2_fc), ", q = ", fmt_p(sig$q_value), ")",
                collapse = "; "), ".")
} else {
  paste0("No gene differs at q < ", FDR_THRESHOLD, " (", n_genes_input, " genes tested).")
}
sens_text <- if (setequal(sig_ls, sig$gene)) {
  "the same set of significant genes"
} else {
  paste0(length(sig_ls), " significant gene(s) (", paste(sig_ls, collapse = ", "), ")")
}

report <- c(
  "# Pseudobulk differential expression: disease vs healthy",
  "",
  paste0("Generated by `pseudobulk-de` v", SKILL_VERSION, " on ", format(Sys.Date()), " with:"),
  "",
  "```bash", command_line, "```",
  "",
  "## Result",
  "",
  sig_text,
  "",
  "## Data",
  "",
  md_table(data.frame(
    Quantity = c("Donors", "Cells", "Genes"),
    Input = c(paste0(n_donors_input, " (", n_by[["healthy"]], " healthy, ", n_by[["disease"]], " disease)"),
              fmt_n(n_cells_input), n_genes_input),
    `Used in model` = c(nrow(donors), fmt_n(ncol(counts)), n_tested),
    check.names = FALSE
  )),
  "",
  paste0("Input: `", counts_path, "` (raw UMI counts, genes x cells) and `", cells_path,
         "` (cell -> donor -> condition). Every cell ID matched between the two files and every ",
         "donor belongs to exactly one condition."),
  "",
  "## Replication unit",
  "",
  paste0("The **donor** is the unit of replication: ", n_by[["healthy"]], " healthy vs ",
         n_by[["disease"]], " disease donors. Condition is assigned per donor, so cells from ",
         "one donor are not independent observations of the condition effect. Testing cells as ",
         "replicates ignores donor-to-donor variation and yields many false positives; the test ",
         "here has ", nrow(donors) - ncol(design), " residual degrees of freedom per gene, not ",
         fmt_n(ncol(counts) - ncol(design)), "."),
  "",
  "## Aggregation",
  "",
  paste0("For each donor and gene, raw UMI counts were **summed over all of that donor's cells** ",
         "(pseudobulk), giving a ", n_tested, "-gene x ", nrow(donors), "-donor count matrix. Sums ",
         "(not means) keep the data as counts, which the negative binomial model requires."),
  "",
  md_table(donor_tbl),
  "",
  "## Filters",
  "",
  "Filters were applied in this order. Each row reports what it removed, including zero.",
  "",
  md_table(filter_tbl),
  "",
  paste0("After filtering, each condition was required to keep at least ", MIN_DONORS_PER_CONDITION,
         " donors; this held."),
  "",
  "## Model",
  "",
  paste0("**edgeR quasi-likelihood (QL) negative binomial GLM** (`glmQLFit` + `glmQLFTest`, edgeR ",
         packageVersion("edgeR"), "), design `~ condition` with `healthy` as reference."),
  "",
  "Why this method:",
  "",
  "- Pseudobulk counts are overdispersed counts; the negative binomial GLM models them directly.",
  paste0("- With only ", n_by[["healthy"]], " vs ", n_by[["disease"]], " donors, per-gene variance ",
         "estimates are unstable. The QL F-test shrinks gene-wise dispersions towards a common trend ",
         "(empirical Bayes) and accounts for the uncertainty in the estimate, which controls the ",
         "type I error better than likelihood-ratio tests at small sample sizes."),
  "- Pseudobulk + edgeR/limma has been shown to control false discoveries in multi-donor scRNA-seq, where cell-level tests do not (Squair et al. 2021; Crowell et al. 2020).",
  "",
  paste0("**Normalization:** TMM (`normLibSizes`) on the pseudobulk library sizes. With a ",
         n_genes_input, "-gene panel, one strongly changed, abundant gene shifts every other gene's ",
         "share of the library, so composition correction matters. TMM assumes most genes are not ",
         "differentially expressed. Factors ranged ", sprintf("%.3f", min(y$samples$norm.factors)),
         "-", sprintf("%.3f", max(y$samples$norm.factors)), "."),
  "",
  paste0("**Dispersion diagnostics:** common biological coefficient of variation (BCV) between donors = ",
         sprintf("%.3f", sqrt(disp$common.dispersion)), "; QL prior degrees of freedom = ",
         format(round(median(fit$df.prior), 1), big.mark = ","),
         ". A large prior df means the genes' QL dispersions are nearly identical, so each gene's test ",
         "borrows almost all of its variance information from the other genes. This is why p-values ",
         "can be very small with only ", n_by[["healthy"]], " vs ", n_by[["disease"]], " donors. ",
         "Read extreme p-values as \"very strong evidence\" under the model, not as exact ",
         "probabilities. For comparison, a distribution-free donor-label permutation test cannot ",
         "go below p = 1/", choose(nrow(donors), n_by[["healthy"]]), " (one-sided)."),
  "",
  "## Contrast",
  "",
  paste0("Coefficient `", coef_name, "`: **disease - healthy** on the log2 scale. `log2_fc` > 0 means ",
         "higher expression in disease. log2 fold changes are edgeR GLM estimates on TMM-normalized ",
         "counts (default prior count 0.125)."),
  "",
  "## Multiple testing",
  "",
  paste0("Benjamini-Hochberg over the **full family of ", n_genes_input, " input genes** (",
         n_tested, " tested; any gene removed by the expression filter counts as untested with p = 1, ",
         "which can only make q-values larger). Significance threshold: q < ", FDR_THRESHOLD, "."),
  "",
  "## Top results",
  "",
  paste0("Top ", min(N_TOP, nrow(de)), " genes by p-value (all genes are in `de.tsv`):"),
  "",
  md_table(top_tbl),
  "",
  "![MDS plot](figures/mds.png)",
  "",
  paste0("**Figure 1. Multidimensional scaling of the ", nrow(donors), " donor pseudobulk samples.** ",
         "Distances are the root-mean-square log2 fold change between samples over the ", n_tested,
         " tested genes (log2 CPM after TMM). Color marks condition (blue = healthy, orange = disease); ",
         "each point is labeled with its donor. Dimension 1 explains ", var_pct[1], "% of the variance. ",
         "Separation along it checks that condition, not a single outlying donor, dominates the data."),
  "",
  "![Volcano plot](figures/volcano.png)",
  "",
  paste0("**Figure 2. Volcano plot of disease vs healthy.** Each point is one gene: x = log2 fold ",
         "change (disease / healthy), y = -log10 QL F-test p-value. Orange, labeled points have ",
         "BH q < ", FDR_THRESHOLD, "; gray points do not."),
  "",
  "## Sensitivity check",
  "",
  paste0("Re-running the same model with library-size normalization only (no TMM) gives ", sens_text,
         " at q < ", FDR_THRESHOLD, "."),
  "",
  "## Limitations",
  "",
  paste0("- ", n_by[["healthy"]], " vs ", n_by[["disease"]], " donors gives low power for modest effects; ",
         "a non-significant gene is not evidence of no difference."),
  "- No covariates (sex, age, batch) were available; any donor-level confounder is aliased with condition.",
  "- All cells were pooled per donor; the analysis assumes the cells are one cell type. A shift in cell-type composition between conditions would appear as differential expression.",
  "",
  "## Outputs",
  "",
  "- `de.tsv`: gene, log2_fc, p_value, q_value (one row per input gene, sorted by p-value)",
  "- `result.json`: parameters, counts, filter log, number significant at q < 0.05",
  "- `figures/mds.png`, `figures/volcano.png`",
  "- `reproducibility/commands.sh`, `reproducibility/sessionInfo.txt`",
  "",
  "## Software",
  "",
  paste0("- ", R.version.string),
  paste0("- ", names(pkg_versions), " ", unlist(pkg_versions))
)
writeLines(report, out("report.md"))

message("Done: ", nrow(sig), " of ", n_genes_input, " genes significant at q < ", FDR_THRESHOLD,
        ". Outputs in ", opt$output)
